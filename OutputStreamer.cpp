#include "OutputStreamer.h"
#include "ThreadPool.h"
#include "Benchmark.h"

#include "OutputFileStream.h"
#include <algorithm>
#include <map>

//**************************** asynchronous output writer ****************************
// One background thread writes full ofbufstream buffers. A stream queues its
// buffers while holding its output_mtx, so each file receives bytes in append
// order. Appenders only wait here for backpressure, never for write() or BGZF
// compression themselves.
class AsyncOutputWriter {
public:
	// nullptr once the writer has shut down (static destruction); callers then write inline.
	static AsyncOutputWriter* instance() {
		static AsyncOutputWriter writer;
		return writer.accepting_.load(std::memory_order_acquire) ? &writer : nullptr;
	}
	// Takes ownership of buffer only when it returns true.
	bool submit(ofbufstream* owner, OutputBuffer&& buffer, size_t size) {
		std::unique_lock<std::mutex> lock(mutex_);
		notFull_.wait(lock, [&] { return stop_ || jobs_.empty() || queuedBytes_ < kMaxQueuedBytes; });
		if (stop_) { return false; }
		queuedBytes_ += size;
		jobs_.push_back(Job{owner, std::move(buffer), size});
		notEmpty_.notify_one();
		return true;
	}
	AsyncOutputWriter(const AsyncOutputWriter&) = delete;
	AsyncOutputWriter& operator=(const AsyncOutputWriter&) = delete;
	~AsyncOutputWriter() {
		accepting_.store(false, std::memory_order_release);
		{
			std::lock_guard<std::mutex> lock(mutex_);
			stop_ = true;
		}
		notEmpty_.notify_all();
		notFull_.notify_all();
		if (thread_.joinable()) { thread_.join(); }
	}

private:
	struct Job {
		ofbufstream* owner;
		OutputBuffer buffer;
		size_t size;
	};
	static constexpr size_t kMaxQueuedBytes = size_t(16) << 20;
	AsyncOutputWriter() : thread_([this] { run(); }) { accepting_.store(true, std::memory_order_release); }
	void run() {
		std::unique_lock<std::mutex> lock(mutex_);
		for (;;) {
			notEmpty_.wait(lock, [&] { return stop_ || !jobs_.empty(); });
			if (jobs_.empty()) { return; }  // stopping and fully drained
			Job job = std::move(jobs_.front());
			jobs_.pop_front();
			lock.unlock();
			job.owner->writeQueued(job.buffer, job.size);
			// The owner may be destroyed as soon as its last buffer is released.
			job.owner->releaseQueued(std::move(job.buffer));
			lock.lock();
			queuedBytes_ -= job.size;
			notFull_.notify_all();
		}
	}
	std::mutex mutex_;
	std::condition_variable notEmpty_, notFull_;
	std::deque<Job> jobs_;
	size_t queuedBytes_ = 0;
	bool stop_ = false;
	std::atomic<bool> accepting_{false};
	std::thread thread_;
};

//**************************** ofbufstream ****************************
ofbufstream::ofbufstream(const string IF, int mif, bool isMC, size_t bufferS)
	: file(IF), keeper(bufferS), modeIO(mif), used(0),
	coutW(false), isGZ(false), doMC(isMC), primary(nullptr), bufS(bufferS) {
    if (bufS == 0) bufS = 5* 1024 * 1024;
	if (keeper.size() != bufS) {
		keeper.resize(bufS);
	}
	if (file == "" || file == "-") {
		coutW = true;
		return;
	}
	activate();
}

ofbufstream::~ofbufstream() {
	try { deactivate(); }
	catch (const std::exception& error) {
		std::cerr << "Failed to finalize " << file << ": " << error.what() << "\n";
	}
}

void ofbufstream::finishWrites() {
	emptyStream();
}

bool ofbufstream::operator! (void) {
	if (coutW) return false;
	if (!primary) return true;
	waitForQueuedWrites();
	return !(*primary);
}

void ofbufstream::operator<< (const string& X) {
	append(X.data(), X.size());
}

void ofbufstream::append(const char* data, size_t size) {
	if (size == 0) return;
	if (coutW) {
		Instr::TimedLockGuard guard(output_mtx); // worker threads may write concurrently
		cout.write(data, static_cast<std::streamsize>(size));
		return;
	}
	if (bufS == 0) {
       Instr::TimedLockGuard guard(output_mtx);
		if (primary) primary->write(data, static_cast<std::streamsize>(size));
		return;
	}
	Instr::TimedLockGuard guard(output_mtx);
	if (!asyncWriterEnabled || !primary) {
		if (used + size > keeper.size()) {
			if (primary && used > 0) {
				primary->write(keeper.data(), (streamsize)used);
			}
			used = 0;
		}

		if (size > keeper.size()) {
			if (primary) primary->write(data, (streamsize)size);
			return;
		}

		memcpy(keeper.data() + used, data, size);
		used += size;
		return;
	}
	// Fill the buffer completely and queue it; the byte stream is unchanged.
	while (size > 0) {
		if (used == keeper.size()) { handOffLocked(); }
		const size_t take = std::min(size, keeper.size() - used);
		memcpy(keeper.data() + used, data, take);
		used += take;
		data += take;
		size -= take;
	}
}

void ofbufstream::handOffLocked() {
	if (used == 0) { return; }
	OutputBuffer next;
	{
		std::lock_guard<std::mutex> lock(asyncMtx_);
		if (!spare_.empty()) {
			next = std::move(spare_.back());
			spare_.pop_back();
		}
		++inflight_;
	}
	if (next.size() != keeper.size()) { next.resize(keeper.size()); }
	std::swap(keeper, next); // next now holds the filled buffer
	const size_t size = used;
	used = 0;
	AsyncOutputWriter* writer = AsyncOutputWriter::instance();
	if (!writer || !writer->submit(this, std::move(next), size)) {
		writeQueued(next, size);
		releaseQueued(std::move(next));
	}
}

void ofbufstream::writeQueued(OutputBuffer& buffer, size_t size) {
	try {
		if (primary && size > 0) {
			primary->write(buffer.data(), static_cast<std::streamsize>(size));
		}
	}
	catch (...) {
		std::lock_guard<std::mutex> lock(asyncMtx_);
		if (!asyncError_) { asyncError_ = std::current_exception(); }
	}
}

void ofbufstream::releaseQueued(OutputBuffer&& buffer) {
	std::lock_guard<std::mutex> lock(asyncMtx_);
	if (spare_.size() < 2) { spare_.push_back(std::move(buffer)); }
	if (--inflight_ == 0) { asyncDone_.notify_all(); }
}

void ofbufstream::waitForQueuedWrites() {
	std::unique_lock<std::mutex> lock(asyncMtx_);
	asyncDone_.wait(lock, [this] { return inflight_ == 0; });
	if (asyncError_) {
		std::exception_ptr error = asyncError_;
		asyncError_ = nullptr;
		lock.unlock();
		std::rethrow_exception(error);
	}
}

void ofbufstream::emptyStream() {
	if (coutW) return;
	{
		Instr::TimedLockGuard guard(output_mtx);
		if (!asyncWriterEnabled || !primary) {
			if (!primary || used == 0) return;
			primary->write(keeper.data(), (streamsize)used);
			used = 0;
			return;
		}
		handOffLocked();
	}
	waitForQueuedWrites();
}

void ofbufstream::activate() {
	if (coutW) return;
	if (primary) return;
    string outFile = file;
	trim(outFile);
	if (outFile.size() >= 2) {
		const char first = outFile.front();
		const char last = outFile.back();
		if ((first == '"' && last == '"') || (first == '\'' && last == '\'')) {
			outFile = outFile.substr(1, outFile.size() - 2);
			trim(outFile);
		}
	}
	file = outFile;
	isGZ = isGZfile(file);
	ios_base::openmode m = ios::out;
	if (modeIO == ios::app) m |= ios::app;
	else m |= ios::trunc;

	primary = openOutputFile(file, m, bufS, doMC);
	asyncWriterEnabled = shouldUseAsyncWriter();
}

void ofbufstream::deactivate() {
	if (coutW) return;
	finishWrites();
	if (retained_) return;
	if (!primary) return;
	try {
		finalizeOutputFile(primary, file);
		primary.reset();
	}
	catch (...) {
		primary.reset();
		throw;
	}
}

bool ofbufstream::shouldUseAsyncWriter() const {
	// -threadIO 1 with more than one worker: overlap file I/O and compression
	// with processing. Single-threaded runs keep writing inline.
	return doMC && !coutW && bufS > 0 && ThreadPool::configured_thread_count() > 1 &&
		AsyncOutputWriter::instance() != nullptr;
}

//**************************** retained output streams ****************************
namespace {
struct RetentionState {
	std::mutex mutex;
	bool enabled = false;
	std::map<std::string, std::shared_ptr<ofbufstream>> streams;
	~RetentionState() {
		for (auto& entry : streams) {
			try {
				entry.second->setRetained(false);
				entry.second->deactivate();
			}
			catch (const std::exception& error) {
				std::cerr << "Failed to finalize " << entry.first << ": " << error.what() << "\n";
			}
		}
	}
};
RetentionState& retentionState() {
	static RetentionState state;
	return state;
}
}

namespace OutputStreamRetention {
void enable(bool on) {
	auto& state = retentionState();
	std::lock_guard<std::mutex> lock(state.mutex);
	state.enabled = on;
}

bool enabled() {
	auto& state = retentionState();
	std::lock_guard<std::mutex> lock(state.mutex);
	return state.enabled;
}

std::shared_ptr<ofbufstream> open(const string& path, int mode, bool isMC, size_t bufferS) {
	auto& state = retentionState();
	std::unique_lock<std::mutex> lock(state.mutex);
	if (!state.enabled || path.empty() || path == "-") {
		lock.unlock();
		return std::make_shared<ofbufstream>(path, mode, isMC, bufferS);
	}
	auto found = state.streams.find(path);
	if (found != state.streams.end()) {
		// Continue the open file; this replaces finalize + append-mode reopen.
		if (mode == ios::app) { return found->second; }
		// A truncating open restarts the file, exactly as without retention.
		std::shared_ptr<ofbufstream> previous = std::move(found->second);
		state.streams.erase(found);
		previous->setRetained(false);
		previous->deactivate();
	}
	auto stream = std::make_shared<ofbufstream>(path, mode, isMC, bufferS);
	stream->setRetained(true);
	state.streams.emplace(path, stream);
	return stream;
}

void closeAll() {
	auto& state = retentionState();
	std::map<std::string, std::shared_ptr<ofbufstream>> streams;
	{
		std::lock_guard<std::mutex> lock(state.mutex);
		streams.swap(state.streams);
	}
	std::exception_ptr failure;
	for (auto& entry : streams) {
		try {
			entry.second->setRetained(false);
			entry.second->deactivate();
		}
		catch (...) {
			if (!failure) { failure = std::current_exception(); }
		}
	}
	streams.clear();
	if (failure) { std::rethrow_exception(failure); }
}
}

dualOfBufStream::dualOfBufStream(void)
   : buf1S(5*1024 * 1024), buf2S(5*1024 * 1024), bufs(2, ""), FileNames(2, ""), dualOutStr(2), opened(2, false), active(false) {
	bufs[0].reserve(buf1S);
	bufs[1].reserve(buf2S);
}

dualOfBufStream::~dualOfBufStream(void) {
	// Destructors are noexcept: a finalize error (e.g. disk full) must not terminate.
	try {
		emptyStreams(true);
		deactivate();
	}
	catch (const std::exception& error) {
		std::cerr << "Failed to finalize " << FileNames[0] << "," << FileNames[1] << ": " << error.what() << "\n";
	}
}

void dualOfBufStream::write(const string& in, int stream) {
   if (stream < 0 || (size_t)stream >= bufs.size()) {
		return;
	}
    {
		Instr::TimedLockGuard lg(dualMtx);
		bufs[stream] += in;
	}
	emptyStreams(false);
}

void dualOfBufStream::write2(const string& in, const string& in2) {
 Instr::TimedLockGuard lg(dualMtx);
	if (dualOutStr.size() > 1 && dualOutStr[1]) {
		(*dualOutStr[1]) << in2;
	}
	if (dualOutStr.size() > 0 && dualOutStr[0]) {
		(*dualOutStr[0]) << in;
	}
}

bool dualOfBufStream::open(const string IF, int mif, int pair, bool isMC, size_t bufferS) {
  if (pair < 0 || (size_t)pair >= dualOutStr.size() || (size_t)pair >= opened.size() || (size_t)pair >= FileNames.size()) {
		return false;
	}
	if (opened[pair]) {
		return false;
	}
	dualOutStr[pair] = OutputStreamRetention::open(IF, mif, isMC, bufferS);
	opened[pair] = true;
	FileNames[pair] = IF;
	if (!dualOutStr[pair]) { return false; }
	return true;
}

bool dualOfBufStream::activate() {
	if (active) { return true; }
	for (size_t i = 0; i < dualOutStr.size(); i++) { if (dualOutStr[i] != nullptr) { dualOutStr[i]->activate(); } }
	active = true;
	cerr << "Activating dual ostreams: " << FileNames[0] << "," << FileNames[1] << endl;
	return true;
}

bool dualOfBufStream::deactivate() {
	Instr::TimedLockGuard lg(dualMtx);
	// Finalize every mate even if one fails, then report the first failure.
	std::exception_ptr failure;
	for (size_t i = 0; i < dualOutStr.size(); i++) {
		if (dualOutStr[i] != nullptr) {
			try {
				dualOutStr[i]->finishWrites();
				dualOutStr[i]->deactivate();
			}
			catch (...) {
				if (!failure) { failure = std::current_exception(); }
			}
		}
	}
	active = false;
	if (failure) { std::rethrow_exception(failure); }
	return true;
}

void dualOfBufStream::emptyStreams(bool force) {
 Instr::TimedLockGuard lg(dualMtx);
	if ((force || bufs[1].size() >= buf2S) && dualOutStr.size() > 1 && dualOutStr[1]) {
		(*dualOutStr[1]) << bufs[1];
		bufs[1].clear();
	}
	if ((force || bufs[0].size() >= buf1S) && dualOutStr.size() > 0 && dualOutStr[0]) {
		(*dualOutStr[0]) << bufs[0];
		bufs[0].clear();
	}
}

//**************************** per-worker output batches ****************************
// Records one worker formats for one input block. Publishing the whole block
// at once keeps R1/R2 (and FASTA/QUAL) files in step with a single lock per
// block instead of several locks per read.
struct OutputStreamer::WriteBatch {
	static constexpr size_t kFlushBytes = size_t(4) << 20;
	struct Dual {
		dualOfBufStream* stream;
		std::string first, second;
	};
	struct Single {
		ostr* stream;
		std::string data;
	};
	struct Demulti {
		ostr* first = nullptr;
		ostr* second = nullptr;
		ostr* merged = nullptr;
		std::string firstData, secondData, mergedData;
		bool touched = false;
	};
	std::vector<Dual> duals;
	std::vector<Single> singles;
	std::vector<Demulti> demulti;      // indexed by sample
	std::vector<size_t> demultiTouched;
	uint64_t records = 0;
	size_t bytes = 0;

	void stageDual(dualOfBufStream* stream, const std::string& first, const std::string& second) {
		Dual* entry = nullptr;
		for (auto& candidate : duals) { if (candidate.stream == stream) { entry = &candidate; break; } }
		if (!entry) { duals.push_back(Dual{stream, {}, {}}); entry = &duals.back(); }
		entry->first += first;
		entry->second += second;
		bytes += first.size() + second.size();
	}
	void stage(ostr* stream, const std::string& data) {
		Single* entry = nullptr;
		for (auto& candidate : singles) { if (candidate.stream == stream) { entry = &candidate; break; } }
		if (!entry) { singles.push_back(Single{stream, {}}); entry = &singles.back(); }
		entry->data += data;
		bytes += data.size();
	}
	void stageDemulti(size_t sample, bool merged, ostr* firstStream, const std::string& first,
		ostr* secondStream, const std::string& second) {
		if (sample >= demulti.size()) { demulti.resize(sample + 1); }
		Demulti& entry = demulti[sample];
		if (!entry.touched) { entry.touched = true; demultiTouched.push_back(sample); }
		if (merged) {
			entry.merged = firstStream;
			entry.mergedData += first;
		} else {
			if (firstStream) { entry.first = firstStream; entry.firstData += first; }
			if (secondStream) { entry.second = secondStream; entry.secondData += second; }
		}
		bytes += first.size() + second.size();
	}
	bool empty() const { return bytes == 0 && records == 0 && demultiTouched.empty(); }
	void clearMain() {
		for (auto& entry : duals) { entry.first.clear(); entry.second.clear(); }
		for (auto& entry : singles) { entry.data.clear(); }
		records = 0;
	}
	void clearDemulti() {
		for (size_t sample : demultiTouched) {
			Demulti& entry = demulti[sample];
			entry.firstData.clear(); entry.secondData.clear(); entry.mergedData.clear();
			entry.touched = false;
		}
		demultiTouched.clear();
	}
};

OutputStreamer::WriteBatch* OutputStreamer::writeBatch(int curThread) {
	if (curThread < 0 || !batchWrites_.load(std::memory_order_acquire) ||
		static_cast<size_t>(curThread) >= writeBatches_.size()) {
		return nullptr;
	}
	return writeBatches_[curThread].get();
}

void OutputStreamer::beginWriteBatching() {
	// Output rotation and read limits count records one at a time; they keep
	// the per-record path so file boundaries and limits stay exact.
	if (maxReadsPerOFile > 0 || maxRdsOut > 0) { return; }
	std::lock_guard<std::mutex> guard(sqfqostrMTX);
	writeBatches_.clear();
	const size_t workers = static_cast<size_t>(std::max(Nthrds, 1));
	for (size_t i = 0; i < workers; ++i) { writeBatches_.push_back(std::make_unique<WriteBatch>()); }
	batchWrites_.store(true, std::memory_order_release);
}

void OutputStreamer::publishBatchLocked(WriteBatch& batch) {
	for (auto& entry : batch.duals) {
		if (!entry.first.empty() || !entry.second.empty()) { entry.stream->write2(entry.first, entry.second); }
	}
	for (auto& entry : batch.singles) {
		if (!entry.data.empty()) { *entry.stream << entry.data; }
	}
	ReadsWritten += batch.records;
	TotalReadsWritten += batch.records;
	batch.clearMain();
}

void OutputStreamer::publishDemultiBatch(WriteBatch& batch) {
	if (batch.demultiTouched.empty()) { return; }
	// R1/R2 of a sample must stay synchronized; one lock per block and worker.
	std::lock_guard<std::mutex> guard(dmltMTX);
	for (size_t sample : batch.demultiTouched) {
		const WriteBatch::Demulti& entry = batch.demulti[sample];
		if (entry.first && !entry.firstData.empty()) { *entry.first << entry.firstData; }
		if (entry.second && !entry.secondData.empty()) { *entry.second << entry.secondData; }
		if (entry.merged && !entry.mergedData.empty()) { *entry.merged << entry.mergedData; }
	}
	batch.clearDemulti();
}

void OutputStreamer::flushWriteBatch(int curThread) {
	WriteBatch* batch = writeBatch(curThread);
	if (!batch || batch->empty()) { return; }
	{
		std::lock_guard<std::mutex> guard(sqfqostrMTX);
		publishBatchLocked(*batch);
	}
	publishDemultiBatch(*batch);
	batch->bytes = 0;
}

void OutputStreamer::flushAllWriteBatchesLocked() {
	for (auto& batch : writeBatches_) {
		if (!batch || batch->empty()) { continue; }
		publishBatchLocked(*batch);
		publishDemultiBatch(*batch);
		batch->bytes = 0;
	}
}

void OutputStreamer::endWriteBatching() {
	if (!batchWrites_.load(std::memory_order_acquire)) { return; }
	// Workers have finished; publish what each staged, then write per record again.
	std::lock_guard<std::mutex> guard(sqfqostrMTX);
	flushAllWriteBatchesLocked();
	batchWrites_.store(false, std::memory_order_release);
	writeBatches_.clear();
}

OutputStreamer::OutputStreamer(Filters* fil, OptContainer* cmdArgs,
	std::ios_base::openmode writeStatus,
	shared_ptr<ReadSubset> RDSset,
	int numThreads, string fileExt, int forceFmt) :
	MFil(nullptr), subFilter(0),
	//locBufGreen(2,""),
	suppressOutWrite(0), write2File(true),
	fastQver(33),
	fastQoutVer(33), BWriteQual(false),
	BWriteFastQ(false), b_multiOutStream(false), pairedSeq(-1),
	b_changeFQheadVer(false),
	b_checkedHeaderChange(false),
	b_oneLinerFasta(false), b_doDereplicate(false), b_writeGreenQual(true), b_writeYellowQual(true),
	maxReadsPerOFile(0),
	demultiBPperSR(0),
	BPwrittenInSR(0), BPwrittenInSRmerg(0),
	ReadsWritten(0), TotalReadsWritten(0),
	maxRdsOut(0), stopAll(false),
	leadingOutf(""), locCmdArgs(cmdArgs), dereplicator(nullptr), cntDerep(0), wrMode(ios::out),
	sFile(0), qFile(0), fqFile(0),
	sPairFile(0), qPairFile(0), fqPairFile(0),
	of_merged_fq(0),
	//sFileStr(0), qFileStr(0), fqFileStr(0), fqNoBCFile(0), 
	totalFileStrms(0), doTIO(true),
	bDoDemultiplexIntoFiles(false), demultiSinglFiles(0), //demultiSinglFilesF(0),
	demultiMergeFiles(0),
	onlyCompletePairsDemulti(false),
	b_merge_pairs_demulti_(false), b_merge_pairs_filter_(false),
	b_merge_pairs_(false), mergers(0),
	Nthrds(numThreads), _benchmark(nullptr)
{
	MFil = fil;

	fastQver = fil->getuserReqFastqVer();
	fastQoutVer = fil->getuserReqFastqOutVer();

	maxReadsPerOFile = fil->maxReadsOutput();
	demultiBPperSR = fil->getDemultiBPperSR();
	ReadsWritten = fil->writtenReads();

	pairedSeq = MFil->isPaired();
	if (cmdArgs->find("-suppressOutput") != cmdArgs->end()) {
		suppressOutWrite = atoi((*cmdArgs)["-suppressOutput"].c_str());
	}
	if (cmdArgs->find("-pairedDemulti") != cmdArgs->end() && (*cmdArgs)["-pairedDemulti"] == "1") { //only demultiplex proper pairs
		onlyCompletePairsDemulti = true;
	}
	if (cmdArgs->find("-oneLineFastaFormat") != cmdArgs->end() && (*cmdArgs)["-oneLineFastaFormat"] == "1") {
		b_oneLinerFasta = true;
	}


	if (suppressOutWrite == 3 || suppressOutWrite == 1) {
		b_writeGreenQual = false;
	}
	if (suppressOutWrite == 3 || suppressOutWrite == 2) {
		b_writeYellowQual = false;
	}
	if (fil->doSubselReads() && RDSset != nullptr) {
		this->openSeveralOutstreams(locCmdArgs, RDSset, writeStatus);
	}
	else {//standard one file output stream
		this->openOutStreams(locCmdArgs, MFil->getFileIncrementor(),
			writeStatus, fileExt, forceFmt);
	}
	if ((*cmdArgs)["-o_fastq_noBC"] != "") {
		openNoBCoutstrean((*cmdArgs)["-o_fastq_noBC"]);
	}
	if (cmdArgs->find("-o_demultiplex") != cmdArgs->end()) {//demulti: always write mode out
		bool gzDemulti = false;
		if (cmdArgs->find("-o_demultiplex_gz") != cmdArgs->end()) {
			if ((*cmdArgs)["-o_demultiplex_gz"] == "1") {
				gzDemulti = true;
				//cerr << "Gzip demulti output\n";
			}
			else {
				//cerr << "\n\nNo Gzip demulti output\n";
			}
		}
		generateDemultiOutFiles((*cmdArgs)["-o_demultiplex"], fil, ios::out, gzDemulti);
	}
	if (cmdArgs->find("-merge_pairs_filter") != cmdArgs->end()
		&& (*cmdArgs)["-merge_pairs_filter"] == "1") {
		b_merge_pairs_filter_ = true;
	}
	if (cmdArgs->find("-merge_pairs_demulti") != cmdArgs->end()
		&& (*cmdArgs)["-merge_pairs_demulti"] == "1") {
		b_merge_pairs_demulti_ = true;
	}
	if (b_merge_pairs_demulti_ || b_merge_pairs_filter_) {
		b_merge_pairs_ = true;
	}
	setSubfilters(Nthrds);

	if (this->pairedSeq != 2) {
		onlyCompletePairsDemulti = false;
	}
	if ((*cmdArgs)["-threadIO"] == "1") { doTIO = true; }
	if ((*cmdArgs)["-threadIO"] == "0") { doTIO = false; }


	//threads = futures(num_threads);
	//if (pairedSeq){cerr<<"paired MFil\n";}
}


OutputStreamer::~OutputStreamer() {
	//delete MFil;
	cdbg("Destr OutputStreamer");
	//delAllDNAvectors();
	//cdbg(".. done\n" );
	try { endWriteBatching(); }
	catch (const std::exception& error) { cerr << "Failed to write buffered output: " << error.what() << "\n"; }
	if (dereplicator) {
		// Queued entries may still publish withheld sets into these streams.
		try { dereplicator->flushPending(); }
		catch (const std::exception& error) { cerr << "Failed to apply queued dereplication entries: " << error.what() << "\n"; }
		dereplicator->setUncountedPublisher(nullptr);
	}
    // Acquire stream-creation/usage mutex first, then demulti mutex to avoid lock-order inversion
    {
		Instr::TimedLockGuard g1(sqfqostrMTX);
		Instr::TimedLockGuard g2(dmltMTX);
		// Call the internal locked variants to avoid double-locking sqfqostrMTX.
		// Callers close explicitly (and get the exception); this is the fallback.
		try { this->closeOutStreams_locked(true); }
		catch (const std::exception& error) { cerr << "Failed to finalize output: " << error.what() << "\n"; }
		try { this->closeOutFilesDemulti_locked(); }
		catch (const std::exception& error) { cerr << "Failed to finalize demultiplexed output: " << error.what() << "\n"; }
		for (uint i = 0; i < subFilter.size(); i++) {
			subFilter[i].reset();
		}
		for (size_t x = 0; x < mergers.size(); x++) {
			mergers[x].reset();
		}
	}
	if (paired_sync_skipped_counter_ > 0) {
		cerr << "Warning: skipped " << paired_sync_skipped_counter_.load()
			<< " paired write attempts to keep R1/R2 synchronized." << endl;
	}
	cdbg("Subfilters deleted, streams closed\n");
	//for (size_t x = 0; x < merger.size(); x++) { delete merger[x]; } merger.clear(); 
	//delete optim;
}
void OutputStreamer::delAllDNAvectors() {
#ifdef DEBUG
	cerr << "cleaning MD..";
#endif
	//locBufGreen[0].clear(); locBufGreen[1].clear();
	//hard delete all saved data, no write to HDD
	for (size_t i = 0; i < sFile.size(); i++) {
		for (size_t j = 0; j < sFile[i].size(); j++) {
			if (sFile[i][j] != nullptr) {
				sFile[i][j]->emptyStream();
			}
		}
	}
	for (size_t i = 0; i < qFile.size(); i++) {
		for (size_t j = 0; j < qFile[i].size(); j++) {
			if (qFile[i][j] != nullptr) {
				qFile[i][j]->emptyStream();
			}
		}
	}
	for (size_t i = 0; i < fqFile.size(); i++) {
		for (size_t j = 0; j < fqFile[i].size(); j++) {
			if (fqFile[i][j] != nullptr) {
				fqFile[i][j]->emptyStream();
			}
		}
	}
	// s/q/fqPairFile destructor: release unique_ptrs
	for (size_t i = 0; i < fqPairFile.size(); i++) {
		if (fqPairFile[i]) {
			fqPairFile[i].reset();
		}
	}
	for (size_t i = 0; i < qPairFile.size(); i++) {
		if (qPairFile[i]) {
			qPairFile[i].reset();
		}
	}
	for (size_t i = 0; i < sPairFile.size(); i++) {
		if (sPairFile[i]) {
			sPairFile[i].reset();
		}
	}

}

void OutputStreamer::analyzeDNA(const shared_ptr<DNA>& d, int FilterUse, int pair, int& idx, int thr) {
	if (!d) {
		return;
	}
	if (d->isFailed()) { return; }
	(void)FilterUse;
	Filters* curFil = this->getFilters(thr);
	const bool isP1 = (pair <= 0);
	const int cut5p = isP1 ? curFil->getcut5PR1() : curFil->getcut5PR2();
	const bool doFilter = curFil->doFilterAtAll();
	const bool dualBC = curFil->doubleBarcodes();

	if (cut5p > 0) {
		d->cutSeq(0, cut5p);
	}

	if (!doFilter) {
		if (idx < 0 && (isP1 || dualBC)) {
			idx = curFil->detectCutBC(d, isP1); //still need to check for BC
		}
		if (idx >= 0) { //prevent second read pair_ from being flagged as true
			d->setPassed(true);
		}
		return;
	}

	//pricipally safe to call from thread
	curFil->checkYellowAndGreen(d, pair, idx, false, Filters::FilteringContext::FinalFiltering);

	//bit redundant.. but better safe than sorry
	if (!curFil->secondaryOutput() && d->isYellowQual() && !d->isGreenQual()) {
		d->failed();
	}
	if (idx == -1) {
		d->setBarcodeDetected(false);
	}
	return;

	/*
	if (curFil->secondaryOutput() ){
		//make this primary routine..

	} else if (FilterUse==-1){
		curFil->check(d, false, pair, idx);
	} else {
		cerr << "Invalid control path in analyzeDNA\n"; exit(55);
		//Filters* f = subFilter[FilterUse];
		//subFilter[FilterUse].check(d, false, pair_, idx);
		//f->check(d, false, pair, idx);
	}
	*/

	//count this as failure if BC was present
	//d->prepareWrite(fastQoutVer);
}

vector<bool> OutputStreamer::analyzeDNA(shared_ptr<DNA> p1, shared_ptr<DNA> p2, shared_ptr<DNA> mid,
	bool changePHead, int FilterUse) {
	cerr << "deprecated analuze DNA"; exit(2323);
	vector<bool> ret(2, true);
	//1st: check if DNA pointer valid
	if (p1 == NULL) {
		ret[0] = false;
	}
	if (p2 == NULL) {
		ret[1] = false;
	}
	if (mid == NULL) {//no MID? kill
		ret[1] = false; ret[0] = false; return ret;
	}
	else {
		mid->setMIDseq(true);
	}

	//2nd: check if basic DNA signatures are valid

	//3rd: check if all quality criteria are met, BC is true etc.
	//if (FilterUse == -1){
	//ret = curFil->check_pairs(p1, p2, mid, ret, changePHead);
	//}	else { //mutithread, to avoid race conditions etc use separate filter
	//	ret = subFilter[FilterUse]->check_pairs(p1, p2, mid, ret, changePHead);
	//}
	//if (ret[0]){ p1->prepareWrite(fastQoutVer); }
	//if (ret[1]){	p2->prepareWrite(fastQoutVer);}

	return ret;
}

bool OutputStreamer::checkFastqHeadVersion(shared_ptr<DNA> d, bool disable) {
    // If we've already checked, return cached result
	if (b_checkedHeaderChange) {
		return b_changeFQheadVer;
	}
	b_checkedHeaderChange = true;
	if (pairedSeq == 1) { b_changeFQheadVer = false; return false; }
    bool ret = false;
	{
		Instr::TimedLockGuard g(drpMTX);
		b_changeFQheadVer = false;
		int fastQheadVer = 0;
		string head = d->getId();
		int shouldVer = MFil->FQheadV();
	if (head.find("/1") != string::npos || head.find("/2") != string::npos) {
		fastQheadVer = 1;
	}
	if (head.find(" 1:") != string::npos || head.find(" 2:") != string::npos) {
		fastQheadVer = 2;
	}
	if (shouldVer != 0 && shouldVer != fastQheadVer) {
		b_changeFQheadVer = true;
	}
        ret = b_changeFQheadVer;
	if (disable) {
		b_changeFQheadVer = false;
	}
    }
	return ret;
}

void OutputStreamer::closeOutFilesDemulti() {
	cdbg("OutputStreamer::closeOutFilesDemulti");
	if (!bDoDemultiplexIntoFiles) { return; }
	endWriteBatching();
	// Acquire locks in the global order: sqfqostrMTX then dmltMTX
    Instr::TimedLockGuard guard1(sqfqostrMTX);
	Instr::TimedLockGuard guard2(dmltMTX);
	closeOutFilesDemulti_locked();
	cdbg(" done closeOutFilesDemulti\n");
}

// Internal variant that assumes caller already holds required locks.
// Finalizes every per-sample file explicitly so that a write failure is
// reported (the destructor only prints it); the first failure is rethrown after
// all files were closed.
void OutputStreamer::closeOutFilesDemulti_locked() {
	if (!bDoDemultiplexIntoFiles) { return; }
	std::exception_ptr failure;
	auto closeStream = [&](std::shared_ptr<ofbufstream>& stream) {
		if (!stream) { return; }
		try { stream->deactivate(); }
		catch (...) { if (!failure) { failure = std::current_exception(); } }
		stream.reset();
	};
	for (auto& sample : demultiSinglFiles) { for (auto& stream : sample) { closeStream(stream); } }
	for (auto& stream : demultiMergeFiles) { closeStream(stream); }
	cdbg(" done closeOutFilesDemulti_locked\n");
	if (failure) { std::rethrow_exception(failure); }
}

bool OutputStreamer::reserveDemultiBases(bool merged, uint64_t bases) {
	auto& counter = merged ? BPwrittenInSRmerg : BPwrittenInSR;
	if (demultiBPperSR == 0) {
		counter += bases;
		return true;
	}
	// Same admission rule as the locked path, as one atomic reservation.
	uint64_t written = counter.load();
	do {
		if (written >= demultiBPperSR || bases > demultiBPperSR - written) { return false; }
	} while (!counter.compare_exchange_weak(written, written + bases));
	return true;
}

bool OutputStreamer::writeDemultiRecords(int idx, const string& first, const string& second,
	uint64_t bases, bool merged, int curThread) {
	if (idx < 0 || (first.empty() && second.empty())) { return false; }
	if (WriteBatch* batch = writeBatch(curThread)) {
		// Demultiplexing streams are opened before processing and stay fixed
		// while batching, so they can be resolved without dmltMTX.
		ostr* firstStream = nullptr;
		ostr* secondStream = nullptr;
		if (merged) {
			if (static_cast<size_t>(idx) >= demultiMergeFiles.size() || !demultiMergeFiles[idx]) { return false; }
			firstStream = demultiMergeFiles[idx].get();
		} else {
			if (static_cast<size_t>(idx) >= demultiSinglFiles.size()) { return false; }
			const auto& streams = demultiSinglFiles[idx];
			if ((!first.empty() && (streams.empty() || !streams[0])) ||
				(!second.empty() && (streams.size() < 2 || !streams[1]))) { return false; }
			if (!first.empty()) { firstStream = streams[0].get(); }
			if (!second.empty()) { secondStream = streams[1].get(); }
		}
		if (!reserveDemultiBases(merged, bases)) { return false; }
		batch->stageDemulti(static_cast<size_t>(idx), merged, firstStream, first, secondStream, second);
		if (batch->bytes >= WriteBatch::kFlushBytes) { flushWriteBatch(curThread); }
		return true;
	}
	std::lock_guard<std::mutex> guard(dmltMTX);
	auto& counter = merged ? BPwrittenInSRmerg : BPwrittenInSR;
	const uint64_t written = counter.load();
	// Admit whole records/pairs and count both mates. Each destination has its
	// own per-run budget; an exhausted merged budget cannot fall back to R1/R2.
	if (demultiBPperSR > 0 && (written >= demultiBPperSR || bases > demultiBPperSR - written)) {
		return false;
	}
	if (merged) {
		if (static_cast<size_t>(idx) >= demultiMergeFiles.size() || !demultiMergeFiles[idx]) { return false; }
		*demultiMergeFiles[idx] << first;
	} else {
		if (static_cast<size_t>(idx) >= demultiSinglFiles.size()) { return false; }
		const auto& streams = demultiSinglFiles[idx];
		if ((!first.empty() && (streams.empty() || !streams[0])) ||
			(!second.empty() && (streams.size() < 2 || !streams[1]))) { return false; }
		if (!first.empty()) { *streams[0] << first; }
		if (!second.empty()) { *streams[1] << second; }
	}
	counter += bases;
	return true;
}

void OutputStreamer::write2Demulti(const shared_ptr<DNA>& d1, const shared_ptr<DNA>& d2, int BCoffset,
	int curThread) {
	if (!Demulti2Fls() || !d1) { return; }
	if (isPEseq() < 2) {
		write2Demulti(d1, 0, BCoffset, curThread);
		return;
	}
	const int idx = d1->getBCnumber() - BCoffset;
	if (idx < 0) { return; }
	const bool green1 = d1->isGreenQual();
	const bool green2 = d2 && d2->isGreenQual();
	if (d2 && (green1 || green2) && b_merge_pairs_demulti_ && d1->merge_seed_pos_ >= 0 &&
		curThread >= 0 && static_cast<size_t>(curThread) < mergers.size() && mergers[curThread]) {
		auto merged = mergers[curThread]->merge(d1, d2);
		if (merged && getFilters(curThread)->checkMergedDNA(merged)) {
			d1->adoptMergeStats(*merged);
			d2->adoptMergeStats(*merged);
			merged->prepareWrite(fastQoutVer);
			writeDemultiRecords(idx, merged->writeFastQ(), "", merged->length(), true, curThread);
			return;
		}
	}
	if (onlyCompletePairsDemulti && (!green1 || !green2)) { return; }
	string first, second;
	uint64_t bases = 0;
	if (green1) {
		d1->prepareWrite(fastQoutVer);
		d1->writeFastQ(first, false);
		if (!first.empty()) { bases += d1->length(); }
	}
	if (green2) {
		d2->prepareWrite(fastQoutVer);
		d2->writeFastQ(second, false);
		if (!second.empty()) { bases += d2->length(); }
	}
	writeDemultiRecords(idx, first, second, bases, false, curThread);
}

void OutputStreamer::write2Demulti(const shared_ptr<DNA>& d, int p, int BCoffset, int curThread) {
	if (!Demulti2Fls() || !d || !d->isGreenQual() || p < 0 || p > 1) { return; }
	d->prepareWrite(fastQoutVer);
	const string record = d->writeFastQ(false);
	writeDemultiRecords(d->getBCnumber() - BCoffset, p == 0 ? record : "",
		p == 1 ? record : "", d->length(), false, curThread);
}

void OutputStreamer::generateDemultiOutFiles(string path, Filters* fil,
	std::ios_base::openmode writeStatus, bool demulti2gz) {
	//fill in demultiSinglFiles vector
	vector<ofbufstream*> empVec(2, NULL);
	//vector<string> empVec2(2, "");

	bool doMC = Nthrds > 1;
	string gzSuffiz = "";
	string gzReport = "";
	if (demulti2gz) {
		gzSuffiz = ".gz";
		gzReport = "gzipped ";
	}

	bool displayInfo(false);
	if (fil->getReadsWritten() == 0) {
		displayInfo = true;
	}

	struct stat info;
	if (stat(path.c_str(), &info) != 0) {
		cerr << "Output path for demultiplexed files does not exist, please create this directory:\n" << path << endl;
		exit(833);
	}
	else if (info.st_mode & S_IFDIR) { // S_ISDIR() doesn't exist on my windows 
		if (displayInfo) {
			cerr << "\nWriting " << gzReport << "demultiplexed files to : " << path << endl;// printf(" % s is a directory\n", path);
		}
	}
	else {
		cerr << path << " is no directory\n"; exit(834);
	}
	path += "/";

	bDoDemultiplexIntoFiles = true;
	bool openOstreams = true; uint ostrCnt(0);
	size_t bufS = OUTPUT_BUFFER_SIZE;
	//save as .gz? -> no for dada2..

	demultiSinglFiles.resize(fil->SampleID.size());
	demultiMergeFiles.resize(fil->SampleID.size());
	//demultiSinglFilesF.resize(fil->SampleID.size(), empVec2);
	for (size_t i = 0; i < fil->SampleID.size(); i++) {
		//actually needs to know if paired files..
		//if (ostrCnt > maxFileStreams) {openOstreams = false;}
		if (pairedSeq == 1 || pairedSeq == -1) {
			string nfile = path + fil->SampleID[i] + ".fq" + gzSuffiz;
			if (openOstreams) {
                demultiSinglFiles[i].resize(1);
				demultiSinglFiles[i][0] = std::make_shared<ofbufstream>(nfile.c_str(), writeStatus, doMC, (size_t)bufS);
			}
			//demultiSinglFilesF[i][0] = nfile;
			ostrCnt++;
		}
		else {
			string nfile = path + fil->SampleID[i] + ".1.fq" + gzSuffiz;
			if (openOstreams) {
                demultiSinglFiles[i].resize(2);
				demultiSinglFiles[i][0] = std::make_shared<ofbufstream>(nfile.c_str(), writeStatus, doMC, size_t(bufS * 0.8));
			}
			//demultiSinglFilesF[i][0] = nfile;
			nfile = path + fil->SampleID[i] + ".2.fq" + gzSuffiz;
			if (openOstreams) {
                demultiSinglFiles[i][1] = std::make_shared<ofbufstream>(nfile.c_str(), writeStatus, doMC, size_t(bufS * 1.2));
			}
			//demultiSinglFilesF[i][1] = nfile;
			nfile = path + fil->SampleID[i] + ".merg.fq" + gzSuffiz;
			if (openOstreams) {
                demultiMergeFiles[i] = std::make_shared<ofbufstream>(nfile.c_str(), writeStatus, doMC, (size_t)bufS);
			}

			ostrCnt += 2;
		}
	}

}
void OutputStreamer::incrementOutputFile() {
	this->closeOutStreams(true);
	MFil->incrementFileIncrementor();
	this->openOutStreams(locCmdArgs, MFil->getFileIncrementor(), ios_base::out);
	ReadsWritten = 0;
}




bool  OutputStreamer::saveForWrite(const shared_ptr<DNA>& d, int Pair, int thr, int& Cstream, bool write) {
	//second most important part: collect stats on DNA passing through here (should be all read)
	//most important part: save DNA to be written later (or discard)
	if (d == NULL || stopAll) {
		Cstream = 100;
		return !stopAll;
	}
	//int Cstream = 0;
	//bool writen(false);

	//threadsafe
	Filters* curFil = this->getFilters(thr);

	curFil->addDNAtoCStats(d, Pair);

	if (suppressOutWrite == 3) {
		Cstream = 100; //set to impossibly high value
		return !stopAll;
	}

	if (curFil->doFilterAtAll()) {
		if ((b_writeGreenQual && d->isGreenQual()) ||
			(b_writeYellowQual && d->isYellowQual())) {
			if (BWriteFastQ && b_changeFQheadVer) {//check if header PE naming needs to be changed
				d->changeHeadPver(curFil->FQheadV());
			}
			if (d->isYellowQual()) {
				Cstream = 1;
			}
			else { Cstream = 0; }
		}
		else {
			Cstream = 100;
			return !stopAll;
		}
	}
	else {
		Cstream = 0;
	}
	d->prepareWrite(fastQoutVer);

	/*if (curFil->passedReads(25000)) {
		if (pairedSeq > 1) {
			cerr << "25k read pairs/" << _benchmark->now_interval_secs(2) << "s\r";
		} else {
			cerr << "25k reads/" << _benchmark->now_interval_secs(2) << "s\r";
		}
		//cerr.flush();
	}*/

	//sqfqostrMTX.lock();
	if (write) {
		writeForWrite(d, Pair, Cstream, nullptr, -1, 100, thr);
	}
	return !stopAll;


}

void OutputStreamer::writeForWrite(const shared_ptr<DNA>& d1, int Pair1, int Cstream1,
	const shared_ptr<DNA>& d2, int Pair2, int Cstream2, int curThread) {
	writeForWriteImpl(d1, Pair1, Cstream1, d2, Pair2, Cstream2, nullptr, curThread);
}

OutputStreamer::StreamSlots OutputStreamer::selectStreams(bool writeFastq, bool isStrictPairedWrite,
	bool mergedOutput, int Pair1, int Cstream1, int Pair2, int Cstream2) const {
	StreamSlots slots;
	const bool p1Valid = Pair1 > 0;
	const bool p2Valid = Pair2 > 0;
	if (writeFastq && isStrictPairedWrite && Cstream1 >= 0 &&
		static_cast<size_t>(Cstream1) < fqPairFile.size()) {
		slots.fqPair = &fqPairFile[Cstream1];
	}
	if (mergedOutput) {
		if (!of_merged_fq.empty()) { slots.fq1 = &of_merged_fq[0]; }
	}
	else if (writeFastq) {
		if (Cstream1 >= 0 && static_cast<size_t>(Cstream1) < fqFile.size() && p1Valid &&
			static_cast<size_t>(Pair1 - 1) < fqFile[Cstream1].size()) {
			slots.fq1 = &fqFile[Cstream1][Pair1 - 1];
		}
		if (Cstream2 >= 0 && static_cast<size_t>(Cstream2) < fqFile.size() && p2Valid &&
			static_cast<size_t>(Pair2 - 1) < fqFile[Cstream2].size()) {
			slots.fq2 = &fqFile[Cstream2][Pair2 - 1];
		}
	}
	else {
		if (Cstream1 >= 0 && static_cast<size_t>(Cstream1) < sFile.size() && p1Valid &&
			static_cast<size_t>(Pair1 - 1) < sFile[Cstream1].size()) {
			slots.seq1 = &sFile[Cstream1][Pair1 - 1];
			if (BWriteQual && static_cast<size_t>(Cstream1) < qFile.size() &&
				static_cast<size_t>(Pair1 - 1) < qFile[Cstream1].size()) {
				slots.qual1 = &qFile[Cstream1][Pair1 - 1];
			}
		}
		if (Cstream2 >= 0 && static_cast<size_t>(Cstream2) < sFile.size() && p2Valid &&
			static_cast<size_t>(Pair2 - 1) < sFile[Cstream2].size()) {
			slots.seq2 = &sFile[Cstream2][Pair2 - 1];
			if (BWriteQual && static_cast<size_t>(Cstream2) < qFile.size() &&
				static_cast<size_t>(Pair2 - 1) < qFile[Cstream2].size()) {
				slots.qual2 = &qFile[Cstream2][Pair2 - 1];
			}
		}
	}
	return slots;
}

bool OutputStreamer::writeForWriteImpl(const shared_ptr<DNA>& d1, int Pair1, int Cstream1,
	const shared_ptr<DNA>& d2, int Pair2, int Cstream2, const string* mergedFastq, int curThread) {
	const bool mergedOutput = mergedFastq != nullptr;
	if (mergedOutput && (!d1 || mergedFastq->empty())) { return false; }
	// Suppressed output is common in dereplication: avoid formatting and locks.
	if (Cstream1 >= 100 && Cstream2 >= 100) { return false; }
	const bool isStrictPairedWrite = (!mergedOutput && pairedSeq > 1 && Pair1 > 0 && Pair1 < 3 && Pair2 > 0 && Pair2 < 3);
	if (isStrictPairedWrite) {
		// Primary suppression must not discard a green mate when its partner
		// qualifies for secondary output: the whole pair also passes that tier.
		if (Cstream1 == 100 && Cstream2 == 1 && d1 && d1->isGreenQual()) { Cstream1 = 1; }
		if (Cstream2 == 100 && Cstream1 == 1 && d2 && d2->isGreenQual()) { Cstream2 = 1; }
		if (!d1 || !d2 || d1->length() == 0 || d2->length() == 0 ||
			Cstream1 < 0 || Cstream2 < 0 || Cstream1 >= 100 || Cstream2 >= 100) {
			paired_sync_skipped_counter_++;
			return false;
		}
		// A pair belongs to the less stringent of its two quality tiers.
		Cstream1 = Cstream2 = std::max(Cstream1, Cstream2);
	}
	WriteBatch* batch = writeBatch(curThread);
	// Batching excludes rotation, so the epoch lock is only needed per record.
	std::shared_lock<std::shared_mutex> epochGuard(outputEpochMTX, std::defer_lock);
	if (!batch) { epochGuard.lock(); }
	// Output rotation re-opens the streams (and sets the format) under the
	// exclusive epoch lock, so read the format only after taking it.
	const bool writeFastq = BWriteFastQ || mergedOutput;
	thread_local string fastqTmp1;
	thread_local string fastqTmp2;
	const bool isStrictPairedFastqWrite = writeFastq && isStrictPairedWrite;
	if (isStrictPairedFastqWrite) {
		fastqTmp1.clear(); fastqTmp2.clear();
		d1->writeFastQ(fastqTmp1);
		d2->writeFastQ(fastqTmp2);
		if (fastqTmp1.empty() || fastqTmp2.empty()) {
			paired_sync_skipped_counter_++;
			return false;
		}
	}

	thread_local string seqTmp1, seqTmp2, qualTmp1, qualTmp2;
	if (!writeFastq) {
		seqTmp1.clear(); seqTmp2.clear(); qualTmp1.clear(); qualTmp2.clear();
		if (d1) {
			d1->writeSeq(seqTmp1, b_oneLinerFasta);
			if (BWriteQual) { d1->writeQual(qualTmp1, b_oneLinerFasta); }
		}
		if (d2) {
			d2->writeSeq(seqTmp2, b_oneLinerFasta);
			if (BWriteQual) { d2->writeQual(qualTmp2, b_oneLinerFasta); }
		}
		if (isStrictPairedWrite && (seqTmp1.empty() || seqTmp2.empty())) {
			paired_sync_skipped_counter_++;
			return false;
		}
	}

	if (batch) {
		// Stream vectors do not change while batching (no rotation or read limit),
		// so records are staged without sqfqostrMTX or shared_ptr copies.
		const StreamSlots slots = selectStreams(writeFastq, isStrictPairedWrite, mergedOutput,
			Pair1, Cstream1, Pair2, Cstream2);
		auto raw = [](const auto* slot) { return slot ? slot->get() : nullptr; };
		dualOfBufStream* fqPair = raw(slots.fqPair);
		ostr* fq1 = raw(slots.fq1);
		ostr* fq2 = raw(slots.fq2);
		ostr* seq1 = raw(slots.seq1);
		ostr* seq2 = raw(slots.seq2);
		ostr* qual1 = raw(slots.qual1);
		ostr* qual2 = raw(slots.qual2);
		if (isStrictPairedWrite && ((writeFastq && !fqPair) ||
			(!writeFastq && (!seq1 || !seq2 ||
				((qual1 || qual2) && (!qual1 || !qual2 ||
					qualTmp1.empty() || qualTmp2.empty())))))) {
			paired_sync_skipped_counter_++;
			return false;
		}
		const bool willWrite = (fqPair && d1 && d2) || (fq1 && d1) ||
			(fq2 && d2) || (seq1 && d1) || (seq2 && d2);
		if (!willWrite) {
			if (isStrictPairedFastqWrite) { paired_sync_skipped_counter_++; }
			return false;
		}
		if (fqPair) {
			batch->stageDual(fqPair, fastqTmp1, fastqTmp2);
		}
		else if (mergedOutput) {
			batch->stage(fq1, *mergedFastq);
		}
		else {
			if (fq1 && d1) { fastqTmp1.clear(); d1->writeFastQ(fastqTmp1); batch->stage(fq1, fastqTmp1); }
			if (fq2 && d2) { fastqTmp2.clear(); d2->writeFastQ(fastqTmp2); batch->stage(fq2, fastqTmp2); }
			if (!writeFastq) {
				if (seq1 && d1) {
					batch->stage(seq1, seqTmp1);
					if (qual1) { batch->stage(qual1, qualTmp1); }
				}
				if (seq2 && d2) {
					batch->stage(seq2, seqTmp2);
					if (qual2) { batch->stage(qual2, qualTmp2); }
				}
			}
		}
		++batch->records;
		if (batch->bytes >= WriteBatch::kFlushBytes) { flushWriteBatch(curThread); }
		return true;
	}

	for (;;) {
		std::shared_ptr<dualOfBufStream> localFqPair;
		std::shared_ptr<ostr> localFq1, localFq2, localSeq1, localSeq2, localQual1, localQual2;
		bool rotate = false;
		{
			std::lock_guard<std::mutex> streamGuard(sqfqostrMTX);
			if (Cstream1 >= 100 && Cstream2 >= 100) { return false; }
			if (maxRdsOut > 0 && TotalReadsWritten >= maxRdsOut) { return false; }
			{
				const StreamSlots slots = selectStreams(writeFastq, isStrictPairedWrite, mergedOutput,
					Pair1, Cstream1, Pair2, Cstream2);
				if (slots.fqPair) { localFqPair = *slots.fqPair; }
				if (slots.fq1) { localFq1 = *slots.fq1; }
				if (slots.fq2) { localFq2 = *slots.fq2; }
				if (slots.seq1) { localSeq1 = *slots.seq1; }
				if (slots.seq2) { localSeq2 = *slots.seq2; }
				if (slots.qual1) { localQual1 = *slots.qual1; }
				if (slots.qual2) { localQual2 = *slots.qual2; }

				if (isStrictPairedWrite && ((writeFastq && !localFqPair) ||
					(!writeFastq && (!localSeq1 || !localSeq2 ||
						((localQual1 || localQual2) && (!localQual1 || !localQual2 ||
							qualTmp1.empty() || qualTmp2.empty())))))) {
					paired_sync_skipped_counter_++;
					return false;
				}
				const bool willWrite = (localFqPair && d1 && d2) || (localFq1 && d1) ||
					(localFq2 && d2) || (localSeq1 && d1) || (localSeq2 && d2);
				if (!willWrite) {
					if (isStrictPairedFastqWrite) { paired_sync_skipped_counter_++; }
					return false;
				}
				if (maxReadsPerOFile > 0 && ReadsWritten >= maxReadsPerOFile) {
					rotate = true;
				} else {
					++ReadsWritten;
					++TotalReadsWritten;
				}
			}
		}

		if (rotate) {
			epochGuard.unlock();
			{
				std::unique_lock<std::shared_mutex> rotationGuard(outputEpochMTX);
				std::lock_guard<std::mutex> incrementGuard(fileIncMTX);
				if (maxReadsPerOFile > 0 && ReadsWritten >= maxReadsPerOFile &&
					(maxRdsOut == 0 || TotalReadsWritten < maxRdsOut)) {
					incrementOutputFile();
				}
			}
			epochGuard.lock();
			continue;
		}

		if (localFqPair) {
			localFqPair->write2(fastqTmp1, fastqTmp2);
			return true;
		}
		if (mergedOutput) {
			*localFq1 << *mergedFastq;
			return true;
		}
		if (localFq1 && d1) {
			fastqTmp1.clear(); d1->writeFastQ(fastqTmp1); *localFq1 << fastqTmp1;
		}
		if (localFq2 && d2) {
			fastqTmp2.clear(); d2->writeFastQ(fastqTmp2); *localFq2 << fastqTmp2;
		}
		if (!writeFastq) {
			// Each stream's own lock protects its bytes, but all related records
			// must share this lock to keep R1/R2 and FASTA/QUAL in the same order.
			std::lock_guard<std::mutex> streamGuard(sqfqostrMTX);
			if (localSeq1 && d1) {
				*localSeq1 << seqTmp1;
				if (localQual1) { *localQual1 << qualTmp1; }
			}
			if (localSeq2 && d2) {
				*localSeq2 << seqTmp2;
				if (localQual2) { *localQual2 << qualTmp2; }
			}
		}
		return true;
	}
}

bool OutputStreamer::saveForWrite_merge(shared_ptr<DNAunique> d,
	string newHeader, int curThread, bool elseWriteD1) {
	if (!d) { return false; }
	shared_ptr<DNA> d2 = d->getPair();
	shared_ptr<DerepRead> dna_merged = d->getMerge();

	if (!dna_merged && d2 != nullptr) {
		findSeedForMerge(d, d2, curThread);
		if (d->merge_seed_pos_ >= 0) {
			shared_ptr<DNA> tdM = mergers[curThread]->merge(d, d2);
			if (tdM) { dna_merged = make_shared<DerepRead>(tdM); }
		}
	}
	// Prepare both cached and newly constructed seeds at their full length.
	if (dna_merged && !getFilters(curThread)->checkMergedDNA(dna_merged, true)) { dna_merged.reset(); }
	d->attachMerge(dna_merged);
	if (!dna_merged) {
		if (!elseWriteD1 || !getFilters(curThread)->checkMergedLength(d)) { return false; }
		dna_merged = d;
	}
	if (!newHeader.empty()) { dna_merged->setNewID(newHeader); }
	dna_merged->prepareWrite(fastQoutVer);
	if (!BWriteFastQ) {
		// The selected merged OTU belongs in the requested main FASTA/QUAL.
		// Reuse normal output accounting, rotation, and serialization.
		return writeForWriteImpl(dna_merged, 1, 0, nullptr, -1, 100, nullptr, curThread);
	}
	thread_local string mergedFastqTmp;
	mergedFastqTmp.clear();
	dna_merged->writeFastQ(mergedFastqTmp);
	return writeForWriteImpl(dna_merged, 1, 0, nullptr, -1, 100, &mergedFastqTmp, curThread);
}


void OutputStreamer::findSeedForMerge(shared_ptr<DNA> dna1, shared_ptr<DNA> dna2, int thrPos) {
    if (dna1) { dna1->clearMergeSeed(); }
    if (dna2) { dna2->clearMergeSeed(); }
    ++total_read_preMerge_;
    if (!dna1 || !dna2 || thrPos < 0 || static_cast<size_t>(thrPos) >= mergers.size() || !mergers[thrPos]) {
        return;
    }
    if (mergers[thrPos]->findSeedForMerge(dna1, dna2)) { ++merged_counter_; }
}

void OutputStreamer::writeSelectiveStream(shared_ptr<DNA> d, int Pair, int FS) {
	//ofstream tmpS, tmpQ, tmpFQ;
	if (d == 0) {
		return;
	}
	d->prepareWrite(fastQoutVer);
	int PairC = Pair;
	if (Pair > 1) {
		PairC = Pair - 2;
	}
	if (d->isGreenQual()) {
		MFil->DNAstatLQ(d, PairC, false);
	}
	else if (d->isYellowQual()) {
		MFil->DNAstatLQ(d, PairC, true);
	}
	if (BWriteFastQ && b_changeFQheadVer) {//check if header PE naming needs to be changed
		d->changeHeadPver(MFil->FQheadV());
	}
	if (BWriteFastQ) {
		*(fqFile[FS][Pair]) << d->writeFastQ();
	}
	else {
		*sFile[FS][Pair] << d->writeSeq(b_oneLinerFasta);
		if (BWriteQual) {
			*qFile[FS][Pair] << d->writeQual(b_oneLinerFasta);
		}
	}


	//	delete d;

}
void OutputStreamer::writeNonBCReads(const shared_ptr<DNA>& d, const shared_ptr<DNA>& d2) {
	if (fqNoBCFile.empty() || !d || d->getBarcodeDetected()) { return; }
	if (pairedSeq > 1 && (!d2 || d2->getBarcodeDetected())) { return; }
	string d1s, d2s;
	d->prepareWrite(fastQoutVer);
	d->writeFastQ(d1s);
	if (pairedSeq > 1) {
		d2->prepareWrite(fastQoutVer);
		d2->writeFastQ(d2s);
		if (d1s.empty() || d2s.empty()) { return; }
	}
	std::lock_guard<std::mutex> guard(dmltMTX);
	*(fqNoBCFile[0]) << d1s;
	if (!d2s.empty()) {
		*(fqNoBCFile[fqNoBCFile.size() == 1 ? 0 : 1]) << d2s;
	}
}

void OutputStreamer::setSubfilters(int num) {
	if (num < 1) { return; }
	subFilter.resize(num);
	for (uint i = 0; i < subFilter.size(); i++) {
		subFilter[i].reset();
		subFilter[i] = std::make_unique<Filters>(MFil, MFil->currentBCnumber(), true);
	}
}



void OutputStreamer::mergeSubFilters() {
	vector<int> idx(MFil->Barcode.size(), 0);
	for (uint i = 0; i < idx.size(); i++) { idx[i] = i; }
	for (uint i = 0; i < subFilter.size(); i++) {
		MFil->addStats(subFilter[i].get(), idx);
	}
}
/*
void OutputStreamer::mergeSubFiltersMT() {
	vector<int> idx (MFil->Barcode.size(),0);

	for (uint i=0;i<idx.size();i++){idx[i]=i;}
	for (uint i=0; i<subFilter.size();i++){
//        MFil->addStats(subFilter[i],idx);
		StatisticsMultithreaded::mergeStatisticsInto(MFil->collectStatistics,
			MFil->statAddition, subFilter[i]->statistics_, idx);
	}
}
*/

void OutputStreamer::attachDereplicator(shared_ptr<Dereplicate> de) {
	if (de != NULL) {
		dereplicator = de; b_doDereplicate = true;
		dereplicator->setPaired(pairedSeq > 1);
		// Withheld mid-quality sets the dereplicator does not count are written
		// here, through the per-record (locked) path since any thread may apply them.
		dereplicator->setUncountedPublisher([this](const shared_ptr<DNA>& d1, const shared_ptr<DNA>& d2,
			const Dereplicate::DerepTierRoute& route) {
			writeForWrite(d1, route.idx1, route.cstr1, d2, route.idx2, route.cstr2, -1);
		});
		//insert code here to fix BC offset in filter & add to derep info on sample names from the current filter
		int curBCOffset = dereplicator->getHighestBCoffset();

     cdbg("BARCODE INCREMENT: " + std::to_string(curBCOffset) + "\n");

		MFil->setBCoffset(curBCOffset);
		for (size_t i = 0; i < subFilter.size(); i++) {
			subFilter[i]->setBCoffset(curBCOffset);
		}
		dereplicator->BCnamesAdding(MFil);

		if (dereplicator->mergeDereRead()) {
			b_merge_pairs_ = true;
		}

	}
}

/*void OutputStreamer::dereplicateDNA(shared_ptr<DNA> tdn) {

	if (b_doDereplicate && tdn->isGreenQual()) {
		cntDerep++;
		dereplicator->addDNA(tdn);
	}
}*/
int OutputStreamer::outputTier(const shared_ptr<DNA>& d, int thr) {
	if (!d || suppressOutWrite == 3) { return 100; }
	Filters* curFil = this->getFilters(thr);
	if (!curFil->doFilterAtAll()) { return 0; }
	if (b_writeGreenQual && d->isGreenQual()) { return 0; }
	if (b_writeYellowQual && d->isYellowQual()) { return 1; }
	return 100;
}

bool OutputStreamer::dereplicateDNA(const shared_ptr<DNA>& tdn, const shared_ptr<DNA>& tdn2, int curThread,
	int idx1, int idx2) {
	if (!b_doDereplicate || !tdn) {
		return false;
	}
	// Sets bound for the mid-quality tier are counted or written, never both:
	// the dereplicator publishes them through the attached callback if it
	// does not count them.
	Dereplicate::DerepTierRoute route{idx1, outputTier(tdn, curThread), idx2, outputTier(tdn2, curThread)};
	const bool midTier = route.cstr1 == 1 || route.cstr2 == 1;
	bool queued = false;
	dereplicator->addDNA(tdn, tdn2, getFilters(curThread), midTier ? &route : nullptr, &queued);
	if (tdn->isDereplicated()) {
		cntDerep++;
	}
	return midTier && queued;
}

void OutputStreamer::flushDereplication() {
	if (dereplicator) { dereplicator->flushPending(); }
}

void OutputStreamer::resetOutFilesAndFilter() {
#ifdef DEBUG
	cerr << "resetOutFilesAndFilter" << endl;
#endif
	//reset count stats
	MFil->resetStats();
	// The worker statistics are merged into MFil later; reset them too.
	for (auto& sub : subFilter) { if (sub) { sub->resetStats(); } }
	//reset all stored DNA
	delAllDNAvectors();
	dereplicator.reset();
	//close streams -- why? no reason, since nothing has been writen
	//this->closeOutStreams();
	totalFileStrms = 0;
}

// Internal implementation: caller must already hold sqfqostrMTX.
// Every stream is finalized explicitly (buffers are written in any case, as the
// destructors did before): a delayed write failure such as a full disk is then
// reported to the caller instead of being swallowed by ~ofbufstream. All
// streams are closed before the first failure is rethrown. Retained streams are
// only flushed here and finalized by OutputStreamRetention::closeAll().
void OutputStreamer::closeOutStreams_locked(bool wr) {
	(void)wr;
	write2File = true;
	std::exception_ptr failure;
	auto remember = [&]() { if (!failure) { failure = std::current_exception(); } };
	auto closeStream = [&](std::shared_ptr<ostr>& stream) {
		if (!stream) { return; }
		try { stream->deactivate(); }
		catch (...) { remember(); }
		stream.reset();
	};
	auto closePair = [&](std::shared_ptr<dualOfBufStream>& stream) {
		if (!stream) { return; }
		try {
			stream->emptyStreams(true);
			stream->deactivate();
		}
		catch (...) { remember(); }
		stream.reset();
	};

	for (auto& tier : sFile) { for (auto& stream : tier) { closeStream(stream); } }
	for (auto& tier : qFile) { for (auto& stream : tier) { closeStream(stream); } }
	for (auto& tier : fqFile) { for (auto& stream : tier) { closeStream(stream); } }
	for (auto& stream : of_merged_fq) { closeStream(stream); }
	for (auto& stream : fqNoBCFile) { closeStream(stream); }
	for (auto& stream : fqPairFile) { closePair(stream); }
	for (auto& stream : sPairFile) { closePair(stream); }
	for (auto& stream : qPairFile) { closePair(stream); }

	if (MFil != nullptr) {
		MFil->setWrittenReads(ReadsWritten.load());
	}
#ifdef DEBUG
	cerr << ".. closed\n";
#endif
	if (failure) { std::rethrow_exception(failure); }
}

// Public wrapper that acquires sqfqostrMTX then calls the locked implementation
void OutputStreamer::closeOutStreams(bool wr) {
	endWriteBatching();
	std::lock_guard<std::mutex> guard(sqfqostrMTX);
	closeOutStreams_locked(wr);
}

void OutputStreamer::openOFstream(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep, int wh) {
	switch (wh) {
	case 1:
		openOFstreamFNA(opOF, wrMode, p1, p2, errMsg, onlyPrep);
		break;
	case 0:
		openOFstreamFQ(opOF, wrMode, p1, p2, errMsg, onlyPrep);
		break;
	case 2:
		openOFstreamQL(opOF, wrMode, p1, p2, errMsg, onlyPrep);
		break;
	default:
		cerr << "Wrong wh specified"; exit(1002);
	}
}
void OutputStreamer::openNoBCoutstrean(const string inS) {
	const vector<string> paths = splitByCommas(inS);
	if (paths.empty() || paths.size() > (pairedSeq > 1 ? 2u : 1u) ||
		std::any_of(paths.begin(), paths.end(), [](const string& path) { return path.empty(); })) {
		throw std::runtime_error("No-barcode output requires one path for single reads, or one or two paths for paired reads.");
	}
	std::lock_guard<std::mutex> guard(sqfqostrMTX);
	fqNoBCFile.clear();
	for (const string& path : paths) {
		fqNoBCFile.push_back(OutputStreamRetention::open(path, wrMode, doTIO));
	}
}

void OutputStreamer::openOFstreamFQpair(const string opOF, const string opOF2, std::ios_base::openmode wrMode,
	int p1, string errMsg, bool onlyPrep) {
	//bool doMC = Nthrds >= 1;
	assert(p1 < 2);
	while (fqPairFile.size() <= p1) {
		// create unique_ptr entries
		fqPairFile.emplace_back();
	}
    {
		std::lock_guard<std::mutex> guard(sqfqostrMTX);
		if (!fqPairFile[p1]) {
			fqPairFile[p1] = std::make_shared<dualOfBufStream>();
		}
		if (!fqPairFile[p1]->open(opOF, wrMode, 0, doTIO, 250000)) {
		cerr << "Could not open " << errMsg << " pair 1fastq output file " << opOF << endl << p1 << " " << totalFileStrms << endl;
		exit(5);
	}
        if (!fqPairFile[p1]->open(opOF2, wrMode, 1, doTIO, 350000)) {
		cerr << "Could not open " << errMsg << " pair 2 fastq output file " << opOF2 << endl << p1 << " " << totalFileStrms << endl;
		exit(5);
	}
        fqPairFile[p1]->activate();
	}
}

void OutputStreamer::openOFstreamFQ(const string opOF, std::ios_base::openmode wrMode,
	int p1, int p2, string errMsg, bool onlyPrep) {
	//p1: 0:green, 1:yellow
	//p2: 0:pair1, 1:pair2, 2:single1, 3:single2
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1 + 1 >= (int)fqFile.size()) {
		fqFile.resize(p1 + 1);
	}
	//bool doMC = Nthrds > 1;

    {
		std::lock_guard<std::mutex> guard(sqfqostrMTX);
		if (fqFile[p1].size() < 4) fqFile[p1].resize(4);
		fqFile[p1][p2].reset();
		fqFile[p1][p2] = OutputStreamRetention::open(opOF, wrMode, doTIO);
		if (!onlyPrep) {
			fqFile[p1][p2]->activate();
		}
		if (!*fqFile[p1][p2]) {
			cerr << "Could not open " << errMsg << " fastq output file " << opOF << endl << p1 << " " << p2 << " " << totalFileStrms << endl;
			exit(4);
		}
	}
}
void OutputStreamer::openOFstreamFQ_mrg(const string opOF, std::ios_base::openmode wrMode,
	int p1, string errMsg, bool onlyPrep) {
	//also set up potential stream for merge out
	if (p1 + 1 >= (int)of_merged_fq.size()) {
		of_merged_fq.resize(p1 + 1);
	}
	//bool doMC = Nthrds > 1;
    {
		std::lock_guard<std::mutex> guard(sqfqostrMTX);
		of_merged_fq.resize(max((size_t)1, (size_t)p1 + 1));
		of_merged_fq[p1] = OutputStreamRetention::open(opOF, wrMode, doTIO);
		if (!*of_merged_fq[p1]) {
			cerr << "Could not open " << errMsg << " fastq merged output file " << opOF << endl << p1 << " " << totalFileStrms << endl;
			exit(4);
		}
	}
}
void OutputStreamer::openOFstreamFNA(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep) {
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1 + 1 >= (int)sFile.size()) {
		sFile.resize(p1 + 1);
	}
	/*if ((int)sFileStr.size() - 1 <= p1) {
		sFileStr.resize(p1 + 1, vector<string>(4, ""));
	}
	sFileStr[p1][p2] = opOF;*/
	//bool doMC = Nthrds > 1;

    {
		std::lock_guard<std::mutex> guard(sqfqostrMTX);
		if (sFile[p1].size() < 4) sFile[p1].resize(4);
		sFile[p1][p2].reset();
		sFile[p1][p2] = OutputStreamRetention::open(opOF, wrMode, doTIO);
		if (!*sFile[p1][p2]) {
			cerr << "Could not open " << errMsg << " fasta output file " << opOF << endl << p1 << " " << p2 << " " << totalFileStrms << endl;
			exit(4);
		}
	}
	//	if (onlyPrep ) { return; }
		//if (p1 == 1 && !b_writeYellowQual){ return; }//p1==1: mid passed suppressOutWrite >= 2
		//if (p1 == 0 && !b_writeGreenQual){ return; }//suppressOutWrite == 1

}
void OutputStreamer::openOFstreamQL(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool onlyPrep) {
	if (p2 > 3) { cerr << "internal error: can't have more than 4 entries in output file stream\n"; exit(1001); }
	if (p1 + 1 >= (int)qFile.size()) {
		qFile.resize(p1 + 1);
	}
	//bool doMC = Nthrds > 1;

	//if (onlyPrep) { return; }
	//if (p1 == 1 && !b_writeYellowQual){ return; }//p1==1: mid passed suppressOutWrite >= 2
	//if (p1 == 0 && !b_writeGreenQual){ return; }//suppressOutWrite == 1
    {
		std::lock_guard<std::mutex> guard(sqfqostrMTX);
		if (qFile[p1].size() < 4) qFile[p1].resize(4);
		qFile[p1][p2].reset();
		qFile[p1][p2] = OutputStreamRetention::open(opOF, wrMode, doTIO);
		if (!*qFile[p1][p2]) {
			cerr << "Could not open " << errMsg << " quality output file " << opOF << endl;
			exit(4);
		}
	}

}

void OutputStreamer::openSeveralOutstreams(OptContainer* cmdArgs, shared_ptr<ReadSubset> RDS, std::ios_base::openmode wrMode) {
#ifdef DEBUG
	cerr << " openining multiple out streams" << endl;
#endif
	string path = "", fileEnd(".fna");
	vector<string> outFile = RDS->getOFiles();
	bool openStrms = true;
	int omode(1);
	if (cmdArgs->find("-o_fastq") != cmdArgs->end() && (*cmdArgs)["-o_fastq"] != "" && (*cmdArgs)["-o_fna"] == "") { //write fastq
		BWriteFastQ = true;
		path = (*cmdArgs)["-o_fastq"];
		omode = 0;
		fileEnd = ".fastq";
	}
	else  if ((*cmdArgs)["-o_fna"] != "") {
		BWriteFastQ = false;
		path = (*cmdArgs)["-o_fna"];
	}
	else {
		cerr << "Could not find valid sequence output file path. Exiting\n";
		exit(55);
	}
	if (RDS->multiFile()) {
		b_multiOutStream = true;
	}
	int i(0);
	for (i = 0; i < (int)outFile.size(); i++) {
		if (totalFileStrms >= MAX_FILE_STREAMS && openStrms) {
			cerr << "Too many output file streams\nSwitching to dynamical file appending\n";
			openStrms = false;
			//this->closeOutStreams();
		}
		string baseFile = path + removeFileEnding(outFile[i]);
		if (pairedSeq > 1) {
			if (i < 4) {
				cerr << "Outfile " << i << ": " << baseFile + ".1" << fileEnd << "\n";
			}
			openOFstream(baseFile + ".1" + fileEnd, wrMode, i, 0, "paired 1st " + itos(i), !openStrms, omode);
			//2nd pair_
			openOFstream(baseFile + ".2" + fileEnd, wrMode, i, 1, "paired 2nd " + itos(i), !openStrms, omode);
			//1st singleton
			//if (openStrms) { openOFstream(baseFile + ".1" + fileEnd + SingletonFileDescr, wrMode, i, 2, "Singleton 1 " + itos(i), omode); }
			//2nd singleton
			//if (openStrms) { openOFstream(baseFile + ".2" + fileEnd + SingletonFileDescr, wrMode, i, 3, "Singleton 2 " + itos(i), omode); }
			totalFileStrms += 2;
		}
		else {
			openOFstream(baseFile + fileEnd, wrMode, i, 0, "main " + itos(i), !openStrms, omode);
			totalFileStrms++;
		}

	}
	//pipe remaining reads into this file
	if ((*cmdArgs)["-excludeFile"] != "") {
		//set a flag
		RDS->setRemainingFilepipe(i);
		if (pairedSeq == 1) {
			openOFstream((*cmdArgs)["-excludeFile"], wrMode, i, 0, "excludeFile file ", !openStrms, omode);
		}
		else {
			string baseFile = path + removeFileEnding((*cmdArgs)["-excludeFile"]);
			if (i < 4) {
				cerr << "out file " << i << baseFile + ".1" << fileEnd << "\n";
			}
			openOFstream(baseFile + ".1" + fileEnd, wrMode, i, 0, "paired 1st " + itos(i), !openStrms, omode);
			//2nd pair_
			openOFstream(baseFile + ".2" + fileEnd, wrMode, i, 1, "paired 2nd " + itos(i), !openStrms, omode);
			//1st singleton
			//if (openStrms) { openOFstream(baseFile + ".1" + fileEnd + SingletonFileDescr, wrMode, i, 2, "Singleton 1 " + itos(i), omode); }
			//2nd singleton
			//if (openStrms) { openOFstream(baseFile + ".2" + fileEnd + SingletonFileDescr, wrMode, i, 3, "Singleton 2 " + itos(i), omode); }
			totalFileStrms += 2;
		}
	}
}

bool OutputStreamer::singleSeedFasta() const {
	const auto path = locCmdArgs->find("-o_fna");
	return MFil->doOptimalClusterSeq() && !BWriteFastQ && path != locCmdArgs->end() &&
		!path->second.empty() && splitByCommas(path->second).size() == 1;
}

void OutputStreamer::openOutStreams(OptContainer* cmdArgs, int fileIt,
	std::ios_base::openmode wrMode_i,
	string fileExt, int forceFmt) {
	this->setwriteMode(wrMode_i);
	if (suppressOutWrite == 3 || ((*cmdArgs)["-o_fastq"] == "" && (*cmdArgs)["-o_fna"] == "" && (*cmdArgs)["-o_qual"] == "")) {
		suppressOutWrite = 3; b_writeGreenQual = false;
		b_writeYellowQual = false; return;
	}
	if (forceFmt != -1) {
		if (forceFmt == 1 && (cmdArgs->find("-o_fna") == cmdArgs->end() || (*cmdArgs)["-o_fna"] == "")) {//force fna, no qual, required for seed ref fastas
			(*cmdArgs)["-o_fna"] = (*cmdArgs)["-o_fastq"];
			(*cmdArgs)["-o_fastq"] = "";
		}
	}


	if (cmdArgs->find("-o_fastq") != cmdArgs->end() && (*cmdArgs)["-o_fastq"] != "" && (*cmdArgs)["-o_fna"] == "") { //write fastq
		this->setFastQWrite(true);
		if (pairedSeq > 1) { //open second pair_ + singleton
#ifdef DEBUG
			cerr << " paired fastq out " << endl;
#endif
			vector<string> tfnaout = splitByCommas((*cmdArgs)["-o_fastq"]);
			if (tfnaout.size() != 2) {
				cerr << "Paired sequences given as input, requires paired output file (2 files separated by \",\"). Given output file = " << (*cmdArgs)["-o_fastq"] << endl; exit(57);
			}
			leadingOutf = applyFileIT(tfnaout[0] + fileExt, fileIt);
			//now save pairs in dedicated object:
			openOFstreamFQpair(leadingOutf.c_str(), applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 0, "paired out");

			/*
			//1st pair (or main) file out
			openOFstreamFQ(leadingOutf.c_str(), wrMode, 0, 0, "paired 1st");
			//2nd pair
			openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 0, 1, "paired 2nd");
			//1st singleton
			*/
			openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 0, 2, "Singleton 1", true);
			//2nd singleton
			openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 0, 3, "Singleton 2", true);

			//merged reads
			int prefL = (int)len_common_prefix_base(tfnaout[0].c_str(), tfnaout[1].c_str());
			string prefix = tfnaout[0].substr(0, prefL);
			openOFstreamFQ_mrg(applyFileIT(prefix + "merg.fq" + fileExt, fileIt).c_str(), wrMode, 0, "merged");

			//additional file
			if (cmdArgs->find("-o_fastq2") != cmdArgs->end() && (*cmdArgs)["-o_fastq2"].length() > 1) { //write fastq
				tfnaout = splitByCommas((*cmdArgs)["-o_fastq2"]);
				if (tfnaout.size() != 2) {
					cerr << "Paired sequences given as input, requires paired additional output file (2 files separated by \",\"). Given -o_fastq2 = " << (*cmdArgs)["-o_fastq2"] << endl; exit(57);
				}
				openOFstreamFQpair(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 1, "additional paired out");
				/*openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional paired 1st");
				openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt).c_str(), wrMode, 1, 1, "additional paired 2nd");
				*/
				openOFstreamFQ(applyFileIT(tfnaout[0] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 1, 2, "additional Singleton 1", true);
				openOFstreamFQ(applyFileIT(tfnaout[1] + fileExt, fileIt, SingletonFileDescr).c_str(), wrMode, 1, 3, "additional Singleton 2", true);
			}
			else { b_writeYellowQual = false; }
		}
		else {
#ifdef DEBUG
			cerr << " single fastq out " << endl;
#endif
			openOFstreamFQ(applyFileIT((*cmdArgs)["-o_fastq"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "the main");
			leadingOutf = applyFileIT((*cmdArgs)["-o_fastq"] + fileExt, fileIt);
			//additional file
			if (cmdArgs->find("-o_fastq2") != cmdArgs->end() && (*cmdArgs)["-o_fastq2"].length() > 1) { //write fastq
				openOFstreamFQ(applyFileIT((*cmdArgs)["-o_fastq2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "the additional");
			}
			else { b_writeYellowQual = false; }
		}

		return;
	}
	BWriteFastQ = false;
	if (pairedSeq == 1 || singleSeedFasta()) {
#ifdef DEBUG
		cerr << " fasta singleton output " << endl;
#endif
		openOFstreamFNA(applyFileIT((*cmdArgs)["-o_fna"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "main");
		leadingOutf = applyFileIT((*cmdArgs)["-o_fna"] + fileExt, fileIt);
		if ((*cmdArgs)["-o_qual"] != "") {
			openOFstreamQL(applyFileIT((*cmdArgs)["-o_qual"] + fileExt, fileIt).c_str(), wrMode, 0, 0, "main");
			this->setQualWrite(true);
		}
		else {
			this->setQualWrite(false);
		}
		//additional file (secondary filter)
		if (cmdArgs->find("-o_fna2") != cmdArgs->end() && (*cmdArgs)["-o_fna2"].length() > 1) { //add file
			openOFstreamFNA(applyFileIT((*cmdArgs)["-o_fna2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional");
			if (cmdArgs->find("-o_qual2") != cmdArgs->end() && (*cmdArgs)["-o_qual2"].length() > 1) { //add file
				openOFstreamQL(applyFileIT((*cmdArgs)["-o_qual2"] + fileExt, fileIt).c_str(), wrMode, 1, 0, "additional");
				this->setQualWrite(true);
			}
		}
		else { b_writeYellowQual = false; }

	}
	else {
#ifdef DEBUG
		cerr << " fasta paired output " << endl;
#endif

		// Open each mate and its singleton through the existing stream helpers,
		// including matching QUAL files when requested.
		auto openPairedFiles = [&](const string& paths, int tier, int format) {
			const vector<string> files = splitByCommas(paths);
			if (files.size() != 2 || files[0].empty() || files[1].empty()) {
				throw std::runtime_error("Paired FASTA/QUAL output requires two comma-separated paths: " + paths);
			}
			for (int mate = 0; mate < 2; ++mate) {
				openOFstream(applyFileIT(files[mate] + fileExt, fileIt), wrMode,
					tier, mate, "paired", false, format);
				openOFstream(applyFileIT(files[mate] + fileExt + SingletonFileDescr, fileIt), wrMode,
					tier, mate + 2, "singleton", true, format);
			}
			if (tier == 0 && format == 1) { leadingOutf = applyFileIT(files[0] + fileExt, fileIt); }
		};
		openPairedFiles((*cmdArgs)["-o_fna"], 0, 1);
		setQualWrite(false);
		if (!(*cmdArgs)["-o_qual"].empty()) {
			openPairedFiles((*cmdArgs)["-o_qual"], 0, 2);
			setQualWrite(true);
		}
		if (!(*cmdArgs)["-o_fna2"].empty()) {
			openPairedFiles((*cmdArgs)["-o_fna2"], 1, 1);
			if (!(*cmdArgs)["-o_qual2"].empty()) {
				openPairedFiles((*cmdArgs)["-o_qual2"], 1, 2);
				setQualWrite(true);
			}
		}
		else { b_writeYellowQual = false; }
	}
}
