#pragma once


#include "DNAconsts.h"
#include "InputStream.h"
#include "Filters.h"
#include "ReadMerger.h"
#include <cstring>
//#include "ThreadPool.h"
#include <fstream>
#include <mutex>
#include <future>
#include <memory>
#include <vector>
#include <atomic>
#include <condition_variable>
#include <deque>
#include <thread>
#include <shared_mutex>
#include <type_traits>
#include <limits>
#include <map>

class AsyncOutputWriter;

// Output buffer bytes are always written before they are read, so creating or
// growing a buffer need not zero it (megabytes per output stream otherwise).
template <class T>
struct UninitializedAllocator : std::allocator<T> {
	UninitializedAllocator() = default;
	template <class U> UninitializedAllocator(const UninitializedAllocator<U>&) noexcept {}
	template <class U> void construct(U* p) noexcept(std::is_nothrow_default_constructible_v<U>) { ::new (static_cast<void*>(p)) U; }
	template <class U, class... Args> void construct(U* p, Args&&... args) { ::new (static_cast<void*>(p)) U(std::forward<Args>(args)...); }
};
using OutputBuffer = std::vector<char, UninitializedAllocator<char>>;

//static mutex output_mtx;
// Buffered output file. Appends only copy into the buffer under output_mtx.
// With multithreaded I/O, full buffers are handed in order to one shared
// writer thread (AsyncOutputWriter), so no appender ever waits for write()
// or BGZF compression while it holds this stream's (or a caller's) lock.
class ofbufstream {//: private std::streambuf, public std::ostream {
public:
    ofbufstream() :file("T"), keeper(0), modeIO(ios::app), used(0),
     coutW(true), isGZ(false), doMC(false), primary(nullptr), bufS(5* 1024 * 1024) {}
	ofbufstream(size_t bufferS) :file("T"), keeper(0), modeIO(ios::app), used(0),
		coutW(true), isGZ(false), doMC(false), primary(nullptr), bufS(bufferS) {

    }
    ofbufstream(const string IF, int mif, bool isMC = false, size_t bufferS = 5* 1024 * 1024);
    ~ofbufstream();
    void finishWrites();
    bool operator! (void);
    void operator<< (const string& X);
    void append(const char* data, size_t size);
    // Hands buffered bytes to the file and waits until they were written.
    void emptyStream();
    void activate();
    // Flushes and finalizes the file; a retained stream is only flushed.
    void deactivate();
    // Retained streams stay open when their OutputStreamer closes them; they
    // are finalized by OutputStreamRetention::closeAll().
    void setRetained(bool retained) { retained_ = retained; }
    const string& fileName() const { return file; }
    // end

private:
    friend class AsyncOutputWriter;
    mutex output_mtx;
	bool shouldUseAsyncWriter() const;
	// Caller holds output_mtx: queues the filled keeper and installs a spare buffer.
	void handOffLocked();
	// Writer thread: writes one queued buffer and recycles it.
	void writeQueued(OutputBuffer& buffer, size_t size);
	void releaseQueued(OutputBuffer&& buffer);
	void waitForQueuedWrites();

    //functions
    string file;
    OutputBuffer keeper;
    int modeIO;
    size_t used;
    //write to cout? gzip input? do multicore?
    bool coutW, isGZ, doMC;
    std::unique_ptr<std::ostream> primary;


    size_t bufS;

    bool asyncWriterEnabled = false;
    bool retained_ = false;
    // Queued-buffer bookkeeping shared with the writer thread.
    std::mutex asyncMtx_;
    std::condition_variable asyncDone_;
    size_t inflight_ = 0;
    std::vector<OutputBuffer> spare_;
    std::exception_ptr asyncError_;
};

typedef ofbufstream ostr;

// Keeps main output streams open for the whole run instead of finalizing and
// re-opening them (in append mode) for every input file. Streams opened while
// retention is disabled behave exactly as before.
namespace OutputStreamRetention {
	void enable(bool on);
	bool enabled();
	std::shared_ptr<ofbufstream> open(const string& path, int mode, bool isMC, size_t bufferS = 5 * 1024 * 1024);
	// Finalizes every retained stream; rethrows the first failure after closing all.
	void closeAll();
}

class Filters;
class Dereplicate;
class ReadSubset;
class Benchmark;




class dualOfBufStream {
public:
    dualOfBufStream(void);
    ~dualOfBufStream(void);
    void write(const string& in, int stream);
    void write2(const string& in, const string& in2);
    bool open(const string IF, int mif, int pair, bool isMC = false, size_t bufferS = 1024 * 1024);
    bool activate();
    bool deactivate();
    void emptyStreams(bool force);
private:
    size_t buf1S, buf2S;
    vector<string> bufs;
    vector<string> FileNames;
	vector<std::shared_ptr<ofbufstream>> dualOutStr;
    vector<bool> opened;
    bool active;
    mutex dualMtx;

};




//writes successful demultis and stores unsuccessful matches of fna/qual_ for later matching
class OutputStreamer {
public:
	//wrStatus controls if this appends or overwrites output
	OutputStreamer(Filters* filters, OptContainer* cmdArgs,
		std::ios_base::openmode wrStatus, shared_ptr<ReadSubset>,
		int numThreads, string fileExt = "", int = -1);
	~OutputStreamer();
	//clean up streams
	void closeOutFilesDemulti();
	void closeOutStreams(bool wr = true);


	//	void threadAnalyzeDNA(shared_ptr<DNA>);
	void setFastQWrite(bool x) { BWriteFastQ = x; BWriteQual = !x; }
	void setQualWrite(bool x) { BWriteQual = x; }
	//void addNoHeadDNA(shared_ptr<DNA> d) { DNAsNoHead.push_back(d); }
	//-1,-1,-2
	void analyzeDNA(const shared_ptr<DNA>& d, int FilterUse, int pair, int& idx, int thr);
	vector<bool> analyzeDNA(shared_ptr<DNA> p1, shared_ptr<DNA> p2, shared_ptr<DNA> mid, bool changePHead, int = -1);
	void findSeedForMerge(shared_ptr<DNA> dna1, shared_ptr<DNA> dna2, int thrPos);

	//void writeAndDel(shared_ptr<DNA> d, int p=1) { writeAndDel(d.get(), p); }
	//void writeAndDel(shared_ptr<DNA> d, int Pair = 1);//1=pair1;2=pair2;3=singleton1,4=singl2
	//Function specifically if several output files are required
	void writeSelectiveStream(shared_ptr<DNA> d, int Pair, int FS);//1=pair1;2=pair2;3=singl1,4=singl2  ;; FS: different multi FileStreams to be used

	 void writeForWrite(const shared_ptr<DNA>& d1, int Pair1, int Cstream1,
		const shared_ptr<DNA>& d2, int Pair2, int Cstream2, int curThread = -1);
	// Per-worker output batching. While enabled, records written by worker
	// `curThread` are staged in its own buffers and published once per block
	// by flushWriteBatch(), instead of taking the output locks for every read.
	// Without ordering, disabled when output rotation or read limits need
	// per-record accounting (with it, records are counted at publication).
	// With orderWindow > 0 the blocks are published in input order: the reader
	// numbers them (reserveOutputBlock) and the worker hands its block over
	// (finishOutputBlock); a block's output is its records, then the mid-quality
	// sets the dereplicator withheld and did not count, in record order. So the
	// files do not depend on the number of threads or on scheduling. A worker
	// whose block is orderWindow blocks or more ahead of the next one to publish
	// waits, which bounds the output held back.
	void beginWriteBatching(size_t orderWindow = 0);
	void flushWriteBatch(int curThread);
	void endWriteBatching();
	static constexpr uint64_t kNoOutputBlock = std::numeric_limits<uint64_t>::max();
	uint64_t reserveOutputBlock();
	// Publishes (or queues) the block's output; without a block number as flushWriteBatch.
	void finishOutputBlock(int curThread, uint64_t block);
	// After an error: nothing more is published, and waiting workers return.
	void failOutputBlocks();
	//pretty final bool, aborts all, so careful with this
	//collects stats on read, writes then. Keep read pairs together by using "writeForWrite"in second mutext step
	bool saveForWrite(const shared_ptr<DNA>& d, int Pair, int thr, int& Cstr, bool = true);//Pair:1=pair1;2=pair2;3=singleton
	bool saveForWrite_merge(shared_ptr<DNAunique> d,
		string newHeader = "", int curThread = 0, bool elseWriteD1 = false);
	//bool saveForWriteMT(shared_ptr<DNA> dna, int thread, int pair = 1);
	Filters* getFilters(size_t w = (size_t)-1) {
		if (w == (size_t)-1) { return MFil; }
		if (w < subFilter.size()) { return subFilter[w].get(); }
		return MFil;
	}
	int isPEseq() { return pairedSeq; }
	bool singleSeedFasta() const;
	//ofstream::app, ios_base::out
	void openOutStreams(OptContainer* cmdArgs, int, std::ios_base::openmode, string = "", int = -1);
	void openSeveralOutstreams(OptContainer* cmdArgs, shared_ptr<ReadSubset>, std::ios_base::openmode);
	string leadOutFile() { return leadingOutf; }
	//void setfastQver(int x){fastQver = x;}
	//void setfastQoutVer(int x){fastQoutVer = x;}

	bool checkFastqHeadVersion(shared_ptr<DNA> d, bool = false);
	//int getFastqMod(){return MFil->FastqModifier();}
	int getFastqVer() { return fastQver; }
	int getfastQoutVer() { return fastQoutVer; }
	void resetOutFilesAndFilter();//MD->closeOutStreams();
	void setBCfixed(bool b, bool fwd) { MFil->setBCfixed(b, fwd); write2File = b; }
	void setSubfilters(int num);
	void mergeSubFilters();
	void activateWrite2File() { write2File = true; }
	void setOneLinerFastaFmt(bool b) { b_oneLinerFasta = b; }
	//void revConstellationCnts(int x) { MFil->revConstellationCnts(x); }
	//dereplication of DNA seqs
	void attachDereplicator(shared_ptr<Dereplicate> de);
	//void dereplicateDNA(shared_ptr<DNA>);
	//dereplicate DNA by looking in dereplicator for 100% id seq. is treadsafe
	// idx1/idx2 are the output slots the set would be written to. Returns true
	// when the set's mid-quality write is withheld: the dereplicator writes it
	// later if it does not count it, so the caller must not write it.
	bool dereplicateDNA(const shared_ptr<DNA>&, const shared_ptr<DNA>&, int curThread, int idx1, int idx2);
	// Output tier of a record as saveForWrite selects it: 0 primary, 1 mid-quality, 100 none.
	int outputTier(const shared_ptr<DNA>& d, int thr);
	// Applies queued dereplication entries, publishing withheld sets, while
	// the output streams are still open.
	void flushDereplication();
	//debug function to look closer at nonBC reads, threadsafe
    void writeNonBCReads(const shared_ptr<DNA>& d, const shared_ptr<DNA>& d2, int curThread = -1);
	void setReadLimit(uint64_t x) { maxRdsOut = x; }

	//demultiplex related
 void write2Demulti(const shared_ptr<DNA>&, int, int BCoffset, int curThread = -1);
	//write DNA to demultiplexed files, is threadsafe
  void write2Demulti(const shared_ptr<DNA>&, const shared_ptr<DNA>&, int BCoffset, int curThread);
	void generateDemultiOutFiles(string, Filters*, std::ios_base::openmode = ios::out, bool = false);

	bool mergeReads() {
		return b_merge_pairs_;
	}

	atomic_size_t merged_counter_ = 0;
	atomic_size_t total_read_preMerge_ = 0;
	atomic_size_t paired_sync_skipped_counter_ = 0;

	uint getDemultiBPperSR() { return demultiBPperSR; }
	void resetDemultiBPperSR() { demultiBPperSR = 0; }
	void setBPwrittenInSR(uint64_t x) { BPwrittenInSR = x; }
	uint64_t getBPwrittenInSR(void) { return BPwrittenInSR.load(); }
	void setBPwrittenInSRmerg(uint64_t x) { BPwrittenInSRmerg = x; }
	uint64_t getBPwrittenInSRmerg(void) { return BPwrittenInSRmerg.load(); }
	uint64_t getReadsWritten() { return TotalReadsWritten.load(); }

	void attachBenchmark(Benchmark* bench) { _benchmark = bench; }
	void addMergeStatsTo(ReadMerger& target) {
		for (const auto& merger : mergers) { target.addRMstats(merger.get()); }
	}
	void activateReadMerger(int sz, bool collectStats = false) {
       if (sz > (int)mergers.size()) { mergers.resize(sz); }
		for (size_t x = 0; x < (size_t)sz; x++) {
			mergers[x].reset();
			mergers[x] = std::make_shared<ReadMerger>(collectStats);
		}
	}
	//void attachReadMerger(ReadMerger* merg) { if (merger.size() > 0 && merger.back() == nullptr) { merger.back() = merg; } else { merger.push_back(merg); } }
	bool Demulti2Fls() { return bDoDemultiplexIntoFiles; }
	bool doDeriplicate() { return	b_doDereplicate; }
	bool doWriteNonBCrds() { return !fqNoBCFile.empty(); } // one path (single/interleaved) or R1,R2
	bool doDerepMrgSrch() { if (dereplicator != nullptr) { return dereplicator->doSearchWithMerge(); } else { return false; } }
	// Ordered dereplication admission: blocks are numbered by the reader in
	// input order and bracketed by the worker that processes them.
	uint64_t reserveDerepBlock() {
		return (b_doDereplicate && dereplicator) ? dereplicator->reserveBlock() : Dereplicate::kNoBlock;
	}
	void beginDerepBlock(uint64_t block, uint64_t outputBlock = kNoOutputBlock) {
		if (block != Dereplicate::kNoBlock && dereplicator) { dereplicator->beginBlock(block, outputBlock); }
	}
	void endDerepBlock(uint64_t block) {
		if (block != Dereplicate::kNoBlock && dereplicator) { dereplicator->endBlock(); }
	}


private:
	bool writeDemultiRecords(int idx, const string& first, const string& second, uint64_t bases, bool merged,
		int curThread = -1);
	bool writeForWriteImpl(const shared_ptr<DNA>& d1, int Pair1, int Cstream1,
		const shared_ptr<DNA>& d2, int Pair2, int Cstream2, const string* mergedFastq = nullptr,
		int curThread = -1);
	struct WriteBatch;
	// Stream slots selected for one record; callers either copy the shared_ptrs
	// (per-record path, which may race with rotation) or use them in place.
	struct StreamSlots {
		const std::shared_ptr<dualOfBufStream>* fqPair = nullptr;
		const std::shared_ptr<ostr>* fq1 = nullptr;
		const std::shared_ptr<ostr>* fq2 = nullptr;
		const std::shared_ptr<ostr>* seq1 = nullptr;
		const std::shared_ptr<ostr>* seq2 = nullptr;
		const std::shared_ptr<ostr>* qual1 = nullptr;
		const std::shared_ptr<ostr>* qual2 = nullptr;
	};
	StreamSlots selectStreams(bool writeFastq, bool isStrictPairedWrite, bool mergedOutput,
		int Pair1, int Cstream1, int Pair2, int Cstream2) const;
	WriteBatch* writeBatch(int curThread);
	void publishBatchLocked(WriteBatch& batch);
	void publishDemultiBatch(WriteBatch& batch);
	void flushAllWriteBatchesLocked();
	bool reserveDemultiBases(bool merged, uint64_t bases);
	std::vector<std::unique_ptr<WriteBatch>> writeBatches_;
	std::atomic<bool> batchWrites_{false};
	// Ordered output with -maxReadsPerOutput or -XfirstReadsWritten: records are staged with their
	// text and written (counted, limited, rotated) as their block is published, in input order.
	struct StagedRecord;
	void replayStagedRecords(WriteBatch& batch);
	std::atomic<bool> replayWrites_{false};
	std::mutex replayMTX;                     // before outputEpochMTX and sqfqostrMTX
	std::unique_ptr<WriteBatch> replayOut_;   // under replayMTX: a block's text per stream
	// Ordered output (beginWriteBatching with an order window), under orderMTX.
	struct OrderedBlock;
	OrderedBlock& orderedBlockLocked(uint64_t block);
	// Publishes the blocks that are ready, in order (all: every block held, at the end).
	void publishOrderedBlocks(std::unique_lock<std::mutex>& lock, bool all);
	void publishOrderedBlock(OrderedBlock& block);
	// A withheld mid-quality set the dereplicator does not count: held for its block's
	// output; false if it is to be written now (no ordered block).
	bool holdWithheld(uint64_t block, uint32_t index, const shared_ptr<DNA>& d1, const shared_ptr<DNA>& d2,
		const Dereplicate::DerepTierRoute& route);
	void derepBlockApplied(uint64_t block);
	std::mutex orderMTX;
	std::condition_variable orderCV;
	std::map<uint64_t, std::unique_ptr<OrderedBlock>> orderedBlocks_;
	std::vector<std::vector<std::unique_ptr<WriteBatch>>> freeBatches_; // per worker: its published buffers
	std::atomic<uint64_t> nextOutputBlock_{0};
	std::atomic<size_t> orderWindow_{0};
	uint64_t nextPublish_ = 0;     // the next block to publish
	bool orderDerep_ = false;      // a block waits for the dereplicator to apply it
	bool orderPublishing_ = false; // a thread is publishing (others leave their blocks to it)
	bool orderFailed_ = false;
	void setwriteMode(std::ios_base::openmode wm) { wrMode = wm; }
	inline void setFilePos(ofstream& str, streamoff& pos) {
		str.seekp(0, ios_base::end);	pos = str.tellp();
	}
	//wh: 0=fastq; 1=fna; 2=qual_
	inline void openOFstream(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool, int);
	inline void openOFstreamFQ(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool = false);
	inline void openOFstreamFQpair(const string opOF, const string opOF2, std::ios_base::openmode wrMode, int p1, string errMsg, bool = false);
	inline void openOFstreamFQ_mrg(const string opOF, std::ios_base::openmode wrMode, int p1, string errMsg, bool = false);
	inline void openOFstreamFNA(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool = false);
	inline void openOFstreamQL(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool = false);
	void openNoBCoutstrean(const string);


	//	void resetOutStreams();
	void delAllDNAvectors();
	void incrementOutputFile();
	// internal variant of closeOutStreams that assumes sqfqostrMTX already held
	// rotating: an output file switch (-maxReadsPerOutput), which keeps the no-barcode files open.
	void closeOutStreams_locked(bool wr, bool rotating = false);
	// internal variant of closeOutFilesDemulti that assumes dmltMTX/sqfqostrMTX already held
	void closeOutFilesDemulti_locked();

	//contains min seq pars & Barcodes etc.
	Filters* MFil;
	//UClinks* optim;
	//for threaded statisitics counting
	//vector<Filters*> subFilter;
	vector<std::unique_ptr<Filters>> subFilter;


	uint totalFileStrms;
	//p1: 0:green, 1:yellow
	//p2: 0:pair1, 1:pair2, 2:single1, 3:single2
    vector<vector<std::shared_ptr<ostr>>> sFile, qFile, fqFile;
	//0:green, 1:yellow..
	vector<std::shared_ptr<dualOfBufStream>> sPairFile, qPairFile, fqPairFile;
	mutex sqfqostrMTX;
	vector<std::shared_ptr<ostr>> fqNoBCFile;
	//mutex nobcostrmMTX;
    vector<std::shared_ptr<ostr>> of_merged_fq;//1D vec, since no read pairs


	//vector<string> locBufGreen;

	vector<string> IDs;
	//controls how memory DNA is written to out file
	int suppressOutWrite;//0=all normal, 1=skip mainfile, 2=skip addfile, 3=skip both
	bool write2File;

	int fastQver; //33, 62 or 59
	int fastQoutVer; //33, 62 or 59
	//write out quality file, is the input paired End sequenced
	// Atomic: output rotation re-applies these (unchanged) while workers read them.
	std::atomic<bool> BWriteQual, BWriteFastQ;
	bool b_multiOutStream;
	int pairedSeq; //1=single, 2=PE, 3=PE+MID
	bool b_changeFQheadVer; // linked to fastQheadVer
	bool b_checkedHeaderChange;
	bool b_oneLinerFasta; // write one line per sequence?
	std::atomic<bool> b_writeGreenQual;
	std::atomic<bool> b_writeYellowQual;

	//asynchronous threads
	//std::vector<std::future<ulong>> threads; 
	int Nthrds;  int thrdsCnt; bool thrdsActive;
	//controls output file size
	uint64_t maxReadsPerOFile;
	std::atomic<uint64_t> ReadsWritten;
	std::atomic<uint64_t> TotalReadsWritten;
	uint demultiBPperSR;
	std::atomic<uint64_t> BPwrittenInSR;
	std::atomic<uint64_t> BPwrittenInSRmerg;
	uint64_t maxRdsOut;
	bool stopAll;//red button, just stop all
	string leadingOutf;
	OptContainer* locCmdArgs;
	shared_ptr<Dereplicate> dereplicator;
	std::atomic<uint64_t> cntDerep;
	mutex drpMTX;

	//abstraction to real file type
	//0,1,2,3 refers to pairs (0,1) & singletons (2,3)
	//0=high qual_, 1=mid qual_
	std::ios_base::openmode wrMode;
	bool doTIO;



	//demultiplexing into files
	bool bDoDemultiplexIntoFiles;
	bool b_doDereplicate;

    //demultiplex files into these (use smart pointers for ownership)
	vector<vector<std::shared_ptr<ofbufstream>>> demultiSinglFiles;
	//vector<vector<string>> demultiSinglFilesF;
	vector<std::shared_ptr<ofbufstream>> demultiMergeFiles;
	bool onlyCompletePairsDemulti; //demultiplex even if other pair_ is not passing?
	mutex dmltMTX;

	//merge the filtered output?
	bool b_merge_pairs_ = false;
	bool b_merge_pairs_filter_ = false;
	bool b_merge_pairs_demulti_ = false;
    vector<std::shared_ptr<ReadMerger>> mergers;
	// serialize file increment (rotate) operations
	mutex fileIncMTX;
	std::shared_mutex outputEpochMTX;

	Benchmark* _benchmark;
};
