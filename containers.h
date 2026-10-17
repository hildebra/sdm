/* sdm: simple demultiplexer
Copyright (C) 2013  Falk Hildebrand

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef _containers_h
#define _containers_h

#include "InputStream.h"
// #include "Filters.h" // avoid circular include; include in .cpp files needing full Filters type
#include "ReadMerger.h"
// #include "OutputStreamer.h" // avoid circular include; forward declare below

#include <memory>
#include <array>
#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <map>
#include <deque>

#include <shared_mutex>

#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include<cstdio>

//#include "include/robin_map.h"

//#include <math.h>


//fwd declaration
#include <algorithm>
#include <functional>
#include <limits>
#include <stdexcept>
class PartitionWorkPool;
class ThreadPool;
class UClinks;
class Dereplicate;
class Filters;
class OutputStreamer;
class ofbufstream;

//definitions



std::ptrdiff_t len_common_prefix_base(char const a[], char const b[]);

struct ltstr
{
  bool operator()( std::string s1,  std::string s2)
  {
    return strcmp(s1.c_str(), s2.c_str()) < 0;
  }
};

// (HashDNA typedef is defined after DNAuniqSet so the type is available)

template<typename K, typename V>
vector<pair<K, V>> mapToVector(const unordered_map<K, V>& map) {
	return std::vector<std::pair<K, V>>(map.begin(), map.end());
}

//output stream for main DNA object
typedef ofbufstream ostr;

class OptContainer {
public:
	using container_type = std::map<std::string, std::string>;
	using iterator = container_type::iterator;
	using const_iterator = container_type::const_iterator;

	struct RuntimeOptions {
		int iniBlockSize = 120;
		int threads = 1;
		bool threadIO = true;
		std::string ignoreIOErrors = "0";
		std::string pairedRDHDOut = "1";
		std::string onlyPair;
		bool onlyCanonicalNTs = false;
		bool strictPositiveQualsOut = false;

		bool noSearchWithMerge = false;
		size_t derepReserve = 1024 * 512;
		bool disableUsearchSizeFmt = false;
		bool derepPerSR = false;
		std::string minDerepCopies;
		std::string derepFormat;
		bool mergePairsDerep = false;
		int derepSrchLen = 150;

		bool hasDerepReserve = false;
		bool hasDereSizeFmt = false;
		bool hasMinDerepCopies = false;
		bool hasDerepFormat = false;
		bool hasMergePairsDerep = false;
		bool hasDerepSrchLen = false;
	};

	std::string& operator[](const std::string& key) {
		runtimeDirty_ = true;
		return options_[key];
	}
	std::string& operator[](std::string&& key) {
		runtimeDirty_ = true;
		return options_[std::move(key)];
	}

	iterator find(const std::string& key) { return options_.find(key); }
	const_iterator find(const std::string& key) const { return options_.find(key); }

	iterator begin() { return options_.begin(); }
	const_iterator begin() const { return options_.begin(); }
	iterator end() { return options_.end(); }
	const_iterator end() const { return options_.end(); }

	bool contains(const std::string& key) const { return options_.find(key) != options_.end(); }

	void set(const std::string& key, const std::string& value) {
		runtimeDirty_ = true;
		options_[key] = value;
	}

	void setDefault(const std::string& key, const std::string& defaultValue) {
		if (!contains(key)) {
			runtimeDirty_ = true;
			options_[key] = defaultValue;
		}
	}

	std::string get(const std::string& key, const std::string& fallback = "") const {
		const_iterator it = options_.find(key);
		if (it != options_.end()) {
			return it->second;
		}
		return fallback;
	}

	container_type& data() {
		runtimeDirty_ = true;
		return options_;
	}
	const container_type& data() const { return options_; }

	const RuntimeOptions& runtime() const {
		refreshRuntimeOptions();
		return runtime_;
	}

	void invalidateRuntime() { runtimeDirty_ = true; }

private:
	void refreshRuntimeOptions() const {
		if (!runtimeDirty_) {
			return;
		}

		RuntimeOptions opts;
		auto getValue = [&](const std::string& key) -> const std::string* {
			const_iterator it = options_.find(key);
			if (it == options_.end()) {
				return nullptr;
			}
			return &it->second;
		};
		auto parseInt = [](const std::string& value, int fallback) {
			size_t parsedChars = 0;
			try {
				int parsed = std::stoi(value, &parsedChars);
				if (parsedChars == value.size()) {
					return parsed;
				}
			}
			catch (...) {}
			return fallback;
		};
		auto parseSizeT = [](const std::string& value, size_t fallback) {
			size_t parsedChars = 0;
			// stoull wraps "-5" to a huge value; a negative size is never valid here.
			if (value.find('-') != std::string::npos) { return fallback; }
			try {
				size_t parsed = std::stoull(value, &parsedChars);
				if (parsedChars == value.size()) {
					return parsed;
				}
			}
			catch (...) {}
			return fallback;
		};

		if (const std::string* value = getValue("-iniBlockSize")) {
			opts.iniBlockSize = parseInt(*value, opts.iniBlockSize);
		}
		if (opts.iniBlockSize < 1) {
			opts.iniBlockSize = 1;
		}

		if (const std::string* value = getValue("-threads")) {
			opts.threads = parseInt(*value, opts.threads);
		}
		if (opts.threads < 1) {
			opts.threads = 1;
		}

		if (const std::string* value = getValue("-threadIO")) {
			opts.threadIO = (*value != "0");
		}

		if (const std::string* value = getValue("-ignore_IO_errors")) {
			opts.ignoreIOErrors = *value;
		}
		if (const std::string* value = getValue("-pairedRD_HD_out")) {
			opts.pairedRDHDOut = *value;
		}
		if (const std::string* value = getValue("-onlyPair")) {
			opts.onlyPair = *value;
		}
		if (const std::string* value = getValue("-onlyCanonicalNTs")) {
			opts.onlyCanonicalNTs = (*value == "1");
		}
		if (const std::string* value = getValue("-strictPositiveQualsOut")) {
			opts.strictPositiveQualsOut = (*value == "1");
		}

		opts.noSearchWithMerge = getValue("-noSearchWithMerge") != nullptr;

		if (const std::string* value = getValue("-derep_reserve")) {
			opts.hasDerepReserve = true;
			opts.derepReserve = parseSizeT(*value, opts.derepReserve);
		}

		if (const std::string* value = getValue("-dere_size_fmt")) {
			opts.hasDereSizeFmt = true;
			opts.disableUsearchSizeFmt = (*value == "1");
		}

		if (const std::string* value = getValue("-derepPerSR")) {
			opts.derepPerSR = (*value == "1");
		}

		if (const std::string* value = getValue("-min_derep_copies")) {
			opts.hasMinDerepCopies = true;
			opts.minDerepCopies = *value;
		}

		if (const std::string* value = getValue("-derep_format")) {
			opts.hasDerepFormat = true;
			opts.derepFormat = *value;
		}

		if (const std::string* value = getValue("-merge_pairs_derep")) {
			opts.hasMergePairsDerep = true;
			opts.mergePairsDerep = (*value == "1");
		}

		if (const std::string* value = getValue("-derepSrchLen")) {
			opts.hasDerepSrchLen = true;
			opts.derepSrchLen = parseInt(*value, opts.derepSrchLen);
		}

		runtime_ = opts;
		runtimeDirty_ = false;
	}

	container_type options_;
	mutable RuntimeOptions runtime_;
	mutable bool runtimeDirty_ = true;
};

void setGlobalOptContainer(OptContainer* options);
OptContainer* getGlobalOptContainer();
OptContainer& globalOptContainer();

typedef std::unordered_map<std::string, int> ClusterIdx;
typedef robin_hood::unordered_map<std::string, int> BarcodeMap;//links directly to entry number in Barcode vector
//typedef std::base_map<std::string, int, ltstr> ClusterIdx;
//used in UCF file
//typedef std::unordered_map<string, int>::iterator DNAidmapsIT;
typedef std::unordered_map<string, int> DNAidmaps;

#ifdef KHASH 
typedef khset_t<const char*> HashDNA;
typedef khset_t<const char*>::iterator HashDNAIT;
#else

/*size_t DNAHasher2(shared_ptr<DNA> k) {
	return ((hash<string>()(k->getSeqPseudo())) >> 1);
}*/

class DNAHasher3 {
public:
	size_t operator() (const shared_ptr<DNA> k) const {
		return ((hash<string>()(k->getSeqPseudo())) >> 1);
	}
};

class DNAequal {
public:
	bool operator()(const shared_ptr<DNA> val1, const shared_ptr<DNA> val2) const {
		return val1->getSeqPseudo() == val2->getSeqPseudo();
	}
};

//typedef std::unordered_set<shared_ptr<DNAunique>, function<decltype(DNAHasher2)>> HashDNA;

//typedef std::unordered_set<shared_ptr<DNAunique>, DNAHasher3, DNAequal> HashDNA;
//typedef tsl::robin_set<shared_ptr<DNAunique>, DNAHasher3, DNAequal> HashDNA;
//typedef robin_hood::unordered_set<shared_ptr<DNAunique>, DNAHasher3, DNAequal> HashDNA;

//typedef std::unordered_map<string, int>::iterator HashDNAIT;
#endif


void trim(std::string& s);//removes white spaces from string
bool is_digits(const std::string &str);

const std::string SingletonFileDescr = ".singl";
const std::string DEFAULT_BarcodeNameSep = "__";
const std::string DEFAULT_output_qual_offset = "33"; //61 or 33
const std::string DEFAULT_ignore_IO_errors = "0"; //61 or 33
const std::string DEFAULT_pairedRD_HD_out = "1"; //1=write /1 or 1:N:00 etc out, has to be surpressed for some applications (dada2)
const std::string DEFAULT_5PR1cut = "0";
const std::string DEFAULT_5PR2cut = "0";




	//functions
string additionalFileName(const string& in);
string additionalFileName2(const string& in);
inline string getFileNoPath(string & s) {
	// Either separator; the old "pos2 > pos" test also fired for pos2 == npos on POSIX paths.
	size_t pos = s.find_last_of("/\\");
	if (pos == string::npos) {
		return s;
	}
	return s.substr(pos+1);
}
inline string removeFileEnding(string & s) {
	size_t pos = s.find_last_of(".");
	if (pos == string::npos) {
		return s;
	}
	return s.substr(0,pos);
}
//function to insert string into filename, e.g. filt.txt -> filt.ax.txt
string subfile(string x, string y);


inline int FastqVerMod(int x){
	if (x==1){
		return 33;
	} else if (x==2){
		return 64;
	} 
	return 59;
}


//filters a fasta file for certain reads
class ReadSubset{
public:
	ReadSubset(const string,const string);
	~ReadSubset(){}
	bool multiFile() {		if (outFiles.size() > 1) { return true; } return false;	}
	vector<string> getOFiles() { return outFiles; }
	void findMatches(shared_ptr<InputStreamer>, OutputStreamer*,bool mocatFix);
	void setRemainingFilepipe(int j) { RemainderStrPos = j; }
private:
	//-1 deactivates
	int RemainderStrPos;
	unordered_map <string, int> Targets;
	vector<string> newHD, outFiles;
	vector<uint> outFilesIdx;
};


bool DNAuPointerCompare(shared_ptr<DNAunique> l, shared_ptr<DNAunique> r);


#include "CoarseDereplication.h"

class Dereplicate{
public:
	Dereplicate(OptContainer*, Filters* mf);//
	~Dereplicate();
	int getHighestBCoffset() { return (int)barcode_number_to_sample_id_.size(); }
	// Where a record set would be written if the dereplicator does not count it
	// (the mid-quality tier). A set is either counted or written, never both, so
	// pipelines that back-map the mid-quality output cannot count it twice.
	struct DerepTierRoute { int idx1 = -1, cstr1 = 100, idx2 = -1, cstr2 = 100; };
	// A set's output block (the caller's number of the read block it came from, beginBlock; kNoBlock
	// outside blocks) and its position among the block's entries, so that its output can be written
	// in input order however the block is applied.
	using UncountedPublisher = std::function<void(const shared_ptr<DNA>&, const shared_ptr<DNA>&, const DerepTierRoute&,
		uint64_t outputBlock, uint32_t indexInBlock)>;
	void setUncountedPublisher(UncountedPublisher publisher) { publisher_ = std::move(publisher); }
	// Called once every entry of an output block has been applied (so every withheld set of it has been
	// published), from the thread that applied it.
	using BlockAppliedListener = std::function<void(uint64_t outputBlock)>;
	void setBlockAppliedListener(BlockAppliedListener listener) { blockApplied_ = std::move(listener); }
	// Whether reserveBlock numbers blocks (ordered admission).
	bool ordersBlocks() const { return orderedAdmission_; }
	// Applies every queued entry (publishing withheld sets); call when the
	// output streams that the publisher writes to are about to close.
	void flushPending();
	// Returns the R1 pass state. With a route, the caller withholds the set's
	// mid-quality write and this object publishes it if it is not counted;
	// queued tells whether the entry was taken at all.
	bool addDNA(shared_ptr<DNA> dna, shared_ptr<DNA> dna2, Filters* filter = nullptr,
		const DerepTierRoute* route = nullptr, bool* queued = nullptr);
	string writeDereplDNA(Filters* fil, string SRblock);
	void writeLog(string logF, string rep) {
		ofstream logx;
		string logPS = logF.substr(0, logF.length() - 4) + "dereplication.log";
		logx.open(logPS.c_str(), ios_base::out);
		logx << "Dereplication log:\n"<<rep;
		logx.close();
	}
	void setPaired(bool b) { b_pairedInput = b; }
	void BCnamesAdding(Filters*fil);
	void reset();
	bool DerepPerSR() { return b_derepPerSR; }
	bool mergeDereRead() {return b_merge_pairs_derep_;	}
	//void attachMerger(ReadMerger* m) {if (merger != nullptr) { delete merger; } merger = m;}
	void activateMerger(bool collectStats = false) {
		if (merger == nullptr) { merger = DBG_NEW ReadMerger(collectStats); }
		else { merger->setStatsTaken(collectStats); }
	}
	void addMergeStatsTo(ReadMerger& target) { target.addRMstats(merger); }
	void printMergeStats(string R1, string R2) {
		if (merger != nullptr) {
			if (R1 != "") { merger->printMergeHisto(R1); }
			if (R2 != "") { merger->printQualHisto(R2); }
		}
	}

	void finishMap();
	bool doSearchWithMerge() { return searchWithMerg; }
	bool usesPrefixDereplication() const { return prefixMode_; }
    bool usesCoarseDereplication() const { return coarse_ != nullptr; }
    size_t coarsePartitionCount() const { return coarse_ ? shards_.size() : 0; }

	// Ordered admission (several workers). The reader numbers read blocks in
	// input order (reserveBlock); a worker brackets its processing of a block
	// with beginBlock/endBlock. Entries are applied per shard in block order,
	// so admission, representatives and counts equal a single-threaded run no
	// matter how the blocks are scheduled. Shards still apply in parallel.
	static constexpr uint64_t kNoBlock = std::numeric_limits<uint64_t>::max();
	uint64_t reserveBlock();
	// outputBlock: the caller's number of the block for its output (publisher, block-applied listener).
	void beginBlock(uint64_t block, uint64_t outputBlock = kNoBlock);
	void endBlock();


private:
	struct PendingDerepEntry {
		string srchSeq;
		shared_ptr<DNA> dna;
		shared_ptr<DNA> dna2;
		shared_ptr<DNA> dna_merged;
		int mergePos;
		int sample_id;
		bool pass;
		size_t shard;
		// Global coarse index, ordered admission: search done before sealing.
		CoarseDereplication::Probe probe{};
		// Mid-quality output withheld by the caller: written by publish_uncounted
		// if this set is not counted.
		bool withheld = false;
		DerepTierRoute route{};
		uint64_t outputBlock = kNoBlock; // the producer's output block, and the entry's place in it
		uint32_t indexInBlock = 0;
	};

	struct alignas(64) ThreadLocalBatch {
		ThreadLocalBatch() : owner(nullptr) {}
		std::vector<PendingDerepEntry> pending;
		Dereplicate* owner;
		std::atomic<bool> active{false};
		uint64_t block = kNoBlock; // block being processed (ordered admission)
		uint64_t outputBlock = kNoBlock; // its output block (beginBlock)
		// Global coarse index: reads this producer allocated, handed back by the
		// threads that committed or applied them and freed here, in its own arena.
		std::mutex returnedMutex;
		std::vector<std::vector<PendingDerepEntry>> returnedEntries;
		std::vector<std::vector<CoarseDereplication::Update>> returnedUpdates;
	};
	// A sealed block and the producer whose thread allocated its reads.
	struct SealedBlock {
		std::vector<PendingDerepEntry> entries;
		ThreadLocalBatch* owner = nullptr;
		uint64_t outputBlock = kNoBlock;
	};
	// Tells the block-applied listener that these blocks are applied.
	void notify_applied(const std::vector<SealedBlock>& blocks) const;

	// deferring: the entry goes into a shard's deferred (rejected) table, which
	// also accepts non-passing reads whose key is unknown.
	bool apply_pending_entry_locked(const PendingDerepEntry& pending, HashDNA& tracker, bool deferring);
    struct Shard;
    void apply_pending_entry(const PendingDerepEntry& pending, Shard& shard,
        const CoarseDereplication::Probe* probe = nullptr, std::vector<CoarseDereplication::Update>* deferred = nullptr);
    // Non-passing reads never create a dereplicate. Whenever exact keys are the
    // output unit they are deferred instead of dropped when their key is still
    // unknown, so recovery does not depend on which copy of a sequence arrives
    // first. Coarse-cluster output admits them to clusters directly.
    // Variant-retaining modes (-derepStoreQuals/-derepStoreDiffs) store every
    // counted read as an exact variant, so exact-key mode does not defer there
    // (it would hold every failed read); prefix mode keeps its full evidence.
    bool defersNonPassing() const { return !coarse_ || coarse_->outputsExact(); }
    // Every counted read is also stored as an exact variant (-derepStoreQuals/-derepStoreDiffs).
    bool retainsVariants() const { return coarse_ && (coarse_->storesQualities() || coarse_->storesDiffs()); }
    bool exact_key_known(const PendingDerepEntry& pending, const Shard& shard) const;
    void publish_uncounted(const PendingDerepEntry& pending);
    UncountedPublisher publisher_;
    BlockAppliedListener blockApplied_;
	void flush_batch(ThreadLocalBatch& batch);
    // Global coarse index with several workers: candidate searches run
    // concurrently against a snapshot, commits are exclusive and in batch
    // order, and per-target updates run afterwards under striped locks.
    void flush_global_coarse(ThreadLocalBatch& batch);
    // Clears every shard (in parallel when idle workers exist). Needs LifecyclePause.
    void clear_shards();
	void flush_all_batches();
	ThreadLocalBatch& get_or_create_thread_batch();

	//is the exact derep string fullfilled?
	inline bool pass_deprep_conditions(shared_ptr<DNAunique>);

	//vector<shared_ptr<DNAunique>> Dnas;
	//vector<shared_ptr<DNAunique>> DNApair;
	vector<string> barcode_number_to_sample_id_;
    // Per-sample counts of deferred non-passing reads of one key: no read data,
    // only what output-time recovery needs (FASTA output collects no qualities).
    struct SampleTally {
        std::vector<std::pair<int, uint64_t>> counts; // sorted by sample id
        uint64_t total = 0;
        void add(int sample, uint64_t copies = 1) {
            if (sample < 0 || copies == 0) { return; }
            auto it = std::lower_bound(counts.begin(), counts.end(), sample,
                [](const std::pair<int, uint64_t>& entry, int s) { return entry.first < s; });
            if (it == counts.end() || it->first != sample) { it = counts.emplace(it, sample, 0); }
            if (copies > std::numeric_limits<uint64_t>::max() - it->second || copies > std::numeric_limits<uint64_t>::max() - total) {
                throw std::overflow_error("Deferred dereplication count overflow");
            }
            it->second += copies;
            total += copies;
        }
    };
    typedef robin_hood::unordered_node_map<string, SampleTally> RejectedTallies;
    // Exact-key mode: non-passing reads whose key has no passing dereplicate yet,
    // kept as sample counts under a 64-bit hash of the key (no read, no key
    // string: about 21 bytes per key). Nearly every key has one sample; keys
    // seen in several samples move to `multi`. The table is capped (see
    // kDeferredKeyCap), so its memory is bounded independently of the input.
    struct DeferredCounts {
        struct Single { int32_t sample; uint32_t copies; };
        robin_hood::unordered_flat_map<uint64_t, Single> single;
        robin_hood::unordered_map<uint64_t, std::vector<std::pair<int, uint64_t>>> multi;
        size_t size() const { return single.size() + multi.size(); }
        // false when the key is new and the table is full (the read is not kept).
        bool add(uint64_t key, int sample, size_t cap);
        void clear() { single.clear(); multi.clear(); }
    };
    static constexpr size_t kDeferredKeyCap = 4000000; // all shards together, ~85 MB
    // Variant-retaining modes must record a recovered read as an exact variant,
    // so they keep the (already allocated) reads, up to a fixed byte budget.
    struct DeferredEvidence { std::shared_ptr<DNA> dna, dna2, merged; int sample; };
    static constexpr uint64_t kDeferredEvidenceBytes = uint64_t(64) << 20;
    std::atomic<uint64_t> deferredEvidenceBytes_{0};
    static uint64_t deferredKeyHash(std::string_view key) { return robin_hood::hash<std::string_view>{}(key); }
    struct Shard {
        std::mutex mutex;
        // Global coarse mode: shared for concurrent probes, exclusive for commits.
        std::shared_mutex snapshot;
        HashDNA tracker, rejected;
        // Prefix mode: deferred non-passing reads when only their counts are needed.
        RejectedTallies rejectedTallies;
        DeferredCounts deferredCounts;
        robin_hood::unordered_node_map<uint64_t, std::vector<DeferredEvidence>> deferredEvidence;
        std::unique_ptr<CoarseDereplication> coarse;
        struct AdmissionStats {
            uint64_t high = 0, mid = 0, failed = 0, unmatched = 0, deferred = 0, published = 0, mergedKeys = 0;
            uint64_t deferredDropped = 0; // deferred-key table full: read not kept
        } admission;
    };
    struct DerepOutputGroup {
        uint64_t count;
        std::string_view key;
        DNAuniqSet* representative;
        std::vector<DNAuniqSet*> members;
        std::shared_ptr<DNAunique> refined;
        std::vector<std::pair<int, uint64_t>> recovered; // deferred sample counts attached to this group
        std::vector<const DeferredEvidence*> evidence;   // deferred reads recorded as variants at write time
    };
    // With a pool, shards are collected and the prefix-mode key sort runs in parallel.
    std::vector<DerepOutputGroup> outputGroups(std::vector<uint64_t>& recoveredCounts, ThreadPool* pool);
	std::vector<std::unique_ptr<Shard>> shards_;
	//vector<int> counts;
	string outfile;
	bool b_usearch_fmt, b_singleLine;
	bool b_pairedInput;
	// -min_derep_copies: a sequence passes if any rule holds.
	// "N": at least N copies in total across all samples (samples == 0).
	// "X:Y": at least X copies in each of at least Y samples.
	struct DerepCopyRule { int64_t copies; int samples; };
	vector<DerepCopyRule> minCopyRules;
	int64_t minCopiesFastFailThreshold;
	string minCopiesStr;
	//int passedSize;
	int tmpCnt;
	int curBCoffset;
	// Output as fasta or fastq
	bool b_derep_as_fasta_;
	bool b_derepPerSR;
	bool b_wroteMapHD;
	bool b_merge_pairs_derep_;
	ReadMerger* merger;

	int searchLength;
	bool prefixMode_ = false;
	// Deferred non-passing reads keep only a sample tally (FASTA output, no
	// variant/quality retention); otherwise the DNAunique-based table is used.
	bool deferCountsOnly_ = false; // prefix mode only; exact-key mode always uses DeferredCounts
	size_t deferredKeyCapPerShard_ = kDeferredKeyCap;
    // Non-owning first index: configuration is shared; mutable data belongs to each shard.
    CoarseDereplication* coarse_ = nullptr;
    std::unique_ptr<PartitionWorkPool> coarseWorkPool_;
    bool parallelGlobalCoarse_ = false;
    std::array<std::mutex, 256> coarseUpdateStripes_;
    // Ordered admission state.
    bool orderedAdmission_ = false;
    std::unique_ptr<PartitionWorkPool> orderedPool_;
    std::atomic<uint64_t> nextBlock_{0};
    std::mutex sequenceMutex_;   // sealed_, nextSubmit_
    std::mutex submitMutex_;     // serializes enqueues into orderedPool_
    std::map<uint64_t, SealedBlock> sealed_;
    uint64_t nextSubmit_ = 0;
    // Applies whole blocks, in the given order, one partition job per shard and block.
    // Afterwards each block's entries go back to their producer (see return_reads),
    // unless that is `self`, the calling producer.
    void apply_blocks_ordered(std::vector<SealedBlock>& blocks,
        std::unique_lock<std::mutex>& submitLock, ThreadLocalBatch* self);
    // Global coarse index, ordered admission: sealed blocks queue here in block
    // order. The worker that finds no commit running becomes the committer and
    // commits the queue; the others return to their next block.
    // The committer posts each block's representative/count/variant updates to
    // per-target lanes, in commit order. Any thread may apply a ready lane;
    // one executor per lane keeps every target's updates in commit order.
    // Producers help with the lanes (instead of reading on) while the queued
    // blocks or updates exceed their backlog, since both hold reads in memory.
    // The committer returns only once the queue and all lanes are empty.
    std::mutex commitMutex_;     // everything below up to commitError_
    std::condition_variable commitSpace_;
    std::deque<SealedBlock> commitQueue_;
    size_t commitBacklog_ = 0;   // queued blocks before producers help
    bool committing_ = false;
    static constexpr size_t kUpdateLanes = 64;
    struct UpdateGroup {
        std::vector<CoarseDereplication::Update> updates;
        ThreadLocalBatch* owner = nullptr;
    };
    struct UpdateLane {
        std::deque<UpdateGroup> queue;
        bool running = false;
    };
    std::array<UpdateLane, kUpdateLanes> updateLanes_;
    std::deque<size_t> readyLanes_;  // lanes with queued updates and no executor
    size_t pendingUpdates_ = 0;      // queued or being applied
    size_t updateBacklog_ = 0;       // pending updates before producers help
    std::exception_ptr commitError_;
    // Set while the committer waits for the exclusive snapshot lock; probes
    // pause between chunks, since the lock itself prefers readers.
    std::atomic<bool> commitWaiting_{false};
    static constexpr size_t kProbeChunk = 8;
    // Caller holds sequenceMutex_. Returns true if the caller must commit.
    bool queue_global_commits(std::vector<SealedBlock>& blocks);
    void commit_global_blocks();
    void commit_global_block(SealedBlock& block);
    void post_updates(std::vector<CoarseDereplication::Update>& updates, ThreadLocalBatch* owner);
    // Hand reads back to the producer that allocated them (see ThreadLocalBatch).
    void return_reads(ThreadLocalBatch* owner, std::vector<PendingDerepEntry>&& entries);
    void return_reads(ThreadLocalBatch* owner, std::vector<CoarseDereplication::Update>&& updates);
    void free_returned_reads(ThreadLocalBatch& batch);
    // Applies one ready lane's queued updates; false if no lane is ready.
    bool apply_ready_lane(std::unique_lock<std::mutex>& lock);
    void fail_commits_locked(std::exception_ptr error);
    void wait_for_commit_backlog();
    // Lifecycle pause: submits sealed blocks left behind a missing block.
    void drain_sealed_blocks();
    void account_search_keys(const std::vector<PendingDerepEntry>& entries);
    bool writeSubclusters_ = false;


	//output files, that need appending
	string mapF;// = baseOF + ".map";
	string outHQf;// = baseOF + ".1.hq.fq";
	string outHQf_p2;// = baseOF + ".2.hq.fq";
	string outRest;// = outfile + ".rest";

	Filters* mainFilter;

	mutable std::mutex batch_registry_mtx_;
	std::vector<std::shared_ptr<ThreadLocalBatch>> batch_registry_;
	// Lifecycle operations serialize with each other; producers announce activity
	// on their own cache line so unrelated reads do not contend on one counter.
	std::mutex lifecycle_operation_mtx_;
	struct LifecyclePause {
		explicit LifecyclePause(Dereplicate& owner);
		~LifecyclePause();
		Dereplicate& owner;
		std::unique_lock<std::mutex> operationLock;
	};
	struct ActiveAddGuard {
		ActiveAddGuard(Dereplicate& owner, ThreadLocalBatch& batch);
		~ActiveAddGuard();
		Dereplicate& owner;
		ThreadLocalBatch& batch;
	};
	void release_batch_activity(ThreadLocalBatch& batch);
	std::atomic<uint64_t> srch_seq_total_len_{ 0 };
	std::atomic<uint64_t> srch_seq_count_{ 0 };
    std::atomic<uint64_t> skippedEmptyMates_{0};
    std::atomic<uint64_t> skippedEmptyKeys_{0};
	alignas(64) std::atomic<bool> lifecycle_pause_{ false };
	mutable std::mutex lifecycle_wait_mtx_;
	mutable std::condition_variable lifecycle_wait_cv_;
	static constexpr size_t kThreadBatchFlushThreshold = 256;

    size_t tracker_size_total() const;


	bool searchWithMerg;// merge pairs during dereplication (subgroups, representative choice)
	// Key on the merged consensus only where it is also written (-merge_pairs_derep);
	// otherwise the written R1 is the key, so no two output records are identical.
	bool mergedSearchKeys_ = false;

};


class betterSeedStruct {
public:
	betterSeedStruct() : tdn1(nullptr), bestPID(0.f), len(0), mergLen(0),
		mergeCnt(0), notMergeCnt(0), cummMergeLen(0){}
	betterSeedStruct(shared_ptr<DNAunique> t) : 
		tdn1(t), bestPID(t->getTempFloat()), len(t->length()), mergLen(0),
		mergeCnt(0), notMergeCnt(0), cummMergeLen(0){
		if (t->getMerge() != nullptr) { 
			mergeCnt++; mergLen = t->getMerge()->length(); 
			cummMergeLen += mergLen;
		}else { notMergeCnt++; }
	}
	bool whoIsBetter2(shared_ptr<DNAunique>, uint64_t candidateAbundance = 0);


	shared_ptr<DNAunique> tdn1;
	//shared_ptr<DNA> tdn2;
	float bestPID;
	uint len, mergLen;
	uint mergeCnt, notMergeCnt;
	long long cummMergeLen;
    // Zero uses the DNAunique sample total (legacy HQ seeds). Reconstructed
    // variants have multiplicities, but no variant-specific sample assignments.
    uint64_t abundance = 0;
};

class UClinks{
public:
	UClinks(OptContainer* );
	~UClinks();
	void findSeq2UCinstruction(shared_ptr<InputStreamer>,bool, Filters* fil);
	void writeNewSeeds(OutputStreamer*, Filters* fil, bool, bool=false);
	void printStats(ostream&);
	void finishUCfile(Filters* fil, string, bool);
	void finishMAPfile();
	void setupDefSeeds(shared_ptr<InputStreamer> FA, const vector<string>& smpls);
	//to add "high qual_" ref sequences
	void addDefSeeds(shared_ptr<InputStreamer> FA, Filters* fil);
	void pairedSeqsMerged(){ pairsMerge = true; }
	void writeOTUmatrix(string outfile);
	void resetInputUcUp(){ UpUcFnd = false; }
	void set2UC(){ UPARSE8up = false; }
	void activateMerger(bool collectStats = false) {
		if (merger == nullptr) { merger = DBG_NEW ReadMerger(collectStats); }
		else { merger->setStatsTaken(collectStats); }
	}
	void addMergeStatsTo(ReadMerger& target) { target.addRMstats(merger); }
	//void attachMerger(ReadMerger* merg) {if (merger != nullptr) { delete merger; } merger = merg;}
	void setRefMode(){ RefDBmode = true; RefDBotuStart = (int)oriKey.size(); }//from now on only count adds or ref DB seqs
private:
	void addUCdo(string,bool );
	matrixUnit OTUmatSum() {
		matrixUnit mcnt = 0;
		for (uint i = 0; i < OTUmat.size(); i++) { for (size_t j = 0; j < OTUmat[i].size(); j++) { mcnt += OTUmat[i][j]; } } 
		return mcnt;
	}
	void add2OTUmat(const string&, int, matrixUnit);
	void add2OTUmat(shared_ptr<DNAunique>, int, matrixUnit);
	bool uclInOldDNA(const string&, const vector<int>&, float, Filters* fil);
	bool uclInOldDNA_simple(const string&, const vector<int>&, int&);
	bool getMAPPERline(string&, string&, float&, vector<int>&, bool addFromHDstring = false);
    void addSeedCounts(const vector<int>&, shared_ptr<DNAunique>);
    void besterDNA(const vector<int>&, shared_ptr<DNAunique>, Filters*,
        bool addCounts = true, uint64_t candidateAbundance = 0);
    shared_ptr<DNAunique> prepareSeedCandidate(shared_ptr<DNA>, shared_ptr<DNA>, Filters*);
    void findSubclusterSeeds(shared_ptr<InputStreamer>, bool readFQ, Filters*);
    string seedParentKey(string);
	void setOTUnms();
	void resetMarks();//reset the uc file format, in case of cd-hit/vsearch input

	inline void removeSizeStr(string&);
	inline void removeSeqStr(string&);
	inline void removeCentrStr(string&);
	inline void removeSampleID(string&, const string &);
	inline void removeSampleID(string&, const string &, string&);
	void readDerepInfo(const string);
	uint64_t oneDerepLine(shared_ptr<DNAunique>);
    bool nextDerepLine(string&);
    uint64_t parseDerepLine(const string&, shared_ptr<DNAunique>, bool strictCounts);

	inline bool getTMPmapperLine(string&);
	bool readMapperFileLine(string&);
	string clusterKey(string);

	//pair_: important to keep track whether to remove BC etc.: -1 to remove BC (454); 0 not to (MID miSeq)
	int CurSetPair;
    //store not matched DNA and keep track
	// Combined mapping from sequence id to DNA object for unmatched DNA
	unordered_map<string, shared_ptr<DNAunique>> unusedDNA;
    struct SubclusterParent {
        uint64_t expected = 0, observed = 0, lastVariant = 0;
        vector<int> clusters;
        float identity = 0.f;
    };
    unordered_map<string, SubclusterParent> subclusterParents;
    bool seedSubclusters = false, subclusterMappingsLoaded = false;
	//map<int,shared_ptr<DNA>> oldDNA;
	//map<int,shared_ptr<DNA>> oldDNA2;
	//DNAidmaps unusedID;
	//std::list<string> oldestID;
	//uint DNAunusedPos;
	string derepMapFile;

	//search terms:  "otu" "chimera" "chimera"
	//string otu_term, chimera_term, chimera_term_noise;

	ClusterIdx seq2CI;
	//vector<shared_ptr<DNAunique>> bestDNA;
	vector<shared_ptr<betterSeedStruct>> bestDNA;

	//vector<shared_ptr<DNA>> bestDNA2;
	vector<string> oriKey;
	list<string> mapLines;
	//vector<float> bestPID;
	//vector<uint> bestLEN;
	int clusCnt, uclines;
	string SEP;
	ifstream ucf, mapdere;
	bool UCread,pairsMerge,MAPread;
	bool b_derepAvailable;//has sdm been run in demultiplexer mode?
	bool UPARSE8up, UPARSE9up, UPARSE11up, UpUcFnd;
	bool cdhit, repFound;
	bool vsearch;
	bool ucispaf;
	string otuTerm,otuOUTterm;
	bool RefDBmode;
	int RefDBotuStart;
	bool SeedsAreWritten;
	uint64_t seedMergeAttempts_ = 0;
	uint64_t seedMergeSuccesses_ = 0;
	//count matrix related
	vector < vector <matrixUnit>> OTUmat;
	unordered_map<string, int> SmplIDs;
	unordered_set<string> perfectChims;
	bool unregistered_samples;
	bool doChimeraCnt;
	bool OTUnumFixed; // can new OTUs be added, after inital reading of DNA OTU.fna?
	uint64_t totalDerepCnt;

	float qCovThr, perIDmatch;

	bool b_merge_pairs_optiSeed_;
	ReadMerger* merger;
};



#endif
