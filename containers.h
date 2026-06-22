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

#include <shared_mutex>

#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif
#include<cstdio>

//#include "include/robin_map.h"

//#include <math.h>


//fwd declaration
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
	size_t pos = s.find_last_of("/");
	size_t pos2 = s.find_last_of("\\");
	if (pos == string::npos || pos2 > pos) {
		pos = pos2;
	}
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


class Dereplicate{
public:
	Dereplicate(OptContainer*, Filters* mf);//
	~Dereplicate();
	int getHighestBCoffset() { return (int)barcode_number_to_sample_id_.size(); }
	//bool addDNA(shared_ptr<DNA> dna);
	bool addDNA(shared_ptr<DNA> dna, shared_ptr<DNA> dna2);
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
	void activateMerger() { if (merger == nullptr) { merger = DBG_NEW ReadMerger(); } }
	void printMergeStats(string R1, string R2) {
		if (merger != nullptr) {
			if (R1 != "") { merger->printMergeHisto(R1); }
			if (R2 != "") { merger->printQualHisto(R2); }
		}
	}

	void finishMap();
	bool doSearchWithMerge() { return searchWithMerg; }


private:
	struct PendingDerepEntry {
		string srchSeq;
		shared_ptr<DNA> dna;
		shared_ptr<DNA> dna2;
		shared_ptr<DNA> dna_merged;
		int mergePos;
		int sample_id;
		bool pass;
	};

	struct ThreadLocalBatch {
		ThreadLocalBatch() : owner(nullptr), registered(false) {}
		std::vector<PendingDerepEntry> pending;
		Dereplicate* owner;
		bool registered;
	};

	bool apply_pending_entry_locked(const PendingDerepEntry& pending);
	void flush_batch(ThreadLocalBatch& batch);
	void flush_all_batches();
	std::shared_ptr<ThreadLocalBatch> get_or_create_thread_batch();

	//is the exact derep string fullfilled?
	inline bool pass_deprep_conditions(shared_ptr<DNAunique>);

	//vector<shared_ptr<DNAunique>> Dnas;
	//vector<shared_ptr<DNAunique>> DNApair;
	vector<string> barcode_number_to_sample_id_;
	HashDNA tracker_;
	//vector<int> counts;
	string outfile;
	bool b_usearch_fmt, b_singleLine;
	bool b_pairedInput;
	vector<int> minCopies;
	size_t minCopiesSiz;
	int minCopiesFastFailThreshold;
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


	//output files, that need appending
	string mapF;// = baseOF + ".map";
	string outHQf;// = baseOF + ".1.hq.fq";
	string outHQf_p2;// = baseOF + ".2.hq.fq";
	string outRest;// = outfile + ".rest";

	Filters* mainFilter;

	mutable std::shared_mutex drpMTX_;
	mutable std::mutex batch_registry_mtx_;
	std::vector<std::weak_ptr<ThreadLocalBatch>> batch_registry_;
	std::atomic<uint32_t> active_adddna_{ 0 };
	std::atomic<uint64_t> srch_seq_total_len_{ 0 };
	std::atomic<uint64_t> srch_seq_count_{ 0 };
	std::atomic<bool> lifecycle_pause_{ false };
	mutable std::mutex lifecycle_wait_mtx_;
	mutable std::condition_variable lifecycle_wait_cv_;
	static constexpr size_t kThreadBatchFlushThreshold = 256;

	inline size_t tracker_size_total() const {
		return tracker_.size();
	}

	bool searchWithMerg;

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
	bool whoIsBetter2(shared_ptr<DNAunique>);


	shared_ptr<DNAunique> tdn1;
	//shared_ptr<DNA> tdn2;
	float bestPID;
	uint len, mergLen;
	uint mergeCnt, notMergeCnt;
	long long cummMergeLen;
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
	void activateMerger() { if (merger == nullptr) { merger = DBG_NEW ReadMerger(); } }
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
	void besterDNA(const vector<int>& curCLIDpre, shared_ptr<DNAunique> tdn1,  Filters* fil);
	void setOTUnms();
	void resetMarks();//reset the uc file format, in case of cd-hit/vsearch input

	inline void removeSizeStr(string&);
	inline void removeSeqStr(string&);
	inline void removeCentrStr(string&);
	inline void removeSampleID(string&, const string &);
	inline void removeSampleID(string&, const string &, string&);
	void readDerepInfo(const string);
	int oneDerepLine(shared_ptr<DNAunique>);

	inline bool getTMPmapperLine(string&);

	//pair_: important to keep track whether to remove BC etc.: -1 to remove BC (454); 0 not to (MID miSeq)
	int CurSetPair;
    //store not matched DNA and keep track
	// Combined mapping from sequence id to DNA object for unmatched DNA
	unordered_map<string, shared_ptr<DNAunique>> unusedDNA;
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
	//count matrix related
	vector < vector <matrixUnit>> OTUmat;
	unordered_map<string, int> SmplIDs;
	unordered_set<string> perfectChims;
	bool unregistered_samples;
	bool doChimeraCnt;
	bool OTUnumFixed; // can new OTUs be added, after inital reading of DNA OTU.fna?
	long totalDerepCnt;

	float qCovThr, perIDmatch;

	bool b_merge_pairs_optiSeed_;
	ReadMerger* merger;
};



#endif