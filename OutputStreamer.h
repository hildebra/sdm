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


#ifdef _isa1gzip
#include "include/GZipStr.h"
#endif

#ifdef _izlib
#include "include/izlib.h"
#else

#ifdef _gzipread
#include "include/gzstream.h"
#include "include/zstr.h"
#endif
/*typedef struct {
    FILE* fp;char* mode;int is_plain;
//    struct isal_gzip_header* gzip_header;
//    struct inflate_state* state;
//    struct isal_zstream* zstream;
    uint8_t* buf_in;size_t buf_in_size;uint8_t* buf_out;size_t buf_out_size;
} gzFile_t;
int gzeof(gzFile_t* fp) { return 0; }
typedef gzFile_t* gzFile;
*/
#endif

//static mutex output_mtx;
class ofbufstream {//: private std::streambuf, public std::ostream {
public:
    ofbufstream() :file("T"), keeper(0), keeperW(0), modeIO(ios::app), used(0), usedW(0),
		coutW(true), isGZ(false), doMC(false), primary(nullptr), bufS(0) {}
	ofbufstream(size_t bufferS) :file("T"), keeper(0), keeperW(0), modeIO(ios::app), used(0), usedW(0),
		coutW(true), isGZ(false), doMC(false), primary(nullptr), bufS(bufferS) {

    }
    ofbufstream(const string IF, int mif, bool isMC = false, size_t bufferS = 20000);
    ~ofbufstream();
    void finishWrites();
    bool operator! (void);
    void operator<< (const string& X);
    void emptyStream();
    void activate();
    void deactivate();
    // end

private:
    mutex append_mtx_;
    mutex output_mtx;
    bool internalWrite(bool closeThis);
    bool internalWriteBuffer(std::vector<char>&& buf, bool closeThis);
    void write(std::string s, std::string file);
    void writeStream(bool doKickoff = true);

    //functions
    string file;
    std::vector<char> keeper;
    std::vector<char> keeperW;
    int modeIO;
    size_t used, usedW;
    //write to cout? gzip input? do multicore?
    bool coutW, isGZ, doMC;
    std::unique_ptr<std::ostream> primary;


    size_t bufS;

    // Multithreading throw threadpool
	//ThreadPool *pool = nullptr;//currently just used to check if MC or not
	// Use a persistent thread-pool for task submission
	bool use_thread_pool = false;
    // end
    // track outstanding async write tasks
	std::vector<std::future<bool>> writeKickoffs;
	std::mutex writeKickoffs_mtx;

};

typedef ofbufstream ostr;
typedef std::map<std::string, std::string> OptContainer;

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
    bool open(const string IF, int mif, int pair, bool isMC = false, size_t bufferS = 20000);
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

	//void writeAllStoredDNA();

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
		const shared_ptr<DNA>& d2, int Pair2, int Cstream2);
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
	bool haveToRestartSet() { return MFil->haveToRestartSet(); }
	void resetOutFilesAndFilter();//MD->closeOutStreams();
	void setBCfixed(bool b, bool fwd) { MFil->setBCfixed(b, fwd); write2File = b; }
	void setSubfilters(int num);
	void mergeSubFilters();
	void activateWrite2File() { write2File = true; }
	void createWriteThread() { writeThreadStatus = 1; }
	void setOneLinerFastaFmt(bool b) { b_oneLinerFasta = b; }
	//void printStorage() { cerr << "Size of MD DNA P1:" << DNAsP1.size() << " P2: " << DNAsP2.size() << endl; }
	//void revConstellationCnts(int x) { MFil->revConstellationCnts(x); }
	//dereplication of DNA seqs
	void attachDereplicator(shared_ptr<Dereplicate> de);
	//void dereplicateDNA(shared_ptr<DNA>);
	//dereplicate DNA by looking in dereplicator for 100% id seq. is treadsafe
  void dereplicateDNA(const shared_ptr<DNA>&, const shared_ptr<DNA>&);
	//debug function to look closer at nonBC reads, threadsafe
    void writeNonBCReads(const shared_ptr<DNA>& d, const shared_ptr<DNA>& d2);
	void setReadLimit(int x) { maxRdsOut = x; }

	//demultiplex related
 void write2Demulti(const shared_ptr<DNA>&, int, int BCoffset);
	//write DNA to demultiplexed files, is threadsafe
  void write2Demulti(const shared_ptr<DNA>&, const shared_ptr<DNA>&, int BCoffset, int curThread);
	void generateDemultiOutFiles(string, Filters*, std::ios_base::openmode = ios::out, bool = false);

	bool mergeReads() {
		return b_merge_pairs_;
	}

	atomic_size_t merged_counter_ = 0;
	atomic_size_t total_read_preMerge_ = 0;

	uint getDemultiBPperSR() { return demultiBPperSR; }
	void resetDemultiBPperSR() { demultiBPperSR = 0; }
	void setBPwrittenInSR(uint x) { BPwrittenInSR = x; }
	uint getBPwrittenInSR(void) { return BPwrittenInSR; }
	void setBPwrittenInSRmerg(uint x) { BPwrittenInSRmerg = x; }
	uint getBPwrittenInSRmerg(void) { return BPwrittenInSRmerg; }
	int getReadsWritten() { return ReadsWritten; }

	void attachBenchmark(Benchmark* bench) { _benchmark = bench; }
	void activateReadMerger(int sz) {
       if (sz > (int)mergers.size()) { mergers.resize(sz); }
		for (size_t x = 0; x < (size_t)sz; x++) {
			mergers[x].reset();
			mergers[x] = std::make_shared<ReadMerger>();
		}
	}
	//void attachReadMerger(ReadMerger* merg) { if (merger.size() > 0 && merger.back() == nullptr) { merger.back() = merg; } else { merger.push_back(merg); } }
	bool Demulti2Fls() { return bDoDemultiplexIntoFiles; }
	bool doDeriplicate() { return	b_doDereplicate; }
	bool doWriteNonBCrds() { return fqNoBCFile.size() == 2; }
	bool doDerepMrgSrch() { if (dereplicator != nullptr) { return dereplicator->doSearchWithMerge(); } else { return false; } }


private:
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
	//void writeAllStoredDNA2();
	//void writeAllStoredDNA2t();
	void incrementOutputFile();
	// internal variant of closeOutStreams that assumes sqfqostrMTX already held
	void closeOutStreams_locked(bool wr);
	// internal variant of closeOutFilesDemulti that assumes dmltMTX/sqfqostrMTX already held
	void closeOutFilesDemulti_locked();

	//contains min seq pars & Barcodes etc.
	Filters* MFil;
	//UClinks* optim;
	//for threaded statisitics counting
	//vector<Filters*> subFilter;
	vector<std::unique_ptr<Filters>> subFilter;


	//contains DNA sequences (that failed to have matching Q and vice versa
	//TODO rausschmeissen andreplace with ofbufstrea,
//	vector<shared_ptr<DNA>> DNAsP1;
//	vector<shared_ptr<DNA>> DNAsP2,DNAsS1,DNAsS2,DNAsNoHead;
//	vector<shared_ptr<DNA>> DNAsP1_alt,DNAsP2_alt,DNAsS1_alt,DNAsS2_alt;
	// bis hier
	//Replace with ofbufstream
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
	//bool mem_used;
	atomic_int DNAinMem, writeThreadStatus;

	int fastQver; //33, 62 or 59
	int fastQoutVer; //33, 62 or 59
	//write out quality file, is the input paired End sequenced
	bool BWriteQual, BWriteFastQ;
	bool b_multiOutStream;
	int pairedSeq; //1=single, 2=PE, 3=PE+MID
	bool b_changeFQheadVer; // linked to fastQheadVer
	bool b_checkedHeaderChange;
	bool b_oneLinerFasta; // write one line per sequence?
	bool b_writeGreenQual;
	bool b_writeYellowQual;

	//asynchronous threads
	//std::vector<std::future<ulong>> threads; 
	int Nthrds;  int thrdsCnt; bool thrdsActive;
	//controls output file size
	int maxReadsPerOFile;
	atomic_int ReadsWritten;
	uint demultiBPperSR;
	atomic_int BPwrittenInSR;
	atomic_int BPwrittenInSRmerg;
	int maxRdsOut;//not used currently?
	bool stopAll;//red button, just stop all
	string leadingOutf;
	OptContainer* locCmdArgs;
	shared_ptr<Dereplicate> dereplicator;
	atomic_int cntDerep;
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
	mutex mergStatMTX;

	Benchmark* _benchmark;
};
