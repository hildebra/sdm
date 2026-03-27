#pragma once

//#include "InputStream.h"
#include "containers.h"
#include "Statistics.h"

class Filters;

inline std::vector<std::string> get_stats_invariant_warnings(const collectstats& stats, const std::string& label) {
	std::vector<std::string> warnings;
	const long long total = (long long)stats.total;
	const long long success = (long long)stats.totalSuccess;
	const long long rejected = (long long)stats.totalRejected;
	const long long mid = (long long)stats.totalMid;
	if (success + rejected != total) {
		warnings.push_back(label + ": totalSuccess + totalRejected != total");
	}
	if (mid > success) {
		warnings.push_back(label + ": totalMid > totalSuccess");
	}
	return warnings;
}


//class Filters does the main demultiplexing of raw DNA/QUAL data
class Filters : public std::enable_shared_from_this<Filters> {
public:
	Filters(OptContainer*);
	Filters(Filters* of, int, bool = false, size_t threads = 1);
	//Filters(Filters* of, int, bool = false, size_t threads=1);
	~Filters();
	//was previously in separateByFile()
	//void ini_filestruct(OptContainer& cmdArgs);

	UClinks* ini_SeedsReadsDerep(UClinks* ucl, shared_ptr<ReadSubset>& RDSset,
		shared_ptr<Dereplicate>& Dere);// , ReadMerger* merg);
	//checks some basic parameters of DNA quality in UClinks seed finding..
	void miniCheckDNA(shared_ptr<DNA> d, shared_ptr<DNA> d2);


	Filters* newFilterPerBCgroup(const vector<int>);

	//pair_:-1: no Pair-sequence_, 0,1=pair_ 1/2 (assumes MID BC)
	//doSeeding: extract longes Seed //false, -1, -2
	//bool check(shared_ptr<DNA> in, bool doSeeding, int pair, int &tagIdx);// , bool checkBC = true);
  bool checkYellowAndGreen(const shared_ptr<DNA>& d, int pairPre, int& tagIdx, bool doSeeding);
	//vector<bool> check_pairs(shared_ptr<DNA>p1, shared_ptr<DNA>p2, shared_ptr<DNA>mid, vector<bool>, bool changePHead);
	void setSeqLength(float minL, int maxL);
	void setMaxAmb(int x) { MaxAmb = x; };
	void setAvgMinQual(float x) { min_q = x; };
	bool readMap();
	void setPrimerErrs(int x) { PrimerErrs = x; }
	void setTagErrs(int x) { barcodeErrors_ = x; }
	void removePrimer(bool x) { BcutPrimer = x; }
	void removeTag(bool x) { BcutTag = x; }
	void setMaxHomo(int x) { maxHomonucleotide = x; }
	void setTrimHomo(int x) { trimHomonucleotide = x; }
	void checkDoubleBarcode();
	void checDoubleSampleID();
	void checkDoubleSampleIDHead();

	//complete: filter whole sequence if any window below threshhold
	void setFloatingQWin(int width, float aveQ) { FQWwidth = width; FQWthr = aveQ; };
	//partial: cut end of sequence_ that is below the window threshold
	void setFloatingEWin(int width, float aveQ) { EWwidth = width; EWthr = aveQ; };
	bool setcmdArgsFiles();
   bool remove_adapter(const shared_ptr<DNA>&);
	vector<string> getFastaFiles() { return FastaF; }
	vector<string> getQualFiles() { return QualF; }
	vector<string> getFastqFiles() { return FastqF; }
	vector<string> getMIDfqFiles() { return MIDfqF; }
	void restartFileSet(bool b) { restartSet = b; }
	void setBCfixed(bool b, bool fwd) {
		if (fwd) { BCdFWDREV[fwd].b_BCdirFix = b; }
		else { BCdFWDREV[fwd].b_BCdirFix = b; }
	}
	bool eval_reversingBC(bool);
	bool haveToRestartSet() { if (restartSet) { restartSet = false; return true; }return false; }

	bool doOptimalClusterSeq() { return b_optiClusterSeq; }
	bool doSubselReads() { return b_subselectionReads; }
	void statAddDerepBadSeq(int BC) { //seq did not pass qual_ filter, but could be dereplicated
		collectStatistics[0]->BarcodeDetected[BC - BCoffset]++;
		//collectStatistics[0]->BarcodeDetectedFail[BC - BCoffset]--;
		collectStatistics[0]->DerepAddBadSeq++;
	}
	void countBCdetected(int BC, int Pair, bool MidQ);

	void allResize(unsigned int x);
	void addPrimerL(string, int);
	void addPrimerR(string, int);
	void BarcodePreStats(void);
	void sanityCheckFilesSRs(void);
	void removeSRs(void);
	void resetStats();
	//idxG needs to be BCoffset free, BC from shared_ptr<DNA> needs to have BCoffset added
  void failedStats2(const shared_ptr<DNA>& d, int);

   void BCintoHead(int idx, const shared_ptr<DNA>& d, std::string_view presentBC, const int, bool, bool = false);
	void setBCdna(int idx, shared_ptr<DNA> d) { d->setBCnumber(idx, BCoffset); }
	void SampleIntoHead(const int idx, shared_ptr<DNA> d, const size_t pos);
	void setMultiDNA(OutputStreamer* m) { lMD = m; }
	//stats... probably mutexed functions
	bool doReversePrimers() { return bPrimerR; }
	//routine checks, and reverses/swaps DNA objects
	bool swapReverseDNApairs(vector< shared_ptr<DNA>>&);
	//for single PB reads: check if they are reversed?
	bool isReversedAmplicon(shared_ptr<DNA>);
 void preFilterSeqStat(const shared_ptr<DNA>& d, int pair);
	//    void preFilterSeqStatMT(shared_ptr<DNA> d, data_MT *data, int pair_);
	inline void updateMaxSeqL(int x);
	bool betterSeed(shared_ptr<DNAunique>, shared_ptr<DNAunique>, float, int, bool);
	bool secondaryOutput() { return bAdditionalOutput; }

	void setGoldAxe(bool b, int a, int i) {
		doGoldAxe = b; GoldAxeMinAmpli = i, GoldAxeMaxAmpli = a;
		if (GoldAxeMinAmpli <= 0) { GoldAxeMinAmpli = -1; }// means not to filter at all
		if (GoldAxeMaxAmpli <= 0) { GoldAxeMaxAmpli = -1; }
	}
	bool isGoldAxe() { return doGoldAxe; }// reads are GoldenAxe PacBio?
	vector<shared_ptr<DNA>>  GoldenAxe(vector< shared_ptr<DNA>>& tdn); //GoldenAxe deconcat
	inline bool checkSwitchedRdPairs() { return b2ndRDBcPrimCk; }
	inline bool checkRevRd() { return bRevRdCk; }

	bool synRdPairs() { return bChkRdPrs; }
	int writtenReads() { return ReadsWritten; }
	int maxReadsOutput() { return maxReadsPerOFile; }
	uint getDemultiBPperSR() { return demultiBPperSR; }
	void setWrittenReads(int x) { ReadsWritten = x; }
	int getFileIncrementor() { return OFileIncre; }
	int getReadsWritten() { return ReadsWritten; }
	void incrementFileIncrementor() { OFileIncre++; ReadsWritten = 0; }//
	void setBCoffset(int x) {
		BCoffset = x;
	}
	inline int getBCoffset() { return BCoffset; }
	//if no qual_ file present, than deactivate qual_ filter
	void deactivateQualFilter() { b_doQualFilter = false; }
	//output file
	int getuserReqFastqOutVer(void) { return userReqFastqOutVer; }
	//input file
	int getuserReqFastqVer(void) { return userReqFastqVer; }
	int isPaired() { return  pairedSeq; }
	int& setPaired() { return  pairedSeq; }
	int FQheadV() { return PEheaderVerWr; }
	inline bool consistentPairs() { return bCompletePairs; }
	bool doDemultiplex() { return bDoMultiplexing; }
	bool doDereplicate() { return bDoDereplicate; }
	bool doFilterAtAll() { return b_doFilter; }
	bool doSeedExtension() { return bDoSeedExtension; }
	uint getcut5PR2() { return cut5PR2; }
	uint getcut5PR1() { return cut5PR1; }


	//*************************
	//DNA statistic collection
	void prepStats();
	void addMergeStats(OutputStreamer* out);
	//void revConstellationCnts(int x) { revConstellationN += x; }//number of read pairs, where pair1/2 are changed (mo)
    void addDNAtoCStats(const shared_ptr<DNA>& d, int);
	void sTotalPlus(int pair) {
		//csMTX[pair]->lock();
		collectStatistics[pair]->total++; //collectStatistics[pair_].totalRejected++;
		//csMTX[pair]->unlock();
	}
	void addStats(Filters* fil, vector<int>& idx);

  void DNAstatLQ(const shared_ptr<DNA>& d, int pair, bool Additional) {
		if (Additional) {
			statAddition[pair]->addPostFilt(d);//PostFilt.addDNAStats(d);
		}
		else {
			collectStatistics[pair]->addPostFilt(d); //->PostFilt.addDNAStats(d);
		}
	}
	int getLocalRdsRead() {
		if (pairedSeq > 1) {
			return collectStatistics[0]->total + collectStatistics[1]->total;
		}
		return collectStatistics[0]->total;
	}
	void printStats(ostream&, string, string, bool);
	void printGC(ostream&, int);
	string shortStats(const string&);
	void SmplSpecStats(ostream&);
	void printHisto(ostream&, int which, int set = 1);//which: 1=qual_ //set:0 only filter, 1 all available
	void printLenVsQual(ostream& give);
	bool combineSamples() { return bDoCombiSamples; }

	//handles setting up file paths, file order, file types.. (input only)
	void FileEssentials(filesStr& files, OptContainer* cmdArgs);// 
	//return a vector that says entry x (from invec) corresponds to group y
	vector<int> combiSmplConvergeVec(const vector<string>&);
	//public version of BC finder..
		//-1= no HIT; -5=reverse hits
   int detectCutBC(const shared_ptr<DNA>& d, bool isPair1);//returns id_, important for cutPrimer()
	//int detectCutBC(shared_ptr<DNA> d, string&, int&,bool);//returns id_, important for cutPrimer()
 int findTag(const shared_ptr<DNA>& d, string&, int&, bool,
		int revChecks, bool cutBC, bool endCheck);//returns id_, important for cutPrimer()
	//2nd BC on same DNA sequence (from the 3' end)
	//int findTag2(shared_ptr<DNA> d, string&, int&, bool,int revChecks);
	inline bool doubleBarcodes() { return bDoBarcode2; }
	inline bool doBarcodes() { return bDoBarcode; }

   void dblBCeval(int& tagIdx, int& tagIdx2, string& presentBC, shared_ptr<DNA> tdn, shared_ptr<DNA> tdn2);
	vector<int> getDrerepSampleSpecifity() { return derepMinNum; }
	bool findPrimer(shared_ptr<DNA> d, int primerID, bool, int);


	// Multithreading
	/*void setThreads(size_t threads) {
		this->threads = threads;
		this->statistics_.resize(this->threads);
	}*/

	//public vars *************************************************
	//check for heterogenity primers (can be useful for chimera estimation)
	bool doHetPrimerExplicit;
	vector<string> PrimerL; vector<string> PrimerR;
	vector<string> PrimerL_RC; vector<string> PrimerR_RC;
	vector<int> PrimerIdx; //one entry per barcode / links to PrimerL 
	vector<int> PrimerIdxRev; //one entry per barcode / links to PrimerR 
	vector<string> Barcode, revBarcode, Barcode2, revBarcode2;
	vector<string> SampleID, SampleID_Combi, HeadSmplID;
	vector<vector<string>> hetPrimer;
	vector<string> SequencingRun;
	map<string, vector<int>> SequencingRun2id;

	// statistics for multithreaded variant
	//std::vector<StatisticsMultithreaded> statistics_;
   // std::size_t threads = 0;

	// Central statistics

	// mostly collect statistics of filter (green filtered)
	vector<shared_ptr<collectstats>> collectStatistics; // top quality statistics
	mutex csMTX1, csMTX2;
	// stats for additional reads to be output (yellow)
	vector<shared_ptr<collectstats>> statAddition; // mid quality statistics
	shared_ptr<GAstats> GAstatistics;//GoldenAxe statistics
	shared_ptr<MEstats> mergeStats;





	//combiner of samples base_map to collect the group number
	unordered_map<string, int> combiMapCollectGrp;
	int getXreadsWr() { return firstXreadsW; }
	int getXreadsRd() { return firstXreadsR; }
	int getLocalAcceptRds() { //just plain number of successes..
		if (pairedSeq > 1) {
			return collectStatistics[0]->totalSuccess + collectStatistics[1]->totalSuccess;
		}
		return collectStatistics[0]->totalSuccess;
	}
	void singReadBC2() {
		if (Barcode2.size() > 0 && doubleBarcodes() && isPaired() == 1) {
			bDoBarcode2Rd1 = true;
		}
	}

	//    void addToStatistics(shared_ptr<DNA> d, Statistics &statistics);
	//    void addDNAtoCStatsMT(shared_ptr<DNA> d, int pair, int thread_id);

	 //   void preFilterSeqStatMT(shared_ptr<DNA> d, int pair, uint thread);

		//void addStatsMT(Filters* fil, vector<int> &idx);
	int currentBCnumber() { return curBCnumber; }//only used in "one sample per file" cases

	//quick check if a rev Primer seq matches correct position -> reverse this seq
	bool checkIfRevPrimerHits(shared_ptr<DNA> d, int primerID, int pair = 0, bool = false);
	bool checkIfPrimerHits(shared_ptr<DNA> d, int primerID, int pair = 0);

	bool passedReads(int n);



protected:
	bool check_lengthXtra(shared_ptr<DNA> d, int hindrance = 0, int leng = -1) {
		if (min_l > 0) {
			if (leng == -1) {
				leng = d->length();
			}
			if (leng - hindrance < min_l) {
				d->QualCtrl.minL = true;
				if (leng - hindrance >= alt_min_l) {
					d->setYellowQual(true);
					return false;
				}
				//statAddition.minL++;
				return true;
			}
		}
		if (max_l > 0 && leng - hindrance > (int)max_l) {
			d->QualCtrl.maxL = true; //sMaxLength(pair_);
			return true;
		}

		return false;
	}
	bool check_length(int leng, int hindrance = 0) {
		if (min_l == 0) { return false; }
		return leng - hindrance < min_l;
	}
	bool cutPrimer(shared_ptr<DNA> d, int primerID, bool&, int, bool = false);
	bool cutPrimerRev(shared_ptr<DNA> d, int primerID, bool&, bool = false);



	inline void scanBC(shared_ptr<DNA> d, int& start, int& stop, int& idx, int c_err, int scanRegion,
		string& presentBC, bool fwdStrand, bool revBC = false, bool endScan = false);
	//reverse BC scan on end of read
/*	inline void scanBC_rev(shared_ptr<DNA> d, int& start, int& stop, int& idx, int c_err, int scanRegion,
		string& presentBC, bool fwdStrand);*/
		//just scan the back of read with normal BCs
	//	inline void scanBC_back(shared_ptr<DNA> d, int& start, int& stop, int& idx, int c_err, int scanRegion,
	//		string& presentBC, bool useBC1, bool revBC);

	void extractMap(int k, int cnt, int tbcnt, string& segments, bool);
	void fakeEssentials(bool all);
	void noMapMode();
	void reverseTS_all_BC();
	void reverseTS_all_BC2();
	void debugVerifyStats(const char* context) const;

	//decides if BC is read from fasta header or looked for in MID seq/ DNA seq
	void decideHeadBC();



	vector<string> FastaF, QualF, FastqF, MIDfqF;
	vector<int> derepMinNum;
	OutputStreamer* lMD;

	//technical adapter removal
	string tAdapter;
	unsigned int tAdapterLength;
	//do adapter removal? Do Barcode checking?
	bool removeAdapter, bDoMultiplexing;
	//which kind of barcoding?
	bool bDoBarcode, bDoBarcode2, bDoBarcode2Rd1, bDoHeadSmplID;
	bool bBarcodeSameSize;

	//related to "one sample per file"
	bool bOneFileSample;
	int curBCnumber, BCoffset;

	//do additional 2nd output file using different filter options
	bool bAdditionalOutput;
	//check if reverse primer + rev BC are present (on 2nd read)
	bool b2ndRDBcPrimCk;
	//check if reads have been reversed_
	bool bRevRdCk;
	//check if read pairs are correctly synced
	bool bChkRdPrs;
	//specialized function for LotuS, which doesn't need all the huge output files..
	//BCs are in mid file
	//bool bHasMidSeq;

	//reverse all BCs and save for later use
	void reverse_all_BC();

	//filter related
	int min_l, alt_min_l;
	float min_l_p, alt_min_l_p;
	int maxReadLength;//stats
	bool norm2fiveNTs; //change IUPAC code to 5 bases (ACTGN)
	uint max_l;
	float min_q, alt_min_q;
	bool BcutPrimer, alt_BcutPrimer, bPrimerR;//cut Primers from seq?
	bool bRequireRevPrim, alt_bRequireRevPrim; // reject seq if reverse primer not found
	bool BextensivePrimerChecks;
	bool bRequireFwdPrim, alt_bRequireFwdPrim;
	bool BcutTag;//cut Tag from seq?
	bool bCompletePairs;//if paired seq, only accept complete pairs
	bool bShortAmplicons;//checks for reverse primer on 1st read
	//minBCLength1_ is Barcode length
	unsigned int minBCLength1_, minBCLength2_, maxBCLength1_, maxBCLength2_, minPrimerLength_, maxHomonucleotide, trimHomonucleotide;
	uint cut5PR1, cut5PR2;
	int PrimerErrs, alt_PrimerErrs, barcodeErrors_, MaxAmb, alt_MaxAmb;//allowed max errs per Primer, Tag; max Ambigous Chars(not ACGT)
	int FQWwidth, EWwidth; //Floating window width for avg quality
	int RevPrimSeedL; // seed length of primer that will be searched for
	bool b_BinFilBothPairs;
	float BinFilErr, BinFilP; //binomal filter parameters
	float FQWthr, EWthr, alt_FQWthr, alt_EWthr; //Floating window avg quality under which seq is kicked
	int PEheaderVerWr;//correct PE header format (0/1/2) this is to accomodate the illumina miSeq paired end annotations 2="@XXX 1:0:4" insteand of 1="@XXX/1". 0=don't change or no PE seq.
	int TrimStartNTs;//remove start NT. -1 indicates auto check for GC infrequencies
	int TruncSeq;//remove trailing NT's after this seq length (length after removal of adapters, primers, Barcodes)
	string iniSpacer; // spacer in fasta file name after barcoding
	int userReqFastqVer;//either 1 (33), 2(59) or 3 (62)
	int userReqFastqOutVer;
	double maxAccumQP, alt_maxAccumQP;

	//paired end sequencing related
	int pairedSeq; //1= single read, 2= PE, 3= PE + 1 file with barcodes
	//int revConstellationN;//number of read pairs, where pair1/2 are changed (mo)


	//flow control bools
	struct BCdecide
	{
		int BChit, BCrevhit;
		bool b_BCdirFix, reversedBCs;
		BCdecide() : BChit(0), BCrevhit(0), b_BCdirFix(false), reversedBCs(false) {}
		void reset() { BChit = 0; BCrevhit = 0; b_BCdirFix = false; reversedBCs = false; }
		void fix() { BChit = 0; BCrevhit = 0; b_BCdirFix = true; reversedBCs = false; }
	};
	vector<BCdecide> BCdFWDREV;
	int firstXreadsW;//just prints the first X reads for experiment (read pairs being counted as 2)
	int firstXreadsR;//just reads the first X reads for experiment (read pairs being counted as 2)
	bool restartSet;//start from beginning, i.e. wrong BC direction
	bool b_optiClusterSeq;//SEED extension
	bool b_subselectionReads;//filter out a specific set of reads
	bool b_doQualFilter;//qulity file provided? Then no qual_ filter
	bool b_doFilter; //option file not provided? just crunch files through, but careful about demultiplexing..

	bool bDoDereplicate;
	bool bDoSeedExtension;
	bool bDoCombiSamples;

	//controls output file size
	int maxReadsPerOFile;
	int OFileIncre;
	atomic_int ReadsWritten;
	uint demultiBPperSR;
	//needed to pass by ref
	BarcodeMap emptyBarcodes;
	//base_map with barcodes.. faster matching (?)
	BarcodeMap barcodes1_, barcodes2_, revBarcodes1_, revBarcodes2_;
	vector<int> barcodeLengths1_;
	vector<int> barcodeLengths2_;

	string illuPEfwd, illuPErev, illuSEuni, illuSEidx;
	bool Bcheck4illuAdapts;

	bool doGoldAxe; // reads are GoldenAxe PacBio?h
	int GoldAxeMinAmpli, GoldAxeMaxAmpli;


	OptContainer* cmdArgs;

	uint passed_interval_reads;


};
