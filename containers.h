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


//#include <math.h>




//definitions

struct ltstr
{
  bool operator()( std::string s1,  std::string s2)
  {
    return strcmp(s1.c_str(), s2.c_str()) < 0;
  }
};
typedef std::map<std::string, std::string> OptContainer;
typedef std::unordered_map<std::string, int> ClusterIdx;
typedef std::unordered_map<std::string, int> BarcodeList;//links directly to entry number in Barcode vector
//typedef std::map<std::string, int, ltstr> ClusterIdx;
//used in UCF file
typedef std::unordered_map<string, int>::iterator DNAidmapsIT;
typedef std::unordered_map<string, int> DNAidmaps;

#ifdef KHASH 
typedef khset_t<const char*> HashDNA;
typedef khset_t<const char*>::iterator HashDNAIT;
#else
typedef std::unordered_map<string, int> HashDNA;
typedef std::unordered_map<string, int>::iterator HashDNAIT;
#endif


void trim(std::string& s);//removes white spaces from string
bool is_digits(const std::string &str);

const std::string SingletonFileDescr = ".singl";
const std::string DEFAULT_BarcodeNameSep = "__";
const std::string DEFAULT_output_qual_offset = "33"; //61 or 33
const std::string DEFAULT_ignore_IO_errors = "0"; //61 or 33



bool betterPreSeed(DNA* d1, DNA* d2, DNAunique* ref);

	//functions
string additionalFileName(const string& in);
string additionalFileName2(const string& in);
inline string getFileNoPath(string & s) {
	size_t pos = s.find_last_of("/");
	if (pos == string::npos) {
		pos = s.find_last_of("\\");
	}
	if (pos == string::npos) {
		return s;
	}
	return s.substr(pos);
}
inline string removeFileEnding(string & s) {
	size_t pos = s.find_last_of(".");
	if (pos == string::npos) {
		return s;
	}
	return s.substr(0,pos);
}

inline int FastqVerMod(int x){
	if (x==1){
		return 33;
	} else if (x==2){
		return 64;
	} 
	return 59;
}


class MultiDNA;




struct collectstats{
	unsigned int maxL, PrimerFail,AvgQual, HomoNT;
	unsigned int PrimerRevFail; //Number of sequences, where RevPrimer was detected (and removed)
	unsigned int minL,minLqualTrim, TagFail, MaxAmb, QualWin;
	unsigned int Trimmed, AccErrTrimmed, QWinTrimmed, total, totalRejected;
	unsigned int fail_correct_BC, suc_correct_BC,failedDNAread;
	unsigned int adapterRem, RevPrimFound;
	uint DerepAddBadSeq;
	//binomial error model
	unsigned int BinomialErr;
	uint dblTagFail;
	//recovered singletons within pairs
	unsigned int singleton;
	vector<int> BarcodeDetected;
	vector<int> BarcodeDetectedFail;
	collectstats() : maxL(0), PrimerFail(0), AvgQual(0), HomoNT(0),
		PrimerRevFail(0),minL(0),minLqualTrim(0),
		TagFail(0), MaxAmb(0), QualWin(0),
		Trimmed(0),AccErrTrimmed(0), QWinTrimmed(0),
		total(0), totalRejected(0),
		fail_correct_BC(0), suc_correct_BC(0),
		failedDNAread(0), adapterRem(0), RevPrimFound(0), 
		DerepAddBadSeq(0),
		BinomialErr(0),
		dblTagFail(0),
		singleton(0), BarcodeDetected(0), BarcodeDetectedFail(0) {}
	void addStats(collectstats&, vector<int>& idx);
	void reset();
	
};


//reported stats on sequence properties
class ReportStats{
public:
	ReportStats(bool MedianDo);
	void reset();
	void addDNAStats(DNA* d);
	void calcSummaryStats(float remSeqs, unsigned int min_l, float min_q);
	float calc_median(vector<unsigned int>& in, float perc);
	void add_median2histo(vector<unsigned int>& in, vector<unsigned int>& histo);
	void add_median2histo(unsigned int in, vector<unsigned int>& histo);
	void addMeanStats(unsigned int NT, unsigned int Qsum, float AccErr){
		rstat_NTs+=NT;rstat_totReads++;rstat_qualSum +=Qsum;rstat_accumError+=AccErr;
	}
	unsigned int lowest(const vector<uint>& in);
	unsigned int highest(const vector<uint>& in);
	void printStats2(ostream& give, float remSeqs,int pair);
	void printGCstats(ostream& give);
	void addStats(ReportStats*);
	bool bMedianCalcs;
	const vector<unsigned int> &get_rstat_Vmed(int x) {
		if (x == 1) { return rstat_VQmed; } else { return rstat_VSmed; }
	}
	//const vector<unsigned int> &get_rstat_VSmed(){return rstat_VSmed;}
	vector<size_t> getVrange(int which);
protected:
	
	//median
	vector<size_t> medVrange(const vector<uint>);
	unsigned int rstat_totReads,rstat_NTs, rstat_qualSum, rstat_Qmed, rstat_Smed;
	//means, Relative Sample Quality Score (RSQS), Unifying Sample Quality Score (USQS)
	float RSQS, USQS;
	float rstat_accumError;
	vector<long> QperNT, NTcounts;
	float GCcontent() { return float(NTcounts[2] + NTcounts[3]) / float(NTcounts[0] + NTcounts[1] + NTcounts[2] + NTcounts[3]); }

	//bin based median calculation's
	vector<unsigned int> rstat_VQmed, rstat_VSmed;
};

class dualPrimerDistrStats{
public:
	dualPrimerDistrStats(){}
	dualPrimerDistrStats(const vector<string>&, const vector<string>&);
	~dualPrimerDistrStats(){}
	void reset(){}
};

/*class dualHetSpacDistrStats{
public:
	dualHetSpacDistrStats(){}
	dualHetSpacDistrStats(const vector<string>&, const vector<string>&);
	~dualHetSpacDistrStats(){}
	void reset(){}
};
*/

//filters a fasta file for certain reads
class ReadSubset{
public:
	ReadSubset(const string,const string);
	~ReadSubset(){}
	bool multiFile() {		if (outFiles.size() > 1) { return true; } return false;	}
	vector<string> getOFiles() { return outFiles; }
	void findMatches(InputStreamer*, MultiDNA*,bool mocatFix);
	void setRemainingFilepipe(int j) { RemainderStrPos = j; }
private:
	//-1 deactivates
	int RemainderStrPos;
	unordered_map <string, int> Targets;
	vector<string> newHD, outFiles;
	vector<uint> outFilesIdx;
};


//class Filters does the main demultiplexing of raw DNA/QUAL data
class Filters{
public:
	Filters(OptContainer&);
	Filters(Filters* of, int, bool = false);
	~Filters() {
		for (size_t i = 0; i < 2; i++) { delete RepStat[i]; delete RepStatAddition[i]; } 
		for (size_t i = 0; i < demultiSinglFiles.size(); i++) { 
			if (demultiSinglFiles[i][0] != NULL) { (*demultiSinglFiles[i][0]).close(); } if (demultiSinglFiles[i][1] != NULL) { (*demultiSinglFiles[i][1]).close(); }
		}
		delete PreFiltP1; delete PreFiltP2;
		delete dPDS; delete dHDS;
	}
	//pair:-1: no Pair-Seq, 0,1=pair 1/2 (assumes MID BC)
	//doSeeding: extract longes Seed //false, -1, -2
	bool check(DNA* in, bool doSeeding, int pair, int &tagIdx);// , bool checkBC = true);
	bool checkXtra(DNA* in, int pair, int &tagIdx );
	vector<bool> check_pairs(DNA*p1, DNA*p2, DNA*mid, vector<bool>, bool changePHead);
	void setSeqLength(float minL, int maxL);
	void setMaxAmb(int x) { MaxAmb = x; };
	void setAvgMinQual(float x) { min_q = x; };
	bool readMap(OptContainer&);
	void setPrimerErrs(int x) { PrimerErrs = x; }
	void setTagErrs(int x) { TagErrs = x; }
	void removePrimer(bool x) { BcutPrimer = x; }
	void removeTag(bool x) { BcutTag = x; }
	void setMaxHomo(int x) { maxHomonucleotide = x; }
	void checkDoubleBarcode();
	void checDoubleSampleID();
	void checkDoubleSampleIDHead();

	//complete: filter whole sequence if any window below threshhold
	void setFloatingQWin(int width, float aveQ) { FQWwidth = width; FQWthr = aveQ; };
	//partial: cut end of Seq that is below the window threshold
	void setFloatingEWin(int width, float aveQ) { EWwidth = width; EWthr = aveQ; };
	bool setcmdArgsFiles(OptContainer&);
	bool remove_adapter(DNA*);
	vector<string> getFastaFiles() { return FastaF; }
	vector<string> getQualFiles() { return QualF; }
	vector<string> getFastqFiles() { return FastqF; }
	vector<string> getMIDfqFiles() { return MIDfqF; }
	void restartFileSet(bool b) { restartSet = b; }
	void setBCfixed(bool b, bool fwd) {
		if ( fwd ) { BCdFWDREV[fwd].b_BCdirFix = b; } else { BCdFWDREV[fwd].b_BCdirFix = b; }
	}
	bool eval_reversingBC(bool);
	bool haveToRestartSet() { if (restartSet) { restartSet = false; return true; }return false; }

	bool doOptimalClusterSeq() { return b_optiClusterSeq; }
	bool doSubselReads() { return b_subselectionReads; }
	void statAddDerepBadSeq(int BC){ //seq did not pass qual filter, but could be dereplicated
		colStats[0].BarcodeDetected[BC - BCoffset]++;
		//colStats[0].BarcodeDetectedFail[BC - BCoffset]--;
		colStats[0].DerepAddBadSeq++;
	}
	void countBCdetected(int BC, int Pair, bool MidQ) {
		if (!bDoMultiplexing) { return; }
		if (!MidQ) {
			if (Pair < 0) { Pair = 0; }
			colStats[Pair].BarcodeDetected[BC - BCoffset]++;
		}
		else {
			if (Pair <= 0) {
				statAddition.BarcodeDetected[BC - BCoffset]++;
			}
		}


	}

	void allResize(unsigned int x);
	void addPrimerL(string, int);
	void addPrimerR(string, int);
	void BarcodePreStats();
	void resetStats();
	//idxG needs to be BCoffset free, BC from DNA* needs to have BCoffset added
	void failedStats2(DNA* d,int);

	void BCintoHead(int idx, DNA* d, const string, const int, bool, bool = false);
	void setBCdna(int idx, DNA* d){ d->setBCnumber(idx + BCoffset); }
	void SampleIntoHead(const int idx, DNA* d, const size_t pos);
	void setMultiDNA(MultiDNA* m) { lMD = m; }
	//stats... probably mutexed functions
	bool doReversePrimers() { return bPrimerR; }
	void preFilterSeqStat(DNA* d,int pair);
	inline void updateMaxSeqL(int x);
	bool betterSeed(DNA*, DNA*, DNA*, DNA*, float, uint, int,bool);
	bool secondaryOutput(){return bAdditionalOutput;}
	inline bool checkBC2ndRd() { return b2ndRDBcPrimCk; }
	inline bool checkRevRd() { return bRevRdCk; }
	bool synRdPairs() { return bChkRdPrs; }
	int writtenReads(){return ReadsWritten;}
	int maxReadsOutput(){return maxReadsPerOFile;}
	void setWrittenReads(int x){ReadsWritten=x;}
	int getFileIncrementor(){return OFileIncre;}
	void incrementFileIncrementor(){ OFileIncre++; ReadsWritten = 0; }//
	void setBCoffset(int x) { BCoffset = x; }
	inline int getBCoffset() { return BCoffset; }
	//if no qual file present, than deactivate qual filter
	void deactivateQualFilter() { b_doQualFilter = false; }
	//output file
	int getuserReqFastqOutVer(void){ return userReqFastqOutVer; }
	//input file
	int getuserReqFastqVer(void){ return userReqFastqVer; }
	int & isPaired(){ return  pairedSeq; }
	int FQheadV(){ return PEheaderVerWr; }
	inline bool consistentPairs(){ return bCompletePairs; }
	bool doDemultiplex(){ return bDoMultiplexing; }
	bool doDereplicate() { return bDoDereplicate; }
	bool doFilterAtAll() { return b_doFilter; }

	//*************************
	//DNA statistic collection
	void prepStats(){
		float remSeqs = float(colStats[0].total - colStats[0].totalRejected);
		RepStat[0]->calcSummaryStats(remSeqs, min_l, min_q);
		RepStat[1]->calcSummaryStats(remSeqs, min_l, min_q);
		if (bAdditionalOutput){
			remSeqs = float(statAddition.total - statAddition.totalRejected);
			RepStatAddition[0]->calcSummaryStats(remSeqs, min_l, min_q);
			RepStatAddition[1]->calcSummaryStats(remSeqs, min_l, min_q);
		}
		PreFiltP1->calcSummaryStats(1, min_l, min_q);
		PreFiltP2->calcSummaryStats(1, min_l, min_q);
	}
	void revConstellationCnts(int x) { revConstellationN += x; }//number of read pairs, where pair1/2 are changed (mo)
	void sPrimerFail(int pair) { colStats[pair].PrimerFail++; }
	void sAvgQual(int pair) {  colStats[ pair].AvgQual++; } 
	void sQualWin(int pair) { colStats[pair].QualWin++; }
	void sBinomError(int pair,float Err) { colStats[pair].BinomialErr++; }//maybe collect later info on expected error?
	void sMaxAmbig(int pair) { colStats[pair].MaxAmb++; }
	void sHomoNT(int pair) {colStats[pair].HomoNT++; } 
	void sReversePrimerFnd(int pair) { colStats[ pair].RevPrimFound++; } 
	void sRevPrimerFail(int pair) { colStats[ pair].PrimerRevFail++; } 
	void sMinLength(int pair) {  colStats[pair].minL++; } 
	void sMinQTrim(int pair) {colStats[pair].minLqualTrim++; } 
	void sMaxLength(int pair) { colStats[ pair].maxL++; } 
	void sTagFail(int pair) { colStats[pair].TagFail++; } 
	void sTagCorrected(int pair) { colStats[ pair].suc_correct_BC++; }
	void sTagNotCorrected(int pair) { colStats[pair].fail_correct_BC++; }
	void sTotalPlus(int pair) { colStats[pair].total++; colStats[pair].totalRejected++; }
	void sTotalMinus(int pair) { colStats[pair].total--; colStats[pair].totalRejected--; }
	void sTotalPlusXtra(int pair) {
		if (pair > 1) { return; } statAddition.total++; statAddition.totalRejected++;
	}
	void addStats(Filters*, vector<int>& idx);
	void DNAstatLQ(DNA* d, int pair,bool Add) { 
		if (Add) {
			RepStatAddition[pair]->addDNAStats(d);
		} else {
			RepStat[pair]->addDNAStats(d); 
		}
	}

	void sTrimmed(int pair) { colStats[pair].Trimmed++; } 
	void printStats(ostream&, string, string, bool);
	void printGC(ostream&,int);
	string shortStats(const string &);
	void SmplSpecStats(ostream&);
	void printHisto(ostream&, int which, int set = 1);//which: 1=qual //set:0 only filter, 1 all available
	bool combineSamples(){ return bDoCombiSamples; }
	//return a vector that says entry x (from invec) corresponds to group y
	vector<int> combiSmplConvergeVec(const vector<string>&);
//public version of BC finder..
	int cutTag(DNA* d, string&, int&,bool);//returns id, important for cutPrimer()
	inline bool doubleBarcodes() { return bDoBarcode2; }
	inline bool doBarcodes() { return bDoBarcode; }

	void dblBCeval(int& tagIdx, int tagIdx2, string presentBC, DNA* tdn, DNA* tdn2);
	vector<int> getDrerepSampleSpecifity() {return derepMinNum;	}
	bool Demulti2Fls() {return bDoDemultiplexIntoFiles;}
	void write2Demulti(DNA*,int,int idx, int fqOvr);

	//public vars *************************************************
	//check for heterogenity primers (can be useful for chimera estimation)
	bool doHetPrimerExplicit;
	vector<string> PrimerL; vector<string> PrimerR;
	vector<string> PrimerL_RC; vector<string> PrimerR_RC;
	vector<int> PrimerIdx; //one entry per barcode / links to PrimerL 
	vector<int> PrimerIdxRev; //one entry per barcode / links to PrimerR 
	vector<string> Barcode, revBarcode, Barcode2, revBarcode2;
	vector<string> SampleID, SampleID_Combi,HeadSmplID;
	vector<vector<string>> hetPrimer;
	//demultiplex files into these
	vector<vector<ofstream*>> demultiSinglFiles;
	//collect statistics of filter reasons
	vector<collectstats> colStats;
	collectstats statAddition; // stats for additional reads to be output
	//combiner of samples map to collect the group number
	unordered_map<string, int> combiMapCollectGrp;
protected:
	bool check_lengthXtra(DNA* d, int hindrance=0, int leng=-1){
		if (min_l==0 &&alt_min_l ==0 ){return false;}
		if (leng==-1){
			leng = d->length();
		}
		if (leng-hindrance < min_l){
			if (leng-hindrance >= alt_min_l){
				d->setMidQual(true);
				return false;
			}
			statAddition.minL++;
			return true;
		}
		return false;
	}
	bool check_length( int leng, int hindrance=0){
		if (min_l==0){return false;}
		return leng-hindrance < min_l;
	}
	bool cutPrimer(DNA* d,int primerID,bool,int);
	bool cutPrimerRev(DNA* d,int primerID,bool);
	//-1= no HIT; -5=reverse hits
	int cutTag(DNA* d, bool);//returns id, important for cutPrimer()

	inline void scanBC(DNA* d,int& start,int& stop,int& idx,int c_err, int scanRegion,
		string & presentBC, bool fwdStrand);
	inline void scanBC_rev(DNA* d,int& start,int& stop,int& idx,int c_err, int scanRegion,
		string & presentBC, bool fwdStrand);

	void extractMap(int k, int cnt, int tbcnt, string & segments, bool);
	void fakeEssentials(void);
	void noMapMode(OptContainer& cmdArgs);
	void reverseTS_all_BC();	
	void reverseTS_all_BC2();

	void decideHeadBC();
	int currentBCnumber() { return curBCnumber; }//only used in "one sample per file" cases

	void generateDemultiOutFiles(const string);


	vector<string> FastaF, QualF, FastqF, MIDfqF;
	vector<int> derepMinNum;
	MultiDNA* lMD;
	vector<ReportStats*> RepStat;//2 entries (2 pairs)
	vector<ReportStats*> RepStatAddition;
	ReportStats* PreFiltP1; ReportStats* PreFiltP2;
	
	//technical adapter removal
	string tAdapter;
	unsigned int tAdapterLength;
	//do adapter removal? Do Barcode checking?
	bool bDoAdapter, bDoMultiplexing;
	//which kind of barcoding?
	bool bDoBarcode, bDoBarcode2, bDoHeadSmplID;
	bool bBarcodeSameSize;
	
	//related to "one sample per file"
	bool bOneFileSample;
	int curBCnumber,BCoffset;
	
	//do additional 2nd output file using different filter options
	bool bAdditionalOutput;
	//check if reverse primer + rev BC are present (on 2nd read)
	bool b2ndRDBcPrimCk;
	//check if reads have been reversed
	bool bRevRdCk;
	//check if read pairs are correctly synced
	bool bChkRdPrs;
	//specialized function for LotuS, which doesn't need all the huge output files..
	//BCs are in mid file
	//bool bHasMidSeq;

	//filter related
	int min_l, alt_min_l;
	float min_l_p, alt_min_l_p;
	int maxReadLength;
	uint max_l;
	float min_q,alt_min_q;
	bool BcutPrimer,alt_BcutPrimer,bPrimerR;//cut Primers from seq?
	bool bRequireRevPrim,alt_bRequireRevPrim; // reject seq if reverse primer not found
	bool bRequireFwdPrim,alt_bRequireFwdPrim;
	bool BcutTag;//cut Tag from seq?
	bool bCompletePairs;//if paired seq, only accept complete pairs
	bool bShortAmplicons;
	unsigned int MinTagLen, MinTagLen2, MaxTagLen, MaxTagLen2, MinPrimLen, maxHomonucleotide;
	int PrimerErrs,alt_PrimerErrs,TagErrs,MaxAmb, alt_MaxAmb;//allowed max errs per Primer, Tag; max Ambigous Chars(not ACGT)
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
	double maxAccumQP,alt_maxAccumQP;
	
	//paired end sequencing related
	int pairedSeq; //1= single read, 2= PE, 3= PE + 1 file with barcodes
	int revConstellationN;//number of read pairs, where pair1/2 are changed (mo)


	//flow control bools
	struct BCdecide
	{
		int BChit, BCrevhit;
		bool b_BCdirFix, reversedBCs;
		BCdecide(): BChit(0), BCrevhit(0), b_BCdirFix(false), reversedBCs(false){}
		void reset() { BChit = 0; BCrevhit = 0; b_BCdirFix = false; reversedBCs = false; }
		void fix() { BChit = 0; BCrevhit = 0; b_BCdirFix = true; reversedBCs = false; }
	};
	vector<BCdecide> BCdFWDREV;
	
	bool restartSet;//start from beginning, i.e. wrong BC direction
	bool b_optiClusterSeq;//SEED extension
	bool b_subselectionReads;//filter out a specific set of reads
	bool b_doQualFilter;//qulity file provided? Then no qual filter
	bool b_doFilter; //option file not provided? just crunch files through
	bool bDoDereplicate;
	bool bDoCombiSamples;
	//demultiplexing into files
	bool bDoDemultiplexIntoFiles;

	//controls output file size
	int maxReadsPerOFile,ReadsWritten,OFileIncre;
	//needed to pass by ref
	BarcodeList emptyBCs;
	//map with barcodes.. faster matching (?)
	BarcodeList BCList, BCList2;
	vector<int> Barcode_len, Barcode2_len;
	
	//double BC / het spacer collect stats
	dualPrimerDistrStats* dPDS;
	dualPrimerDistrStats* dHDS;
};


bool DNAuPointerCompare(DNAunique* l, DNAunique* r);

class Dereplicate{
public:
	Dereplicate(OptContainer&);
	~Dereplicate() {
		for (size_t i = 0; i < Dnas.size(); i++) { delete Dnas[i]; }
	}
	int getHighestBCoffset() { return (int)BCN2SmplID.size(); }
	//bool addDNA(DNA* d);
	bool addDNA(DNA* d,DNA* d2,bool& added);
	string writeDereplDNA(Filters*);
	void writeLog(string logF, string rep) {
		ofstream logx;
		string logPS = logF.substr(0, logF.length() - 3) + "dereplication.log";
		logx.open(logPS.c_str(), ios_base::out);
		logx << "Dereplication log:\n"<<rep;
		logx.close();
	}
	void setPaired(bool b) { b_pairedInput = b; }
	void BCnamesAdding(Filters* fil);
	void reset();
private:
	//is the exact derep string fullfilled?
	inline bool pass_deprep_conditions(DNAunique*);

	vector<DNAunique*> Dnas;
	//vector<DNAunique*> DNApair;
	vector<string> BCN2SmplID;
	HashDNA Tracker;
	//vector<int> Counts;
	string outfile;
	bool b_usearch_fmt, b_singleLine;
	bool b_pairedInput;
	vector<int> minCopies;
	size_t minCopiesSiz;
	string minCopiesStr;
	int totSize;
	int tmpCnt;
	int curBCoffset;
};


class UClinks{
public:
	UClinks(OptContainer& );
	~UClinks();
	void findSeq2UCinstruction(InputStreamer*,bool, Filters* );
	void writeNewSeeds(MultiDNA*, Filters*,bool, bool=false);
	void printStats(ostream&);
	void finishUCfile(Filters*, string, bool);
	void finishMAPfile();
	void setupDefSeeds(InputStreamer* FA, Filters* fil);
	//to add "high qual" ref sequences
	void addDefSeeds(InputStreamer* FA, Filters* fil);
	void pairedSeqsMerged(Filters* fil){ pairsMerge = true; fil->setFloatingEWin(0, 0.f); }
	void writeOTUmatrix(string, Filters* fil);
	void resetInputUcUp(){ UpUcFnd = false; }
	void set2UC(){ UPARSE8up = false; }
	void setRefMode(){ RefDBmode = true; RefDBotuStart = (int)oriKey.size(); }//from now on only count adds or ref DB seqs
private:
	void addUCdo(string,bool );
	void add2OTUmat(const string&, int, matrixUnit);
	void add2OTUmat(DNAunique*, int, matrixUnit);
	bool uclInOldDNA(const string&, const vector<int>&, float, Filters*);
	bool uclInOldDNA_simple(const string&, const vector<int>&);
	bool getUCFlineInfo(string&, string&, float&, vector<int>&, bool addFromHDstring = false);
	void besterDNA(const vector<int> curCLID, DNAunique* tdn1, DNA* tdn2, Filters*);
	void setOTUnms();

	inline void removeSizeStr(string&);
	inline void removeSampleID(string&, const string &);
	inline void removeSampleID(string&, const string &, string&);
	void readDerepInfo(string);
	void oneDerepLine(DNAunique*);

	//pair: important to keep track whether to remove BC etc.: -1 to remove BC (454); 0 not to (MID miSeq)
	int CurSetPair;
	//store not matched DNA and keep track
	uint maxOldDNAvec;
	vector<string> oldDNAid; 
	vector<DNAunique*> oldDNA;
	vector<DNA*> oldDNA2;
	DNAidmaps unusedID;
	//std::list<string> oldestID;
	uint DNAunusedPos;
	string derepMapFile;

	//search terms:  "otu" "chimera" "chimera"
	//string otu_term, chimera_term, chimera_term_noise;

	ClusterIdx seq2CI;
	vector<DNAunique*> bestDNA;
	vector<DNA*> bestDNA2;
	vector<string> oriKey;
	vector<float> bestPID;
	vector<uint> bestLEN;
	int clusCnt, uclines;
	string SEP;
	ifstream ucf, mapdere;
	bool UCread,pairsMerge,MAPread;
	bool b_derepAvailable;//has sdm been run in demultiplexer mode?
	bool UPARSE8up, UPARSE9up, UpUcFnd;
	bool RefDBmode;
	int RefDBotuStart;
	bool SeedsAreWritten;
	//count matrix related
	vector < vector <matrixUnit>> OTUmat;
	unordered_map<string, int> SmplIDs;
	bool unregistered_samples;
	bool doChimeraCnt;
	bool OTUnumFixed; // can new OTUs be added, after inital reading of DNA OTU.fna?
};


//writes successful demultis and stores unsuccessful matches of fna/qual for later matching
class MultiDNA{
public:
	//wrStatus controls if this appends or overwrites output
	MultiDNA(Filters* filter, OptContainer& cmdArgs, 
		std::ios_base::openmode wrStatus, ReadSubset* = NULL, 
		string fileExt = "",int=-1);
	~MultiDNA();
	//	void threadAnalyzeDNA(DNA*);
	void setFastQWrite(bool x) { BWriteFastQ = x; BWriteQual = !x; }
	void setQualWrite(bool x) { BWriteQual = x; }
	void addNoHeadDNA(DNA* d) { DNAsNoHead.push_back(d); }
	//-1,-1,-2
	void analyzeDNA(DNA* d, int FilterUse, int pair, int &idx);
	void writeAllStoredDNA();
	vector<bool> analyzeDNA(DNA* p1, DNA* p2, DNA* mid, bool changePHead, int = -1);

	void writeAndDel(DNA* d, int Pair = 1);//1=pair1;2=pair2;3=singleton1,4=singl2
	//Function specifically if several output files are required
	void writeSelectiveStream(DNA* d, int Pair, int FS);//1=pair1;2=pair2;3=singl1,4=singl2  ;; FS: different multi FileStreams to be used
	void saveForWrite(DNA* d, int Pair = 1);//1=pair1;2=pair2;3=singleton
	Filters* getFilters(int w = -1) { if (w == -1) { return MFil; } else { return subFilter[w]; } }
	int isPEseq() { return pairedSeq; }
	void closeOutStreams(bool wr = true);
	//ofstream::app, ios_base::out
	void openOutStreams(OptContainer& cmdArgs, int, std::ios_base::openmode, string = "",int=-1);
	void openSeveralOutstreams(OptContainer& cmdArgs, ReadSubset*, std::ios_base::openmode);
	string leadOutFile() { return leadingOutf; }
	//void setfastQver(int x){fastQver = x;}
	//void setfastQoutVer(int x){fastQoutVer = x;}

	bool checkFastqHeadVersion(DNA* d, bool = false);
	//int getFastqMod(){return MFil->FastqModifier();}
	int getFastqVer() { return fastQver; }
	int getfastQoutVer() { return fastQoutVer; }
	bool haveToRestartSet() { return MFil->haveToRestartSet(); }
	void resetOutFilesAndFilter();//MD->closeOutStreams();
	void setBCfixed(bool b,bool fwd) { MFil->setBCfixed(b,fwd); write2File = b; }
	void setSubfilters(int num);
	void mergeSubFilters();
	void activateWrite2File() { write2File = true; }
	void createWriteThread() { writeThreadStatus = 1; }
	void setOneLinerFastaFmt(bool b) { b_oneLinerFasta = b; }
	void printStorage() { cout << "Size of MD DNA P1:" << DNAsP1.size() << " P2: " << DNAsP2.size() << endl; }
	void revConstellationCnts(int x) { MFil->revConstellationCnts(x); }
	//dereplication of DNA seqs
	void attachDereplicator(Dereplicate* de);
	//void depPrep(DNA*);
	void depPrep(DNA*,DNA*);
private:
	void setwriteMode(std::ios_base::openmode wm) {	wrMode = wm;}
	inline void setFilePos(ofstream& str,streamoff& pos){
			str.seekp(0,ios_base::end);	pos = str.tellp();
	}
	//wh: 0=fastq; 1=fna; 2=qual
	inline void openOFstream(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg, bool, int);
	inline void openOFstreamFQ(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg,bool=false );
	inline void openOFstreamFNA(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg,bool=false);
	inline void openOFstreamQL(const string opOF, std::ios_base::openmode wrMode, int p1, int p2, string errMsg,bool=false);

//	void resetOutStreams();
	void delAllDNAvectors();
	void writeAllStoredDNA2();
	void writeAllStoredDNA2t();
	void incrementOutputFile();

	//contains min seq pars & Barcodes etc.
	Filters* MFil;
	//UClinks* optim;
	//for threaded statisitics counting
	vector<Filters*> subFilter;
	//contains DNA sequences (that failed to have matching Q and vice versa
	vector<DNA*> DNAsP1,DNAsP2,DNAsS1,DNAsS2,DNAsNoHead;
	vector<DNA*> DNAsP1_alt,DNAsP2_alt,DNAsS1_alt,DNAsS2_alt;
	vector<string> IDs;
	//controls how memory DNA is written to out file
	int suppressOutWrite;//0=all normal, 1=skip mainfile, 2=skip addfile, 3=skip both
	bool write2File;
	bool mem_used;
	int DNAinMem, writeThreadStatus;
#ifdef _THREADED
	std::thread wrThread;
	std::mutex mutex;
	vector<std::thread> threads;
#endif
	int fastQver; //33, 62 or 59
	int fastQoutVer; //33, 62 or 59
	//write out quality file, is the input paired End sequenced
	bool BWriteQual, BWriteFastQ;
	bool b_multiOutStream;
	int pairedSeq; //1=single, 2=PE, 3=PE+MID
	bool b_changeFQheadVer; // T/F 0=no PE, 1= XX/1, 2=XX 1:0:3
	bool b_oneLinerFasta; // write one line per sequence?
	bool b_doDereplicate;
	bool b_writePassed;
	bool b_writeMidPass;

	//asynchronous threads
	//std::vector<std::future<ulong>> threads; 
	int Nthrds,thrdsCnt; bool thrdsActive;
	//controls output file size
	int maxReadsPerOFile,ReadsWritten;
	string leadingOutf;
	OptContainer locCmdArgs;
	Dereplicate* Derepl;
	int cntDerep;

	//abstraction to real file type
	//0,1,2,3 refers to pairs (0,1) & singletons (2,3)
	//0=high qual, 1=mid qual
	std::ios_base::openmode wrMode;
	vector<vector<ostream*>> sFile, qFile, fqFile;
	vector<vector<string>> sFileStr, qFileStr, fqFileStr;
	uint totalFileStrms;

	//future<void> derepThread;
	//vector<ostream*> sFile_alt, qFile_alt, fqFile_alt;
	/*ofstream qFile, sFile, fqFile;
	ofstream qFile_alt, sFile_alt, fqFile_alt;
	ofstream qFile2, sFile2, fqFile2;//second pair
	ofstream qFile2_alt, sFile2_alt, fqFile2_alt;
	ofstream qFileS, sFileS, fqFileS;//singleton
	ofstream qFileS_alt, sFileS_alt, fqFileS_alt;
	ofstream qFileS2, sFileS2, fqFileS2;//singleton
	ofstream qFileS2_alt, sFileS2_alt, fqFileS2_alt;
	*/
	
//	streamoff qFilePos, sFilePos, fqFilePos;
//	streamoff qFile2Pos, sFile2Pos, fqFile2Pos;//second pair
//	streamoff qFileSPos, sFileSPos, fqFileSPos;//singleton
//	streamoff qFileS2Pos, sFileS2Pos, fqFileS2Pos;//singleton

};

//fwd declarations
//bool read_fasta_entry(ifstream&fna,ifstream&qual,DNA* in,DNA*,int&);
//DNA* read_fastq_entry(ifstream & fna,int fastQver, int &minQScore,
//					  long& pos);

#endif