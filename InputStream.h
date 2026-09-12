/* sdm: simple demultiplexer
Copyright (C) 2026  Falk Hildebrand

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

#pragma once
#include "DNA.h"
#include "Common.h"
#include "InputFileStream.h"
#include "HtsInput.h"
#include <chrono>
#include <thread>


extern char DNA_trans[256];
extern short DNA_amb[256];
extern short NT_POS[256];
extern short DNA_IUPAC[256 * 256];
typedef double matrixUnit;
// read_occ is defined in DNA.h to avoid circular header dependencies



static void cdbg(const string& x) {
#ifdef DEBUG
    cerr << x;
#endif
}

struct filesStr {//used in separateByFile
    string cramRef;
    bool htsBarcode = false;
    vector<string> FastaF;
    vector<string> QualF;
    vector<string> FastqF;
    vector<string> MIDfq;
    vector<string> fastXtar;
    // Indicates if FASTQ files were submitted
    bool isFastq = true;
    //input path
    string path = "";
    //set up some log structures
    string deLog = "";//dereplication main log
    string logF = "";// (*cmdArgs)["-log"];
    string logFA = "";// (*cmdArgs)["-log"].substr(0, (*cmdArgs)["-log"].length() - 3) + "add.log";

    // Unique Fas initialized with first element of tar (can be b_derep_as_fasta_ and fastq)
    // Contains all unique b_derep_as_fasta_ or fastq files from the mapping file
    unordered_map<string, int> uniqueFastxFiles;

    // idx content: [ [0] ]
    // idx contains one row (vector) for each unique string in tar
    // This vector then contains the indices at which this string occurs in tar
    vector < vector<int> > idx = vector<vector<int>>(0);

    //rewrite uniqueFastxFiles to get it sorted after seqRun..
    vector<pair<string, int>> uniqFxFls;


};




struct multi_tmp_lines {
    multi_tmp_lines() :tmp(0) {}
    multi_tmp_lines(int s) :tmp(0)
    {
        string empty(""); empty.reserve(151); vector<string>tmpLines2(4, empty);
        vector<vector<string>> tmpLines(3, tmpLines2);
        tmp.resize(s, tmpLines);
    }


    size_t size() { return tmp.size(); }
    void setSize(size_t X) { tmp.resize(X); }
    vector< vector< vector< string>>> tmp;
    bool lastInFile = false;
    // Timestamp for instrumentation: set by producer (submit) in microseconds since epoch
    long long submit_ts = 0;
};

//static mutex input_mtx;

class ifbufstream {
public:
    ifbufstream(const string& inF, size_t buf1 = 20000, bool isMC = false, bool test = false)
        : file(inF), modeIO(ios::in), at(0), isGZ(false), atEnd(false), doMC(isMC),
        primary(nullptr), bufS(buf1), bufSW(0), bufferCapacity(buf1),
        totalLines(0), total4Lines(0) {
        if (bufS < 10) {
            cerr << "Buffer size chosen too small: " << bufS << endl << "class ifbufstream\n";
            exit(236);
        }
        iniBufStrm();

        if (isGZfile(file)) {
            isGZ = true;
        }
        openFstream();
        if (!primary) {
            atEnd = true;
            return;
        }
        if (test) return;

        {
            Instr::TimedSharedLockGuard input_lock(input_mtx);
            primary->read(
                keeper.data(),
                static_cast<std::streamsize>(bufferCapacity));
            const size_t read = static_cast<size_t>(primary->gcount());
            bufS = read;
            if (read == 0) {
                atEnd = true;
                keeperW.clear();
            }
            else if (primary->eof() || !primary->good()) {
                atEnd = true;
            }
        }
    }

    ~ifbufstream() {
        cdbg(("destroy ibufstream, total lines: " + to_string(totalLines) + " 4Lines: " + to_string(total4Lines) + "\n"));
        {
            Instr::TimedSharedLockGuard input_lock_del(input_mtx);
            keeper.clear();
            keeperW.clear();
            primary.reset();
        }
    }

    void reset() {
        totalLines = 0;
        total4Lines = 0;

        Instr::TimedSharedLockGuard lock(input_mtx);
        at = 0;
        atEnd = false;
        bufS = bufferCapacity;
        bufSW = 0;
        if (primary) primary->clear();
        primary.reset();
        openFstream();
        iniBufStrm();
        primary->read(
            keeper.data(),
            static_cast<std::streamsize>(bufferCapacity));
        const size_t read = static_cast<size_t>(primary->gcount());
        bufS = read;
        if (read == 0) {
            atEnd = true;
            keeperW.clear();
        }
        else if (primary->eof() || !primary->good()) {
            atEnd = true;
        }
    }

    void setMC(bool b) { doMC = b; }
    bool eof() { return atEnd && at >= bufS; }
    bool operator! (void) { return eof(); }
    void jumpLines(int x = 1) { string mpt; for (int y = 0; y < x; y++) this->getlines(mpt, false); }
    int peekFirstChar() const {
        if (bufS == 0 || keeper.empty()) {
            return std::char_traits<char>::eof();
        }
        return static_cast<unsigned char>(keeper[0]);
    }
    int peerFirstChar() const { return peekFirstChar(); }


    bool getlines(string& ret, bool untilNextFasta = false, bool nwlRspace = false) { // int& linesRead,
        if (atEnd && at >= bufS) return false;
        ret.clear();
        for (;;) {
            // Cache bufS locally to avoid mid-iteration changes
            size_t currentBufS = bufS;
            if (at >= currentBufS && !readChunk()) return !ret.empty();
            // find newline
            currentBufS = bufS; // Re-read after potential readChunk
            size_t rem = (currentBufS > at) ? (currentBufS - at) : 0;
            if (rem == 0) {
                if (!readChunk()) return !ret.empty();
                currentBufS = bufS;
                rem = (currentBufS > at) ? (currentBufS - at) : 0;
                if (rem == 0) return !ret.empty();
            }
            char* start = keeper.data() + at;
            void* found = memchr(start, '\n', rem);
            if (found) {
                size_t pos = (char*)found - keeper.data();
                size_t len = pos - at;
                if (len > 0) {
                    if (keeper[pos - 1] == '\r') ret.append(keeper.data() + at, len - 1);
                    else ret.append(keeper.data() + at, len);
                }
                else {
                    // Handle CRLF split across chunks: previous chunk ended with '\r', current starts with '\n'
                    if (!ret.empty() && ret.back() == '\r') ret.pop_back();
                }
                at = pos + 1;
#ifdef _DEBUG
                totalLines++;
#endif
                if (!untilNextFasta) return true;
                //                linesRead++;
                if (nwlRspace) ret.push_back(' ');
                currentBufS = bufS; // Re-read before next check
                if (at >= currentBufS && !readChunk()) return !ret.empty();
                currentBufS = bufS;
                if (at < currentBufS && keeper[at] == '>') return true;
                continue;
            }
            else {
                ret.append(start, rem);
                at = currentBufS;
                continue;
            }
        }
        return true;
    }

    //specific function that gets fasta line looking for pattern '\n>'
    bool get4lines(vector<string>& in) {
        if (in.size() < 4) in.resize(4);
        in[0].clear(); in[1].clear(); in[2].clear(); in[3].clear();

        const bool l0 = this->getlines(in[0], false);
        const bool l1 = this->getlines(in[1], false);
        const bool l2 = this->getlines(in[2], false);
        const bool l3 = this->getlines(in[3], false);
#ifdef _DEBUG
        total4Lines++;
#endif
        return l0 && l1 && l2 && l3;
    }



    string getInFile() { return file; }

private:
    bool internalReadChunk() {
        const size_t localBufS = bufferCapacity;
        if (localBufS == 0 || !primary || !(*primary)) {
            bufSW = 0;
            atEnd = true;
            keeperW.clear();
            return false;
        }
        if (keeperW.size() < localBufS) {
            keeperW.resize(localBufS);
        }

        primary->read(
            keeperW.data(),
            static_cast<std::streamsize>(localBufS));
        const size_t read = static_cast<size_t>(primary->gcount());
        bufSW = read;
        if (read == 0) {
            atEnd = true;
            keeperW.clear();
        }
        else if (primary->eof() || !primary->good()) {
            atEnd = true;
        }
        return bufSW > 0;
    }

    void openFstream() {
        primary = openInputFile(file, modeIO, bufferCapacity, doMC);
    }

    bool readChunk() {
        {
            std::lock_guard<std::mutex> lg(localMTX);
            bool ok = internalReadChunk();
            if (!ok || bufSW == 0) { return false; }
            keeper.swap(keeperW);
            size_t newBufS = bufSW;
            at = 0;
            bufS = newBufS; // Update bufS after at is reset
        }
        return true;
    }

    void iniBufStrm() {
        cdbg("Ini ibufstream");
        keeper.resize(bufferCapacity);
        keeperW.resize(bufferCapacity);
    }

    string file;
    std::vector<char> keeper;
    std::vector<char> keeperW;
    string locBuffer;
    ios_base::openmode modeIO;
    size_t at;
    bool isGZ, atEnd, doMC;
    std::unique_ptr<std::istream> primary;

    size_t bufS, bufSW, bufferCapacity;
    mutex localMTX;
    shared_mutex input_mtx;
    long long totalLines, total4Lines;
};



class InputStreamer{
public:
    InputStreamer(bool fnRd, qual_score fq, string ignore_IO_errors,
        string pairedRD_HD_out, int threads) :
            _fileLength(10), _max(60), _last(0),
			_globalRdsRead(0), _globalMaxRdsRead(0), _localRdsRead(0),
            fasta_istreams(3, NULL), quality_istreams(3, NULL), fastq_istreams(3, NULL),
            fastaFilepathTemp(3, ""), qualityFilepathTemp(3, ""), 
            dnaTemp1(3, NULL), dnaTemp2(3, NULL),
            isFasta(fnRd), hasMIDs(false),
            keepPairHD(true),
            lnCnt(3, 0), fastQver(fq),
            minQScore(SCHAR_MAX), maxQScore((qual_score) -1),
            QverSet(true), numPairs(1),
            pairs_read(3, 0), opos(3,0),
            currentFile(0), totalFileNumber(0), BCnumber(0),
            qualAbsent(false),
            fqReadSafe(true), fqPassedFQsdt(true),
            fqSolexaFmt(false), 
            ErrorLog(0), DieOnError(true), at_thread(0), num_threads(1), 
		    doTIO(true), verbose(1), N_DNApairsReturned(0), N_DNAreturned(0)
	{
		cdbg("Ini inputstreamer");
		opos[0] = 1; if (fastQver == 0) { QverSet = false; }
		if (ignore_IO_errors =="1") { DieOnError = false; }
		if (pairedRD_HD_out == "0") { keepPairHD = false; }
		num_threads = threads;
		slots.resize(1);//never more than one, otherwise threads read at same time from IO
	}
	~InputStreamer();
//most used routine to get a new DNA entry		// stillMore = is fastq file empty? pos = read pair to return [0/1/-1]
	shared_ptr<DNA> getDNA(int pos);
	bool getDNAlines(vector<string>&,int pos);
	bool getDNAlines(multi_tmp_lines*, int blocks, bool,bool=false);
	vector<shared_ptr<DNA>> getDNApairs();
	//mutli core version
	//vector < shared_ptr<DNA>> getDNAMC();

	//path, fasta, qual_, pairNum
	string setupInput(string path, int tarID, const string& uniqueFastxFile, 
		filesStr& files,
		 int &paired, string onlyPair,
		string& mainFilename, bool simulate = false);
	bool setupFastaQual(string,string, string, int&, string,bool=false);
	void setupFna(string);
	//path, fastq, fastqVer, pairNum
	bool setupFastq(string,string, int&,string,bool simu= false, bool verbose=false);
    void configureHts(const string& reference, bool barcodeTags) {
        cramReference = reference; useHtsBarcode = barcodeTags;
    }
	//shared_ptr<DNA> getDNA(bool &has_next, bool &repair_input_stream, int pos);
	void jumpToNextDNA(bool&, int);
	//shared_ptr<DNA> getDNA2(bool&);
	//shared_ptr<DNA> getDNA_MID(bool&);
	bool hasMIDseqs(){return hasMIDs;}
	void allStreamClose();
	void allStreamReset();//go back to line 1
	void setTIO(bool x) { doTIO = x; }
	void setGlobalRdsRead(uint64_t x) { _globalRdsRead = x; }
	void setMaxRdsRead(uint64_t x) { _globalMaxRdsRead = x; }

	void openMIDseqs(string,string);
	int pairNum() { return numPairs; }
	bool qualityPresent() { return !qualAbsent; }
	bool checkInFileStatus();
	void atFileYofX(uint cF, uint tF, uint BCn) { currentFile = cF; totalFileNumber = tF; BCnumber = BCn; }
	uint getCurFileN() { return currentFile; }

	bool keepPairedHD() { return keepPairHD; }
	qual_score fastQscore() { return fastQver; }
    
	
    vector<ifbufstream*> fasta_istreams, quality_istreams;
	vector<ifbufstream*> fastq_istreams;
    
    

private:
    std::unique_ptr<HtsInput> htsInput;
    string cramReference;
    bool useHtsBarcode = false;
	string current_infiles();
	shared_ptr<DNA> getDNAUnlocked(int pos);
	inline qual_score minmaxQscore(qual_score t, bool& print);// , int lnCnt);
	void minmaxQscore(shared_ptr<DNA> t, bool&);// , int lnCnt);
	bool setupFastq_2(string, string, string);
	bool setupFastaQual2(string, string, string = "fasta file");
	shared_ptr<DNA> read_fastq_entry(istream & fna, qual_score &minQScore,
		int&,bool&,bool);
	shared_ptr<DNA> read_fastq_entry_fast(istream & fna, int&,bool&);
	void jmp_fastq(istream &, int&);
	bool read_fasta_entry(ifbufstream*fasta_is, ifbufstream*quality_is, shared_ptr<DNA> in, shared_ptr<DNA>, int&);
	//static bool getFastaQualLine(istream&fna,  string&);
	void maxminQualWarns_fq();
	qual_score auto_fq_version();
	qual_score auto_fq_version(qual_score minQScore, qual_score maxQScore=0);
	qual_score auto_fq_version(const vector<string>& ret);
	bool detectAndSwitchFileFormat(ifbufstream* stream, const std::string& filepath);
	void resetStats() {
		for (size_t i = 0; i < lnCnt.size(); i++) { lnCnt[i] = 0; }
		for (size_t i = 0; i < pairs_read.size(); i++) { pairs_read[i] = 0; }
		for (size_t i = 0; i < opos.size(); i++) { opos[i] = 0; }
		currentFile = 0; totalFileNumber = 0; BCnumber = 0;
		_localRdsRead = 0;
	}
	bool desync(int pos) { if ( abs(pairs_read[pos] - pairs_read[opos[pos]]) > 1 ) {return true; } return false; }
	void IO_Error(string x);
	//bar on file read progress
	void _measure(istream &);
	inline bool _drawbar(istream &);
	inline void _print(int cur, float prog);
	
	int _fileLength, _max, _last;
	uint64_t _globalRdsRead, _globalMaxRdsRead, _localRdsRead; //processed read counts and optional limit
	
	//abstraction to real file type

	//required for Fasta in term storage
	vector<shared_ptr<DNA>> dnaTemp1;
	vector<shared_ptr<DNA>> dnaTemp2;

	vector<string> fastaFilepathTemp, qualityFilepathTemp;
	//fastqFilepathTemp;
	//0,1,2 refers to pairs / MID fasta files
	//vector<ifstream> fna, qual_, fastq;
	//ifstream qual_, fastq,
		//second pair_
		//ifstream fna2, qual2, fastq2,
		//usually used for MID
		//fna3, qual3, fastq3;

	shared_ptr<DNA> tdn21; shared_ptr<DNA> tdn22;
	shared_ptr<DNA> tdn31; shared_ptr<DNA> tdn32;
	bool isFasta, hasMIDs;
	bool keepPairHD;
	vector<int> lnCnt;// , lnCnt2, lnCnt3;//line count
	qual_score fastQver,minQScore,maxQScore;//which version of Fastq? minima encountered Qscore
	mutex fqverMTX;
	bool QverSet;
	//1 or 2?
	int numPairs;
	//keep track of sequences read for each pair_; other position (1=2,2=1)
	vector<int> pairs_read, opos;
	//some stats to print, nothing really relevant
	uint currentFile, totalFileNumber, BCnumber;
	//is quality information even available?
	bool qualAbsent;
	//fq format not checked for completeness
	bool fqReadSafe, fqPassedFQsdt, fqSolexaFmt;

	//collects errors, handles errors
	vector<string> ErrorLog;
	bool DieOnError;

    shared_mutex protect;
	int num_threads;
	bool doTIO;//multi thread IO?
	int at_thread; 
#ifdef togRead
	vector <jobC> slots;
#else
	vector <job3> slots;
#endif
	int verbose;
    uint64_t N_DNApairsReturned, N_DNAreturned;
};




//outdated, don't use
bool whoIsBetter_old(shared_ptr<DNA> d1, shared_ptr<DNA> d2, shared_ptr<DNA> dM,
    shared_ptr<DNA> r1, shared_ptr<DNA> r2, shared_ptr<DNA> rM,
    float& ever_best, bool forSeed);

void setDNAseedSelectionClassic(bool enabled);
bool getDNAseedSelectionClassic();

//true: d is better, false: ref is better
bool whoIsBetter(shared_ptr<DNA> d1, shared_ptr<DNA> d2, shared_ptr<DNA> dM,
	shared_ptr<DNA> r1, shared_ptr<DNA> r2, shared_ptr<DNA> rM,
    float& ever_best, bool forSeed, DNAunique* merge_stats_owner = nullptr, 
    uint64_t dSiz = 2, uint64_t rSiz=2);


bool whoIsBetter(shared_ptr<DNAunique> d1, shared_ptr<DNAunique> r1, float& ever_best, bool forSeed);

double logRatio(double v1, double v2);
double loglogRatio(double v1, double v2);
double pointFiveRatio(double v1, double v2, double scale = 0.5);






