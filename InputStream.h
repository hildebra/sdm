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


#ifndef _InputStr_h
#define _InputStr_h

//input through combined getDNApairs
#define togRe//ad
//input through getDNAlines
#define IOsep

#include "DNAconsts.h"
#include "FastxReader.h"
#include <functional> 
#include <cctype>
#include <locale>
#include <climits>
#include "ThreadPool.h"
#include <fstream>

extern char DNA_trans[256];
extern short DNA_amb[256];
extern short NT_POS[256];
extern short DNA_IUPAC[256 * 256];
typedef double matrixUnit;
typedef robin_hood::unordered_flat_map<int, long> read_occ;


string spaceX(uint k);
int digitsInt(int x);
int digitsFlt(float x);
string intwithcommas(int value);
std::string itos(int number);
std::string ftos(float number);
bool isGZfile(const string fileS);//test if file is gzipped input



static void rtrim(std::string& s) {
	//s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	//return s;
	s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
		return !std::isspace(ch);
		}).base(), s.end());
}
inline int parseInt(const char** p1);

static void cdbg(const string& x) {
#ifdef DEBUG
	cerr << x;
#endif
}

//MOCAT header fix
std::vector<std::string> header_string_split(const std::string str, const std::string sep);
void remove_paired_info(string&, short = -1);
//MOCAT header fix
std::string header_stem(string& header);
std::istream& safeGetline(std::istream& is, std::string& t);
string reverseTS2(const std::string & Seq);
void reverseTS(std::string & Seq);


bool any_lowered(const string& is);
//this function changes input string (file location) to have consistent file names
string applyFileIT(string x, int it, const string xtr = "");
bool fileExists(const std::string& name, int i=-1,bool extiffail=true);
//vector<int> orderOfVec(vector<int>&);

static std::mutex input_mtx;

class ifbufstream {//: private std::streambuf, public std::ostream {
public:
	ifbufstream(const string& inF, size_t buf1=20000,bool test=false) :
		file(inF),modeIO(ios::in),at(0),isGZ(false), atEnd(false), hasKickoff(false), bufS(buf1)
	{
		if (bufS < 10) {
			cerr << "Buffer size chosen too small: " << bufS << endl << "class ifbufstream\n";
			exit(236);
		}

		keeper = new char[bufS];
		keeperW = new char[bufS];
		for (int x = 0; x < bufS; x++) {
			keeper[x] = EOF; keeperW[x] = EOF;
		}
		
		if (isGZfile(file)) { //write a gzip out??
			isGZ = true;
#ifndef _gzipread
			cerr << "ifbufstream::gzip input not supported in your sdm build\n" << file << endl;
			exit(51);
#endif
		}
		if (isGZ) {
#ifdef _gzipread
			primary = new zstr::ifstream(file.c_str());
#else
			cerr << "gzip not supported in your sdm build\n (ifbufstream) " << file; exit(50);
#endif
		}
		else {
			primary = new ifstream(file, modeIO);
		}
		//primary->set_rdbuf(0);

		if (!primary){
			atEnd = true;
			return;
		}
		if (test) {
			return;
		}
		/*primary.seekg(0, primary.end);
		int length = primary.tellg();
		primary.seekg(0, is.beg);*/
		//first round read..
		primary->read(keeper, bufS);
		if (!(*primary) || bufS > primary->gcount()) {
			bufS = primary->gcount()+1;
			atEnd = true;
			delete[] keeperW;
			keeperW = nullptr;
		} else {
			kickOff();
		}
	}
	~ifbufstream() {
		if (hasKickoff) { hasKickoff = false; readKickoff.get(); }
		delete[] keeper;
		if (keeperW != nullptr) { delete[] keeperW; }
		delete primary;
	}
	void clear() {
		at = 0;
		primary->clear();
		delete primary;
		primary = new ifstream(file, modeIO);
	}
	bool eof() {
		return atEnd && at >= bufS;
	}
	bool operator! (void) {
		return !atEnd;
	}
	void jumpLines(int x=1) {
		string mpt;
		for (int y = 0; y < x; y++) {
			this->getline(mpt);
		}
	}
	//specific function that gets fasta line looking for pattern '\n>'
	bool getlines(string& ret,int & linesRead) {
		ret.clear();
		if (atEnd && at >= bufS) {
			return false;
		}
		for (;;) {
			char c = keeper[at];
			at++;
			if (at >= bufS && !readChunk()) {
				return false;// read next chunk already
			}

			switch (c) {
			case '\n':
				linesRead++;
				if (keeper[at] == '>') {
					return true;
				}
				break;
			case '\r':
				if (keeper[at] == '\n') {
					linesRead++;
					at += 1;
					if (at < bufS && keeper[at] == '>') {
						return true;
					}
					break;
				}
			case EOF:
				// Also handle the case when the last line has no line ending
				atEnd = true;
				//primary->setstate(std::ios::eofbit);
				return false;
			default:
				ret += c;
			}
		}
		if (atEnd && at >= bufS) {
			return false;
		}
		return true;
	}
	bool getline(string& ret) {
		ret.clear();
		if (atEnd && at >= bufS) {
			return false;
		}
		for (;;) {
			char c = keeper[at];
			at++;
			if (at >= bufS && !readChunk()) {
				return false;// read next chunk already
			}

			switch (c) {
			case '\n':				
				return true;
			case '\r':
				if (keeper[at] == '\n') {
					at+=1; 
					return true;
				}
			case EOF:
				// Also handle the case when the last line has no line ending
				atEnd = true;
				//primary->setstate(std::ios::eofbit);
				return false;
			default:
				ret += c;
			}
		}
		if (atEnd && at >= bufS) {
			return false;
		}
		return true;
	}
private:
	bool internalReadChunk() {
		if (!*(primary)) {
			if (keeperW != nullptr) { delete[] keeperW; }
			keeperW = nullptr;
			atEnd = true;
			bufS = primary->gcount()+1;
			//keeperW[XX] = char_traits<char>::eof();
			return false;
		}
		input_mtx.lock();
		primary->read(keeperW, bufS);
		
		input_mtx.unlock();
		return true;
	}
	bool kickOff() {
		if (true ) {//do MC?
			if (!atEnd) {
				assert(!hasKickoff);
				localMTX.lock();
				hasKickoff = true;
				readKickoff = async(std::launch::async, &ifbufstream::internalReadChunk,
					this);
				localMTX.unlock();
				return true;
			}	else {
				return false;
			}
		}
		else {
			// Without multithreading
			return internalReadChunk();
		}
	}
	bool readChunk() {
		if (hasKickoff) { hasKickoff = false; readKickoff.get(); }
		if (atEnd && keeperW == nullptr) {
			return false;
		}
		
		//do a swap, no memcpy needed
		char* tmp = keeper;
		keeper = keeperW;
		keeperW = tmp;
		//memcpy(keeper, keeperW, bufS);

		at = 0;
		kickOff();
		

		//primary->read(keeperW, bufS);
		return true;
	}
	string file;
	char* keeper; char* keeperW;
	ios_base::openmode modeIO;
	size_t at;
	bool  isGZ, atEnd, hasKickoff;
	istream* primary;
	future<bool> readKickoff;
	size_t bufS ;
	mutex localMTX;
};
std::ptrdiff_t len_common_prefix_base(char const a[], char const b[]);
static std::mutex output_mtx;
class ofbufstream {//: private std::streambuf, public std::ostream {
public:
	ofbufstream(size_t bufferS):file("T"), modeIO(ios::app), used(0),
		coutW(true), isGZ(false){
		int x = 0;
	}
	ofbufstream(const string IF, int mif, size_t bufferS=20000) :
		file(IF), modeIO(mif), used(0),
		coutW(false), isGZ(false),primary(nullptr), bufS(bufferS){
		
		if (modeIO == ios::out) {
			remove(file.c_str());
		}
		if (file == "T") {//write to ostream
			coutW = true;
			keeper = new char[0];
			keeperW = new char[0];
			return;
		}
		keeper = new char[bufS];
		//second keeper swappable with keeper to always have one added to, one written out (kickoff system)
		keeperW = new char[bufS];
		if (isGZfile(IF)) { //write a gzip out??
			isGZ = true; 
#ifndef _gzipread
			cerr << "ofbufstream::gzip outpout not supported in your sdm build\n" << file<<endl; 
			exit(51);
#endif
			// ostream* os = new zstr::ofstream("output.txt", std::ios::app);
		}
	}

	~ofbufstream() {
		if (hasKickoff) { hasKickoff = false; writeKickoff.get(); }
		writeStream(false);
		deactivate();
		delete[] keeper;
		delete[] keeperW;
	}
	bool operator! (void) {
		return false;
	}
	void operator<< (const string& X) {
        {
			if (coutW) {
				cout << X;
				return;
			}
			
			//std::unique_lock<std::mutex> lock(append_mtx_);
            size_t lX(X.length());
            //writeStream();
            if (lX + used > bufS) {
                writeStream();
            }
			append_mtx_.lock();
            memcpy(keeper + used, X.c_str(), lX);
            used += lX;
			append_mtx_.unlock();
        }
    }

    // Multithreading through Threadpool
    void setThreadPool(ThreadPool *pool) {this->pool = pool;}
	void emptyStream() {
		used = 0;
	}
	void activate() {
		if (primary != nullptr) {
			return;
		}
		if (isGZ) {
#ifdef _gzipread
			primary = new zstr::ofstream(file.c_str(), std::ios::app);
#else
			cerr << "ofbufstream::gzip outpout not supported in your sdm build\n" << file << endl;
			exit(51);
#endif
		}
		else {
			primary = new ofstream(file.c_str(), ios::app);
		}
	}
	void deactivate() {
		if (primary == nullptr) {
			return;
		}
		//primary->close();
		delete primary;// ->close();
		primary = nullptr;

	}
    // end

private:
    std::mutex append_mtx_;
	bool internalWrite(bool closeThis){
		output_mtx.lock();
		primary->write(keeperW, used);
		output_mtx.unlock();
		//delete[] keeperW;//clean up
		return true;
		if (closeThis) {
			deactivate();
		}
	}
	void writeStream(bool doKickoff=true) {
        static int counter = 0;
        //out << omp_get_thread_num << ": " << (counter++) << endl;
        if (used == 0) {
            return;
        }
		//do this first to prevent keeper/keeperW getting corrupted
		if (hasKickoff) { hasKickoff = false; writeKickoff.get(); }
		bool closeThis = true;
		if (primary == nullptr) {
			this->activate();
			closeThis = false;
		}
		//swap opeartion
		char* keeperTmp = keeperW;
		keeperW = keeper;
		keeper = keeperTmp;
		//keeper = new char[bufS];
        //if (false) {
        /*if ( pool ) {
            // With multithreading
            std::string keeper_copy = string(keeper, used);
            std::string file_copy = file;
            pool->enqueue([keeper_copy, file_copy] {
                write(std::move(keeper_copy), file_copy);
            });
        /**/
		//multithreading via kickoff
		if (doKickoff){
			writeKickoff = async(std::launch::async, &ofbufstream::internalWrite, 
				this, closeThis);
			hasKickoff = true;
		} else {
            // Without multithreading
			internalWrite(closeThis);
        }
		used = 0;
    }
	string file;
	char *keeper;
	char* keeperW;
	int modeIO;
	size_t used;
	bool coutW, isGZ;
	ostream* primary;
	
    size_t bufS;

    // Multithreading throw threadpool
    ThreadPool *pool = nullptr;//currently just used to check if MC or not
    // end
	future<bool> writeKickoff;
	bool hasKickoff=false;

    void write(std::string s, std::string file) {
        {
            // Lock the input area
            //std::lock_guard<std::mutex> guard(output_mtx);
            ofstream of(file.c_str(), ios::app);
            of.write(s.c_str(), s.size());
            of.close();
        }
    }
};



inline vector<string> splitByComma(const string& fileS,bool requireTwo, char SrchStr=','){
	string::size_type pos = fileS.find(SrchStr);
	if (pos == string::npos){
		if (requireTwo){
			cerr<<fileS<<endl;
			cerr << "Could not find \"" << SrchStr<<"\" in input file (required for paired end sequences, as these are two files separated by \",\")\n";
			exit(14);
		}else {
			return vector<string> (1,fileS);
		}
	}
	vector<string> tfas(2,"");
	tfas[0] = fileS.substr(0,pos);
	tfas[1] = fileS.substr(pos+1);
	return tfas;
} 

inline vector<string> splitByCommas(const string& fileS, char SrchStr = ',') {
	if (fileS.find(SrchStr) == string::npos) { return vector<string>(1, fileS); }
	vector<string> res = splitByComma(fileS, true, SrchStr);
	vector<string> ret(0); ret.push_back(res[0]);
	while (res[1].find(SrchStr) != string::npos) {
		res = splitByComma(res[1], true, SrchStr);
		ret.push_back(res[0]);
	}
	ret.push_back(res[1]);
	return ret;
}


//requires sorted vector with the entries being actual datapoints
template<class TYPE>
TYPE calc_median2(vector<TYPE>& in, float perc){
	size_t sum = in.size();
	size_t tar = (size_t)(((float)sum) * perc);
	return in[tar];
}


//returns "i_fna" or "i_fastq"
string detectSeqFmt(const string);


class Filters;

class DNA{
	friend class DNAunique;
public:
	DNA(string seq, string names) : sequence_(seq), sequence_length_(sequence_.length()),
                                    id_(names), new_id_(names),
                                    qual_(0), qual_traf_(""), sample_id_(-1), avg_qual_(-1.f),
                                    quality_sum_(0), accumulated_error_(0.), good_quality_(false), mid_quality_(false),
                                    reversed_(false),
                                    read_position_(-1),
                                    FtsDetected(),
                                    id_fixed_(false), tempFloat(0.f),merge_offset_(-1){}
	DNA(): sequence_(""), sequence_length_(0), id_(""), new_id_(""), qual_(0), qual_traf_(""),
           sample_id_(-1), avg_qual_(-1.f),
           quality_sum_(0), accumulated_error_(0.), good_quality_(false), mid_quality_(false),
           reversed_(false),
           read_position_(-1),
           FtsDetected(),
           id_fixed_(false), tempFloat(0.f), merge_offset_(-1) {
	}
	//starts with fastx record
	DNA(FastxRecord*, qual_score & minQScore, qual_score & maxQScore, qual_score fastQver);
	//works directly with 4 lines in fastq format
	DNA(vector<string> fq, qual_score fastQver);
	DNA(vector<string> fq);

	~DNA() {
	    //cout << "destruct" << endl;
	}
	
	bool operator==(DNA i) {
		if (i.getSeqPseudo() == this->getSeqPseudo()) {
			return true;
		}	else {
			return false;
		}
	}
	bool operator==(shared_ptr<DNA> i) {
		if (i->getSeqPseudo() == this->getSeqPseudo()) {
			return true;
		} else {
			return false;
		}
	}
	//something wrong with DNA object, just del all info
	void delself() {
		sequence_ = ""; sequence_length_ = 0; id_ = ""; new_id_ = ""; qual_.resize(0);
		good_quality_ = false; mid_quality_ = false;
	}
	//~DNA(){}
	void appendSequence(const string &s) { sequence_ += s; sequence_length_ = sequence_.length(); }
	void appendQuality(const vector<qual_score> &q);
	void fixQ0(void);
    void setSequence(string & s) {
        sequence_ = s;
        sequence_length_ = sequence_.length();
    }	void setSequence(string && s) {
        sequence_ = s;
        sequence_length_ = sequence_.length();
    }

	string &getSequence() { return sequence_; }

	string getSeqPseudo() {
	    return sequence_.substr(0, sequence_length_);
	}

	void setQual(vector<qual_score>& Q) { qual_ = Q; avg_qual_ = -1.f; }
    void setQual(vector<qual_score>&& Q) { qual_ = Q; avg_qual_ = -1.f; }

	const string& getId() {
	    if (id_fixed_) {
	        return new_id_;
	    }
	    return id_;
	}

	string getPositionFreeId(); // remove /1 /2 #1:0 etc
	const string& getOldId() { return id_; }
	string getShortId() { return id_.substr(0, getShorterHeadPos(id_)); }
	string getNewIDshort() { return new_id_.substr(0, getShorterHeadPos(new_id_)); }
	bool seal(bool isFasta=false);

	bool isEmpty() {
	    if (id_.empty() && sequence_.empty()) {
	        this->setPassed(false);
	        return true;
	    }
	    return false;
	}

	int numACGT();
	void stripLeadEndN();
	int numNonCanonicalDNA(bool);
	float getAvgQual();
	unsigned int getQsum(){return quality_sum_;}
	float qualWinfloat(unsigned int,float,int&);
	
	
	float binomialFilter(int, float);
	//float qualWinfloat_hybr(int,float,int,float,int&);

	bool qualWinPos(unsigned int,float);
	bool qualAccumTrim(double d);
	int qualAccumulate(double d);

	double getAccumError(){
		if (accumulated_error_ == 0.f) {
		    for (uint i = 0; i < qual_.size(); i++) {
		        if (qual_[i] >= 0) {
                    accumulated_error_ += SAqualP[qual_[i]];
		        }
		    }
		}
		if (std::isinf((double) accumulated_error_)) {
            accumulated_error_ = 5.f;
		}
		return accumulated_error_;}
	int minQual(){int mq=50; for (uint i=0; i < qual_.size(); i++){if (qual_[i] < mq){ mq=qual_[i];}}return mq;}
	void ntSpecQualScores(vector<long>&, vector<long>&);
    void ntSpecQualScoresMT(vector<long>&, vector<long>&);


	inline uint length() { return (uint) sequence_length_; }
	inline uint mem_length() { return (uint) sequence_.length(); }
	bool cutSeq(int start, int stop=-1, bool = false);
	bool cutSeqPseudo(int start) { return cutSeq(start, -1, true); }
	bool HomoNTRuns(int);
	int matchSeq(string, int, int, int);
	void reverse_transcribe();
	int matchSeqRev(string, int, int, int=0);
	int matchSeq_tot(string, int, int, int&);
	void writeSeq(ostream&, bool singleLine = false);
	string writeSeq( bool singleLine = false);
	void writeQual(ostream&, bool singleLine = false);
	string writeQual(bool singleLine = false);
	void writeFastQ(ostream&, bool = true);
	string writeFastQ(bool = true);
	//void writeFastQ(ofbufstream&, bool = true);
	void writeFastQEmpty(ostream&);
	void setNewID(string x) { new_id_ = x; }
	void setHeader(string x){ new_id_=x;id_=x;}
	void changeHeadPver(int ver);
	void setTA_cut(bool x) { FtsDetected.TA_cut = x; }
	bool getTA_cut() { return FtsDetected.TA_cut; }
	void setBarcodeCut() { FtsDetected.barcode_cut = true; FtsDetected.barcode_detected = true; }
	bool getBarcodeCut() { return FtsDetected.barcode_cut; }
	void setBarcodeDetected(bool x){ FtsDetected.barcode_detected = x; }
	bool getBarcodeDetected() { return FtsDetected.barcode_detected; }
	bool isMIDseq() { if (read_position_ == 3) { return true; } return false; }
	void setMIDseq(bool b){ if (b){ read_position_ = 3; } }
	void setpairFWD(){ read_position_ = 0; }
	void setpairREV(){ read_position_ = 1; }
	int getReadMatePos() { return (int) read_position_; }
	bool sameHead(shared_ptr<DNA>);
	bool sameHead(const string&);
	//inline void reverseTranscribe();
	void setTempFloat(float i){tempFloat = i;}
	float getTempFloat(){return tempFloat;}
	//void adaptHead(shared_ptr<DNA>,const int,const int);
	void failed(){ good_quality_=false;mid_quality_=false;}
	bool control(){ if (qual_.size() == 0){return false;}return true;}
	
	void setBCnumber(int i, int BCoff) {
	    if (i < 0) {
            sample_id_ = i ;
	        FtsDetected.barcode_detected = false;
	    } else {
            sample_id_ = i + BCoff;
	        FtsDetected.barcode_detected = true;
	    }
	}
	
	
	int getBarcodeNumber() const;//always return BC tag IDX global (no local filter idx accounted for, use getBCoffset() to correct)

	void prepareWrite(int fastQver);
	void reset();
	void resetTruncation() { sequence_length_ = sequence_.length(); }
	void setPassed(bool b);
	void setMidQual(bool b) { mid_quality_ = b; }
	bool isGreenQual(void){return good_quality_;}
	bool isYellowQual(void){return mid_quality_;}
	string getSubSeq(int sta, int sto){return sequence_.substr(sta, sto);}
	void resetQualOffset(int off, bool solexaFmt);
	
	//control & check what happened to any primers (if)
	bool has2PrimersDetected() { return (FtsDetected.reverse && FtsDetected.forward); }
	bool getRevPrimCut() { return FtsDetected.reverse; }
	bool getFwdPrimCut() { return FtsDetected.forward; }
	void setRevPrimCut() { FtsDetected.reverse = true; }
	void setFwdPrimCut() { FtsDetected.forward = true; }
	void setDereplicated() { FtsDetected.dereplicated = true; }
	bool isDereplicated() { return FtsDetected.dereplicated ; }
	void constellationPairRev(bool b) { FtsDetected.revPairConstellation = b; }
	bool isConstellationPairRev() { return FtsDetected.revPairConstellation ; }
	//only used in pre best seed step
	//float getSeedScore() { return tempFloat; }
	//void setSeedScore(float i) { tempFloat = (float)i; }

	struct QualStats {
		bool maxL; bool PrimerFail; bool AvgQual; //sAvgQual
			bool HomoNT; bool PrimerRevFail; bool minL; 
			bool minLqualTrim; //<-sMinQTrim trimmed due to quality
			bool TagFail; bool MaxAmb; bool QualWin;//sQualWin 
			bool AccErrTrimmed; bool QWinTrimmed;  // either of these makes bool Trimmed; 
			bool fail_correct_BC; bool suc_correct_BC; bool
			failedDNAread; 
			//bool adapterRem; -> setTA_cut
			bool RevPrimFound; 
			bool BinomialErr; bool dblTagFail;
		QualStats() :
			maxL(false), PrimerFail(false), AvgQual(false), HomoNT(false),
			PrimerRevFail(true), minL(false), minLqualTrim(false),
			TagFail(false), MaxAmb(false), QualWin(false),
			AccErrTrimmed(false), QWinTrimmed(false),
			fail_correct_BC(false), suc_correct_BC(false),
			failedDNAread(false), RevPrimFound(false),
			BinomialErr(false),
			dblTagFail(false)
		{}
	} QualCtrl;

protected:
	size_t getShorterHeadPos(const string & x, int fastQheadVer = -1);
	//mainly used to mark if rev/Fwd primer was detected
	string xtraHdStr();
	size_t getSpaceHeadPos(const string & x);
	//binomial accumulated error calc
	inline float interpolate(int errors1, float prob1, int errors2, float prob2, float alpha);
	float sum_of_binomials(const float j, int k, float n, int qual_length, const vector<float>& error_probs, const vector< vector<float>> & per_position_accum_probs);
	inline float prob_j_errors(float p, float j, float n);	
	inline bool matchDNA(char,char);

	std::string sequence_;
	size_t sequence_length_;
	string id_, new_id_; //original and newly constructed id_
	vector<qual_score> qual_;

public:
    const vector<qual_score> &getQual() const;

    int merge_seed_pos_ = -1;
    int seed_length_ = -1;
    int merge_offset_ = -1;
    bool reversed_merge_ = false;

	void getEssentialsFts(shared_ptr<DNA> oD) {
		FtsDetected.forward = oD->FtsDetected.forward;
		FtsDetected.reverse = oD->FtsDetected.reverse;
		FtsDetected.TA_cut = oD->FtsDetected.TA_cut;
		FtsDetected.barcode_detected = oD->FtsDetected.barcode_detected;
		FtsDetected.barcode_cut = oD->FtsDetected.barcode_cut;
		FtsDetected.dereplicated = oD->FtsDetected.dereplicated;

		read_position_ = oD->read_position_;
		good_quality_ = oD->good_quality_;
		mid_quality_ = oD-> mid_quality_;
		sample_id_ = oD->sample_id_;

	}

protected:
    std::string qual_traf_;
	int sample_id_;

	//const char* DN;
	float avg_qual_;
	unsigned int quality_sum_;
	double accumulated_error_;
	bool good_quality_, mid_quality_;
	bool reversed_;


	short read_position_;//-1=unkown; 0=pair1 (fwd primer); 1=pair2 (rev primer); 3=MID seq ;

	struct ElementsDetection {
		bool forward; bool reverse;//primers detected
		bool TA_cut; bool barcode_detected;  bool barcode_cut; bool dereplicated;
		bool revPairConstellation;

		ElementsDetection() : forward(false), reverse(false), TA_cut(false), barcode_detected(false),
                              barcode_cut(false), dereplicated(false), revPairConstellation(false) {}
		void transferFrom(ElementsDetection& o) { 
			forward = o.forward;
			reverse = o.reverse;
			TA_cut = o.TA_cut;
			barcode_detected = o.barcode_detected;
			barcode_cut = o.barcode_cut;
			dereplicated = o.dereplicated;
			revPairConstellation = o.revPairConstellation;
		}
		void reset() {
		    forward = false;
		    reverse = false;
		    TA_cut = false;
		    barcode_detected = false;
            barcode_cut = false;
            dereplicated = false;
			revPairConstellation = false;
		}

	} FtsDetected;



	bool id_fixed_;
	float tempFloat;
};

struct DNAHasher
{
	size_t operator()(shared_ptr<DNA> k) const
	{
		// Compute individual hash values for two data members and combine them using XOR and bit shifting
		return ((hash<string>()(k->getSeqPseudo())) >> 1);
	}
};





class DNAunique : public DNA {//used for dereplication
	friend class DNA;
public:
    // Constructors and Destructors
	DNAunique() : DNA(), count_(0), pair_(0) {}
	DNAunique(string s, string x) : DNA(s, x), count_(1) {}

	// Mostly used constructor
	DNAunique(shared_ptr<DNA>d, int BC) : DNA(*d), count_(0), best_seed_length_((uint)sequence_.size()), pair_(0) {
        incrementSampleCounter(BC);
	}
	~DNAunique() {
		/*if (quality_sum_per_base_) {
			std::cout << "delete quality_per_sum_base" << std::endl;
		}*/
//        std::cout << (quality_sum_per_base_ == nullptr) << std::endl;
        delete[] quality_sum_per_base_;
    }

	//string sequence_;	string id_;
	void Count2Head(bool);
	void incrementSampleCounter(int sample_id);
	void writeMap(ofstream & os, const string&, vector<int>&, const vector<int>&);
	inline int getCount() { return count_; }
	uint getBestSeedLength() { return best_seed_length_; }
	void setBestSeedLength(uint i) { DNAuniMTX.lock(); best_seed_length_ = i; DNAuniMTX.unlock();}
	void incrementSampleCounterBy(int sample_id, long count);
	void transferOccurence(shared_ptr<DNAunique>);
	const read_occ& getDerepMap() { return occurence; }
	vector<int> getDerepMapSort(size_t);
	//vector<pair_<int, int>> getDerepMapSort2(size_t wh);
	//void getDerepMapSort(vector<int>&, vector<int>&);

	void prepSumQuals() {
		if (!quality_sum_per_base_) {
			quality_sum_per_base_ = new uint64_t[length()];
			for (uint i = 0; i < length(); i++) {
				quality_sum_per_base_[i] = qual_[i];
			}
		}
	}
	void sumQualities(shared_ptr<DNAunique> dna) {
		if (dna == nullptr) return;
		prepSumQuals();

		if (dna->quality_sum_per_base_) {
			for (uint i = 0; i < length(); i++) {
				quality_sum_per_base_[i] += dna->quality_sum_per_base_[i];
			}
		}
		else {
			for (uint i = 0; i < length(); i++) {
				quality_sum_per_base_[i] += dna->qual_[i];
			}
		}
	}
	void sumQualities(shared_ptr<DNA> dna) {
		if (dna == nullptr) return;
		DNAuniMTX.lock();
		prepSumQuals();
		const vector<qual_score> quals = dna->getQual();
		for (uint i = 0; i < length(); i++) {
			quality_sum_per_base_[i] += quals[i];
		}
		DNAuniMTX.unlock();
	}

	void prepareDerepQualities(int ofastQver);
    void writeDerepFastQ(ofstream &, bool = true);

	void saveMem() {
        qual_traf_ = "";
        new_id_ = id_.substr(0, getSpaceHeadPos(id_));
        id_ = "";
	}


	void attachPair(shared_ptr<DNAunique> dna_unique) {
	    pair_ = dna_unique;
	    pair_->saveMem();
	}

	shared_ptr<DNAunique> getPair(void) {
	    return pair_;
	}

	//estimates if one sample occurence covers the unique counts required for sample specific derep min counts
	bool pass_deprep_smplSpc( const vector<int>&);
	int counts() const { return count_; }

	void takeOver(shared_ptr<DNAunique> dna_unique_old, shared_ptr<DNA> dna2);
	void takeOverDNA(shared_ptr<DNA> dna1, shared_ptr<DNA> dna2);
	uint64_t * transferPerBaseQualitySum();

    uint64_t* quality_sum_per_base_ = nullptr;
private:
	int count_;
	//matrixUnit chimeraCnt;
	int best_seed_length_;
	
	read_occ occurence;
	mutex DNAuniMTX;
	shared_ptr<DNAunique> pair_;

	// to calculate mean
    std::string qualities_avg_;

	//threadsafe
};


//multi threading Input Streamer
struct job3 {
	job3() :inUse(false) {}
	std::future <shared_ptr<DNA>> fut;
	std::future <shared_ptr<DNA>> fut2;
	std::future <shared_ptr<DNA>> fut3;
	bool inUse = false;
};
//multi threading Input Streamer
struct jobC {
	jobC() :inUse(false) {}
	std::future <vector<shared_ptr<DNA>>> Cfut;
	bool inUse = false;
};

shared_ptr<DNA> str2DNA(vector<string> in, bool keepPairHD, int fastQver, int readpos);

class InputStreamer{
public:
	InputStreamer(bool fnRd, qual_score fq, string ignore_IO_errors,string pairedRD_HD_out, int threads) :
            _fileLength(10), _max(60), _last(0),
            fasta_istreams(3, NULL), quality_istreams(3, NULL), fastq_istreams(3, NULL),
            fastaFilepathTemp(3, ""), qualityFilepathTemp(3, ""), 
		//fna(3,NULL), qual_(3,NULL), fastq(3,NULL),
		
#ifdef _gzipread	
		//gzfna(3,NULL), gzqual(3,NULL), gzfastq(3,NULL),
#endif
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
            ErrorLog(0), DieOnError(true), at_thread(0), num_threads(1)
	{
		opos[0] = 1; if (fastQver == 0) { QverSet = false; }
		if (ignore_IO_errors =="1") { DieOnError = false; }
		if (pairedRD_HD_out == "0") { keepPairHD = false; }
		num_threads = threads;
		slots.resize(1);//never more than one, otherwise threads read at same time from IO
	}
	~InputStreamer();
	
//most used routine to get a new DNA entry		// stillMore = is fastq file empty? pos = read pair to return [0/1/-1]
	shared_ptr<DNA> getDNA(int pos);
	void getDNAlines(vector<string>&,int pos);
	vector<shared_ptr<DNA>> getDNApairs();
	//mutli core version
	vector < shared_ptr<DNA>> getDNAMC();

	//path, fasta, qual_, pairNum
	string setupInput(string path, int tarID, const string& uniqueFastxFile, 
		const vector<string>& fastqFiles, 
		const vector<string>& fastaFiles, const vector<string>& qualityFiles,
		const vector<string>& midFiles, int &paired, string onlyPair,
		string& mainFilename, bool simulate = false);
	bool setupFastaQual(string,string, string, int&, string,bool=false);
	void setupFna(string);
	//path, fastq, fastqVer, pairNum
	bool setupFastq(string,string, int&,string,bool simu= false, bool verbose=false);
	//shared_ptr<DNA> getDNA(bool &has_next, bool &repair_input_stream, int pos);
	//void getDNA(FastxRecord* FR1, FastxRecord* FR2, FastxRecord* FRM,shared_ptr<DNA>* dna1, shared_ptr<DNA>* dna2, shared_ptr<DNA>* mid);
	void jumpToNextDNA(bool&, int);
	//shared_ptr<DNA> getDNA2(bool&);
	//shared_ptr<DNA> getDNA_MID(bool&);
	bool hasMIDseqs(){return hasMIDs;}
	void allStreamClose();
	void allStreamReset();
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
    
    
  //  void getDNA(FastxRecord *read1, FastxRecord *read2, FastxRecord *quality1, FastxRecord *quality2, FastxRecord *FRM,
   //             shared_ptr<DNA> *dna1, shared_ptr<DNA> *dna2, shared_ptr<DNA> *mid);

private:
	inline qual_score minmaxQscore(qual_score t);// , int lnCnt);
	void minmaxQscore(shared_ptr<DNA> t);// , int lnCnt);
	bool setupFastq_2(string, string, string);
	bool setupFastaQual2(string, string, string = "fasta file");
	shared_ptr<DNA> read_fastq_entry(istream & fna, qual_score &minQScore,
		int&,bool&,bool);
	shared_ptr<DNA> read_fastq_entry_fast(istream & fna, int&,bool&);
	void jmp_fastq(istream &, int&);
	bool read_fasta_entry(ifbufstream*fasta_is, ifbufstream*quality_is, shared_ptr<DNA> in, shared_ptr<DNA>, int&);
	//static bool getFastaQualLine(istream&fna,  string&);
	void maxminQualWarns_fq();
	int auto_fq_version();
	int auto_fq_version(qual_score minQScore, qual_score maxQScore=0);
	void resetLineCounts(){ lnCnt[0] = 0;	lnCnt[1] = 0;	lnCnt[2] = 0; }
	bool desync(int pos) { if ( abs(pairs_read[pos] - pairs_read[opos[pos]]) > 1 ) {return true; } return false; }
	void IO_Error(string x);
	//bar on file read progress
	void _measure(istream &);
	inline bool _drawbar(istream &);
	inline void _print(int cur, float prog);
	int _fileLength, _max, _last;
	
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

	mutex protect;
	int num_threads;
	int at_thread; 
#ifdef togRead
	vector <jobC> slots;
#else
	vector <job3> slots;
#endif
    
    /*shared_ptr<DNA> getDNA(bool &has_next);
    
    shared_ptr<DNA> getDNAFromRecord(int &pos);
    
    shared_ptr<DNA> getDNAFromRecord(FastxRecord &record, int &pos);
    
    shared_ptr<DNA> getDNAFromRecord(FastxRecord &record);
    
    void addQuality(shared_ptr<DNA> fasta, FastxRecord *quality);
    */
};






#ifdef _gzipread2
std::vector< char > readline(gzFile f);
#endif 

#endif