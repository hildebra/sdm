/* sdm: simple demultiplexer
Copyright (C) 2013  Falk Hildebrand
email: Falk.Hildebrand@gmail.com

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


#define _CRT_SECURE_NO_DEPRECATE

#include <mutex>
#include <thread>
#include "InputStream.h"
#include "FastxReader.h"

#include <iostream>
#include <vector>
#include <functional>
#include <chrono>
#include <string>

#ifdef _gzipread
#include "include/gzstream.h"
//#include <include/zstr.h>
#endif
string spaceX(uint k){
	string ret = "";
	for (uint i = 0; i < k; i++){
		ret += " ";
	}
	return ret;
}

int digitsInt(int x){
	int length = 1;
	while (x /= 10)
		length++;
	return length;
}
int digitsFlt(float x){ 
	std::stringstream s;
	s << x;
	return (int)s.str().length(); 
}
string intwithcommas(int value) {
	string numWithCommas = std::to_string((long long)value);
	int insertPosition = (int)numWithCommas.length() - 3;
	while (insertPosition > 0) {
		numWithCommas.insert(insertPosition, ",");
		insertPosition -= 3;
	}
	return (numWithCommas);
}

std::string itos(int number) {
	std::stringstream ss;
	ss << number;
	return ss.str();
}
std::string ftos(float number,int digits) {
	std::stringstream ss;
	ss << std::setprecision(digits)<< number;
	return ss.str();
}
bool isGZfile(const string fi) {
	string subst = fi.substr(fi.length() - 3);
	if (subst == ".gz") {
		return true;
	}
	return false;
}

std::istream& safeGetline(std::istream& is, std::string& t) {
	t.clear();
	//from http://stackoverflow.com/questions/6089231/getting-std-ifstream-to-handle-lf-cr-and-crlf
	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	std::istream::sentry se(is, true);
	std::streambuf* sb = is.rdbuf();

	for (;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return is;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return is;
		case EOF:
			// Also handle the case when the last line has no line ending
			if (t.empty())
				is.setstate(std::ios::eofbit);
			return is;
		default:
			t += (char)c;
		}
	}
}


/*
void rtrim(std::string& s) {
	//s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
	//return s;
	s.erase(std::find_if(s.rbegin(), s.rend(), [](unsigned char ch) {
		return !std::isspace(ch);
		}).base(), s.end());
}

*/
int parseInt(const char** p1) {//,int& pos){//, const char ** curPos) {
	//from http://stackoverflow.com/questions/5830868/c-stringstream-is-too-slow-how-to-speed-up
	//size_t nxtPos = input.find_first_of(' ', curPos);
	//const char *p = input.substr(curPos,nxtPos).c_str();
	//if (!*p || *p == '?')		return 0;
	//int s = 1;
	//p = (const char*)curPos;
	const char* p = *p1;
	while (*p == ' ') p++;

	int acc = 0;
	while (*p >= '0' && *p <= '9')
		acc = acc * 10 + *p++ - '0';


	*p1 = p;
	//curPos = (size_t) p;
	return acc;
}


//MOCAT
std::vector<std::string> header_string_split(const std::string str, const std::string sep) {
	std::vector<std::string> tokens;
	tokens.reserve(13);
	size_t start = 0;
	size_t pos = 0;
	while ((pos = str.find_first_of(sep, start)) != std::string::npos) {
		tokens.push_back(str.substr(start, pos - start));
		start = pos + 1;
	}
	if (start < str.length()) {
		tokens.push_back(str.substr(start));
	} else if (start == str.length()) {
		tokens.push_back("0");
	}
	return tokens;
}
void remove_paired_info(string &s, short RP) {
	size_t f1 = s.find(" ");
	if (f1 != string::npos) {
		s = s.substr(0, f1);
		f1 = string::npos;
	}
	size_t  tarPos = s.size() - 2;
	/*if (RP == 0) {
		f1 = s.rfind("/1");
		if (f1 == string::npos) { 
			f1 = s.rfind(".1");
		}
	}
	else if (RP == 1) {
		f1 = s.rfind("/2");
		if (f1 == string::npos) {
			f1 = s.rfind(".2");
		}
	}
	else {*/
	f1 = s.rfind("/1");
	if (f1 == string::npos ||  f1 != tarPos) { f1 = s.rfind("/2"); }
	if (f1 == string::npos || f1 != tarPos) { f1 = s.rfind(".1"); }
	if (f1 == string::npos || f1 != tarPos) { f1 = s.rfind(".2"); }
	//}
	if (f1 != string::npos && f1 == tarPos) {
		s = s.substr(0, f1);
	}
}
std::string header_stem(string& header) {
	const size_t slash = header.find('/');
	if (slash != std::string::npos) {
		return header.substr(0, slash);
	}
	
	std::vector<std::string> tokens = header_string_split(header, ": ");

	if (tokens.size() == 11) {
		//if (tokens[8] == "Y") seq.clear();
		return tokens[0] + ":" + tokens[3] + ":" + tokens[4]
			+ ":" + tokens[5] + ":" + tokens[6];// +"#" + tokens[10];
	} else if (tokens.size() == 6 || tokens.size() == 7) {
		return tokens[0] + ":" + tokens[2] + ":" + tokens[3]
			+ ":" + tokens[4] + ":" + tokens[5];// + "#0";
	} else {
		cerr << "fastq header format error\n"; exit(98);
	}
		//throw std::runtime_error("fastq header format error\n" + header); }
}
void reverseTS(string & Seq) {
	int qs = (int)Seq.length() - 1;
	string S2 = Seq.c_str();
	for (int i = qs; i >= 0; i--) {
		Seq[i] = DNA_trans[(int)S2[qs - i]];
	}
}
string reverseTS2(const string & Seq) {
	int qs = (int)Seq.length() - 1;
	string S2 = Seq.c_str();
	for (int i = qs; i >= 0; i--) {
		S2[i] = DNA_trans[(int)Seq[qs - i]];
	}
	return S2;
}
bool any_lowered(const string& is) {
	for (uint i = 0; i < is.length(); i++) {
		char c = is[i];
		if (islower(c)) { return true; }
	}
	return false;
}
std::ptrdiff_t len_common_prefix_base(char const a[], char const b[])
{
	if (std::strlen(b) > std::strlen(a)) {
		return std::distance(a, std::mismatch(a, a + std::strlen(a), b).first);
	}
	return std::distance(b, std::mismatch(b, b + std::strlen(b), a).first);
}

string applyFileIT(string x, int it,const string xtr){

	size_t pos = x.find_last_of(".");
	if (pos != string::npos && isGZfile(x)) {
		pos = x.find_last_of(".", pos-1);
	}
	if (it == 0) {
		if (pos == string::npos) {
			return x + xtr;
		}
		return x.substr(0, pos) + xtr + x.substr(pos);
	} 
	ostringstream ss;
	ss << it;
	if (pos == string::npos) {
		return x + "." + ss.str() + xtr ;
	}
	return x.substr(0,pos) + "."+ss.str() + xtr + x.substr(pos);
	
	
}
bool fileExists(const std::string& name, int i, bool extiffail) {
	if (name == "") { return true; }
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);		return true;
	} else {
		if (extiffail) {
			cerr << "ERROR: Could not find file " << name << endl;
			if (i >= 0) {
				cerr << "on mapping file line " << i << endl;
			}
			exit(92);
		}
		return false; 
	}
}

string detectSeqFmt(const string inF) {

	vector<string> tfasP = splitByCommas(inF, ';');
	for (size_t i = 0; i < tfasP.size(); i++) {
		vector<string> tfas = splitByCommas(tfasP[i]);
		string fileS = tfas[0];
		istream* fnax(NULL);
		string file_type = "test file";
		string tmp("");
		string ret = "";
		if (fileS != "") {
            // auto detect gzip
#ifdef _gzipread
            fnax = DBG_NEW zstr::ifstream(fileS.c_str(), ios::in);
#else
            fnax = DBG_NEW ifstream(fileS.c_str(), ios::in);
#endif
		    

			if (!*(fnax)) { cerr << "\nCouldn't find " << file_type << " file \"" << fileS << "\"!\n Aborting..\n"; exit(4); }
			//char Buffer[RDBUFFER];
			//fasta_istreams[0]->rdbuf()->pubsetbuf(Buffer, RDBUFFER);
			while (safeGetline(*fnax, tmp)) {
				if (tmp[0] == '>') {
					ret = "-i_fna"; break;
				}
				else if (tmp[0] == '@') {
					ret = "-i_fastq"; break;
				}
				else if (tmp.length() == 0) {//do nothing
					;
				}
				else {
					cerr << " Could not auto detect input format. First non-empty line of your file looked like:\n" << tmp << endl;
					exit(888);
				}
			}
	}
		
		delete fnax;
		if (ret == "") {
			cerr << "Empty input file detected:\n" << fileS << endl;
		} else {
			return ret;
		}

	}

	return "empty";
}



	//compares two DNA entries, decides which one has overall better stats
bool whoIsBetter(shared_ptr<DNA> d1, shared_ptr<DNA> d2, shared_ptr<DNA> dM, 
	shared_ptr<DNA> r1, shared_ptr<DNA> r2, shared_ptr<DNA> rM, 
	float& ever_best, bool forSeed) {

	//if (forSeed) {		return false;	}
	//check if two primers present
	if (d2 == nullptr) {//hard reason .. only for PacBio etc reads
		if (d1->has2PrimersDetected() && !r1->has2PrimersDetected()) { return true; }
		if (!d1->has2PrimersDetected() && r1->has2PrimersDetected()) { return false; }
	} else {
		if (d2->getRevPrimCut() && !r1->getFwdPrimCut() && !r2->getRevPrimCut()) {	return true;}
	}
	//check if at least 1 primers present
	if (d1->getFwdPrimCut() && !r1->getFwdPrimCut()) { return true; }//hard reason
	if (!d1->getFwdPrimCut() && r1->getFwdPrimCut()) { return false; }

	double d1pid(1.), refpid(1.);
	if (ever_best >=0) {
		d1pid = (double) d1->getTempFloat(); refpid = (double)r1->getTempFloat();
		if (d1pid > ever_best) {
			ever_best = (float) d1pid;
		}
		//everbest is likely 100.f (ref OTUs)
		if (d1pid < (refpid - 0.3f) || d1pid < (ever_best - 0.4f ) ) { return false; }
	}

	bool dMerge(true), rMerg(true);
	double curL = (double)d1->getMergeLength();
	if (curL < 0) {
		curL = (double) d1->length();
		//if (d2 != NULL) { curL += (double) d2->length(); }
		dMerge = false;
	}
	double refL = (double)r1->getMergeLength();
	if (refL < 0) {
		refL = (double) r1->length();
		//if (r2 != NULL) { refL += (double) r2->length(); }
		rMerg = false;
	}

	//first check if d1 has merged, but ref did not.. clearly go for d, hard filter
	//if (d1->getMergeLength() != -1 && r1->getMergeLength() == -1) { return true; }

	//hard check on length ratios.. too drastically small , don't use d1
	if ((curL) / (refL) < BestLengthRatio ) { return false; }
	//at least 90% length of "good" hit
	//if (r1->getMergeErrors() < 0) {//no merge, can look at read1 only
	//	if (d1->length() / r1->length() < RefLengthRatio) { return false; }
	//}


	float dmergErrSco = 0.f;	float refMergErrSco = 0.f;
	//only check further if both comparisons did merge
	if (r1->getMergeLength() >=0 && d1->getMergeLength() >=0) {
		if (d1->getMergeErrors() > 0) {
			dmergErrSco = (float)d1->getMergeErrors();// *log10((maxQErr - (float)d1->getMergeErrorsQual()));
			dmergErrSco += d1->getMergeErrorsQual()/30;
		}
		if (r1->getMergeErrors() > 0) {
			refMergErrSco = (float)r1->getMergeErrors() + r1->getMergeErrorsQual() / 30;// *log10((maxQErr - (float)r1->getMergeErrorsQual()));
		}
		//scale to 1
		float maxMerr = max(refMergErrSco, dmergErrSco);
		dmergErrSco /= maxMerr;	refMergErrSco /= maxMerr;
		dmergErrSco = 1.f - dmergErrSco; refMergErrSco = 1.f - refMergErrSco;
	}

	//choose merged DNA if possible; d2 no longer needed then
	shared_ptr<DNA> dx, rx;
	bool allowMergeGuide = false;
	if (allowMergeGuide && dM != nullptr) { dx = dM;	d2 = nullptr; }	else { dx = d1; }
	if (allowMergeGuide && rM != nullptr) { rx = rM; r2 = nullptr;}else { rx = r1; }
	float thScore =  dx->getAvgQual(); //*(d1pid / 100)* log((float)curL);
	float rScore =  rx->getAvgQual();// *(refpid / 100)* log((float)refL);//r1->length()
	//if (thScore > rScore) {
		//also check for stable lowest score
		// if (d1->minQual() > r1->minQual() - MinQualDiff) { return true; }
	//}
	float maxScore = max(thScore, rScore);
	thScore /= maxScore;	rScore /= maxScore;


	//checks if the new DNA has a better overall quality
	double dAcSc = dx->getAccumError();	double dLen = dx->mem_length();
	double tAcSc = rx->getAccumError();	double tLen = rx->mem_length();
	if (d2 != nullptr) { dAcSc += d2->getAccumError(); dLen += d2->mem_length(); }
	if (r2 != nullptr) { tAcSc += r2->getAccumError(); tLen += r2->mem_length(); }
	dAcSc /= dLen;	tAcSc /= tLen;
	
	//calculate ratios to compare more easily among metrics
	double ratAccErr =  log(dAcSc)/log(tAcSc); //smaller better, log changes terms around
	double ratLength = double(curL) / (double)refL; //higher better
	double ratId = d1pid / refpid * 10 - 9.; //higher better //10-fold weighting
	if ((ratAccErr * ratLength * ratId) > 1.) {
		return true;
	}
	//normalize and invert
	//normalize to gene length, and invert to convert to positive score system
	/*maxScore = max(dAcSc, tAcSc);
	dAcSc /= maxScore;	tAcSc /= maxScore;
	dAcSc = 1.f - dAcSc; tAcSc = 1.f - tAcSc;
	//norm to 1, to compare to other terms
	double maxEr = max(dAcSc, tAcSc); 
	dAcSc /= maxEr; tAcSc /= maxEr;
	//dmergErrSco refMergErrSco dAcSc + tAcSc +   thScore  rScore
	if (  (dAcSc) *(d1pid / 100) * log((float)curL)
		> (tAcSc)  *(refpid / 100) * log((float)refL )) {
		
		if (dx->minQual() > rx->minQual() - MinQualDiff) { return true; }
//		return true;

	}
	*/

	return false;
}

void DNA::fixQ0(void) {
	return;//still depends on sequencer what is actually returned..
	for (uint i = 0; i < qual_.size(); i++) {
		if (qual_[i] <= 0) {
			sequence_.replace(i, 1, "N");
		}
	}
}
bool DNA::seal(bool isFasta) {//DN = sequence_.c_str();
	size_t QSi = qual_.size();
	if (QSi == 0 && isFasta) {
		return true;//nothing to be done, just empty DNA
	} else if (QSi == 0 && sequence_ == "" ) {
		this->setPassed(false);
		return true;
	} else if (QSi != sequence_.length()) {
		cerr << "Unequal length of seq and quality for name " << this->getId() << "\n";
		this->setPassed(false);
		return false;
	}
	//uppercase DNA
	std::transform(sequence_.begin(), sequence_.end(), sequence_.begin(), ::toupper);
    sequence_length_ = sequence_.length();
	this->fixQ0();
	return true;
}
DNA::DNA(vector<string> fas):DNA() {
	// Abort if input stream is at the end
	if (fas[0].length() == 0 ) {
		return;
	}
	if (fas[0][0] != '>') {
		cerr << "ERROR: Line 1 in fasta file does not start with \">\" :\n"<<fas[0]<<endl;
		exit(23);
	}
	// Storage objects
	string tqual;
	this->setHeader(fas[0].substr(1));
	this->setSequence( fas[1]);
	uint lsize = (uint) fas[1].length();
	vector<qual_score> Iqual(lsize, 11);
	if (fas[2].length() > 0) {
		tqual = fas[2];
		rtrim(tqual);
		const char* lQ = tqual.c_str();
		uint ii(0);
		qual_score nn(0);

		for (; ii < lsize; ii++) {
			nn = (qual_score)parseInt(&lQ);// , posStr);
			Iqual[ii] = nn;
			if (*lQ == '\0') {
				break;
			}
			//issQ >>  Iqual[ii];
		}

		if (Iqual.size() != fas[1].length()) {
			cerr << "Unequal fasta (" << fas[1].length() << ")/qual (" << fas[2].length() << ") length:\n";
			cerr << fas[0] << endl << fas[1] << endl << fas[2] << endl;
			exit(923);
		}
	}
	this->setQual(move(Iqual));
	
}

DNA::DNA(vector<string> fq, qual_score fastQver):DNA(){
	//string line;
	if (fq.size() != 4  || fq[0].length() == 0) { delself(); return; }
	while (fq[0][0] != '@') {
		cerr<< "ERROR on line " + fq[0] + ": Could not find \'@\' when expected (file likely corrupt, trying to recover):\n" ;// << endl;
		delself();
		return;

	}
	this->setHeader(fq[0].substr(1));
	sequence_ = fq[1];
	sequence_length_ = sequence_.size();
	if (fq[2][0] != '+') {
		cerr<<"Error input line " + fq[2] + ": Could not find \'+\' when expected (file likely corrupt, aborting):\n" + fq[0];// << endl;
		delself();
		return;
		//if (!safeGetline(fna, line)) { delete tdn; return NULL; }
	}

	//qual_ score
	string& line = fq[3];


	if (line.length()  != this->mem_length()) {
		string sizdif = "More";
		if (line.length() > this->mem_length()) {
			sizdif = "Less";
		}
		//check that quality gets not more length than DNA
		delself();
		cerr << "Error input line " + fq[3] + "'\n" + fq[1] + "'\n: "+ sizdif+" quality positions than nucleotides detected for sequence\n " +
			fq[0];// << endl;
		return;
	}

	uint qcnt(0); uint lline = (uint)this->mem_length();// line.length();
	vector<qual_score> Iqual(this->mem_length(), 0);
	
	for (; qcnt < lline; qcnt++) {
		qual_score q = (qual_score)line[qcnt] - fastQver;
		Iqual[qcnt] = q;// minmaxQscore();
	}
	this->setQual(move(Iqual));
}

/*
DNA::DNA(FastxRecord* FR, qual_score & minQScore, qual_score & maxQScore, qual_score fastQver) :DNA(){
    id_ = FR->header.substr(1);
    sequence_=FR->sequence;
	

	qual_.resize(FR->quality.length(), 0);
	uint qcnt(0); uint lline = (uint)FR->quality.length();
	for (; qcnt < lline; qcnt++) {
		qual_score t((qual_score)FR->quality[qcnt] - fastQver);
		if (minQScore > t) {
			minQScore = t;
			if (minQScore < 0) {
			}
		}
		else if (maxQScore < t) {
			maxQScore = t;
		}


        qual_[qcnt] = t;

	}
}
*/

string DNA::getPositionFreeId() { // remove /1 /2 #1:0 etc
	string s = this->getShortId();
	remove_paired_info(s, read_position_);
	return s;
}

int DNA::numNonCanonicalDNA(bool all) {
	int DNAch = 0;
	size_t DNAl = length();
	if (all) {
		DNAl=sequence_.size();
	}
	for (size_t i = 0; i < DNAl; i++) {
		DNAch += DNA_amb[sequence_[i]];
	}
	return DNAch;
}
int DNA::numACGT(){
	int DNAch = 0;
	uint maxL = length();
	for (unsigned int i = 0; i < length(); i++){
		DNAch += DNA_amb[(int) sequence_[i]];
	}
	return DNAch;
}
void DNA::stripLeadEndN() {
	int i = 0;
	while (DNA_amb[(int)sequence_[i]]) {
		i++;
	}
	if (i) {
		this->cutSeq(0, i);
	}
//cut at end N's
	i = (int)sequence_.length();
	while (DNA_amb[(int)sequence_[i]]) {
		i--;
	}
	if (i != (int)sequence_.length()) {
		this->cutSeq(i+1, (int)sequence_.length());
	}

}
void DNA::appendQuality(const vector<qual_score> &q)
{
	qual_.insert(qual_.end(), q.begin(), q.end());
	avg_qual_ = -1.f;
}

float DNA::getAvgQual(){
	if (avg_qual_ < 0.f){
		if (quality_sum_ == 0){
			for (unsigned int i = 0; i < this->length(); i++){
                quality_sum_ += qual_[i];
			}
		}
        avg_qual_ = static_cast<float> (quality_sum_) / static_cast<float> (this->length());

	}
	return avg_qual_;
}
int DNA::getMedianQual() {
	/*vector<qual_score> meq = qual_;
	sort(meq.begin(), meq.end());
	return calc_median2(meq, 0.5);*/
	return median(qual_);
}

/*float DNA::qualWinfloat_hybr(int W, float T, int W2, float T2, int& reason){//not used
	//if (T==0.f){return true;}
	int AQS=0, AQL;
	int TotQ = 0;
	int upTs = static_cast<int>(T * W);
	int upTl = static_cast<int>(T2 * W2);
	int QS = int (qual_.size());
	if (W>=QS){W = QS;} // too short
	int smallerW = W, largerW = W2;

	bool W1IsSmall = true;
	if (smallerW > W2){ largerW=W; smallerW = W2; W1IsSmall=false;
	std::swap(upTl,upTs); }

	for (unsigned int i=0; i<(unsigned int) smallerW; i++){
	AQS += qual_[i];
	}
	AQL = AQS;

	//hybrid schleife
	for (unsigned int i=smallerW; i<(unsigned int) largerW; i++){
	AQL += qual_[i];
	AQS += qual_[i]; AQS -=  qual_[i-smallerW];
	if (AQS < upTs){
	if (W1IsSmall){		reason=0;	return 0.f;
	} else {reason=1; this->cutSeq(i/2,this->length()); return float(AQL)/float(QS);}
	}
	}

	TotQ = AQL;
	for (unsigned int i=W; i<(unsigned int) QS; i++){
	AQS += qual_[i]; AQL += qual_[i];
	TotQ += qual_[i];
	AQS -= qual_[i-W];AQL -= qual_[i-W];
	if (AQS < upTs || AQL < upTl){
	if (W1IsSmall){		reason=0;	return 0.f;
	} else {reason=1; this->cutSeq(i/2,this->length()); return float(AQL)/float(QS);}
	}
	}
	//if (averageQ > static_cast<float> (TotQ) /static_cast<float> ( qual_.size())){return false;}
	return static_cast<float> (TotQ) /static_cast<float> ( QS);
	}
	*/

// modified from https://github.com/fpusan/moira/blob/master/moira/bernoullimodule.c
float DNA::interpolate(int errors1, float prob1, int errors2, float prob2, float alpha)
{
	float result = errors1 + ((errors2 - errors1) * ((1 - alpha) - prob1) / (prob2 - prob1));
	if (result < 0) //Happens only for very-short high qual_ sequences in which the probability of having 0 errors is higher than 1 - alpha.
	{
		result = 0;
	}
	return result;
}
float DNA::prob_j_errors(float p, float j, float n) //Where p is the error probability, j is the number of errors and n the number of observations.
{
	float per_position_accum_probs;
	if (j > n)	{
		return 0.0f; //The formula below would also return 0.
	}
	per_position_accum_probs = pow((1 - p), n);	//For j == 0.
	float i(1);
	for (; i <= j; i += 1.f)	{//For j > 0.
		per_position_accum_probs = ((n - i + 1.f) / (1.0f*i)) * (p / (1.f - p)) * per_position_accum_probs;
	}
	return per_position_accum_probs;
	
}
float DNA::sum_of_binomials(const float j, int k, float n, int qual_length, const vector<float>& error_probs, const vector< vector<float>> & per_position_accum_probs)
//#Where j is the number of errors and k is the position in the sequence.
{
	float probability = 0;
	float i(0); int k1 = (int) k - 1;

	for (; i <= j; i+=1.f)
	{
		probability += DNA::prob_j_errors(error_probs[k], i, n) *  per_position_accum_probs[int(j - i)][k1];
		//Where error_probs[k] is the error probability of the k-th position.
		//Where per_position_accum_probs[j-i][k-1] is the probability that all the bases from position 0 to k-1 had a total of j-i errors.
	}

	return probability;
}

float DNA::binomialFilter(int maxErr, float alpha){

	if (alpha == -1.f|| this->length()<3){ return 0.f; }//deactivated
	if (maxErr > 1e5) {return 0.f;	} // impossibly large..

	///Initialize some variables.
	
	int SeqLengthI = (int)sequence_length_;
	int n = 1; //Since we have a Bernoulli distribution.
	float n_f = (float)n;
	float alpha1 = 1.f - alpha;

	///Translate quality scores into error probabilities.
	vector<float> error_probs (sequence_length_, 1.f);
	for (size_t i = 0; i < sequence_length_; i++){
		if (qual_[i]<0 || qual_[i]>maxSAqualP) { continue; }
		error_probs[i] = (float)SAqualP[qual_[i]];//pow(10, (contig_quals[i] / -10.0)); //Since we want a continuous list of non-N error_probs.
	}

	///Actually get the job done.
	int max_expected_errors = maxErr + 3;// (int)sequence_length_ + 1;
	int expected_errors = 0;
	float probability;
	vector<float> accumulated_probs (max_expected_errors,0.f);
	//int j;
	int k;
	vector<float> empty(sequence_length_, 0.f);
	vector< vector<float>> per_position_accum_probs(max_expected_errors, empty);

	while (1)
	{
		float expected_errors_f = (float)expected_errors;
		//vector<float > per_position_accum_probs(sequence_length_, 0.f);
		for (k = 0; k < (int)sequence_length_; k++) {
			if (k == 0)	{
				per_position_accum_probs[expected_errors][k] = DNA::prob_j_errors(error_probs[k], expected_errors_f, n_f);
			} else {
				
				per_position_accum_probs[expected_errors][k] = DNA::sum_of_binomials((float)expected_errors, k, n_f, SeqLengthI, error_probs, per_position_accum_probs);
			}
		}
		probability = per_position_accum_probs[expected_errors][SeqLengthI - 1];

		if (expected_errors == 0){
			accumulated_probs[expected_errors] = probability;
		}else{
			accumulated_probs[expected_errors] = accumulated_probs[expected_errors - 1] + probability;
		}

		if (accumulated_probs[expected_errors] > (alpha1) || expected_errors >= (max_expected_errors-1)){
			break;
		}else{
			expected_errors++;
		}
	}
	if (expected_errors == 0){
		return 0.f;
	}
	float EXE = interpolate(expected_errors - 1, accumulated_probs[expected_errors - 1], expected_errors, accumulated_probs[expected_errors], alpha);
	return EXE;
}

float DNA::qualWinfloat(uint W, float T, int& reason){
	//if (T==0.f){return true;}
	int AQS=0;
	int TotQ = 0;
	int upTs = int(T) * int(W);
	uint QS = this->length();//static_cast<unsigned int> (qual_.size());
	if (W>=QS){W = QS;} // too short

//1st loop to ini window
	for (size_t i = 0; i < size_t (W); i++) {
		AQS += qual_[i];
	}
	TotQ = AQS;

	for (size_t i=W; i<(size_t) QS; i++){
		AQS += qual_[i] - qual_[i - W];
		TotQ += qual_[i];
		if (AQS < upTs ){
			reason=1;	return 0.f;
		}
	}
	//if (averageQ > static_cast<float> (TotQ) /static_cast<float> ( qual_.size())){return false;}
	return float(TotQ) /float( QS);
}

int DNA::qualAccumulate(double d){
	unsigned int i(0);double accErr(0.0);
	for (; i<this->length(); i++) {
		accErr+=SAqualP[qual_[i]];
		if (accErr >= d){break;}
	}

	this->accumulated_error_ = accErr;

	return i;
}

/*std::mutex qsc_mutex;
std::mutex ntcnt_mutex;
std::vector<std::mutex> qsc_mutexes;
std::vector<std::mutex> ntcnt_mutexes;*/

void DNA::ntSpecQualScores(vector<long>& qsc, vector<long>& ntcnt) {
	size_t sql = sequence_.length();

	if (qsc.size() < sql) {
		qsc.resize(sql,0);
	}
	if (ntcnt.size() < sql) {
		ntcnt.resize(sql,0);
	}
	for (uint i = 0; i < sql; i++ ) {
		short p = NT_POS[(int) sequence_[i]];
		qsc[p] += qual_[i];
		ntcnt[p]++;
	}
}

void DNA::ntSpecQualScoresMT(vector<long>& qsc, vector<long>& ntcnt) {
    size_t sql = sequence_.length();

    if (qsc.size() < sql) {
        qsc.resize(sql,0);
    }
    if (ntcnt.size() < sql) {
        ntcnt.resize(sql,0);
    }
    for (uint i = 0; i < sql; i++ ) {
        short p = NT_POS[(int) sequence_[i]];
        qsc[p] += qual_[i];
        ntcnt[p]++;
    }
}

bool DNA::qualAccumTrim(double d){
	if (d == -1.) {
		return true;
	}
	unsigned int i(qualAccumulate(d));
	if (i != this->length()){
		//cut 3' end
		this->cutSeqPseudo(i);
		this->QualCtrl.AccErrTrimmed = true;
		return false;
	}
	//did not cut this sequence:
	return true;
}
bool DNA::qualWinPos(unsigned int W, float T){
	if (T==0.f){return true;}
	int AQ=0;
	int unT = static_cast<int>((float)W*T);
	unsigned int QS = this->length();// (unsigned int)qual_.size();
	unsigned int QSh = QS >> 2;
	QSh = max(QSh,W);
	if (W>=QS){return true;} // too short

	vector<float> WQ((int) W);
	//TODO: check that the right num of nucs is taken..
	int cnt=0;
	if (QS < W) {return true;}
	for (unsigned int i=QS-1; i> QS-(unsigned int) W-1; i--){
		AQ += qual_[i];
		cnt++;
	}
	if (AQ > unT){return true;}
	int curW = QS-(unsigned int) W;
	for (uint i=QS-(unsigned int) W-1; i > QSh; i--){
		AQ += qual_[i];
		AQ -= qual_[i + W];

		if (AQ < unT){ //min Window qual_ was broken.. kick seq
			curW = i;
		} else {
			break;
		}
	}
	
	//partial seq  removal
	int pos = curW - (W>>1);
	if (pos < 0) {pos = 0;}
	this->cutSeqPseudo(pos);
	this->QualCtrl.QWinTrimmed = true;
	return false;
}


shared_ptr<DNA> DNA::getDNAsubseq(int start, int end, string id) {
	shared_ptr<DNA> rdn = make_shared<DNA>();
	rdn->setSequence(this->getSequence().substr(start,(end-start)) );
	vector<qual_score> tq = this->getQual(start, end) ;
	rdn->setQual(tq);
	rdn->setHeader(id);
	if (this->getBCnumber() >= 0) {
		rdn->setBCnumber(this->getBCnumber(),0);
		if (getBarcodeCut()) {
			rdn->setBarcodeCut();
		}
		rdn->setBarcodeDetected(getBarcodeDetected());

	}
	return rdn;
}

//removes part of seq and qual_ indexes between start & stop
//set stop to -2 if cut till end in pseudo == true mode
bool DNA::cutSeq(int start, int stop, bool pseudo){

	 if (stop == -1) {
		if (start >= (int) sequence_length_ || start < 0) { return false; }
	} else if (start >= stop ||  start >= (int)qual_.size()) {
		return false;
	}
	if (stop > (int)qual_.size()) {
		stop = (int)qual_.size();
	}
	
	//pseudo deactivates cutting of 3'
	if (pseudo) {
		if (stop == -1 ){//|| stop <= (int) sequence_length_) {
            sequence_length_ = start;
			return true;
		}
	}

	// this is to account for any already maybe incurred pseudo cuttings
	int seqLDiff = int(sequence_.length() - sequence_length_);

	string se = sequence_;
	if (stop == -1) {
		stop = (int)sequence_.length();
	}
	if (start == 0) {
        sequence_ = se.substr(stop);
	} else {
        sequence_ = se.substr(0, start) + se.substr(stop);
	}

	//DN = sequence_.c_str();
	//Quali
	qual_.erase(qual_.begin() + start, qual_.begin() + stop);

    sequence_length_ = sequence_.length() - seqLDiff;

	return true;
}

int DNA::matchSeq(std::string PrSt,int Err,int searchSpace, int startPos,bool exhaustive){
	//const char* DN = sequence_.c_str();
	//const char* Pr = PrSt.c_str();
	int PrL = (int) PrSt.length();
	int mthL = this->length() - PrL;
	//int wantSc = PrL - Err;
	int endPos(-1),pos(startPos), Prp(0), c_err(0),Prp2(0), c_points(0), point_aim(1);
	point_aim = max(1, int(PrL / 2));
	//bool res(false);
	vector<int> potentialMatches(0);
	for (; pos< searchSpace; pos++){
		if (pos > mthL) {	break;	}
		c_err = 0; Prp = 0; Prp2 = pos; c_points = 0;
		do {
		
#ifdef _NEWMATCH
			//new vector based matching
			c_err += DNA_IUPAC[sequence_[Prp2] + 256 * PrSt[Prp]]; 
			if (c_err > Err){break;}
#else
			//old, direct match
			if (!matchDNA(sequence_[Prp2],PrSt[Prp])){c_err++;if (c_err > Err){break;}}
#endif
			Prp++; Prp2++;
			if (sequence_[Prp2] != 'N') { c_points++; }
		} while ( Prp < PrL);
		/*if (c_points >= point_aim) {
			int x = 0; //DEBUG
		}*/
		if (c_err<=Err && c_points >= point_aim){
			endPos=pos;
			potentialMatches.push_back(endPos);
			if (!exhaustive) {
				break;
			}
			//break;
		}
	}
	//if(!suc){pos=-1;}
	if (potentialMatches.size() > 0) {
		return potentialMatches.back();
	} else {
		return -1;
	}
}
void DNA::reset() {
    accumulated_error_ = 0.; good_quality_ = false; mid_quality_ = false;
	//FtsDetected.reset();
	avg_qual_ = -1.f; quality_sum_ = 0; tempFloat = 0.f;
    qual_traf_ = "";

	this->resetTruncation();
}

void DNA::reverse_compliment(bool reset) {
	reverseTS(sequence_);
	std::reverse(qual_.begin(), qual_.end());
    reversed_ = !reversed_;
	qual_traf_ = "";
	if (!reset) {
		return;
	}
    accumulated_error_ = 0.; good_quality_ = false; mid_quality_ = false;
    avg_qual_ = -1.f; quality_sum_ = 0; tempFloat = 0.f;
	
}

//match from end of sequence_ to find rev primer
int DNA::matchSeqRev(const string& PrSt,int Err, int searchSpace,
				  int start,bool exhaustive){
	//fail::ret -1
	int PrL = (int) PrSt.length();
	if (start ==0){ 
		start = PrL; //default seed set to 5
	}
	int SeL = (int) sequence_.size();
	//int wantSc = PrL - Err;
	int pos(SeL- start), Prp(0), c_err(0), endPos(-1), c_points(0), point_aim(1);
	vector<int> potentialMatches(0);
	for (; pos> searchSpace; pos--){
		c_err = 0; Prp = 0; c_points = 0;
		int PrL2 = min(PrL,SeL-pos);
		point_aim = max(1,int((float)PrL2 *0.9f));
		do {
#ifdef _NEWMATCH
			uint lpos = pos + Prp;
			c_err += DNA_IUPAC[sequence_[lpos] + 256 * PrSt[Prp]]; if (c_err > Err){break;}
#else
			if(!matchDNA(sequence_[pos+Prp],PrSt[Prp])){c_err++;if (c_err > Err){break;}	}
#endif
			Prp++;
			if (sequence_[lpos] != 'N') { c_points++; }
		} while (Prp < PrL2);
		if (c_err<=Err && c_points >=point_aim){
			endPos=pos;
			potentialMatches.push_back(endPos);
			if (!exhaustive) {
				break;
			}
		}
	}
	//secondary check for last few NT's
	if (endPos==-2 && potentialMatches.size()==0){
		pos = (SeL-1);
		for (; pos> (SeL- start); pos--){
			c_err=0;Prp=0;
			int PrL2 = min(PrL,SeL-pos);
			do {
#ifdef _NEWMATCH
				c_err += DNA_IUPAC[sequence_[pos + Prp] + 256 * PrSt[Prp]]; if (c_err > Err){break;}
#else
				if(!matchDNA(sequence_[pos+Prp],PrSt[Prp])){c_err++;	if (c_err > Err){break;}}
#endif
				Prp++;
			} while (Prp < PrL2);
			if (c_err<=Err ){
				endPos=pos;
				potentialMatches.push_back(endPos);
				if (!exhaustive) {
					break;
				}
			}
		}
	}
	if (potentialMatches.size() > 0) {
		return potentialMatches.back();
	} else {
		return -1;
	}
}
// looks through total DNA seq
int DNA::matchSeq_tot(const string& Pr, int error, int maxPos, int& c_err){
	int PrL = (int) Pr.length();
	int pos(0), Prp(0), Prp2(0);
	bool success(false);
	for (pos=0; pos < maxPos; pos++){
		c_err=0;Prp=0;Prp2=pos;
		do {
			if(sequence_[Prp2] != Pr[Prp]){
				c_err++;
			}
			Prp++; Prp2++;
		} while (c_err <= error && Prp < PrL);
		if (c_err <= error){ success=true;break;}
	}
	if(!success){ pos=-1;}
	return pos;
}


bool DNA::matchDNA(char t1,char t2){
	if (t1==t2){
		return true;
	}
	switch (t2){
		case 'N': return true;
		case 'R': if (t1=='A' || t1=='G' ) {return true;}break;
		case 'Y': if (t1=='T' || t1=='C' ) {return true;}break;
		case 'M': if (t1=='C' || t1=='A' ) {return true;}break;
		case 'K': if (t1=='T' || t1=='G' ) {return true;}break;
		case 'W': if (t1=='T' || t1=='A' ) {return true;}break;
		case 'S': if (t1=='C' || t1=='G' ) {return true;}break;
		case 'B': if (t1!='A') {return true;}break;
		case 'D': if (t1!='C' ) {return true;}break;
		case 'H': if (t1!='G' ) {return true;}break;
		case 'V': if (t1!='T' ) {return true;}break;
	}
	return false;
}
bool DNA::HomoNTRuns(int max){
	char lastC = sequence_[0];
	int rowC=1;
	for (unsigned int i=1;i<length();i++){
		if (sequence_[i] == lastC){
			rowC++;
			if (rowC>= max){
				return false;
			}
		} else {
			rowC = 1;
			lastC = sequence_[i];
		}
	}
	return true;
}

int DNA::HomoNTTrim(int max) {
	// Return the sequence length if trimming takes place, else return 0.
	unsigned int i = length() - 1;
	char lastC = sequence_[i];
	int rowC = 0;
	for (;i > 0; i--){
		if (sequence_[i] == lastC) {
			rowC++;
		} else {
			break;
		}
	}

	if (rowC >= max){
		return i + 1;
	} 
	rowC = 0;
	lastC = sequence_[i - 1];
	i = 0;
	return 0;
}

/*
void DNA::writeSeq(ofstream& wr){
	int cnt=0;
	if (sequence_.size()==0){return;}
	wr<<">"<<new_id_<<std::endl;
	for (unsigned int i=0;i<sequence_.size();i++){
		cnt++;
		if (cnt<80){wr<<sequence_[i];
		} else {wr<<sequence_[i] << endl;	cnt=0;
		}
	}
	if (cnt>0){
		wr<<endl;
	}
}
*/
size_t DNA::getSpaceHeadPos(const string & x) {
	size_t pos = x.find_first_of(" \t"); // x.find(' ');
	if (pos == string::npos) {pos = x.length(); }
	return pos;
}
size_t DNA::getShorterHeadPos(const string & x, int fastQheadVer ) {

	size_t pos(string::npos);
	size_t strL = x.length();
	if (fastQheadVer != 0) {
		if (read_position_ == 1) {
			pos = x.find("/2", strL-3);
			//			if (pos == string::npos) {	pos = x.find_first_of(" 1:");}
				//		if (pos == string::npos) { pos = x.find_first_of("/1"); }
		}
		else if (read_position_ == 0) {
			pos = x.find("/1", strL - 3);
		}
		else {
			pos = x.find("/1", strL - 3);
			if (pos == string::npos) { pos = x.find("/2", strL - 3); }
		}
	}
	//if (pos == string::npos){pos=x.length()-min((size_t)5,x.length());}}
	//if(pos<0){pos=0;}

	if (pos == string::npos) {
		pos = this->getSpaceHeadPos(x);
	}
		
		//pos = x.find_first_of(" \t"); }//, x.find('\t');
		//	if (pos == string::npos) { pos = x.length(); }
	return pos;
}
string DNA::xtraHdStr(){
	string xtr = " ";
	if (FtsDetected.forward){ xtr += "F"; }
	if (FtsDetected.reverse){ xtr += "R"; }
	return xtr;
}
void DNA::writeSeq(ostream& wr, bool singleLine) {
	wr << writeSeq(singleLine);
}
void DNA::writeSeq(string& str, bool singleLine ){
	if (sequence_.size() == 0) { return ; }
	str = ">" + new_id_ + "\n";
	if (singleLine) {
		str += sequence_.substr(0, length()) + "\n";
		return ;
	}
	unsigned int SeqS = length();
	unsigned int leftOver = SeqS % 80;
	int last = 0;
	for (unsigned int i = 0; i < (SeqS - leftOver) / 80; i++) {
		str += sequence_.substr(last, 80) + "\n"; last += 80;
	}
	if (leftOver > 0) {
		str += sequence_.substr(last) + "\n";
	}
}
string DNA::writeSeq(bool singleLine ){
	string str = string("");
	this->writeSeq(str, singleLine);
	return str;
	
}
/**/
void DNA::writeQual(ostream& wr, bool singleLine) {
	wr << writeQual(singleLine);
}
string DNA::writeQual(bool singleLine) {
	string str = string("");
	this->writeQual(str, singleLine);
	return str;
}
void DNA::writeQual(string& str, bool singleLine) {
	int cnt=0;
	if (qual_.size() == 0){return ;}
	str += ">" + new_id_ + "\n";
	//wr<< qual_traf_<<endl;
	string endlChr("\n");
	if (singleLine) { endlChr = " "; }
	for (unsigned int i=0;i<length();i++){
		cnt++;
		if (cnt<80){
			str += itos((int)qual_[i]) + " ";
		} else {
			str += itos((int)qual_[i]) + endlChr;
			cnt=0;
		}
	}
	if (cnt>0){
		str += "\n";
	}
}


void DNA::writeFastQ(ostream& wr, bool newHD) {//, int fastQver){
	wr << writeFastQ(newHD);
}
string DNA::writeFastQ(bool newHD) {
	string ret = string("");
	this->writeFastQ(ret, newHD);
	return ret;
}
void DNA::writeFastQ(string& ret, bool newHD) {
	if (qual_.size() == 0 || sequence_.length() == 0 || length() == 0){return ;}
	string xtr = xtraHdStr();
	if (newHD) {
		ret += "@" + new_id_ + "\n";
	} else {
		ret += "@" + id_ + "\n";
	}
	ret += sequence_.substr(0, this->length()) + "\n";
	ret += "+\n" ;//new_id_<<endl;
	ret += qual_traf_ + "\n";
}
/*
void DNA::writeFastQ(ofbufstream& wr, bool newHD) {//, int fastQver){
    if (!(qual_.size() == 0 || sequence_.length() == 0 || length() == 0)) {
        string xtr = xtraHdStr();
        if (FtsDetected.forward) { xtr += "F"; }
        if (FtsDetected.reverse) { xtr += "R"; }
        if (newHD) {
            wr << "@" + new_id_ + "\n";
        } else {
            wr << "@" + id_ + "\n";
        }
        wr << sequence_.substr(0, length()) + "\n";
        wr << "+\n";//new_id_<<endl;
        wr << qual_traf_ + "\n";
    }
}
*/
void DNA::writeFastQEmpty(ostream& wr) {

	wr << "@" << new_id_ << endl;
	wr << "" << endl;
	wr << "+" << endl;//new_id_<<endl;
	wr << "" << endl;
}

void DNA::changeHeadPver(int ver){
	string& oID(id_);
	if (id_fixed_){
		oID = new_id_;
	}
	
	if (ver==1){ //change from XX 1: to XX/1
		size_t pos = oID.find_first_of(" ");
		//unsigned int pos2 = id_.find(" ",pos);
		if (pos != string::npos){
            new_id_ = oID.substr(0, pos) + "/" + oID.substr(pos + 1, 1);
			//if (pos2 != string::npos){new_id_ += id_.substr(pos2);}
		}
	} else if (ver==2){
		cerr<<"Head change Not implemented\n";exit(23);
	}
    id_fixed_=true;
}

/*void DNA::reverseTranscribe(){
	int qs = (int)qual_.size()-1;
	vector<int> Q2(qual_);
	for (int i=qs;i>=0;i--){
		qual_[i] = Q2[qs-i];
	}
	reverseTS(sequence_);
}
*/
bool DNA::sameHeadPosFree(shared_ptr<DNA> d) {
	if (d == NULL) { return false; }
	return (getPositionFreeId() == d->getPositionFreeId());

	//return sameHead(d->getShortId());
}

bool DNA::sameHead(const string& oID) {
	

	size_t pos = getShorterHeadPos(id_);
	if (oID.size() < pos) { return false; }
	if (id_.substr(0, pos) == oID.substr(0, pos)) {
		return true;
	}
	return false;
}



void DNA::setPassed(bool b){
	//once mid qual is active, there is no more good qual possible
	if (failed_) {
		good_quality_ = false; mid_quality_ = false;
	} else if (mid_quality_) {
		good_quality_ = false;
	}
	else {
		good_quality_ = b;
		failed_ = false;
	}
 /*   good_quality_ = b;
	if (good_quality_ && mid_quality_) {
        mid_quality_ = false;
	}
	*/
}


void DNA::prepareWrite(int ofastQver) {
	uint len = length();
	if (qual_traf_.size() == len) {
		return;
	}
	qual_traf_.resize(len);

	unsigned int i = 0;
	for (; i < len; i++) {
        qual_traf_[i] = char(qual_[i] + ofastQver);
	}
    qual_traf_[i] = '\0';
}
void DNA::resetQualOffset(int x, bool fqSol) {
	for (size_t i = 0; i < qual_.size(); i++) {
	    qual_[i] += x;
	}
	if (fqSol) {//quick hack, since q<13 scores are the only deviants and uninteresting in most cases..
		for (size_t i = 0; i < qual_.size(); i++) {
			if (qual_[i] < 0) {
                qual_[i] = 0;
			}
		}
	}
}



///////////////////////////////////////////////////////////////
//INPUT STREAMER

bool DNAunique::betterPreSeed(shared_ptr<DNA> d1, shared_ptr<DNA> d2) {
	//0.2% difference is still ok, but within 0.5% of the best found seed (prevent detoriating sequence match)
	//float blen = (float)ref->length() + (float)d1->length();
	shared_ptr<DNAunique> ref2 = this->getPair();
	//shared_ptr<DNAunique> ref = this;

	

	//uint bestL = this->getBestSeedLength();
	//if (float(curL) / float(bestL) < BestLengthRatio) { return false; }
	cerr << "should not be here DNAunique::betterPreSeed\n"; exit(1243);
	return false;
	//whoIsBetter(d1, d2, shared_from_this(), ref2, this->getBestSeedLength(), -1.f);



	/*float dmergErrSco = (float)d1->getMergeErrors() * (float)d1->getMergeErrorsQual();
	float refMergErrSco = (float)this->getMergeErrors() * (float)this->getMergeErrorsQual();
	if (dmergErrSco) {
		int x = 0;
	}
	*/

	//at least 90% length of "good" hit
	/*
	if (this->getMergeErrors() < 0) {//no merge, can look at read1 only
		if (d1->length() / this->length() < RefLengthRatio) { return false; }
	}
	
	//checks if the new DNA has a better overall quality
	double dAcSc = d1->getAccumError();
	double tAcSc = this->getAccumError();
	if (d2 != nullptr) { dAcSc += d2->getAccumError(); }
	if (ref2 != nullptr) { tAcSc += ref2->getAccumError(); }
	if (dAcSc < tAcSc) {
		return true;
	}
	*/
	
	/* // remove old betterDNA algo, switched to accum error
	//1 added to qual, in case no qual DNA is used
	float thScore = (1 + d1->getAvgQual()) * log((float)d1->mem_length());
	float rScore = (1 + this->getAvgQual()) * log((float)this->mem_length());
	if (thScore > rScore) {
		//also check for stable lowest score
		if (d1->minQual() > this->minQual() - MinQualDiff && (d2 == NULL || ref2 == NULL)) {
			if (curL > bestL) { this->setBestSeedLength(curL); }
			return true;
		}
	}
	*/



	/* whole second read is not useful..
	if (d2 == NULL || ref2 == NULL) {
		return false;
	}
	//if (d2->getRevPrimCut() && !ref2->getRevPrimCut()) { return true; }//hard reason
	//if (!d2->getRevPrimCut() && ref2->getRevPrimCut()) { return false; }//hard reason

	//at least 90% length of "good" hit
	if (d2->mem_length() / ref2->mem_length() < RefLengthRatio) { return false; }

	//checks if the new DNA has a better overall quality
	//weigh with average id_ to OTU seed
	thScore += (1 + d2->getAvgQual()) * log((float)d2->mem_length()) * 97;
	rScore += (1 + ref2->getAvgQual()) * log((float)ref2->mem_length()) * 97;
	if (thScore > rScore) {
		//update best seed length score
		if (curL > bestL) { this->setBestSeedLength(curL); }
		return true;
	}
	*/

	return false;
}

void DNAunique::matchedDNA(shared_ptr<DNA> dna, shared_ptr<DNA> dna2,
		shared_ptr<DNA> dnaM,
		int sample_id, bool b_derep_as_fasta_){
	dna->setDereplicated();
	DNAuniMTX.lock();
	// increment sample counter
	incrementSampleCounter(sample_id);
	if (!b_derep_as_fasta_) { // only improve qualities if we do not replace
		sumQualities(dna);
	}
	// only replace old unique dna with new if it has good quality and its seed is better
	// betterPreSeed makes sense with uparse/
	float tmp(-1.f);
	if (dna->isGreenQual() && 
		whoIsBetter(dna, dna2, dnaM, shared_from_this(), this->getPair(), this->getMerge(),
			 tmp, false)  ){ //this->getBestSeedLength(),
		//betterPreSeed(dna, dna2)) {
		// Uparse does NOT use qualities for clustering, therefore we do not need to calculate average qualities for this dna object
		// Lotus uses qualities for constructing taxonomy (therefore we do replace dna if theres a better quality read)
		takeOverDNA(dna, dna2, dnaM);
	}
	DNAuniMTX.unlock();
}

void DNAunique::incrementSampleCounter(int sample_id) {
	if (sample_id < 0) {
		return;
	}
	//count_++;
#ifdef _MAPDEREPLICATE
	auto sample_counter = occurence.find(sample_id);
	if (sample_counter == occurence.end()) {
		occurence[sample_id] = 1;
	} else {
		sample_counter->second++;
	}
#endif
}
void DNAunique::incrementSampleCounterBy(int sample_id, long count) {
	auto sample_counter = occurence.find(sample_id);
	if (sample_counter == occurence.end()) {
		occurence[sample_id] = count;
	} else {
        sample_counter->second += count;
	}
}
/*vector<pair_<int, int>> DNAunique::getDerepMapSort2(size_t wh ){
	typedef std::pair_<int, int> mypair;
	size_t siz = occurence.size();
	if (wh > siz) { wh = siz; }

	struct IntCmp {
		bool operator()(const mypair &lhs, const mypair &rhs) {
			return lhs.second > rhs.second;
		}
	};


	vector<mypair> myvec(occurence.begin(), occurence.end());
	std::partial_sort(myvec.begin(), myvec.begin() + wh, myvec.end(), IntCmp());

	return myvec;
}*/

bool sortDescending(int i, int j) { return (i>j); }//descending sort
vector<int>  DNAunique::getDerepMapSort(size_t wh) {
	vector<int> vals;
	size_t siz = occurence.size();
	if (wh > siz) {	wh = siz;}
	vals.reserve(siz);
	for (auto kv = occurence.begin(); kv != occurence.end(); kv++) {
		vals.push_back(kv->second);
	}
	//partial sort doesn't make sense, as I want to break border asap
	//partial_sort(vals.begin(), vals.begin() + wh, vals.end(), sortDescending);
	sort(vals.begin(), vals.end(), sortDescending);
	return vals;
}

bool DNAunique::pass_deprep_smplSpc(const vector<int>& cv) {
	//unordered_map<int, int> occ;
	//combined samples will not be considered
	//occ = occurence;
	for (read_occ::iterator iter = occurence.begin(); iter != occurence.end(); ++iter) {
		//int cnts = iter->second;
		int ref = cv[iter->first];
		if (ref != -1 && iter->second >= ref ) {
			return true;
		}
	}
	return false;

}


void DNAunique::transferOccurence(const shared_ptr<DNAunique> dna_unique) {

	if (occurence.size() == 0) {
		occurence = dna_unique->occurence;
        //count_ = dna_unique->count_;
	} else {
		//which sample contains this dna?
		read_occ dereplication_map = dna_unique->getDerepMap();

		for (auto sample_iterator = dereplication_map.begin(); sample_iterator != dereplication_map.end(); ++sample_iterator) {
            auto sample_counter = occurence.find(sample_iterator->first);
			if (sample_counter == occurence.end()) {
				occurence[sample_iterator->first] = sample_iterator->second;
			} else {
                sample_counter->second += sample_iterator->second;
			}
		}
		//size track
		//count_ += dna_unique->count_;
	}
}


void DNAunique::writeMap(ofstream & os, const string &hd, vector<int> &counts_per_sample, const vector<int>& combiID) {
	if (occurence.empty()) { return; }
	int total_count(0);
	read_occ sample_counters;

	if (combiID.size() > 0){//combine all counts on combined categories

		for (read_occ::iterator iter = occurence.begin(); iter != occurence.end(); ++iter) {
			//aim: sample_counters[combiID[iter->first]] += iter->second;
			read_occ::iterator sample_counter = sample_counters.find(combiID[iter->first]);
			if (sample_counter == sample_counters.end()){
			    // if sample_counter exists
                sample_counters[combiID[iter->first]] = iter->second;
			}else{
                sample_counter->second += iter->second;
			}
		}
	}
	else {
        sample_counters = occurence;
	}

	//prints combined sample counts
	os << hd;
	for (read_occ::iterator iter = sample_counters.begin(); iter != sample_counters.end(); ++iter) {
		int count = iter->second;
        total_count += count;
		os << "\t";
		os << iter->first << ":" << count;
	}

	//counts non-combined sample counts
	for (read_occ::iterator iter = occurence.begin(); iter != occurence.end(); ++iter) {
		//int cnts = iter->second;
		counts_per_sample[iter->first] += iter->second;
	}
	os << endl;
	
	/*if (total_count != count_) {
		cerr << "Unequal counts in Map(" << total_count << ") and HeadCnt(" << count_ << "):" << endl << this->getId() << endl;
		exit(82);
	}*/
}
void DNAunique::Count2Head(bool usFmt) {
	//string tmp = this->getShortId();
	if (!usFmt) {
        new_id_ += "_" + itos(totalSum());
	} else {
        new_id_ += ";size=" + itos(totalSum()) + ";";
	}
    id_fixed_ = true;
}

void DNAunique::takeOver(shared_ptr<DNAunique> const dna_unique_old) {
	this->saveMem();
	//this->setBestSeedLength(better_dna->getBestSeedLength());
	this->transferOccurence(dna_unique_old);
	if (dna_unique_old->getPair() != nullptr) {
		this->attachPair(make_shared<DNAunique>(dna_unique_old->getPair(), -1));
	}
	if (dna_unique_old->getMerge() != nullptr) {
		this->attachMerge(make_shared<DNAunique>(dna_unique_old->getMerge(), -1));
	}

	quality_sum_per_base_ = dna_unique_old->transferPerBaseQualitySum();
}
void DNAunique::takeOverDNA(shared_ptr<DNA> const better_dna, shared_ptr<DNA> const dna2, shared_ptr<DNA> const dnaMerge) {

	
	if (dna2 != nullptr) {
		this->attachPair(make_shared<DNAunique>(dna2, -1));
	}
	if (dnaMerge != nullptr) {
		this->attachMerge(make_shared<DNAunique>(dnaMerge, -1));
	}
	
	avg_qual_=better_dna->avg_qual_;
	accumulated_error_ = better_dna-> accumulated_error_;
	FtsDetected=better_dna->FtsDetected;
	good_quality_=better_dna->good_quality_;
	//new_id_ = better_dna->id_;
	id_ = better_dna->id_;
	id_fixed_=better_dna->id_fixed_;
	new_id_ = better_dna-> new_id_;
	mid_quality_=better_dna->mid_quality_;
	merge_offset_=better_dna->merge_offset_;
	merge_seed_pos_=better_dna->merge_seed_pos_;
	reversed_merge_ = better_dna-> reversed_merge_;
	seed_length_ = better_dna-> seed_length_;
	quality_sum_=better_dna->quality_sum_;
	QualCtrl=better_dna->QualCtrl;
	qual_=better_dna->qual_;
	qual_traf_=better_dna->qual_traf_;
	sequence_=better_dna->sequence_;
	sequence_length_ = better_dna-> sequence_length_;
	sample_id_ = better_dna-> sample_id_;
	reversed_ = better_dna-> reversed_;
	read_position_ = better_dna-> read_position_;
	good_quality_=better_dna->good_quality_;
	accumulated_error_ = better_dna-> accumulated_error_;
	tempFloat = better_dna-> tempFloat;
	
	
	this->saveMem();
}

uint64_t * DNAunique::transferPerBaseQualitySum() {
    uint64_t* tmp = quality_sum_per_base_;
    quality_sum_per_base_ = nullptr;
    return tmp;
}

void DNAunique::prepareDerepQualities(int ofastQver) {
	prepSumQuals();
    if (qualities_avg_.size() == length()) {
        return;
    }
    qualities_avg_.resize(length());

    unsigned int i = 0;
    for (; i < length(); i++) {
        int newbase = (int)round((double)quality_sum_per_base_[i]/(double)totalSum()) + ofastQver;
        qualities_avg_[i] = char(newbase);
    }
    qualities_avg_[i] = '\0';

}


void DNAunique::writeDerepFastQ(ofstream & wr, bool newHD) {
    if (!(qual_.size() == 0 || sequence_.length() == 0 || length() == 0)) {
        string xtr = xtraHdStr();
        if (FtsDetected.forward) { xtr += "F"; }
        if (FtsDetected.reverse) { xtr += "R"; }
        if (newHD) {
            wr << "@" + new_id_ + "\n";
        } else {
            wr << "@" + id_ + "\n";
        }
        wr << sequence_.substr(0, length()) + "\n";
        wr << "+\n";//new_id_<<endl;
        wr << qualities_avg_ + "\n";
    }
}



///////////////////////////////////////////////////////////////
//INPUT STREAMER

InputStreamer::~InputStreamer(){
	allStreamClose();
//	for (uint i = 0; i < dnaTemp1.size(); i++) { if (dnaTemp1[i] != NULL) { delete dnaTemp1[i]; } }
//	for (uint i = 0; i < dnaTemp2.size(); i++) { if (dnaTemp2[i] != NULL) { delete dnaTemp2[i]; } }
}
/*
bool InputStreamer::getFastaQualLine(istream& fna,  string& line) {
	
	if (!safeGetline(fna, line))
	    return false;
	while (line[0] == '$') { //$ marks comment
		safeGetline(fna, line);
	}
	return true;
}
*/
bool InputStreamer::read_fasta_entry(ifbufstream* fasta_is, ifbufstream* quality_is, shared_ptr<DNA> tdn1, shared_ptr<DNA> tdn2, int &cnt){
	// Abort if input stream is at the end
	if (fasta_is->eof())
	    return false;
	
	// Storage objects
    string line;
    string tseq;
	string tqual;
	string lineQ;
	
	// READ fastQ - FILE
	fasta_is->getline(line);
	//if (!getFastaQualLine(fasta_is, line) || (line.empty() && fasta_is.eof()))
	//    return false;
	
	if (!qualAbsent)
		quality_is->getline(lineQ);
	    //getFastaQualLine(quality_is, lineQ);
	cnt++;

	if (cnt == 1) {
	    // fasta header has to start with ">", otherwise exit
		if (line[0] != '>' && fasta_is) {
		    cerr << "ERROR: Line 1 in fasta file does not start with \">\" \n";
		    exit(23);
		}
		// init new DNA object in dnaTemp1 ( which is the local dnaTemp1, not the global)
        tdn1->setHeader(line.substr(1));
		
		if (!fasta_is->getline(line))
		    return false;
			
		if (!qualAbsent) {
			quality_is->getline(lineQ);
		    //getFastaQualLine(quality_is, lineQ);
		}
		// increase count by one
		cnt++;
	}
	
	//continous read in until ">" is hit
	while (line[0] != '>') {
		tseq += line;
		if (!fasta_is->getline(line))//getFastaQualLine(fasta_is, line))
		    break;
	}
	if (!qualAbsent) {
		while (lineQ[0] != '>') {
			tqual += " " + lineQ;
			if (!quality_is->getline(lineQ)){//getFastaQualLine(quality_is, lineQ)) { 
				break; }
		}
	}

	// Set sequence of new DNA object
    tdn1->setSequence(tseq);
	
	// length of seq
	size_t lsize = tseq.size();
	
	// store qualities in vector
	vector<qual_score> Iqual(lsize, 11);
	
	// If quality is there
	if (!qualAbsent) {
		rtrim(tqual);
		const char* lQ = tqual.c_str();
		uint ii(0);
		qual_score nn(0);
		
		for (; ii < lsize; ii++) {
			nn = (qual_score)parseInt(&lQ);// , posStr);
			Iqual[ii] = nn;
            if (*lQ == '\0') {
                break;
            }
			//issQ >>  Iqual[ii];
		}

	}
    tdn1->appendQuality(Iqual);


	//since already read in 1 more line, this line needs to be used on new object
	if (fasta_is ) {
		if (!qualAbsent && (line[0] != '>' || lineQ[0] != '>')) { cerr << "ERROR: Desynced fasta reader\n"; exit(23); }
        tdn2->setHeader(line.substr(1));
	} else {
		return false;
	}
	//a new DNA obj was set up, return to process tdni, tdno will be completed next call
	return true;
}

void InputStreamer::jmp_fastq(istream & fna, int &lnCnt) {
	string line;
	string tname = "", tseq = ""; //temporary storage
	size_t cnt = 0, qcnt = 0, DNAlength = 0;
	bool mode = true; //mode=T:fna,F:qual_
	bool needsAT = true, needsPlus = false; // checks if quality was completly read in and a '@' char is now expected in the next line
	//	getline(fna2,line2,'\n');

	while (getline(fna, line, '\n')) {
		lnCnt++;
		if (line == "") { return ; }
		if (line[0] == '@' && needsAT) { //fasta description
			if (cnt != 0) {
				cerr << "Line " << lnCnt << ": Could not find \'@\' when expected: on line \n" << line << endl;
			}
			//tname = line;
			needsAT = false;
			needsPlus = true;
			continue;
		} else if (line[0] == '+' && needsPlus) {
			needsPlus = false;
			mode = false;
			continue;
		} else if (needsAT) {
			cerr << "Line " << lnCnt << " :Could not find \'@\' symbol where expected on line \n " << line << endl;
			exit(6);
		}

		istringstream iss(line);
		iss >> tseq;
		if (mode) {
			DNAlength += tseq.length();
		} else {
			qcnt += tseq.length();
			if (DNAlength == qcnt) { return; }
		}

		cnt++;
	}
}
//reads out single seq + quality entry from fastq file
shared_ptr<DNA> InputStreamer::read_fastq_entry(istream & fna,  qual_score &minQScore, int& lnCnt,
			bool&corrupt,bool repairInStream) {
	string line;
	//string tseq = ""; //temporary storage
	uint qcnt = 0;
	bool mode = true; //mode=T:fna,F:qual_
	shared_ptr<DNA> tdn = make_shared<DNA>("", "");
	vector<qual_score> Iqual(0);
	bool needsAT = true, needsPlus = false; // checks if quality was completly read in and a '@' char is now expected in the next line

	if (!fna) { return NULL;  cerr << "Read Fastq: Input stream does not exist" << lnCnt << endl; exit(53); }


	while (safeGetline(fna, line)) {

		lnCnt++;
		while (repairInStream) {
			if (line[0] == '@')
			    repairInStream = false;
			else
			    safeGetline(fna, line);
		}
		//if (line.length() == 0) { return tdn; }
		if (needsAT) { //fasta description
			if (line[0] == '@') {
				if ((lnCnt-1) % 4 != 0) {
					fqPassedFQsdt = false;
				}
                tdn->setHeader(line.substr(1));
				mode = true; qcnt = 0;
				needsAT = false;
				needsPlus = true;
				continue;
			} else {
				corrupt = true;//try again,could be last empty line in file..
				if (line.length() != 0) {
					IO_Error("Line " + itos(lnCnt) + " :Could not find \'@\' symbol where expected on line \n " +line);// << endl;
				}
				return tdn;
			}
		} else if (needsPlus && line[0] == '+') {
			if ((lnCnt-3) % 4 != 0) {
				fqPassedFQsdt = false;
			}

			
			Iqual.resize(tdn->length(), 0);
			needsPlus = false;
			mode = false;
			continue;
		}

		//istringstream iss(line);
		//iss >> tseq;
		if (mode) {
			//fna
			if ((lnCnt-2) % 4 != 0) {
				fqPassedFQsdt = false;
			}

			//if (std::any_of(line.begin(), line.end(), [](char c) {return (islower(c)); })) {
			if (any_lowered(line)){
				fqPassedFQsdt = false;
				std::transform(line.begin(), line.end(), line.begin(), ::toupper);
			}
            tdn->appendSequence(line);
		} else if (!mode) {
			//qual_
			if ((lnCnt ) % 4 != 0) {
				fqPassedFQsdt = false;
			}
			bool printB(false);
			for (size_t i = 0; i < line.length(); i++) {
				//really 33?
				Iqual[qcnt] = minmaxQscore((qual_score)line[i] - fastQver, printB);
				qcnt++;
			}
			if (qcnt == tdn->length()) {
				//needsAT=true;
                tdn->appendQuality(Iqual);
				
			} else if (line.length() + qcnt > tdn->length()) {
				//check that quality gets not more length than DNA
				IO_Error("ERROR: More quality positions than nucleotides detected for sequence\n " + tdn->getId());
				//tdn->setPassed(false);
				corrupt = true;
				return tdn;
			}
			break;

		}
	}
	//if (tdn!= NULL && !tdn->control()){delete tdn; tdn=NULL;}
	corrupt = false;
	return tdn;

	//to check for fast fastq reader: 1. DNA in uppercase? 2. DNA/QUAL in single line?
}
void InputStreamer::IO_Error(string x) {
	fqverMTX.lock();
	cerr << x << endl;
	if (DieOnError) {
		exit(632);
	}
	ErrorLog.push_back(x);
	fqverMTX.unlock();
}
//reads out single seq + quality entry from fastq file
shared_ptr<DNA> InputStreamer::read_fastq_entry_fast(istream & fna, int& lnCnt, bool& corrupt) {
	string line;
	if (!fna) { return NULL; cerr << "Read Fastq_f: Input stream does not exist " << lnCnt << ".\n"; exit(53); }

	if (!safeGetline(fna, line)) { return NULL; }
	shared_ptr<DNA> tdn = make_shared<DNA>("", "");
	
	if (tdn == nullptr) {
	    cout << "nulli" << endl;
	    exit(0);
	}
	
	if (line.length() == 0) { return tdn; }
	while (line[0] != '@') {
		IO_Error("ERROR on line " + itos(lnCnt) + ": Could not find \'@\' when expected (file likely corrupt, trying to recover):\n" + line);// << endl;
		//exit(55);
		//recover instead and go to next entry..
		corrupt = true;
		if (!safeGetline(fna, line)) { corrupt = true;  return NULL; }//delete tdn;
	}
    tdn->setHeader(line.substr(1));
	//cerr << line.substr(1);
	if (!safeGetline(fna, tdn->getSequence())) { corrupt = true;  return NULL; }//delete tdn;
	
	//if (line.length() == 0) { return NULL; }
	//std::transform(line.begin(), line.end(), line.begin(), ::toupper);
	//tdn->appendSequence(line);
	//"+"
	if (!safeGetline(fna, line)) { corrupt = true; return NULL; }//delete tdn;
	if (line[0] != '+') {
		//recovery is hard, just give up this read
		IO_Error("Error input line " + itos(lnCnt + 2) + ": Could not find \'+\' when expected (file likely corrupt, aborting):\n" + line);// << endl;
		corrupt = true;
		return tdn;
		//if (!safeGetline(fna, line)) { delete tdn; return NULL; }
	}

	//qual_ score
	vector<qual_score> Iqual(tdn->mem_length(), 0);
	if (!safeGetline(fna, line)) { corrupt = true;  return NULL; }//delete tdn;
	uint qcnt(0); uint lline = (uint)line.length();
	for (; qcnt < lline; qcnt++) {
		qual_score q = (qual_score)line[qcnt] - fastQver;
		Iqual[qcnt] = q;// minmaxQscore();

	}

	if (qcnt == tdn->mem_length()) {
		//needsAT=true;
		tdn->setQual(move(Iqual));
	} else if (line.length() + qcnt != tdn->length()) {
		//check that quality gets not more length than DNA
		corrupt = true;
		IO_Error("Error input line " + itos(lnCnt + 3) + ": More quality positions than nucleotides detected for sequence\n " +
                                                         tdn->getId());// << endl;
//		exit(7);
	}
	lnCnt+=4;
	//if (tdn!= NULL && !tdn->control()){delete tdn; tdn=NULL;}
	corrupt = false;
	return tdn;
}
void InputStreamer::minmaxQscore(shared_ptr<DNA> t, bool& print) {
	fqverMTX.lock();
	vector<qual_score> Qs = t->getQual();
	for (size_t i = 0; i < Qs.size(); i++) {
		minmaxQscore(Qs[i], print);
	}
	fqverMTX.unlock();
}

qual_score InputStreamer::minmaxQscore(qual_score t,bool& print) {
	if (t < 0) { ////quick hack, since q<13 scores are the only deviants and uninteresting in most cases..
		if (fqSolexaFmt ){
			if (t < -5 && print) {
				cerr << "Unusually low sloexa quality score (" << t << "); setting to 0.\n";
				print = false;
			}
		} else {
			if (t >= -5) {
				if (print) {
					cerr << "Resetting auto format to Solexa (illumina 1.0-1.3) format.\n";
					print = false;
				}
				//fqSolexaFmt = true; 
			} else {
				if (print) {
					cerr << "Unusually low quality score (" << t << "); setting to 0.\n";
					print = false;
				}
			}
		}
		t = 0;
	}
	if (minQScore > t) {
		minQScore = t;
		if (minQScore < 0) {
		}
	} else if (maxQScore < t) {
		maxQScore = t;
	}
	return t;
}
bool InputStreamer::checkInFileStatus() {
	for (uint i = 0; i < 3; i++) {
		if (isFasta) {
			if (fasta_istreams[i] != NULL && !fasta_istreams[i]->eof()) {
				return true;
			}
		} else {
			if (fastq_istreams[i] != NULL && !fastq_istreams[i]->eof()) {
				return true;
			}
		}
	}
	return false;
}
void InputStreamer::allStreamReset() {
	resetStats();
#ifdef DEBUG
	cerr << "Resetting input streams" << endl;
#endif
	//reopen streams in gz case // sdm 1.01: make default
	
	for (uint i = 0; i < 3; i++) {
		// && !fasta_istreams[i]->eof()
		if (fasta_istreams[i] != NULL) { fasta_istreams[i]->reset(); }//fasta_istreams[i]->seekg(0, ios::beg);
		if (quality_istreams[i] != NULL ) { quality_istreams[i]->reset(); }// quality_istreams[i]->seekg(0, ios::beg);
		if (fastq_istreams[i] != NULL ) { fastq_istreams[i]->reset();  }
	}
	
	//checkInFileStatus();
}
void InputStreamer::allStreamClose(){
	for (uint i = 0; i < 3; i++){
/*		if (*(fna[i])){ fna[i]->close(); }
		if (*(qual_[i])){ qual_[i]->close(); }
		if (*(fastq[i])){ fastq[i]->close(); }
#ifdef _gzipread
		if (gzfna[i]){ gzfna[i].close(); }
		if (gzqual[i]){ gzqual[i].close(); }
		if (gzfastq[i]){ gzfastq[i].close(); }
#endif
		*/
		if (fasta_istreams[i] != NULL) { delete fasta_istreams[i];  fasta_istreams[i] = NULL;}
		if (quality_istreams[i] != NULL) { delete quality_istreams[i];  quality_istreams[i] = NULL; }
		if (fastq_istreams[i] != NULL) { delete fastq_istreams[i];   fastq_istreams[i] = NULL; }
	}

	if (!isFasta && minQScore < SCHAR_MAX) {
			maxminQualWarns_fq( );
	}
}
void InputStreamer::jumpToNextDNA(bool& stillMore, int pos) {
	if (isFasta) {//get DNA from fasta + qual_ files
	    // Read fasta entry
		stillMore = read_fasta_entry(fasta_istreams[pos], quality_istreams[pos], dnaTemp1[pos], dnaTemp2[pos], lnCnt[pos]);
        dnaTemp1[pos] = dnaTemp2[pos];
		dnaTemp2[pos] = make_shared <DNA>("", "");
	} else {
		//jmp_fastq(*fastq_istreams[pos], lnCnt[pos]);
		fastq_istreams[pos]->jumpLines(4);
		if (fastq_istreams[pos]->eof()) {
			stillMore = false;
		}
	}
}
/*
shared_ptr<DNA> InputStreamer::getDNA(bool &has_next, bool &repair_input_stream, int pos) {
    
	if (pos == 1 && numPairs < 2) {
        return nullptr;
    }
    
    shared_ptr<DNA> dna;
    
    // Only for fastq
    bool corrupt(true); // Set true initially
    bool sync(false);
    
    //while (corrupt) {
    if (isFasta) { //get DNA from fasta + qual_ files
        has_next = read_fasta_entry(*(fasta_istreams[pos]), *(quality_istreams[pos]), dnaTemp1[pos], dnaTemp2[pos], lnCnt[pos]);
        corrupt = false;
        
        //dnaTemp1 will be completed, dnaTemp2 will have a header set up
        dna = dnaTemp1[pos];
        dnaTemp1[pos] = dnaTemp2[pos];
        dnaTemp2[pos].reset(new DNA("", ""));
        if (dna->isEmpty()) { // SEAL was here
            dna = nullptr;
        }
    
        // DISPLAY PROGRESS (NO IMPORTANT COMPUTATION)
        if (pos == 0 && pairs_read[pos] % 100 == 0) {
            _drawbar(*(fasta_istreams[pos]));
        } else if (!has_next) {
            _drawbar(*(fasta_istreams[pos]));
        }
    } else { // Read fastq file>
        if (fqReadSafe) {
            dna = read_fastq_entry(*(fastq_istreams[pos]), minQScore, lnCnt[pos], corrupt, repair_input_stream);
            
            if (fastQver == 0 && dna->length() > 5 && !corrupt) {//autodetect
                dna->resetQualOffset(auto_fq_version(), fqSolexaFmt);
                //reset streams
            }
            if (lnCnt[pos] > 100) {//tmp set back to 500
                if (fqPassedFQsdt)
                    fqReadSafe = false;
            }
        } else {
            dna = read_fastq_entry_fast(*(fastq_istreams[pos]), lnCnt[pos], corrupt);
        }
        
        // THIS part can be put somewhere else
        // If has no next or fastq stream is end of file or is nullptr?
        if (!has_next || fastq_istreams[pos]->eof() || (!*(fastq_istreams[pos]))) {
            if (dna != nullptr) {
                if (dna->isEmpty()) // SEAL HERE
                    dna = nullptr;
            } //delete ret;
            has_next = false;
            return nullptr;
            
        } else if (dna == nullptr || dna->isEmpty()) {
            corrupt = true;
        }
        
        // DISPLAY PROGRESS
        if (pos == 0 && pairs_read[pos] % 100 == 0) {
            if (_drawbar(*(fastq_istreams[pos]))) {
                has_next = false;
                return nullptr;
            }
        } else if (!has_next) {
            _drawbar(*(fastq_istreams[pos]));
        }
        
        
        if (corrupt) {
            //delete dna object;
            dna = nullptr;
            sync = true;
            repair_input_stream = true;
        }
    }
    if (!keepPairHD) {//better cut in early
        dna->setHeader(dna->getPositionFreeId());
    }
    pairs_read[pos]++;
    
    return dna;
}
*/

/*
inline void InputStreamer::addQuality(shared_ptr<DNA> fasta, FastxRecord* quality) {
    if (quality != nullptr) {
		rtrim(quality->sequence);
        const char* qualityString = quality->sequence.c_str();
        vector<qual_score> qualityVector(fasta->getSequence().size(), 0);
        qual_score qualScore(0);
        int i(0);
        
        for (; i < qualityVector.size(); i++) {
            qualScore = (qual_score)parseInt(&qualityString);
            qualityVector[i] = qualScore;
            if (*qualityString == '\0') {
                break;
            }
        }
        
        if (fasta->getSequence().size() != qualityVector.size()) {
            cerr << "ERROR: quality counts (" << i << ") not the same length as DNA base counts in Sequence (" << (qualityVector.size() - 1) << ")\n" << fasta->getId() << "\n";
            exit(54);
        }
        fasta->appendQuality(qualityVector);
    }
}

void InputStreamer::getDNA(
        FastxRecord* read1, FastxRecord* read2, FastxRecord* quality1, FastxRecord* quality2, FastxRecord* FRM,
        shared_ptr<DNA>* dna1, shared_ptr<DNA>* dna2, shared_ptr<DNA>* mid) {
	if (read1 != nullptr) {
        *dna1 = make_shared<DNA>(read1, minQScore, maxQScore, fastQver);
        
		if (quality1 != nullptr) {
		    addQuality(*dna1, quality1);
		}
  
		bool seal = (*dna1)->seal();
		bool isempty = (*dna1)->isEmpty();
		
		if (!(*dna1)->seal() || (*dna1)->isEmpty()) {
		    exit(0);
 		    
            dna1 = nullptr;
		} else {
			if (!keepPairHD) {//better cut in early
                (*dna1)->setHeader((*dna1)->getPositionFreeId());
			}
			if (read2 != nullptr) {
                (*dna1)->setpairFWD();
			}
		}
	}
	if (read2 != nullptr) {
        *dna2 = make_shared<DNA>(read1, minQScore, maxQScore, fastQver);
        
        if (quality2 != nullptr) {
        
        }
		
		if (!(*dna2)->seal() || (*dna2)->isEmpty()) {
            (*dna2) = nullptr;
		}
		else {
			if (!keepPairHD) {//better cut in early
                (*dna2)->setHeader((*dna2)->getPositionFreeId());
			}
			if (read1 != nullptr) {
                (*dna2)->setpairREV();
			}
		}

	}
	if (FRM != nullptr) {
        *mid = make_shared<DNA>(read1, minQScore, maxQScore, fastQver);
		if (!(*mid)->seal() || (*mid)->isEmpty()) {
            mid = nullptr;
		}
		else {
            (*mid)->setMIDseq(true);
		}
	}

}
*/

vector<shared_ptr<DNA>> InputStreamer::getDNApairs() {
	vector<shared_ptr<DNA>> ret(3, nullptr);
	ret[0] = getDNA(0);

	//first entry should always exist!
	if (ret[0] == nullptr) {
		return ret;
	}
	if (!ret[0]->seal() || ret[0]->isEmpty()) { ret[0] = NULL; }

	if (numPairs == 2) {
		ret[1] = getDNA(1);
		if (!ret[1]->seal() || ret[1]->isEmpty()) { ret[1] = NULL; }
	}
	if (hasMIDs) {
		ret[2] = getDNA(2);
		if (ret[2] == NULL || !ret[2]->seal() || ret[2]->isEmpty()) {
			ret[2] = NULL; 
		}
		else {
			ret[2]->setMIDseq(true);
		}
	}
	
	return ret;
}


shared_ptr<DNA> str2DNA(vector<string>& in, bool keepPairHD,int fastQver, int readpos) {
	shared_ptr<DNA> ret(nullptr);
	bool isFasta(false);
	if (in.size()==3 || in[0][0] == '>'  ) {//fasta
		ret = make_shared < DNA >(in);
		isFasta = true;
	} else {
		ret = make_shared < DNA >(in, fastQver);
	}
	if (!keepPairHD) {//better cut in early
		ret->setHeader(ret->getPositionFreeId());
	}
	if (readpos == 0) {
		ret->setpairFWD();
	}else if (readpos == 1) {
		ret->setpairREV();
	}else if (readpos == 2) {
		ret->setMIDseq(true);
	}
	if (!ret->seal(isFasta) || ret->isEmpty()) {
		ret = nullptr;
	}
	if (in[2] == "") {
		int x = 0;
	}
	return ret;
}

/*vector<shared_ptr<DNA>> InputStreamer::getDNAMC() {
	
	vector<shared_ptr<DNA>> ret(3,nullptr);
	//there can only be 1 stream with this..
	if (at_thread >= num_threads) { at_thread = 0; }
	if (slots[at_thread].inUse == false) {
		slots[at_thread].inUse = true;
#ifdef togRead
		slots[at_thread].Cfut = async(std::launch::async, &InputStreamer::getDNApairs, this);
#elif defined IOsep
		vector<string>tmpLines(4, "");	getDNAlines(tmpLines, 0);
		slots[at_thread].fut = async(std::launch::async, str2DNA, tmpLines, keepPairHD, fastQver,0);
		if (numPairs == 2) {
			vector<string>tmpLines(4, "");	getDNAlines(tmpLines, 1);
			slots[at_thread].fut2 = async(std::launch::async, str2DNA, tmpLines, keepPairHD, fastQver,1);
		}
		if (hasMIDs) {
			vector<string>tmpLines(4, "");	getDNAlines(tmpLines, 2);
			slots[at_thread].fut3 = async(std::launch::async, str2DNA, tmpLines, keepPairHD, fastQver,2);
		}
#else
		slots[at_thread].fut = async(std::launch::async, &InputStreamer::getDNA, this, 0);
		if (numPairs == 2) { slots[at_thread].fut2 = async(std::launch::async, &InputStreamer::getDNA, this, 1); }
		if (hasMIDs) { slots[at_thread].fut3 = async(std::launch::async, &InputStreamer::getDNA, this, 2); }
#endif// togRead
	}
	//always get DNA in pairs..
	if (slots[at_thread].inUse  ){
		//&& slots[at_thread].Cfut.wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
		slots[at_thread].inUse = false;
#ifdef togRead
		ret = slots[at_thread].Cfut.get();
#else
		ret[0] = slots[at_thread].fut.get();
		if (numPairs == 2) {
			ret[1] = slots[at_thread].fut2.get();
		}
		if (hasMIDs) {
			ret[2] = slots[at_thread].fut3.get();
			if (ret[2] != nullptr) {
				ret[2]->setMIDseq(true);
			}
		}
#endif// togRead


		if (num_threads == 1) {
			return ret;
		}

		//shoot off next job already
		slots[at_thread].inUse = true;
#ifdef togRead
		slots[at_thread].Cfut = async(std::launch::async, &InputStreamer::getDNApairs, this);
#elif defined IOsep
		vector<string>tmpLines(4, "");	getDNAlines(tmpLines, 0);
		slots[at_thread].fut = async(std::launch::async, str2DNA, tmpLines, keepPairHD, fastQver,0);
		if (numPairs == 2) {vector<string>tmpLines(4, "");	getDNAlines(tmpLines, 1);
			slots[at_thread].fut2 = async(std::launch::async, str2DNA, tmpLines, keepPairHD, fastQver,1);	}
		if (hasMIDs) {vector<string>tmpLines(4, "");	getDNAlines(tmpLines, 2);
			slots[at_thread].fut3 = async(std::launch::async, str2DNA, tmpLines, keepPairHD, fastQver,2);	}

#else
		slots[at_thread].fut = async(std::launch::async, &InputStreamer::getDNA, this, 0);
		if (numPairs == 2) {slots[at_thread].fut2 = async(std::launch::async, &InputStreamer::getDNA, this, 1);		}
		if (hasMIDs) {slots[at_thread].fut3 = async(std::launch::async, &InputStreamer::getDNA, this, 2);		}
#endif // togRead
	}
/*		if (slots[at_thread].inUse == false) {
			slots[at_thread].inUse = true;
			slots[at_thread].Cfut = async(std::launch::async, &InputStreamer::getDNApairs, this);


			/*
			//bit more object disentangled, but problems with non-fastq's, so not used
			//1: get new fastq lines
			//ret[0] = str2DNA(tmpLines, keepPairHD, fastQver);
			vector<string>tmpLines(4, "");	getDNAlines(tmpLines,0);
			slots[at_thread].fut = async(std::launch::async, str2DNA, tmpLines, keepPairHD, fastQver);
			
			//ret[0] = slots[at_thread].fut.get();
			//slots[at_thread].inUse = false;
			//return ret;
			
			if (numPairs == 2) {
				vector<string>tmpLines(4, "");	getDNAlines(tmpLines, 1);
				slots[at_thread].fut2 = async(std::launch::async, str2DNA, tmpLines, keepPairHD, fastQver);
//				slots[at_thread].fut2 = async(std::launch::async, &InputStreamer::getDNA, this, 1);
			}
			if (hasMIDs) {
				vector<string>tmpLines(4, "");	getDNAlines(tmpLines, 2);
				slots[at_thread].fut3 = async(std::launch::async, str2DNA, tmpLines, keepPairHD, fastQver);
//				slots[at_thread].fut3 = async(std::launch::async, &InputStreamer::getDNA, this, 2);
			}
			
			
		}
	return ret;
	
}
*/

string InputStreamer::current_infiles() {
	string ret(""); int pos(0);
	if (isFasta) {
		fasta_istreams[pos]->getInFile();
	} else {// fastq format
		fastq_istreams[pos]->getInFile();
	}
	return ret;
}

bool InputStreamer::getDNAlines(multi_tmp_lines* tmpO, int blocks, bool MIDuse,bool safe) {

	


	if (safe) {
		protect.lock();
	}

	assert(tmpO->size() == blocks);
	size_t k(0); bool b1(true), b2(true);
	for (k = 0; k < blocks; k++) {
		if (_globalMaxRdsRead > 0 && _globalRdsRead + _localRdsRead > _globalMaxRdsRead) {
			if (safe) { protect.unlock(); }
			return false;
		}
		b1 = this->getDNAlines(tmpO->tmp[k][0], 0);
		_localRdsRead++;
		if (numPairs > 1) {
			b2 = this->getDNAlines(tmpO->tmp[k][1], 1);
			_localRdsRead++;
		}
		if (MIDuse) {
			this->getDNAlines(tmpO->tmp[k][2], 2);
		}
		if (!b1 || !b2) {
			if ((b1 != b2) && tmpO->tmp[k][1].size() != tmpO->tmp[k][0].size() && numPairs == 2) {
				//cerr << "Currently reading: "<<current_infiles()<<endl;
				cerr << "Problem: in file "<< current_infiles () << " read1 (" << tmpO->tmp[k][0].size() << ") and read2 (" << tmpO->tmp[k][1].size()<< ") appear not be of the same size!\n";
			}
			tmpO->setSize(k+1);
			if (tmpO->tmp[k][0][0].length() == 0) {//extra check if maybe entire attempt was empty
				tmpO->setSize(k);
			}

			if (safe) { protect.unlock(); }
			return false;
		}
	}
	if (safe) {protect.unlock();}

	return true;
}
bool InputStreamer::getDNAlines(vector<string>& ret, int pos) {
	
	//vector<string> ret(4,"");
	bool stillMore = true;
	if (pos == 1 && numPairs <= 1) {
		return false;
	}
	bool repairInStream(false);
	//vector<string>tmpLines(4, "");
	if (isFasta) {//get DNA from fasta + qual_ files
		//cerr << "Can't read fasta in multicore mode!\n";
		//exit(232);
		stillMore = fasta_istreams[pos]->getline(ret[0]);//header
		int lnRd = 1;//read 1 header line
		stillMore = fasta_istreams[pos]->getlines(ret[1], lnRd);//fasta
		lnCnt[pos] += lnRd;
		if (quality_istreams[pos] != nullptr) {
			stillMore = quality_istreams[pos]->getline(ret[2]);//qual
			stillMore = quality_istreams[pos]->getlines(ret[2], lnRd, true);//qual
		}
	} else {// fastq format
		//fqRead
		stillMore = fastq_istreams[pos]->get4lines(ret);
		/*for (int xx = 0; xx < 4; xx++) {
			stillMore = fastq_istreams[pos]->getline(ret[xx]);
		}*/
		lnCnt[pos] += 4;
		if (fastQver == 0 || lnCnt[pos] == 4) { // ok first time just has to be done in function, after this can be unloaded outside
			bool printB(true);
			//cerr << "AutoFQ!!";
			shared_ptr<DNA> tdn = make_shared<DNA>(ret, fastQver);
			minmaxQscore(tdn,printB);
			auto_fq_version();
			//tdn->resetQualOffset(auto_fq_version(), fqSolexaFmt);
		}
	}
	/*
	if (pos == 0 && pairs_read[pos] % 100 == 0) {
		if (_drawbar(*(fastq_istreams[pos]))) { stillMore = false; }
	}
	else 	if (!stillMore) { _drawbar(*(fastq_istreams[pos])); }*/
	if (stillMore) {
		pairs_read[pos]++;
	}
	return stillMore;
}

shared_ptr<DNA> InputStreamer::getDNA(int pos){
	//if (sync) {
	//	while (desync(pos)) {
	//		jumpToNextDNA(stillMore, pos);
	//	}
	//}
	bool stillMore = true;
	if (pos == 1 && numPairs <= 1) {
		return NULL;
	}
	shared_ptr<DNA> ret;
	bool corrupt(false); //corrupt state isn't implemented for fnaread
	bool repairInStream(false);
	if (isFasta) {//get DNA from fasta + qual_ files
		vector<string>tmpStr(3, "");
		stillMore = fasta_istreams[pos]->getline(tmpStr[0]);//header
		int lnRd = 1;//read 1 header line
		stillMore = fasta_istreams[pos]->getlines(tmpStr[1], lnRd);//fasta
		lnCnt[pos] += lnRd;
		if (quality_istreams[pos] != nullptr) {
			stillMore = quality_istreams[pos]->getline(tmpStr[2]);//qual
			stillMore = quality_istreams[pos]->getlines(tmpStr[2], lnRd);//qual
		}

		//ret = make_shared<DNA>(tmpStr);
			
		ret = str2DNA(tmpStr, keepPairHD, fastQver,pos);

		if (!stillMore || fasta_istreams[pos]->eof() || (!*(fasta_istreams[pos]))) {
			if (ret != NULL) { if (!ret->seal() || ret->isEmpty()) { ret = NULL; } } //delete ret;
			stillMore = false; 
		}
		else if (ret == NULL || !ret->seal() || ret->isEmpty()) {
			corrupt = true;
		}

		//old way, don't use any longer...
		/*stillMore = read_fasta_entry(fasta_istreams[pos], quality_istreams[pos], dnaTemp1[pos], dnaTemp2[pos], lnCnt[pos]);
		corrupt = false;
		//dnaTemp1 will be completed, dnaTemp2 will have a header set up
		ret = dnaTemp1[pos];
        dnaTemp1[pos] = dnaTemp2[pos]; dnaTemp2[pos].reset(new DNA("", ""));
		if (!ret->seal() || ret->isEmpty()) { ret = NULL; }
		/*if (pos == 0 && pairs_read[pos] % 100 == 0) {
			_drawbar(*(fasta_istreams[pos]));
		}
		else 	if (!stillMore) { _drawbar(*(fasta_istreams[pos])); }*/
	}
	else { //fqRead
		/*if ( fqReadSafe) {
			ret = read_fastq_entry(*(fastq_istreams[pos]), minQScore, lnCnt[pos], corrupt, repairInStream);
			if (fastQver == 0 && ret->length() > 5 && !corrupt) {//autodetect
				minmaxQscore(ret);
				ret->resetQualOffset(auto_fq_version(), fqSolexaFmt);
				//reset streams
			}
			if (lnCnt[pos] > 100) {//tmp set back to 500
				if (fqPassedFQsdt) {
					fqReadSafe = false;
					cdbg("Switching to fast fastq reader..\n ");
				}
			}
		}
		else {
			ret = read_fastq_entry_fast(*(fastq_istreams[pos]), lnCnt[pos], corrupt);
		}
		/**/
		vector<string>tmpLines(4, "");
		getDNAlines(tmpLines,pos);
		ret = str2DNA(tmpLines, keepPairHD, fastQver,pos);

		if (!stillMore || fastq_istreams[pos]->eof() ) {
			if (ret != NULL) { if (!ret->seal() || ret->isEmpty()) { ret = NULL; } } //delete ret;
			stillMore = false;
		} else if (ret == NULL || !ret->seal() || ret->isEmpty()) {
			corrupt = true;
		}
			
		/*if (pos == 0 && pairs_read[pos] % 100 == 0) {
			if (_drawbar(*(fastq_istreams[pos]))) { stillMore = false; break; }
		}
		else if (!stillMore) { _drawbar(*(fastq_istreams[pos])); }
		*/
			
		if (!keepPairHD) {//better cut in early
			ret->setHeader(ret->getPositionFreeId());
		}
		if (corrupt) {
			//delete ret;
			ret = NULL;
			//sync = true;
			repairInStream = true;
		}
	}
	if (!corrupt) {
		pairs_read[pos]++;
	}
	//last check
	
	//
	return ret;
}

void InputStreamer::openMIDseqs(string p,string in){
	if (in==""){return;}
	
	if (fastq_istreams[2] != NULL){
		cerr << "MID file was already initialized" << endl;
	}
#ifdef DEBUG
	cerr << "Open Mid sequence_ file" << endl;
#endif

	string file_type = "MID specific fastq";
	string tmp = (p + in);
	bool doMC = num_threads > 1;

	fastq_istreams[2] = DBG_NEW ifbufstream(tmp, (size_t)round(INPUT_BUFFER_SIZE*0.6),doTIO);
	if (fastq_istreams[2]->eof()) {
		cerr << "\nCouldn't find " << file_type << " " << tmp << "!\n Aborting..\n";		exit(4);
	}
	hasMIDs=true;
}

bool InputStreamer::setupFastq(string path, string fileS, int& pairs, string subsPairs,
	bool simu, bool verbose) {
	allStreamClose();
	minQScore = SCHAR_MAX;
	maxQScore = -1;
	fqSolexaFmt = false;
	resetStats();
	vector<string> tfas = splitByCommas(fileS);
	if ( pairs == -1 ) {
		pairs = (int) tfas.size();
		if ( pairs > 1 ) { cerr << "Paired input (\",\" separated) detected\n"; }
	}
	numPairs = pairs;
	string p1(""), p2(""), midp("");
	string xtraMsg = "";
	if (getCurFileN() > 0) {
		xtraMsg = " " + itos(getCurFileN()) + " of " + itos(totalFileNumber);
		if (BCnumber > 1) {
			xtraMsg = ", looking for " + itos(BCnumber) + "BCs.\n";
		} else {
			xtraMsg = "";
		}
	}

	if (tfas.size() != (uint)pairs && subsPairs == "") {
		cerr << "Unequal number of files (" << tfas.size() << ") and option-set paired files (" << pairs << ").\n Aborting...\n"; exit(52);
	}
	string file1("");
	if (tfas.size() == 3) {
		if (tfas.size() != 3) { cerr << "Could not detect 3 input files in string\n" << fileS << "\n Aborting.." << endl; exit(76); }
		midp = path + tfas[1];
		p1 = path + tfas[0];
		p2 = path + tfas[2];
		file1 = tfas[0];
//		cerr << p1 << " + " << p2 << " and " << midp << endl;
	} else if (tfas.size() == 2) {
		p1 = path + tfas[0];
		p2 = path + tfas[1];
		file1 = tfas[0];
	} else if (tfas.size() == 1) {
		p1 = path + fileS;
		file1 = fileS;
		//cerr << "Reading fastq " << p1 << endl;
	}
	if (subsPairs == "1") {
		p2 = "";
	} else 	if (subsPairs == "2") {
		p1 = p2; p2 = "";
	}

	if (simu) {

		return fileExists(p1) && fileExists(p2) && fileExists(midp);
	}


	if (verbose && !p1.empty() && !p2.empty()) {
		if (!midp.empty()) {
			//cerr << "Reading paired fastq + MID file" << xtraMsg<<"."<<endl;
		} else {
			//cerr << "Reading paired fastq" << xtraMsg << "." << endl;
			//cerr << p1 << " + " << p2 << endl;
		}
	} else  {
		//cerr << "Reading fastq" << xtraMsg << "." << endl;
		//cerr << p1 << endl;
	}
	if (verbose) {
		cerr << "At " << file1 << ":\n";
	}


	bool suc = setupFastq_2(p1, p2, midp);
	//read progress  bar setup
	//_measure(*fastq_istreams[0]);
	//allStreamClose();
	//setupFastq_2(p1, p2, midp);
	return suc;
}

// Prints out the progress bar
inline void InputStreamer::_print(int cur, float prog) {

	std::cerr << std::fixed << std::setprecision(2)
		<< "\r   [" << std::string(cur, '#')
		<< std::string(_max + 1 - cur, ' ') << "] " << 100 * prog << "%";

	if (prog == 1) std::cerr << std::endl;
	else std::cerr.flush();

}
inline bool InputStreamer::_drawbar(istream & tar) {
    if (_fileLength <= 0) { return false; }
    int pos((int) tar.tellg());
    float prog(pos / float(_fileLength)); // percentage of infile already read
    if (pos == -1 || prog > 1.f) {
        _print(_max + 1, 1);
        _fileLength = 0;
        return false;
    }
    
    // Number of #'s as function of current progress "prog"
    int cur((int) ceil(prog * (float) _max));
    if (_last != cur) _last = cur, _print(cur, prog);
    return prog == 1.f;
}

inline void InputStreamer::_measure(istream& tar) {
	tar.seekg(0, ios_base::end);
	_fileLength = (int) tar.tellg();
	tar.seekg(0, ios_base::beg);
	tar.clear();
}

bool InputStreamer::setupFastq_2(string p1, string p2, string midp) {
	//setupFastq_2(fastqFilepathTemp[0],fastqFilepathTemp[1],fastqFilepathTemp[2]);
#ifdef DEBUG
	cerr << "setupFastq2 " << p1<<endl<<midp << endl;
	//cerr << fastqFilepathTemp.size() << endl;
#endif
	string file_type = "";

	size_t bufS = INPUT_BUFFER_SIZE ;

	//bool doMC = num_threads > 1;

	//        INPUT   files
	if (!p1.empty()){ // first file exists
		file_type = "fastq file 1";
		//setup buffer of different sizes for p1,p2 to avoid simultaneous read
		fastq_istreams[0] = DBG_NEW ifbufstream(p1, (size_t)round(bufS*0.8), doTIO);
		if (fastq_istreams[0]->eof()){ 
			cerr << "\nCouldn't find " << file_type << " " << p1 << "!\n Aborting..\n";		exit(4); 
		}
	}
	//second pair_
	if (!p2.empty()){
		file_type = "fastq file 2";
		fastq_istreams[1] = DBG_NEW ifbufstream(p2, (size_t)round(bufS * 1.2), doTIO);
		if (fastq_istreams[1]->eof()) {
			cerr << "\nCouldn't find " << file_type << " " << p2 << "!\n Aborting..\n";		exit(4);
		}
	}
	//MID file
	if (midp != ""){
		this->openMIDseqs("", midp);
	}
	return true;
}
//mainFile = IS->setupInput(path, uniqueFas[i], FastqF[tarID], FastaF[tarID], QualF[tarID], MIDfq[tarID], fil->isPaired(), (*cmdArgs)["-onlyPair"]);
string InputStreamer::setupInput(string path, int t, const string& uniqueFastxFile, 
		filesStr& files, int &paired, string onlyPair,
          string& mainFilename, bool simulate) {
	string mainFilepath("");
	vector<string> fastqFiles = files.FastqF;
	vector<string> fastaFiles = files.FastaF;
	vector<string> qualityFiles = files.QualF;
	vector<string> midFiles = files.MIDfq;

	
	if (isFasta) { // input is fasta
		if (fastaFiles[t] != uniqueFastxFile) {
			cerr << "Error in matching FASTA target filenames.\n";
			exit(11);
		}
		this->setupFastaQual(path, fastaFiles[t], qualityFiles[t], paired, onlyPair);
        mainFilepath = path + fastaFiles[t];
        mainFilename = fastaFiles[t];
	} else { // Input is fastq
		if (fastqFiles[t] != uniqueFastxFile) {
			cerr << "Error in matching target filenames.\n";
			exit(11);
		}
//        std::cout << "setup fastqFiles with: " << path << ", " << fastqFiles[t] << std::endl;
		this->setupFastq(path, fastqFiles[t], paired, onlyPair, simulate, !simulate);
        mainFilepath = path + fastqFiles[t];
        mainFilename = fastqFiles[t];
	}
	this->openMIDseqs(path, midFiles[t]);
	return mainFilepath;
}


bool InputStreamer::setupFastaQual(string path, string sequenceFile, string qualityFile, int& paired, string onlyPair, bool simulate) {
	resetStats();
	allStreamClose();
	vector<string> splitSequenceFile = splitByCommas(sequenceFile);

    //Debug
    std::cout << "sequence filename: " << sequenceFile << std::endl;
    std::cout << "quality filename: " << qualityFile << std::endl;
	
	if (paired == -1 ) {
        paired = (int) splitSequenceFile.size();
		if (paired > 1 ) { cerr << "Paired input (\",\" separated) detected\n"; }
	}
	if (paired > 1) { cerr << "Paired fasta+qual_ is currently not implemented\n"; exit(72); }
	if (onlyPair == "1") {
		//
	} else 	if (onlyPair == "2") {
		//
	}

	numPairs = paired;
	dnaTemp1[0]= make_shared<DNA>("", "");
	dnaTemp2[0] = make_shared <DNA>("", "");
    
    fastaFilepathTemp[0] = path + sequenceFile;
    qualityFilepathTemp[0] = path + qualityFile;
	
	cout << path << " " << qualityFilepathTemp[0] << endl;
	cout << qualityFile << endl;
	
	if (simulate) { // SET TO FALSE
		if (!fileExists(qualityFilepathTemp[0], -1, false)) {
			cerr << "Warning: corresponding quality file is missing: " << qualityFilepathTemp[0] << endl;
		}
		return fileExists(fastaFilepathTemp[0]) ;
	}

	bool successful = setupFastaQual2(fastaFilepathTemp[0], qualityFilepathTemp[0]);
	
	//activate length measure
	//_measure(*fasta_istreams[0]);
	cout << "success" << endl;
	return successful;
}
bool InputStreamer::setupFastaQual2(string sequenceFile, string qualityFile, string fileType) {
	string file_typeq = "quality file";
	size_t iBufS = INPUT_BUFFER_SIZE;
	//bool doMC = num_threads > 1;
	//        INPUT   file
	if (sequenceFile != "") { // sequence not empty
		fasta_istreams[0] = DBG_NEW ifbufstream(sequenceFile.c_str(), (size_t)round(iBufS * 0.4), doTIO);

		if (fasta_istreams[0]->eof()) {
		    cerr << "\nCouldn't find " << fileType << " file \"" << sequenceFile << "\"!\n Aborting..\n";
		    exit(4);
		}
	}
	//quality file
	if (!qualityFile.empty()) {
		quality_istreams[0] = DBG_NEW ifbufstream(qualityFile, (size_t)round(iBufS*0.5),doTIO);
		if (quality_istreams[0]->eof()) {
			cerr << "\nCouldn't find " << file_typeq << " file \"" << qualityFile << "\"!\n Running in no qual_ filter mode\n";
			qualAbsent = true;
			//exit(55);
		}
		if (getCurFileN() > 0) {
			cerr << "Reading Fasta + Quality file " << getCurFileN() << " of " << totalFileNumber << ".\n";
		} else {
			cerr << "Reading Fasta + Quality file.\n";
		}
		cerr << fileType << " : " << sequenceFile << endl << file_typeq << " : " << qualityFile << endl;
		
	} else if (sequenceFile != "") {
		if (getCurFileN() > 0) {
			cerr << "Reading Fasta file " << getCurFileN() << " of " << totalFileNumber << ".\n";
		} else {
			cerr << "Reading Fasta.\n";
		}
		quality_istreams[0] = nullptr;
		qualAbsent = true;
	}
	return true;
}
void InputStreamer::setupFna(string fileS){
	allStreamClose();
	resetStats();
	numPairs = 1;
	dnaTemp1[0] = make_shared<DNA>("", ""); dnaTemp2[0] = make_shared <DNA>("", "");
	setupFastaQual2(fileS, "","seq");
}

int InputStreamer::auto_fq_version() {
	//if (QverSet) { return 33; }
	int fqDiff(0);
	if (minQScore >= 100 || maxQScore < 2 || fastQver!=0) {
		//don;t know what to do here..
		return fqDiff;
	}
	fqverMTX.lock();
	fqSolexaFmt = false;

	if (minQScore >= 33 && maxQScore <= 88 ){
		fqDiff = (fastQver - 33); fastQver = 33; 
		
	} else if (minQScore >= 59 && maxQScore > 74) {
		fqDiff = (fastQver - 64); fastQver = 64;
		if (minQScore < 64) { //set to illumina1.0 (solexa)
			fqSolexaFmt = true;
			//cerr << "Setting to illumina 1.0-1.3 (solexa) fastq version (q offset = 64, min Q=-5).";
		}
		else {
			//cerr << "Setting to illumina 1.3-1.8 fastq version (q offset = 64).";
		}
		//cerr << "Setting to Sanger fastq version (q offset = 33).";
	} else {
		if (verbose) {
			cerr << "\nUndecided fastq version.. reverting to +33\n"; verbose = 0;
		}

		//fallback to sensible std..
		fqDiff = (fastQver - 33); fastQver = 33;
		//exit(53);
	}
	QverSet = true;
	fqverMTX.unlock();
	return fqDiff;
}

void InputStreamer::maxminQualWarns_fq(){
	if (minQScore >= 59-33 && fastQver==33 ){ //set to sanger version, but no low qual_ over whole dataset -> probably illumina version
		//cerr << "     WARNING :: \nQuality scores in your dataset are unusually high (min Q=" << int(minQScore)<<"). Please check that you have a fastQ file in NCBI SRA, Sanger or Illumina 1.8+ version.\nIf not, set fastqVersion in option file to \"3\" (Illumina 1.0) or \"2\" (Illumina 1.3 < 1.8) .\n\n";
	} 
	if (minQScore < 0 ){ //set to sanger version, but no low qual_ over whole dataset -> probably illumina version
		cerr << "     WARNING :: \nQuality scores in your dataset are unusually low (min Q=" << int(minQScore)<< "). Please check that you have a fastQ file in Illumina 1.0 or Illumina 1.3 < 1.8 version.\nIf not, set fastqVersion in option file to \"1\" (NCBI SRA, Sanger or Illumina 1.8+ version).\n\n";
	} 
}


ofbufstream::ofbufstream(const string IF, int mif, bool isMC , size_t bufferS ) :
	file(IF), modeIO(mif), used(0), usedW(0),
	coutW(false), isGZ(false), doMC(isMC), primary(nullptr), bufS(bufferS) {
	cdbg("Ini obufstream");

	if (modeIO == ios::out) {
		remove(file.c_str());
	}
	if (file == "T") {//write to ostream
		coutW = true;
		keeper = DBG_NEW char[0];
		keeperW = DBG_NEW char[0];
		return;
	}
	keeper = DBG_NEW char[bufS];
	//second keeper swappable with keeper to always have one added to, one written out (kickoff system)
	keeperW = DBG_NEW char[bufS];
	if (isGZfile(IF)) { //write a gzip out??
		isGZ = true;
#ifndef _gzipread
		cerr << "ofbufstream::gzip outpout not supported in your sdm build\n" << file << endl;
		exit(51);
#endif
		// ostream* os = DBG_NEW zstr::ofstream("output.txt", std::ios::app);
	}
}
ofbufstream::~ofbufstream() {
	cdbg("destroy obufstream ");
	finishWrites();
	deactivate();
	delete[] keeper;
	delete[] keeperW;
}
void ofbufstream::finishWrites() {
	append_mtx_.lock();
	if (hasKickoff) { hasKickoff = false; writeKickoff.get(); }
	if (used > 0) {
		writeStream(false);
	}
	used = 0; usedW = 0;
	append_mtx_.unlock();
}
void ofbufstream::emptyStream() {
	finishWrites();
	used = 0; usedW = 0;
}
void ofbufstream::activate() {
	if (primary != nullptr) {
		return;
	}
	delete primary;
	if (isGZ) {
#ifdef _gzipread
		primary = DBG_NEW zstr::ofstream(file.c_str(), std::ios::app);
#else
		cerr << "ofbufstream::gzip outpout not supported in your sdm build\n" << file << endl;
		exit(51);
#endif
	}
	else {
		//cout << "new str" << endl;	
		primary = DBG_NEW ofstream(file.c_str(), ios::app);
	}
}
void ofbufstream::deactivate() {
	if (primary == nullptr) {
		return;
	}
	//primary->close();
	//output_mtx.lock();
	delete primary;// ->close();
	primary = nullptr;
	//output_mtx.unlock();

}

bool ofbufstream::operator! (void) {
	return false;
}
void ofbufstream::operator<< (const string& X) {
	//cerr << "ostr";
	size_t lX(X.length());
	if (lX == 0) {
		return;
	}
	if (coutW) {
		cout << X;
		return;
	}
	size_t at = 0;
	append_mtx_.lock();
	if (used + lX >= bufS) {
		writeStream();
	}
	while (at < lX) {
		keeper[used] = X[at];
		used++; at++;
		if (used >= bufS) {
			//should normally not be called.. but in case of big X and small buff might be needed..
			writeStream();
		}
	}
	append_mtx_.unlock();
}

void ofbufstream::writeStream(bool doKickoff ) {
	//int counter = 0;
	//out << omp_get_thread_num << ": " << (counter++) << endl;
	if (used == 0) {
		return;
	}
	//do this first to prevent keeper/keeperW getting corrupted
	output_mtx.lock();
	if (hasKickoff) { writeKickoff.get(); hasKickoff = false; }
	bool closeThis = false;
	if (primary == nullptr) {
		this->activate();
		closeThis = false;
	}
	//swap operation
	/*char* keeperTmp = keeperW;
	keeperW = keeper;
	keeper = keeperTmp;*/
	swap(keeper, keeperW);
	usedW = used;
	used = 0;
	//rewrite keeperW with keeper, so append can continue..
	//strcpy(keeperW, keeper);
	output_mtx.unlock();
	//keeper = DBG_NEW char[bufS];
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
	if (doMC && doKickoff) {
		writeKickoff = async(std::launch::async, &ofbufstream::internalWrite,
			this, closeThis);
		hasKickoff = true;
	}
	else {
		// Without multithreading
		internalWrite(closeThis);
	}
}

dualOfBufStream::dualOfBufStream(void):buf1S(20000), buf2S(30000), bufs(2, ""), 
FileNames(2,""),ostr(2, nullptr), opened(2, false), active(false)
{
	bufs[0].reserve(buf1S); bufs[1].reserve(buf2S);
}

dualOfBufStream::~dualOfBufStream() {
	emptyStreams(true);
	for (size_t i = 0; i < ostr.size(); i++) { delete ostr[i]; }
}


#ifdef _gzipread2

std::vector< char > readline(gzFile f) {
	//  gzFile fp =gzopen(fname,"r");

	std::vector< char > v(2056);
	int pos = 0;
	for (;;) {
		if (gzgets(f, &v[pos], (int)v.size() - pos) == 0) {
			// end-of-file or error
			int err;
			//const char *msg = gzerror(f, &err);
			if (err != Z_OK) {
				// handle error
			}
			break;
		}
		unsigned read = strlen(&v[pos]);
		if (v[pos + read - 1] == '\n') {
			if (v[pos + read - 2] == '\r') {
				pos = pos + read - 2;
			}
			else {
				pos = pos + read - 1;
			}
			break;
		}
		if (read == 0 || pos + read < v.size() - 1) {
			pos = read + pos;
			break;
		}
		pos = v.size() - 1;
		v.resize(v.size() * 2);
	}
	v.resize(pos);
	return v;
}
#endif