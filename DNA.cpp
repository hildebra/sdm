#include "DNA.h"
#include "InputStream.h" // for helpers and globals (rtrim, parseInt, DNA_amb, NT_POS, DNA_IUPAC, itos, whoIsBetter)




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
	}
	else if (QSi == 0 && sequence_ == "") {
		this->setPassed(false);
		return true;
	}
	else if (QSi != sequence_.length()) {
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
DNA::DNA(const vector<string>& fas) :DNA() {
	// Abort if input stream is at the end
	if (fas[0].length() == 0) {
		return;
	}
	if (fas[0][0] != '>') {
		cerr << "ERROR: Line 1 in fasta file does not start with \">\" :\n" << fas[0] << endl;
		exit(23);
	}
	// Storage objects
	string tqual;
	this->setHeader(fas[0].substr(1));
	sequence_ = fas[1];
	sequence_length_ = sequence_.length();
	uint lsize = (uint)fas[1].length();
	qual_.assign(lsize, 11);
	if (fas[2].length() > 0) {
		tqual = fas[2];
		rtrim(tqual);
		const char* lQ = tqual.c_str();
		uint ii(0);
		qual_score nn(0);

		for (; ii < lsize; ii++) {
			nn = (qual_score)parseInt(&lQ);// , posStr);
			qual_[ii] = nn;
			if (*lQ == '\0') {
				break;
			}
			//issQ >>  Iqual[ii];
		}

		if (qual_.size() != fas[1].length()) {
			cerr << "Unequal fasta (" << fas[1].length() << ")/qual (" << fas[2].length() << ") length:\n";
			cerr << fas[0] << endl << fas[1] << endl << fas[2] << endl;
			exit(923);
		}
	}
	avg_qual_ = -1.f;

}

DNA::DNA(const vector<string>& fq, qual_score fastQver) :DNA() {
	//string line;
	if (fq.size() != 4 || fq[0].length() == 0) { delself(); return; }
	while (fq[0][0] != '@') {
		cerr << "ERROR on line " + fq[0] + ": Could not find \'@\' when expected (file likely corrupt, trying to recover):\n";// << endl;
		delself();
		return;

	}
	this->setHeader(fq[0].substr(1));
	sequence_ = fq[1];
	sequence_length_ = sequence_.size();
	if (fq[2][0] != '+') {
		cerr << "Error input line " + fq[2] + ": Could not find \'+\' when expected (file likely corrupt, aborting):\n" + fq[0];// << endl;
		delself();
		return;
		//if (!safeGetline(fna, line)) { delete tdn; return NULL; }
	}

	//qual_ score
	const string& line = fq[3];


	if (line.length() != this->mem_length()) {
		string sizdif = "More";
		if (line.length() > this->mem_length()) {
			sizdif = "Less";
		}
		//check that quality gets not more length than DNA
		delself();
		cerr << "Error input line " + fq[3] + "'\n" + fq[1] + "'\n: " + sizdif + " quality positions than nucleotides detected for sequence\n " +
			fq[0];// << endl;
		return;
	}

	uint qcnt(0); uint lline = (uint)this->mem_length();// line.length();
	qual_.resize(this->mem_length(), 0);

	for (; qcnt < lline; qcnt++) {
		qual_score q = (qual_score)line[qcnt] - fastQver;
		qual_[qcnt] = q;// minmaxQscore();
	}
	avg_qual_ = -1.f;
}

// Move-aware constructors (definitions)
DNA::DNA(vector<string>&& fq, qual_score fastQver) :DNA() {
	if (fq.size() != 4 || fq[0].length() == 0) { delself(); return; }
	while (fq[0][0] != '@') {
		cerr << "ERROR on line " + fq[0] + ": Could not find '@' when expected (file likely corrupt, trying to recover):\n";
		delself();
		return;

	}
	this->setHeader(std::move(fq[0]).substr(1));
	sequence_ = std::move(fq[1]);
	sequence_length_ = sequence_.size();
	if (fq[2][0] != '+') {
		cerr << "Error input line " + fq[2] + ": Could not find '+' when expected (file likely corrupt, aborting):\n" + fq[0];
		delself();
		return;
	}

	// qual_ score
	const string& line = fq[3];

	if (line.length() != this->mem_length()) {
		string sizdif = "More";
		if (line.length() > this->mem_length()) {
			sizdif = "Less";
		}
		delself();
		cerr << "Error input line " + fq[3] + "'\n" + fq[1] + "'\n: " + sizdif + " quality positions than nucleotides detected for sequence\n " + fq[0];
		return;
	}

	uint qcnt(0); uint lline = (uint)this->mem_length();
	qual_.resize(this->mem_length(), 0);
	for (; qcnt < lline; qcnt++) {
		qual_score q = (qual_score)line[qcnt] - fastQver;
		qual_[qcnt] = q;
	}
	avg_qual_ = -1.f;
}

DNA::DNA(vector<string>&& fas) :DNA() {
	if (fas[0].length() == 0) {
		return;
	}
	if (fas[0][0] != '>') {
		cerr << "ERROR: Line 1 in fasta file does not start with '>' :\n" << fas[0] << endl;
		exit(23);
	}
	this->setHeader(std::move(fas[0]).substr(1));
	sequence_ = std::move(fas[1]);
	sequence_length_ = sequence_.length();
	uint lsize = (uint)sequence_.length();
	qual_.assign(lsize, 11);
	if (fas[2].length() > 0) {
		string tqual = std::move(fas[2]);
		rtrim(tqual);
		const char* lQ = tqual.c_str();
		uint ii(0);
		qual_score nn(0);

		for (; ii < lsize; ii++) {
			nn = (qual_score)parseInt(&lQ);
			qual_[ii] = nn;
			if (*lQ == '\0') {
				break;
			}
		}

		if (qual_.size() != sequence_.length()) {
			cerr << "Unequal fasta (" << sequence_.length() << ")/qual (" << tqual.length() << ") length:\n";
			cerr << fas[0] << endl << sequence_ << endl << tqual << endl;
			exit(923);
		}
	}
	avg_qual_ = -1.f;
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
	string s = id_;
	remove_paired_info(s, read_position_);
	return s;
}

std::string_view DNA::getPositionFreeIdView() const {
	if (id_.empty()) return std::string_view();
	size_t len = id_.size();
	size_t sp = id_.find(' ');
	if (sp != std::string::npos) {
		len = sp;
	}
	if (len >= 2) {
		size_t tarPos = len - 2;
		if ((id_[tarPos] == '/' || id_[tarPos] == '.') && (id_[tarPos + 1] == '1' || id_[tarPos + 1] == '2')) {
			len = tarPos;
		}
	}
	return std::string_view(id_.data(), len);
}

void DNA::normalizePositionFreeId() {
	remove_paired_info(id_, read_position_);
	new_id_ = id_;
}

int DNA::numNonCanonicalDNA(bool all) {
	int DNAch = 0;
	size_t DNAl = length();
	if (all) {
		DNAl = sequence_.size();
	}
	for (size_t i = 0; i < DNAl; i++) {
		DNAch += DNA_amb[sequence_[i]];
	}
	return DNAch;
}
int DNA::numACGT() {
	int DNAch = 0;
	uint maxL = length();
	for (unsigned int i = 0; i < length(); i++) {
		DNAch += DNA_amb[(int)sequence_[i]];
	}
	return DNAch;
}
void DNA::stripLeadEndN() {
	int i = 0;
	while (i >= 0 && DNA_amb[(int)sequence_[i]]) {
		i++;
	}
	if (i) {
		this->cutSeq(0, i);
	}
	//cut at end N's
	i = (int)sequence_.length();
	while (i >= 0 && DNA_amb[(int)sequence_[i]]) {
		i--;
	}
	if (i != (int)sequence_.length()) {
		this->cutSeq(i + 1, (int)sequence_.length());
	}

}
void DNA::appendQuality(const vector<qual_score>& q)
{
	qual_.insert(qual_.end(), q.begin(), q.end());
	avg_qual_ = -1.f;
}

float DNA::getAvgQual() {
	if (avg_qual_ < 0.f) {
		if (quality_sum_ == 0) {
			for (unsigned int i = 0; i < this->length(); i++) {
				quality_sum_ += qual_[i];
			}
		}
		avg_qual_ = static_cast<float>(quality_sum_) / static_cast<float>(this->length());

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
	if (n == 1.0f) {
		if (j <= 0.0f) {
			return 1.0f - p;
		}
		if (j == 1.0f) {
			return p;
		}
		return 0.0f;
	}

	float per_position_accum_probs;
	if (j > n) {
		return 0.0f; //The formula below would also return 0.
	}
	per_position_accum_probs = pow((1 - p), n);	//For j == 0.
	float i(1);
	for (; i <= j; i += 1.f) {//For j > 0.
		per_position_accum_probs = ((n - i + 1.f) / (1.0f * i)) * (p / (1.f - p)) * per_position_accum_probs;
	}
	return per_position_accum_probs;

}
float DNA::sum_of_binomials(int j, int k, float n, int qual_length, const vector<float>& error_probs, const vector<float>& per_position_accum_probs)
//#Where j is the number of errors and k is the position in the sequence.
{
	float probability = 0.f;
	const int k1 = k - 1;
	const float pk = error_probs[k];

	for (int i = 0; i <= j; ++i)
	{
		const int rowOffset = (j - i) * qual_length;
		probability += DNA::prob_j_errors(pk, static_cast<float>(i), n) * per_position_accum_probs[rowOffset + k1];
		//Where error_probs[k] is the error probability of the k-th position.
		//Where per_position_accum_probs[j-i][k-1] is the probability that all the bases from position 0 to k-1 had a total of j-i errors.
	}

	return probability;
}

float DNA::binomialFilter(int maxErr, float alpha) {

	if (alpha == -1.f || this->length() < 3) { return 0.f; }//deactivated
	if (maxErr > 1e5) { return 0.f; } // impossibly large..

	///Initialize some variables.

	int SeqLengthI = (int)sequence_length_;
	int n = 1; //Since we have a Bernoulli distribution.
	float n_f = (float)n;
	float alpha1 = 1.f - alpha;

	///Translate quality scores into error probabilities.
	vector<float> error_probs(sequence_length_, 1.f);
	for (size_t i = 0; i < sequence_length_; i++) {
		if (qual_[i]<0 || qual_[i]>maxSAqualP) { continue; }
		error_probs[i] = (float)SAqualP[qual_[i]];//pow(10, (contig_quals[i] / -10.0)); //Since we want a continuous list of non-N error_probs.
	}

	///Actually get the job done.
	int max_expected_errors = maxErr + 3;// (int)sequence_length_ + 1;
	int expected_errors = 0;
	float probability;
	vector<float> accumulated_probs(max_expected_errors, 0.f);
	//int j;
	int k;
	vector<float> per_position_accum_probs(max_expected_errors * SeqLengthI, 0.f);

	while (1)
	{
		float expected_errors_f = (float)expected_errors;
		const int expectedRowOffset = expected_errors * SeqLengthI;
		//vector<float > per_position_accum_probs(sequence_length_, 0.f);
		for (k = 0; k < (int)sequence_length_; k++) {
			if (k == 0) {
				per_position_accum_probs[expectedRowOffset + k] = DNA::prob_j_errors(error_probs[k], expected_errors_f, n_f);
			}
			else {

				per_position_accum_probs[expectedRowOffset + k] = DNA::sum_of_binomials(expected_errors, k, n_f, SeqLengthI, error_probs, per_position_accum_probs);
			}
		}
		probability = per_position_accum_probs[expectedRowOffset + (SeqLengthI - 1)];

		if (expected_errors == 0) {
			accumulated_probs[expected_errors] = probability;
		}
		else {
			accumulated_probs[expected_errors] = accumulated_probs[expected_errors - 1] + probability;
		}

		if (accumulated_probs[expected_errors] > (alpha1) || expected_errors >= (max_expected_errors - 1)) {
			break;
		}
		else {
			expected_errors++;
		}
	}
	if (expected_errors == 0) {
		return 0.f;
	}
	float EXE = interpolate(expected_errors - 1, accumulated_probs[expected_errors - 1], expected_errors, accumulated_probs[expected_errors], alpha);
	return EXE;
}

float DNA::qualWinfloat(uint W, float T, int& reason) {
	//if (T==0.f){return true;}
	int AQS = 0;
	int TotQ = 0;
	int upTs = int(T) * int(W);
	uint QS = this->length();//static_cast<unsigned int> (qual_.size());
	if (W >= QS) { W = QS; } // too short

	//1st loop to ini window
	for (size_t i = 0; i < size_t(W); i++) {
		AQS += qual_[i];
	}
	TotQ = AQS;

	for (size_t i = W; i < (size_t)QS; i++) {
		AQS += qual_[i] - qual_[i - W];
		TotQ += qual_[i];
		if (AQS < upTs) {
			reason = 1;	return 0.f;
		}
	}
	//if (averageQ > static_cast<float> (TotQ) /static_cast<float> ( qual_.size())){return false;}
	return float(TotQ) / float(QS);
}

int DNA::qualAccumulate(double d) {
	unsigned int i(0); double accErr(0.0);
	for (; i < this->length(); i++) {
		accErr += SAqualP[qual_[i]];
		if (accErr >= d) { break; }
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
		qsc.resize(sql, 0);
	}
	if (ntcnt.size() < sql) {
		ntcnt.resize(sql, 0);
	}
	for (uint i = 0; i < sql; i++) {
		short p = NT_POS[(int)sequence_[i]];
		qsc[p] += qual_[i];
		ntcnt[p]++;
	}
}

void DNA::ntSpecQualScoresMT(vector<long>& qsc, vector<long>& ntcnt) {
	size_t sql = sequence_.length();

	if (qsc.size() < sql) {
		qsc.resize(sql, 0);
	}
	if (ntcnt.size() < sql) {
		ntcnt.resize(sql, 0);
	}
	for (uint i = 0; i < sql; i++) {
		short p = NT_POS[(int)sequence_[i]];
		qsc[p] += qual_[i];
		ntcnt[p]++;
	}
}

bool DNA::qualAccumTrim(double d) {
	if (d == -1.) {
		return true;
	}
	unsigned int i(qualAccumulate(d));
	if (i != this->length()) {
		//cut 3' end
		this->cutSeqPseudo(i);
		this->QualCtrl.AccErrTrimmed = true;
		return false;
	}
	//did not cut this sequence:
	return true;
}
bool DNA::qualWinPos(unsigned int W, float T) {
	if (T == 0.f) { return true; }
	int AQ = 0;
	int unT = static_cast<int>((float)W * T);
	unsigned int QS = this->length();// (unsigned int)qual_.size();
	unsigned int QSh = QS >> 2;
	QSh = max(QSh, W);
	if (W >= QS) { return true; } // too short

	vector<float> WQ((int)W);
	//TODO: check that the right num of nucs is taken..
	int cnt = 0;
	if (QS < W) { return true; }
	for (unsigned int i = QS - 1; i > QS - (unsigned int)W - 1; i--) {
		AQ += qual_[i];
		cnt++;
		if (AQ > unT) { return true; }
	}
	int curW = QS - (unsigned int)W;
	for (uint i = QS - (unsigned int)W - 1; i > QSh; i--) {
		AQ += qual_[i];
		AQ -= qual_[i + W];

		if (AQ < unT) { //min Window qual_ was broken.. kick seq
			curW = i;
		}
		else {
			break;
		}
	}

	//partial seq  removal
	int pos = curW - (W >> 1);
	if (pos < 0) { pos = 0; }
	this->cutSeqPseudo(pos);
	this->QualCtrl.QWinTrimmed = true;
	return false;
}


shared_ptr<DNA> DNA::getDNAsubseq(int start, int end, string id) {
	shared_ptr<DNA> rdn = make_shared<DNA>();
	rdn->setSequence(this->getSequence().substr(start, (end - start)));
	vector<qual_score> tq = this->getQual(start, end);
	rdn->setQual(tq);
	rdn->setHeader(id);
	if (this->getBCnumber() >= 0) {
		rdn->setBCnumber(this->getBCnumber(), 0);
		if (getBarcodeCut()) {
			rdn->setBarcodeCut();
		}
		rdn->setBarcodeDetected(getBarcodeDetected());

	}
	return rdn;
}

//removes part of seq and qual_ indexes between start & stop
//set stop to -2 if cut till end in pseudo == true mode
bool DNA::cutSeq(int start, int stop, bool pseudo) {

	if (stop == -1) {
		if (start >= (int)sequence_length_ || start < 0) { return false; }
	}
	else if (start >= stop || start >= (int)qual_.size()) {
		return false;
	}
	if (stop > (int)qual_.size()) {
		stop = (int)qual_.size();
	}

	//pseudo deactivates cutting of 3'
	if (pseudo) {
		if (stop == -1) {//|| stop <= (int) sequence_length_) {
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
	}
	else {
		sequence_ = se.substr(0, start) + se.substr(stop);
	}

	//DN = sequence_.c_str();
	//Quali
	qual_.erase(qual_.begin() + start, qual_.begin() + stop);

	sequence_length_ = sequence_.length() - seqLDiff;

	return true;
}

int DNA::matchSeq(std::string PrSt, int Err, int searchSpace, int startPos, bool exhaustive) {
	//const char* DN = sequence_.c_str();
	//const char* Pr = PrSt.c_str();
	int PrL = (int)PrSt.length();
	int mthL = this->length() - PrL;
	//int wantSc = PrL - Err;
	int endPos(-1), pos(startPos), Prp(0), c_err(0), Prp2(0), c_points(0), point_aim(1);
	point_aim = max(1, int(PrL / 2));
	//bool res(false);
	vector<int> potentialMatches(0);
	for (; pos < searchSpace; pos++) {
		if (pos > mthL) { break; }
		c_err = 0; Prp = 0; Prp2 = pos; c_points = 0;
		do {

#ifdef _NEWMATCH
			//new vector based matching
			c_err += DNA_IUPAC[sequence_[Prp2] + 256 * PrSt[Prp]];
			if (c_err > Err) { break; }
#else
			//old, direct match
			if (!matchDNA(sequence_[Prp2], PrSt[Prp])) { c_err++; if (c_err > Err) { break; } }
#endif
			Prp++; Prp2++;
			if (sequence_[Prp2] != 'N') { c_points++; }
		} while (Prp < PrL);
		/*if (c_points >= point_aim) {
			int x = 0; //DEBUG
		}*/
		if (c_err <= Err && c_points >= point_aim) {
			endPos = pos;
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
	}
	else {
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
int DNA::matchSeqRev(const string& PrSt, int Err, int searchSpace,
	int start, bool exhaustive) {
	//fail::ret -1
	int PrL = (int)PrSt.length();
	if (start == 0) {
		start = PrL; //default seed set to 5
	}
	int SeL = (int)sequence_.size();
	//int wantSc = PrL - Err;
	int pos(SeL - start), Prp(0), c_err(0), endPos(-1), c_points(0), point_aim(1);
	vector<int> potentialMatches(0);
	for (; pos > searchSpace; pos--) {
		c_err = 0; Prp = 0; c_points = 0;
		int PrL2 = min(PrL, SeL - pos);
		point_aim = max(1, int((float)PrL2 * 0.9f));
		do {
#ifdef _NEWMATCH
			uint lpos = pos + Prp;
			c_err += DNA_IUPAC[sequence_[lpos] + 256 * PrSt[Prp]]; if (c_err > Err) { break; }
#else
			if (!matchDNA(sequence_[pos + Prp], PrSt[Prp])) { c_err++; if (c_err > Err) { break; } }
#endif
			Prp++;
			if (sequence_[lpos] != 'N') { c_points++; }
		} while (Prp < PrL2);
		if (c_err <= Err && c_points >= point_aim) {
			endPos = pos;
			potentialMatches.push_back(endPos);
			if (!exhaustive) {
				break;
			}
		}
	}
	//secondary check for last few NT's
	if (endPos == -2 && potentialMatches.size() == 0) {
		pos = (SeL - 1);
		for (; pos > (SeL - start); pos--) {
			c_err = 0; Prp = 0;
			int PrL2 = min(PrL, SeL - pos);
			do {
#ifdef _NEWMATCH
				c_err += DNA_IUPAC[sequence_[pos + Prp] + 256 * PrSt[Prp]]; if (c_err > Err) { break; }
#else
				if (!matchDNA(sequence_[pos + Prp], PrSt[Prp])) { c_err++;	if (c_err > Err) { break; } }
#endif
				Prp++;
			} while (Prp < PrL2);
			if (c_err <= Err) {
				endPos = pos;
				potentialMatches.push_back(endPos);
				if (!exhaustive) {
					break;
				}
			}
		}
	}
	if (potentialMatches.size() > 0) {
		return potentialMatches.back();
	}
	else {
		return -1;
	}
}
// looks through total DNA seq
int DNA::matchSeq_tot(const string& Pr, int error, int maxPos, int& c_err) {
	int PrL = (int)Pr.length();
	int pos(0), Prp(0), Prp2(0);
	bool success(false);
	for (pos = 0; pos < maxPos; pos++) {
		c_err = 0; Prp = 0; Prp2 = pos;
		do {
			if (sequence_[Prp2] != Pr[Prp]) {
				c_err++;
			}
			Prp++; Prp2++;
		} while (c_err <= error && Prp < PrL);
		if (c_err <= error) { success = true; break; }
	}
	if (!success) { pos = -1; }
	return pos;
}


bool DNA::matchDNA(char t1, char t2) {
	if (t1 == t2) {
		return true;
	}
	switch (t2) {
	case 'N': return true;
	case 'R': if (t1 == 'A' || t1 == 'G') { return true; }break;
	case 'Y': if (t1 == 'T' || t1 == 'C') { return true; }break;
	case 'M': if (t1 == 'C' || t1 == 'A') { return true; }break;
	case 'K': if (t1 == 'T' || t1 == 'G') { return true; }break;
	case 'W': if (t1 == 'T' || t1 == 'A') { return true; }break;
	case 'S': if (t1 == 'C' || t1 == 'G') { return true; }break;
	case 'B': if (t1 != 'A') { return true; }break;
	case 'D': if (t1 != 'C') { return true; }break;
	case 'H': if (t1 != 'G') { return true; }break;
	case 'V': if (t1 != 'T') { return true; }break;
	}
	return false;
}
bool DNA::HomoNTRuns(int max) {
	char lastC = sequence_[0];
	int rowC = 1;
	for (unsigned int i = 1; i < length(); i++) {
		if (sequence_[i] == lastC) {
			rowC++;
			if (rowC >= max) {
				return false;
			}
		}
		else {
			rowC = 1;
			lastC = sequence_[i];
		}
	}
	return true;
}

int DNA::HomoNTTrim(int max) {
	// Return the sequence length if trimming takes place, else return 0.
	if (max >= (int)length()) { return 0; }
	unsigned int i = length() - 1;
	char lastC = sequence_[i];
	int rowC = 0;
	for (; i > 0; i--) {
		if (sequence_[i] == lastC) {
			rowC++;
		}
		else {
			break;
		}
	}

	if (rowC >= max) {
		return i + 1;
	}
	rowC = 0;
	//lastC = sequence_[i - 1];
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
size_t DNA::getSpaceHeadPos(const string& x) {
	size_t pos = x.find_first_of(" \t"); // x.find(' ');
	if (pos == string::npos) { pos = x.length(); }
	return pos;
}
size_t DNA::getShorterHeadPos(const string& x, int fastQheadVer) {

	size_t pos(string::npos);
	size_t strL = x.length();
	if (fastQheadVer != 0) {
		if (read_position_ == 1) {
			pos = x.find("/2", strL - 3);
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
string DNA::xtraHdStr() {
	string xtr = " ";
	if (FtsDetected.forward) { xtr += "F"; }
	if (FtsDetected.reverse) { xtr += "R"; }
	return xtr;
}
void DNA::writeSeq(ostream& wr, bool singleLine) {
	wr << writeSeq(singleLine);
}
void DNA::writeSeq(string& str, bool singleLine) {
	if (sequence_.size() == 0) { return; }
	str = ">" + new_id_ + "\n";
	if (singleLine) {
		str += sequence_.substr(0, length()) + "\n";
		return;
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
string DNA::writeSeq(bool singleLine) {
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
	int cnt = 0;
	if (qual_.size() == 0) { return; }
	str += ">" + new_id_ + "\n";
	//wr<< qual_traf_<<endl;
	string endlChr("\n");
	if (singleLine) { endlChr = " "; }
	for (unsigned int i = 0; i < length(); i++) {
		cnt++;
		if (cnt < 80) {
			str += itos((int)qual_[i]) + " ";
		}
		else {
			str += itos((int)qual_[i]) + endlChr;
			cnt = 0;
		}
	}
	if (cnt > 0) {
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
	if (qual_.size() == 0 || sequence_.length() == 0 || length() == 0) { return; }
	const size_t seqLen = this->length();
	const string& hdr = newHD ? new_id_ : id_;
	ret.reserve(ret.size() + hdr.size() + seqLen + qual_traf_.size() + 5);

	ret.push_back('@');
	ret += hdr;
	ret.push_back('\n');
	ret.append(sequence_.data(), seqLen);
	ret += "\n+\n";
	ret.append(qual_traf_.data(), qual_traf_.size());
	ret.push_back('\n');
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

void DNA::changeHeadPver(int ver) {
	string& oID(id_);
	if (id_fixed_) {
		oID = new_id_;
	}

	if (ver == 1) { //change from XX 1: to XX/1
		size_t pos = oID.find_first_of(" ");
		//unsigned int pos2 = id_.find(" ",pos);
		if (pos != string::npos) {
			new_id_ = oID.substr(0, pos) + "/" + oID.substr(pos + 1, 1);
			//if (pos2 != string::npos){new_id_ += id_.substr(pos2);}
		}
	}
	else if (ver == 2) {
		cerr << "Head change Not implemented\n"; exit(23);
	}
	id_fixed_ = true;
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



void DNA::setPassed(bool b) {
	//once mid qual is active, there is no more good qual possible
	if (failed_) {
		good_quality_ = false; mid_quality_ = false;
	}
	else if (mid_quality_) {
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


//////////////////////////////////////////////////////////////////
//  DNAunique class
//////////////////////////////////////////////////////////////////



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
	int sample_id, bool b_derep_as_fasta_) {
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
			tmp, false)) { //this->getBestSeedLength(),
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
	}
	else {
		sample_counter->second++;
	}
#endif
}
void DNAunique::incrementSampleCounterBy(int sample_id, long count) {
	auto sample_counter = occurence.find(sample_id);
	if (sample_counter == occurence.end()) {
		occurence[sample_id] = count;
	}
	else {
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

bool sortDescending(int i, int j) { return (i > j); }//descending sort
vector<int>  DNAunique::getDerepMapSort(size_t wh) {
	vector<int> vals;
	size_t siz = occurence.size();
	if (wh > siz) { wh = siz; }
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
		if (ref != -1 && iter->second >= ref) {
			return true;
		}
	}
	return false;

}


void DNAunique::transferOccurence(const shared_ptr<DNAunique> dna_unique) {

	if (occurence.size() == 0) {
		occurence = dna_unique->occurence;
		//count_ = dna_unique->count_;
	}
	else {
		//which sample contains this dna?
		read_occ dereplication_map = dna_unique->getDerepMap();

		for (auto sample_iterator = dereplication_map.begin(); sample_iterator != dereplication_map.end(); ++sample_iterator) {
			auto sample_counter = occurence.find(sample_iterator->first);
			if (sample_counter == occurence.end()) {
				occurence[sample_iterator->first] = sample_iterator->second;
			}
			else {
				sample_counter->second += sample_iterator->second;
			}
		}
		//size track
		//count_ += dna_unique->count_;
	}
}


void DNAunique::writeMap(ofstream& os, const string& hd, vector<int>& counts_per_sample, const vector<int>& combiID) {
	if (occurence.empty()) { return; }
	int total_count(0);
	read_occ sample_counters;

	if (combiID.size() > 0) {//combine all counts on combined categories

		for (read_occ::iterator iter = occurence.begin(); iter != occurence.end(); ++iter) {
			//aim: sample_counters[combiID[iter->first]] += iter->second;
			read_occ::iterator sample_counter = sample_counters.find(combiID[iter->first]);
			if (sample_counter == sample_counters.end()) {
				// if sample_counter exists
				sample_counters[combiID[iter->first]] = iter->second;
			}
			else {
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
	//unordered_map<int, int> occurence;
	//occ = occurence;

	if (occurence.empty()) {
		return;
	}
	if (usFmt) {
		new_id_ += ";size=" + itos(totalSum()) + ";";
	}
	else {
		new_id_ += "_" + itos(totalSum());
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

	avg_qual_ = better_dna->avg_qual_;
	accumulated_error_ = better_dna->accumulated_error_;
	FtsDetected = better_dna->FtsDetected;
	good_quality_ = better_dna->good_quality_;
	//new_id_ = better_dna->id_;
	id_ = better_dna->id_;
	id_fixed_ = better_dna->id_fixed_;
	new_id_ = better_dna->new_id_;
	mid_quality_ = better_dna->mid_quality_;
	merge_offset_ = better_dna->merge_offset_;
	merge_seed_pos_ = better_dna->merge_seed_pos_;
	reversed_merge_ = better_dna->reversed_merge_;
	seed_length_ = better_dna->seed_length_;
	quality_sum_ = better_dna->quality_sum_;
	QualCtrl = better_dna->QualCtrl;
	qual_ = better_dna->qual_;
	qual_traf_ = better_dna->qual_traf_;
	sequence_ = better_dna->sequence_;
	sequence_length_ = better_dna->sequence_length_;
	sample_id_ = better_dna->sample_id_;
	reversed_ = better_dna->reversed_;
	read_position_ = better_dna->read_position_;
	good_quality_ = better_dna->good_quality_;
	accumulated_error_ = better_dna->accumulated_error_;
	tempFloat = better_dna->tempFloat;


	this->saveMem();
}

uint64_t* DNAunique::transferPerBaseQualitySum() {
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
		int newbase = (int)round((double)quality_sum_per_base_[i] / (double)totalSum()) + ofastQver;
		qualities_avg_[i] = char(newbase);
	}
	qualities_avg_[i] = '\0';

}


void DNAunique::writeDerepFastQ(ofstream& wr, bool newHD) {
	if (!(qual_.size() == 0 || sequence_.length() == 0 || length() == 0)) {
		string xtr = xtraHdStr();
		if (FtsDetected.forward) { xtr += "F"; }
		if (FtsDetected.reverse) { xtr += "R"; }
		if (newHD) {
			wr << "@" + new_id_ + "\n";
		}
		else {
			wr << "@" + id_ + "\n";
		}
		wr << sequence_.substr(0, length()) + "\n";
		wr << "+\n";//new_id_<<endl;
		wr << qualities_avg_ + "\n";
	}
}


shared_ptr<DNA> str2DNA(vector<string>& in, bool keepPairHD, int fastQver, int readpos) {
	shared_ptr<DNA> ret(nullptr);
	bool isFasta(false);
	if (in.size() == 3) {// || in[0][0] == '>'  ) {//fasta
		ret = make_shared<DNA>(in);
		isFasta = true;
	}
	else {
		ret = make_shared<DNA>(in, fastQver);
	}
	if (!keepPairHD) {//better cut in early
		ret->normalizePositionFreeId();
	}
	if (readpos == 0) {
		ret->setpairFWD();
	}
	else if (readpos == 1) {
		ret->setpairREV();
	}
	else if (readpos == 2) {
		ret->setMIDseq(true);
	}
	if (!ret->seal(isFasta) || ret->isEmpty()) {
		ret = nullptr;
	}
	return ret;
}

// rvalue overload - move strings into DNA to avoid copies when possible
shared_ptr<DNA> str2DNA(vector<string>&& in, bool keepPairHD, int fastQver, int readpos) {
	shared_ptr<DNA> ret(nullptr);
	bool isFasta(false);
	if (in.size() == 3) {
		ret = make_shared<DNA>(std::move(in));
		isFasta = true;
	}
	else {
		ret = make_shared<DNA>(std::move(in), fastQver);
	}
	if (!keepPairHD) {
		ret->normalizePositionFreeId();
	}
	if (readpos == 0) {
		ret->setpairFWD();
	}
	else if (readpos == 1) {
		ret->setpairREV();
	}
	else if (readpos == 2) {
		ret->setMIDseq(true);
	}
	if (!ret->seal(isFasta) || ret->isEmpty()) {
		ret = nullptr;
	}
	return ret;
}



