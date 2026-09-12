#include "DNA.h"
#include "InputStream.h" // for helpers and globals (rtrim, parseInt, DNA_amb, NT_POS, DNA_IUPAC, itos, whoIsBetter)
#include "containers.h"

namespace {
const int MAX_FASTQ_HEADER_WARNINGS = 20;
int g_fastqHeaderWarnings = 0;
}

namespace {
	size_t getFastqRecordLength(size_t logicalLength, size_t sequenceSize, size_t qualitySize) {
		return std::min(logicalLength, std::min(sequenceSize, qualitySize));
	}

	inline char sanitizeFastqQualityChar(int qWithOffset, int offset) {
		const int minFastqChar = std::max(33, offset);
		if (qWithOffset < minFastqChar) {
			qWithOffset = minFastqChar;
		}
		if (qWithOffset > 126) {
			qWithOffset = 126;
		}
		return static_cast<char>(qWithOffset);
	}

	void buildFastqQualityString(std::string& out, const std::vector<qual_score>& quals, size_t len, qual_score fastqOffset) {
		out.resize(len);
		const int offset = static_cast<int>(fastqOffset);
		for (size_t i = 0; i < len; ++i) {
			const int qWithOffset = static_cast<int>(quals[i]) + offset;
			out[i] = sanitizeFastqQualityChar(qWithOffset, offset);
		}
	}

	std::vector<qual_score> parseFastaQualities(const std::string& text, size_t expectedCount,
		const std::string& header) {
		std::istringstream input(text);
		std::vector<qual_score> qualities;
		qualities.reserve(expectedCount);
		int value = 0;
		while (input >> value) {
			if (value < 0 || value > maxSAqualP) {
				throw std::runtime_error("FASTA/QUAL score outside supported range for " + header);
			}
			qualities.push_back(static_cast<qual_score>(value));
		}
		if (!input.eof()) {
			throw std::runtime_error("Invalid token in FASTA/QUAL scores for " + header);
		}
		if (qualities.size() != expectedCount) {
			throw std::runtime_error("FASTA/QUAL length mismatch for " + header + ": sequence has " +
				std::to_string(expectedCount) + " bases but QUAL has " +
				std::to_string(qualities.size()) + " scores");
		}
		return qualities;
	}
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
	// uppercase DNA only when needed (avoids unconditional per-char writes)
	bool needsUppercase = false;
	for (const char c : sequence_) {
		if (c >= 'a' && c <= 'z') {
			needsUppercase = true;
			break;
		}
	}
	if (needsUppercase) {
		for (char& c : sequence_) {
			if (c >= 'a' && c <= 'z') {
				c = static_cast<char>(c - ('a' - 'A'));
			}
		}
	}
	sequence_length_ = sequence_.length();
	this->fixQ0();
	return true;
}
DNA::DNA(const vector<string>& fas) :DNA() {
	// Abort if input stream is at the end
	if (fas.size() < 2 || fas[0].empty()) {
		delself();
		return;
	}
	if (fas[0][0] != '>') {
		cerr << "ERROR: Line 1 in fasta file does not start with \">\" :\n" << fas[0] << endl;
		exit(23);
	}
	// Storage objects
	this->setHeader(fas[0].substr(1));
	sequence_ = fas[1];
	sequence_length_ = sequence_.length();
	const size_t lsize = sequence_.length();
	if (fas.size() > 2 && !fas[2].empty()) {
		qual_ = parseFastaQualities(fas[2], lsize, fas[0]);
	}
	else {
		qual_.assign(lsize, 11);
	}
	avg_qual_ = -1.f;

}

DNA::DNA(const vector<string>& fq, qual_score fastQver) :DNA() {
    if (fq.size() == 5) { sample_barcode_ = fq[4]; fastQver = 33; }
	//string line;
	if ((fq.size() != 4 && fq.size() != 5) || fq[0].length() == 0) { delself(); return; }
	if (fq[0][0] != '@') {
		throw std::runtime_error("Malformed FASTQ header: " + fq[0]);
	}
	this->setHeader(fq[0].substr(1));
	sequence_ = fq[1];
	sequence_length_ = sequence_.size();
	if (fq[2].empty() || fq[2][0] != '+') {
		throw std::runtime_error("Malformed FASTQ plus line for " + fq[0]);
	}

	//qual_ score
	const string& line = fq[3];


	if (line.length() != this->mem_length()) {
		throw std::runtime_error("FASTQ sequence/quality length mismatch for " + fq[0]);
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
    if (fq.size() == 5) { sample_barcode_ = std::move(fq[4]); fastQver = 33; }
	if ((fq.size() != 4 && fq.size() != 5) || fq[0].length() == 0) { delself(); return; }
	if (fq[0][0] != '@') {
		throw std::runtime_error("Malformed FASTQ header: " + fq[0]);
	}
	const string headerForError = fq[0];
	fq[0].erase(0, 1);
	this->setHeader(std::move(fq[0]));
	sequence_ = std::move(fq[1]);
	sequence_length_ = sequence_.size();
	if (fq[2].empty() || fq[2][0] != '+') {
		throw std::runtime_error("Malformed FASTQ plus line for " + headerForError);
	}

	// qual_ score
	const string& line = fq[3];

	if (line.length() != this->mem_length()) {
		throw std::runtime_error("FASTQ sequence/quality length mismatch for " + headerForError);
	}

	const size_t lline = this->mem_length();
	qual_.resize(lline, 0);
	for (size_t qcnt = 0; qcnt < lline; ++qcnt) {
		qual_[qcnt] = static_cast<qual_score>(line[qcnt] - fastQver);
	}
	avg_qual_ = -1.f;
}

DNA::DNA(vector<string>&& fas) :DNA() {
	if (fas.size() < 2 || fas[0].empty()) {
		delself();
		return;
	}
	if (fas[0][0] != '>') {
		cerr << "ERROR: Line 1 in fasta file does not start with '>' :\n" << fas[0] << endl;
		exit(23);
	}
	fas[0].erase(0, 1);
	this->setHeader(std::move(fas[0]));
	sequence_ = std::move(fas[1]);
	sequence_length_ = sequence_.length();
	const size_t lsize = sequence_.length();
	if (fas.size() > 2 && !fas[2].empty()) {
		string tqual = std::move(fas[2]);
		qual_ = parseFastaQualities(tqual, lsize, id_);
	}
	else {
		qual_.assign(lsize, 11);
	}
	avg_qual_ = -1.f;
}


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
 if (len >= 4) {
		size_t tarPos = len - 4;
		if (id_[tarPos] == '#' && (id_[tarPos + 1] == '1' || id_[tarPos + 1] == '2') && id_[tarPos + 2] == ':' && id_[tarPos + 3] == '0') {
			len = tarPos;
		}
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
		DNAch += DNA_amb[static_cast<unsigned char>(sequence_[i])];
	}
	return DNAch;
}
int DNA::numACGT() {
	int DNAch = 0;
	const uint len = length();
	for (unsigned int i = 0; i < len; i++) {
		DNAch += DNA_amb[static_cast<unsigned char>(sequence_[i])];
	}
	return DNAch;
}
bool DNA::scanSequenceChecks(int maxHomoNT, int& ambNTs_out) {
	const uint len = length();
	ambNTs_out = 0;
	if (len == 0) { return true; }
	char lastC = sequence_[0];
	int rowC = 1;
	ambNTs_out += DNA_amb[static_cast<unsigned char>(lastC)];
	for (unsigned int i = 1; i < len; i++) {
		const char c = sequence_[i];
		ambNTs_out += DNA_amb[static_cast<unsigned char>(c)];
		if (c == lastC) {
			rowC++;
			if (maxHomoNT != 0 && rowC >= maxHomoNT) {
				return false;
			}
		}
		else {
			rowC = 1;
			lastC = c;
		}
	}
	return true;
}
void DNA::stripLeadEndN() {
	size_t first = 0;
	while (first < sequence_.size() && DNA_amb[static_cast<unsigned char>(sequence_[first])]) {
		++first;
	}
	size_t last = sequence_.size();
	while (last > first && DNA_amb[static_cast<unsigned char>(sequence_[last - 1])]) {
		--last;
	}
	if (last < sequence_.size()) {
		this->cutSeq(static_cast<int>(last), static_cast<int>(sequence_.size()));
	}
	if (first > 0) {
		this->cutSeq(0, static_cast<int>(first));
	}

}
void DNA::appendQuality(const vector<qual_score>& q)
{
	qual_.insert(qual_.end(), q.begin(), q.end());
	invalidateQualityCaches();
}

float DNA::getAvgQual() {
	if (avg_qual_ < 0.f) {
		quality_sum_ = 0;
		const size_t logicalLength = std::min(sequence_length_, qual_.size());
		for (size_t i = 0; i < logicalLength; ++i) {
			quality_sum_ += static_cast<uint64_t>(std::max<int>(0, qual_[i]));
		}
		if (logicalLength == 0) {
			avg_qual_ = 0.f;
		}
		else {
			avg_qual_ = static_cast<float>(quality_sum_) / static_cast<float>(logicalLength);
		}
	}
	return avg_qual_;
}
int DNA::getMedianQual() {
	const size_t logicalLength = std::min(sequence_length_, qual_.size());
	if (logicalLength == 0) { return 0; }
	vector<qual_score> logicalQualities(qual_.begin(), qual_.begin() + logicalLength);
	return median(logicalQualities);
}

double DNA::getAccumError() {
	if (accumulated_error_ == 0.) {
		const size_t logicalLength = std::min(sequence_length_, qual_.size());
		for (size_t i = 0; i < logicalLength; ++i) {
			const int boundedQuality = std::max(0, std::min(maxSAqualP, static_cast<int>(qual_[i])));
			accumulated_error_ += SAqualP[boundedQuality];
		}
	}
	if (!std::isfinite(accumulated_error_)) {
		accumulated_error_ = 5.;
	}
	return accumulated_error_;
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

	const int SeqLengthI = (int)sequence_length_;
	constexpr float n_f = 1.f; //Since we have a Bernoulli distribution.
	const float alpha1 = 1.f - alpha;

	///Translate quality scores into error probabilities.
	thread_local vector<float> error_probs;
	error_probs.assign(sequence_length_, 1.f);
	for (size_t i = 0; i < sequence_length_; i++) {
		if (qual_[i]<0 || qual_[i]>maxSAqualP) { continue; }
		error_probs[i] = (float)SAqualP[qual_[i]];
	}

	///Actually get the job done.
	const int max_expected_errors = maxErr + 3;
	int expected_errors = 0;
	float probability;
	thread_local vector<float> accumulated_probs;
	accumulated_probs.assign(max_expected_errors, 0.f);
	const size_t ppSize = (size_t)max_expected_errors * SeqLengthI;
	thread_local vector<float> per_position_accum_probs;
	per_position_accum_probs.assign(ppSize, 0.f);

	while (1)
	{
		const float expected_errors_f = (float)expected_errors;
		const int expectedRowOffset = expected_errors * SeqLengthI;
		for (int k = 0; k < (int)sequence_length_; k++) {
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
	int AQS = 0;
	int TotQ = 0;
	const uint QS = static_cast<uint>(std::min(sequence_length_, qual_.size()));
	if (QS == 0 || qual_.empty()) { return 0.f; }
	if (W >= QS) { W = QS; } // too short
	if (W == 0) { return 0.f; }
	const float threshold = T * static_cast<float>(W);

	//1st loop to ini window
	for (size_t i = 0; i < size_t(W); i++) {
		AQS += qual_[i];
	}
	TotQ = AQS;
	if (static_cast<float>(AQS) < threshold) {
		reason = 1;
		return 0.f;
	}

	for (size_t i = W; i < (size_t)QS; i++) {
		AQS += qual_[i] - qual_[i - W];
		TotQ += qual_[i];
		if (static_cast<float>(AQS) < threshold) {
			reason = 1;	return 0.f;
		}
	}
	//if (averageQ > static_cast<float> (TotQ) /static_cast<float> ( qual_.size())){return false;}
	return float(TotQ) / float(QS);
}

int DNA::qualAccumulate(double d) {
	unsigned int i(0); double accErr(0.0);
	const size_t logicalLength = std::min(sequence_length_, qual_.size());
	for (; i < logicalLength; i++) {
		const int boundedQuality = std::max(0, std::min(maxSAqualP, static_cast<int>(qual_[i])));
		accErr += SAqualP[boundedQuality];
		if (accErr >= d) { break; }
	}

	this->accumulated_error_ = accErr;

	return i;
}

/*std::mutex qsc_mutex;
std::mutex ntcnt_mutex;
std::vector<std::mutex> qsc_mutexes;
std::vector<std::mutex> ntcnt_mutexes;*/

void DNA::ntSpecQualScores(vector<uint64_t>& qsc, vector<uint64_t>& ntcnt) {
	const size_t sql = std::min(sequence_length_, std::min(sequence_.length(), qual_.size()));

	if (qsc.size() < 6) {
		qsc.resize(6, 0);
	}
	if (ntcnt.size() < 6) {
		ntcnt.resize(6, 0);
	}
	for (size_t i = 0; i < sql; i++) {
		const short p = NT_POS[static_cast<unsigned char>(sequence_[i])];
		qsc[p] += static_cast<uint64_t>(std::max<int>(0, qual_[i]));
		ntcnt[p]++;
	}
}

void DNA::ntSpecQualScoresMT(vector<uint64_t>& qsc, vector<uint64_t>& ntcnt) {
	const size_t sql = std::min(sequence_length_, std::min(sequence_.length(), qual_.size()));

	if (qsc.size() < 6) {
		qsc.resize(6, 0);
	}
	if (ntcnt.size() < 6) {
		ntcnt.resize(6, 0);
	}
	for (size_t i = 0; i < sql; i++) {
		const short p = NT_POS[static_cast<unsigned char>(sequence_[i])];
		qsc[p] += static_cast<uint64_t>(std::max<int>(0, qual_[i]));
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
	unsigned int QS = this->length();
	unsigned int QSh = QS >> 2;
	QSh = max(QSh, W);
	if (W >= QS) { return true; } // too short

	for (unsigned int i = QS - 1; i > QS - (unsigned int)W - 1; i--) {
		AQ += qual_[i];
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
			invalidateQualityCaches();
			return true;
		}
	}

	const size_t oldLogicalLength = std::min(sequence_length_, sequence_.length());

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

	const size_t logicalStart = std::min(static_cast<size_t>(start), oldLogicalLength);
	const size_t logicalStop = std::min(static_cast<size_t>(stop), oldLogicalLength);
	const size_t removedLogical = logicalStop > logicalStart ? logicalStop - logicalStart : 0;
	sequence_length_ = std::min(sequence_.length(), oldLogicalLength - removedLogical);
	invalidateQualityCaches();

	return true;
}

DNA::PrimerMatchDetails DNA::getPrimerMatchDetails(const string& primer, int start) const {
	PrimerMatchDetails details;
	details.start = start;
	if (start < 0 || static_cast<size_t>(start) + primer.length() > sequence_length_) {
		return details;
	}

	details.fullLength = true;
	details.observed = sequence_.substr(start, primer.length());
	for (size_t i = 0; i < primer.length(); ++i) {
		if (DNA_IUPAC[sequence_[start + i] + 256 * primer[i]] != 0) {
			++details.mismatches;
		}
	}
	return details;
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
			if (sequence_[Prp2] != 'N') { c_points++; }
			Prp++; Prp2++;
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

	this->resetTruncation();
}

void DNA::reverse_compliment(bool reset) {
	const size_t logicalLength = std::min(sequence_length_, std::min(sequence_.size(), qual_.size()));
	sequence_.resize(logicalLength);
	qual_.resize(logicalLength);
	sequence_length_ = logicalLength;
	reverseTS(sequence_);
	std::reverse(qual_.begin(), qual_.end());
	reversed_ = !reversed_;
	invalidateQualityCaches();
	if (!reset) {
		return;
	}
	accumulated_error_ = 0.; good_quality_ = false; mid_quality_ = false;
	tempFloat = 0.f;

}

//match from end of sequence_ to find rev primer
int DNA::matchSeqRev(const string& PrSt, int Err, int searchSpace,
	int start, bool exhaustive) {
	//fail::ret -1
	int PrL = (int)PrSt.length();
	if (start == 0) {
		start = PrL; //default seed set to 5
	}
	int SeL = static_cast<int>(std::min(sequence_length_, sequence_.size()));
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
	if (PrL == 0 || error < 0 || static_cast<size_t>(PrL) > sequence_.size()) {
		c_err = error + 1;
		return -1;
	}
	const int validStarts = static_cast<int>(sequence_.size()) - PrL + 1;
	const int searchEnd = std::max(0, std::min(maxPos, validStarts));
	for (pos = 0; pos < searchEnd; pos++) {
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
	if (length() == 0) { return true; }
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
	if (max == 0 || max >= (int)length()) { return 0; }
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

void DNA::removeNegativeQuals() {
	//remove negative qual_ scores from sequence
	if (qual_.size() == 0) { return; }
	for (unsigned int i = 0; i < qual_.size(); i++) {
		if (qual_[i] < 0) {
			qual_[i] = 0;
		}
	}
}

void DNA::removeNonCanonicalNTs() {
	//remove non canonical NTs from sequence
	if (sequence_.size() == 0) { return; }
	for (unsigned int i = 0; i < sequence_.size(); i++) {
		if (!canonicalDNA(sequence_[i])) {
			sequence_[i] = 'A';
			//qual_[i] = 5;// static_cast<char>(5);
			qual_[i] =  static_cast<char>(0);
		}
	}
}
void DNA::writeSeq_derep(ostream& wr, bool singleLine) {
	//specialized functions called from derep objects to write sequences with dereplication counts in header
	if ( getGlobalOptContainer()->runtime().onlyCanonicalNTs){
		this->removeNonCanonicalNTs();
	}
	
	writeSeq(wr, singleLine);

}


void DNA::writeSeq(ostream & wr, bool singleLine) {
		if (sequence_.size() == 0) { return; }
	wr << ">" << new_id_ << "\n";
	const size_t seqLen = length();
	if (singleLine) {
		wr.write(sequence_.data(), seqLen);
		wr.put('\n');
		return;
	}

	size_t last = 0;
	while ((last + 80) <= seqLen) {
		wr.write(sequence_.data() + last, 80);
		wr.put('\n');
		last += 80;
	}
	if (last < seqLen) {
		wr.write(sequence_.data() + last, seqLen - last);
		wr.put('\n');
	}
}
void DNA::writeSeq(string& str, bool singleLine) {
	if (sequence_.size() == 0) { return; }
 const size_t seqLen = length();
	str.clear();
	str.reserve(new_id_.size() + seqLen + 8 + (singleLine ? 1 : (seqLen / 80) + 1));
	str.push_back('>');
	str += new_id_;
	str.push_back('\n');
	if (singleLine) {
        str.append(sequence_.data(), seqLen);
		str.push_back('\n');
		return;
	}

	size_t last = 0;
	while ((last + 80) <= seqLen) {
		str.append(sequence_.data() + last, 80);
		str.push_back('\n');
		last += 80;
	}
 if (last < seqLen) {
		str.append(sequence_.data() + last, seqLen - last);
		str.push_back('\n');
	}
}
string DNA::writeSeq(bool singleLine) {
	string str = string("");
	this->writeSeq(str, singleLine);
	return str;

}
/**/
void DNA::writeQual(ostream& wr, bool singleLine) {
    if (qual_.size() == 0) { return; }
	int cnt = 0;
	wr << ">" << new_id_ << "\n";
	const char endChr = singleLine ? ' ' : '\n';
	for (unsigned int i = 0; i < length(); i++) {
		cnt++;
		wr << static_cast<int>(qual_[i]);
		if (cnt < 80) {
			wr.put(' ');
		}
		else {
			wr.put(endChr);
			cnt = 0;
		}
	}
	if (cnt > 0) {
		wr.put('\n');
	}
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
void DNA::writeFastQ_derep(ostream& wr, bool newHD) {//
	const OptContainer::RuntimeOptions& rt = globalOptContainer().runtime();
	if (rt.onlyCanonicalNTs) {
		this->removeNonCanonicalNTs();
	}
	if (rt.strictPositiveQualsOut) {
		this->removeNegativeQuals();
	}

	writeFastQ(wr, newHD);
}

void DNA::writeFastQ(ostream& wr, bool newHD) {//
	if (qual_.size() == 0 || sequence_.length() == 0 || length() == 0) { return; }
	const size_t seqLen = getFastqRecordLength(this->length(), sequence_.size(), qual_.size());
	if (seqLen == 0) { return; }
	const string& hdr = newHD ? new_id_ : id_;
	thread_local std::string qual_traf;
	buildFastqQualityString(qual_traf, qual_, seqLen, fastq_write_offset_);

	wr.put('@');
	wr << hdr;
	wr.put('\n');
	wr.write(sequence_.data(), seqLen);
	wr.write("\n+\n", 3);
	wr.write(qual_traf.data(), static_cast<std::streamsize>(qual_traf.size()));
	wr.put('\n');
}
string DNA::writeFastQ(bool newHD) {
	string ret = string("");
	this->writeFastQ(ret, newHD);
	return ret;
}
void DNA::writeFastQ(string& ret, bool newHD) {
	if (qual_.size() == 0 || sequence_.length() == 0 || length() == 0) { return; }
	const size_t seqLen = getFastqRecordLength(this->length(), sequence_.size(), qual_.size());
	if (seqLen == 0) { return; }
	const string& hdr = newHD ? new_id_ : id_;
	thread_local std::string qual_traf;
	buildFastqQualityString(qual_traf, qual_, seqLen, fastq_write_offset_);
	ret.reserve(ret.size() + hdr.size() + seqLen + qual_traf.size() + 5);

	ret.push_back('@');
	ret += hdr;
	ret.push_back('\n');
	ret.append(sequence_.data(), seqLen);
	ret += "\n+\n";
	ret.append(qual_traf.data(), qual_traf.size());
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
 return (getPositionFreeIdView() == d->getPositionFreeIdView());

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
	fastq_write_offset_ = static_cast<qual_score>(ofastQver);
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
	int sample_id, bool b_derep_as_fasta_, bool collectMergedQualities, bool markDereplicated) {
	if (dna == nullptr) {
		return;
	}
	if (markDereplicated) { dna->setDereplicated(); }
   const bool doCounterUpdate = (sample_id >= 0);
	const bool doQualityUpdate = !b_derep_as_fasta_;
	const bool doReplacementCheck = dna->isGreenQual();
	if (!doCounterUpdate && !doQualityUpdate && !doReplacementCheck) {
		return;
	}

	if (doCounterUpdate) {
		incrementSampleCounter(sample_id);
	}
	if (doQualityUpdate) { // only improve qualities if we do not replace
		sumQualities(dna);
        if (collectMergedQualities && dnaM) {
            auto evidence = merge_ ? merge_ : merge_quality_evidence_;
            if (evidence) { evidence->sumQualities(dnaM); }
            else { merge_quality_evidence_ = make_shared<DNAunique>(dnaM, -1); }
        }
	}
	// only replace old unique dna with new if it has good quality and its seed is better
	// betterPreSeed makes sense with uparse/
	float tmp(-1.f);
	// Borrow this representative for the synchronous comparison without allocating a control block.
	const shared_ptr<DNA> self_view(shared_ptr<DNA>{}, this);
	if (doReplacementCheck &&
		whoIsBetter(dna, dna2, dnaM, self_view, this->getPair(), this->getMerge(),
		 tmp, false, this)) { //this->getBestSeedLength(),
		//betterPreSeed(dna, dna2)) {
		// Uparse does NOT use qualities for clustering, therefore we do not need to calculate average qualities for this dna object
		// Lotus uses qualities for constructing taxonomy (therefore we do replace dna if theres a better quality read)
		takeOverDNA(dna, dna2, dnaM, true);
	}
}

read_occ::iterator DNAunique::findOccurrenceEntry(int sample_id) {
	return std::find_if(occurence.begin(), occurence.end(), [sample_id](const read_occ_entry& entry) {
		return entry.first == sample_id;
		});
}

read_occ::const_iterator DNAunique::findOccurrenceEntry(int sample_id) const {
	return std::find_if(occurence.begin(), occurence.end(), [sample_id](const read_occ_entry& entry) {
		return entry.first == sample_id;
		});
}

void DNAunique::incrementSampleCounter(int sample_id) {
	if (sample_id < 0) {
		return;
	}
	//count_++;
#ifdef _MAPDEREPLICATE
	auto it = findOccurrenceEntry(sample_id);
	if (it == occurence.end()) {
		occurence.emplace_back(sample_id, static_cast<sample_count_t>(1));
	}
	else {
		if (it->second == std::numeric_limits<sample_count_t>::max()) {
			throw std::overflow_error("Dereplication sample count overflow");
		}
		++(it->second);
	}
	invalidateTotalSumCache();
#endif
}
void DNAunique::incrementSampleCounterBy(int sample_id, uint64_t count) {
	if (sample_id < 0 || count == 0) {
		return;
	}
	auto it = findOccurrenceEntry(sample_id);
	if (it == occurence.end()) {
		occurence.emplace_back(sample_id, count);
	}	else {
		if (it->second > std::numeric_limits<sample_count_t>::max() - count) {
			throw std::overflow_error("Dereplication sample count overflow");
		}
		it->second += count;
	}
	invalidateTotalSumCache();
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

bool sortDescending(uint64_t i, uint64_t j) { return (i > j); }//descending sort
vector<uint64_t> DNAunique::getDerepMapSort(size_t wh) {
	vector<uint64_t> vals;
	vals.reserve(occurence.size());
	for (const auto& entry : occurence) {
		if (entry.second == 0) { continue; }
		vals.push_back(entry.second);
	}
	if (wh > vals.size()) { wh = vals.size(); }
	if (vals.size() <= 1 || wh == 0) {
		return vals;
	}
	if (wh >= vals.size()) {
		sort(vals.begin(), vals.end(), sortDescending);
	}
	else {
		std::partial_sort(vals.begin(), vals.begin() + wh, vals.end(), sortDescending);
	}
	return vals;
}

bool DNAunique::pass_deprep_smplSpc(const vector<int>& cv) {
	//unordered_map<int, int> occ;
	//combined samples will not be considered
	//occ = occurence;
	for (const auto& entry : occurence) {
		if (entry.second == 0 || entry.first < 0 || static_cast<size_t>(entry.first) >= cv.size()) { continue; }
		int ref = cv[entry.first];
		if (ref != -1 && entry.second >= static_cast<sample_count_t>(ref)) {
			return true;
		}
	}
	return false;

}


void DNAunique::transferOccurence(const shared_ptr<DNAunique> dna_unique) {
	if (dna_unique == nullptr) {
		return;
	}
	const read_occ& dereplication_map = dna_unique->getDerepMap();
	if (dereplication_map.empty()) {
		return;
	}
	for (const auto& entry : dereplication_map) {
		if (entry.second == 0) { continue; }
		auto it = findOccurrenceEntry(entry.first);
		if (it == occurence.end()) {
			occurence.push_back(entry);
		}
		else if (it->second > (std::numeric_limits<sample_count_t>::max() - entry.second)) {
			throw std::overflow_error("Dereplication sample count overflow while merging occurrences");
		}
		else {
			it->second += entry.second;
		}
	}
	invalidateTotalSumCache();
}


void DNAunique::writeMap(ofstream& os, const string& hd, vector<uint64_t>& counts_per_sample, const vector<int>& combiID, vector<uint64_t>& combined_counts, vector<int>& touched_groups) {
	if (occurence.empty()) { return; }
	uint64_t total_count(0);

	//prints combined sample counts
	os << hd;
	if (!combiID.empty()) {
		for (const auto& entry : occurence) {
			if (entry.second == 0 || entry.first < 0 || static_cast<size_t>(entry.first) >= combiID.size()) { continue; }
			const int combined_id = combiID[entry.first];
			if (combined_id < 0 || static_cast<size_t>(combined_id) >= combined_counts.size()) { continue; }
			if (combined_counts[combined_id] == 0) {
				touched_groups.push_back(combined_id);
			}
			combined_counts[combined_id] += entry.second;
		}

		for (size_t i = 0; i < touched_groups.size(); ++i) {
			const int combined_id = touched_groups[i];
			const uint64_t rawCount = combined_counts[combined_id];
			total_count += rawCount;
			os << "\t";
			os << combined_id << ":" << rawCount;
			combined_counts[combined_id] = 0;
		}
		touched_groups.clear();
	}
	else {
		for (const auto& entry : occurence) {
			if (entry.second == 0) { continue; }
			const uint64_t rawCount = entry.second;
			total_count += rawCount;
			os << "\t";
			os << entry.first << ":" << rawCount;
		}
	}

	//counts non-combined sample counts
	for (const auto& entry : occurence) {
		if (entry.second == 0 || entry.first < 0 || static_cast<size_t>(entry.first) >= counts_per_sample.size()) { continue; }
		const uint64_t add = entry.second;
		if (counts_per_sample[entry.first] > (std::numeric_limits<uint64_t>::max() - add)) {
			throw std::overflow_error("Dereplication map count overflow");
		}
		else {
			counts_per_sample[entry.first] += add;
		}
	}
	os << endl;

	/*if (total_count != count_) {
		cerr << "Unequal counts in Map(" << total_count << ") and HeadCnt(" << count_ << "):" << endl << this->getId() << endl;
		exit(82);
	}*/
}

uint64_t DNAunique::totalSum() {
	if (!total_sum_cache_valid_) {
		total_sum_cache_ = 0;
		for (const auto& xx : occurence) {
			if (total_sum_cache_ > std::numeric_limits<uint64_t>::max() - xx.second) {
				throw std::overflow_error("Dereplication total count overflow");
			}
			total_sum_cache_ += xx.second;
		}
		total_sum_cache_valid_ = true;
	}
	if (total_sum_cache_ == 0) {
		return 1;
	}
	return total_sum_cache_;
}


void DNAunique::addToQualSum(qual_accum_t& target, qual_accum_t value) {
	const qual_accum_t maxVal = std::numeric_limits<qual_accum_t>::max();
	if (target > (maxVal - value)) {
		throw std::overflow_error("Dereplication quality accumulation overflow");
	}
	else {
		target += value;
	}
}



shared_ptr<DNAunique> DNAunique::cloneForOutput() const {
    auto copy = make_shared<DNAunique>(*this);
    if (pair_) { copy->pair_ = pair_->cloneForOutput(); }
    if (merge_) { copy->merge_ = merge_->cloneForOutput(); }
    if (merge_quality_evidence_) { copy->merge_quality_evidence_ = merge_quality_evidence_->cloneForOutput(); }
    return copy;
}

void DNAunique::sumMergeQualities(const shared_ptr<DNAunique>& other) {
    if (!merge_ || !other) { return; }
    auto evidence = other->merge_ ? other->merge_ : other->merge_quality_evidence_;
    if (evidence) { merge_->sumQualities(evidence); }
}

void DNAunique::saveMem() {
	if (!id_.empty()) {
		new_id_.assign(id_.data(), getSpaceHeadPos(id_));
	}
	id_.clear();
	id_.shrink_to_fit();
}


void DNAunique::attachPair(shared_ptr<DNAunique> dna_unique) {
	if (dna_unique == nullptr) { return; }
	pair_ = dna_unique;
	pair_->saveMem();
}
shared_ptr<DNAunique> DNAunique::getPair(void) {
	return pair_;
}
shared_ptr<DNAunique> DNAunique::getMerge(void) {
	return merge_;
}

void DNAunique::attachMerge(shared_ptr<DNAunique> dnamerge) {
	if (dnamerge == nullptr) { return; }
	merge_ = dnamerge;
	merge_->saveMem();
	FtsDetected.mergeLength = merge_->length();
}

void DNAunique::promoteQualitiesToSums(bool releaseRawQualities) {
	// Preserve evidence beyond the current representative's logical end: a
	// later representative may cover those positions again.
	if (quality_sum_per_base_.empty() && !qual_.empty()) {
		const size_t n = std::min(static_cast<size_t>(length()), qual_.size());
		quality_sum_per_base_.resize(n, 0);
		quality_coverage_per_base_.assign(n, 1);
		for (size_t i = 0; i < n; ++i) {
			quality_sum_per_base_[i] = static_cast<qual_accum_t>(std::max<int>(0, qual_[i]));
		}
	}
	if (releaseRawQualities && !quality_sum_per_base_.empty()) {
		qual_.clear();
		qual_.shrink_to_fit();
	}
}

void DNAunique::prepSumQuals() {
	promoteQualitiesToSums(false);
}
void DNAunique::sumQualities(const shared_ptr<DNAunique>& dna) {
	if (!dna || dna.get() == this) { return; }
	if (dna->quality_sum_per_base_.empty()) {
		sumQualities(std::static_pointer_cast<DNA>(dna));
		return;
	}
	prepSumQuals();
	const size_t n = dna->quality_sum_per_base_.size();
	quality_sum_per_base_.resize(std::max(n, quality_sum_per_base_.size()), 0);
	quality_coverage_per_base_.resize(quality_sum_per_base_.size(), 0);
	for (size_t i = 0; i < n; ++i) {
		addToQualSum(quality_sum_per_base_[i], dna->quality_sum_per_base_[i]);
		addToQualSum(quality_coverage_per_base_[i], dna->quality_coverage_per_base_[i]);
	}
}
void DNAunique::sumQualities(const shared_ptr<DNA>& dna) {
	if (!dna) { return; }
	const auto& quals = dna->getQual();
	const size_t n = std::min(static_cast<size_t>(dna->length()), quals.size());
	if (n == 0) { return; }
	prepSumQuals();
	quality_sum_per_base_.resize(std::max(n, quality_sum_per_base_.size()), 0);
	quality_coverage_per_base_.resize(quality_sum_per_base_.size(), 0);
	for (size_t i = 0; i < n; ++i) {
		addToQualSum(quality_sum_per_base_[i], static_cast<qual_accum_t>(std::max<int>(0, quals[i])));
		addToQualSum(quality_coverage_per_base_[i], 1);
	}
}

void DNAunique::releaseRawQualitiesIfPromoted() {
	promoteQualitiesToSums(true);
}

void DNAunique::recordWhoIsBetterCompare(int dMergeLength, int rMergeLength) {
	if (dMergeLength >= 0) {
		++who_is_better_merge_stats_.merged_compares;
		who_is_better_merge_stats_.accumulated_merge_length += static_cast<uint64_t>(dMergeLength);
	}
	else {
		++who_is_better_merge_stats_.non_merged_compares;
	}

	if (rMergeLength >= 0) {
		++who_is_better_merge_stats_.merged_compares;
		who_is_better_merge_stats_.accumulated_merge_length += static_cast<uint64_t>(rMergeLength);
	}
	else {
		++who_is_better_merge_stats_.non_merged_compares;
	}
}


void DNAunique::Count2Head(bool usFmt) {
	//unordered_map<int, int> occurence;
	//occ = occurence;

	if (occurence.empty()) {
		return;
	}
	if (usFmt) {
		new_id_ += ";size=" + std::to_string(totalSum()) + ";";
	}
	else {
		new_id_ += "_" + std::to_string(totalSum());
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
	quality_coverage_per_base_ = std::move(dna_unique_old->quality_coverage_per_base_);
}
void DNAunique::takeOverDNA(shared_ptr<DNA> const better_dna, shared_ptr<DNA> const dna2, shared_ptr<DNA> const dnaMerge, bool keepQualities) {

	// Representative metadata may change; accumulated observations must not.
	pair_.reset();
	if (dna2 != nullptr) {
		auto pair_unique = make_shared<DNAunique>(dna2, -1);
		this->attachPair(pair_unique);
	}
	// Keep accumulated consensus evidence separate from the selected attachment.
    auto mergeEvidence = merge_ ? merge_ : merge_quality_evidence_;
	merge_.reset();
	if (dnaMerge != nullptr) {
		auto merge_unique = make_shared<DNAunique>(dnaMerge, -1);
        if (mergeEvidence && !mergeEvidence->quality_sum_per_base_.empty()) {
            merge_unique->quality_sum_per_base_ = std::move(mergeEvidence->quality_sum_per_base_);
            merge_unique->quality_coverage_per_base_ = std::move(mergeEvidence->quality_coverage_per_base_);
        }
		this->attachMerge(merge_unique);
        merge_quality_evidence_.reset();
	} else {
        merge_quality_evidence_ = std::move(mergeEvidence);
    }

	avg_qual_ = keepQualities ? better_dna->avg_qual_ : 0.f;
	accumulated_error_ = keepQualities ? better_dna->accumulated_error_ : 0.;
	FtsDetected = better_dna->FtsDetected;
	FtsDetected.mergeLength = merge_ ? static_cast<int>(merge_->length()) : -1;
  failed_ = better_dna->failed_;
	id_ = better_dna->id_;
	id_fixed_ = better_dna->id_fixed_;
	new_id_ = better_dna->new_id_;
	mid_quality_ = better_dna->mid_quality_;
  good_quality_ = better_dna->good_quality_;
	merge_offset_ = better_dna->merge_offset_;
	merge_seed_pos_ = better_dna->merge_seed_pos_;
	reversed_merge_ = better_dna->reversed_merge_;
	seed_length_ = better_dna->seed_length_;
	quality_sum_ = keepQualities ? better_dna->quality_sum_ : 0;
	QualCtrl = better_dna->QualCtrl;
	qual_ = keepQualities ? better_dna->qual_ : vector<qual_score>();
	sequence_ = better_dna->sequence_;
	sequence_length_ = better_dna->sequence_length_;
	sample_id_ = better_dna->sample_id_;
	reversed_ = better_dna->reversed_;
	read_position_ = better_dna->read_position_;
	tempFloat = better_dna->tempFloat;


	this->saveMem();
}

void DNAunique::compactForFastaDerepStorage() {
	qual_.clear();
	qual_.shrink_to_fit();
	quality_sum_per_base_.clear();
	quality_sum_per_base_.shrink_to_fit();
	quality_coverage_per_base_.clear();
	quality_coverage_per_base_.shrink_to_fit();
	quality_sum_ = 0;
	avg_qual_ = 0.f;
	accumulated_error_ = 0.;
}

std::vector<DNAunique::qual_accum_t> DNAunique::transferPerBaseQualitySum() {
	return std::move(quality_sum_per_base_);
}

void DNAunique::prepareDerepQualities(int ofastQver) {
	prepSumQuals();
	derep_fastq_offset_ = static_cast<qual_score>(ofastQver);
}

void DNAunique::writeFastQ(ostream& wr, bool newHD) {
	if (!qual_.empty()) {
		DNA::writeFastQ(wr, newHD);
		return;
	}
	if (sequence_.length() == 0 || length() == 0) {
		return;
	}
	prepSumQuals();
	if (quality_sum_per_base_.empty()) {
		return;
	}
	const size_t seqLen = this->length();
	const string& hdr = newHD ? new_id_ : id_;
	thread_local std::string qual_traf;
	qual_traf.assign(seqLen, static_cast<char>(fastq_write_offset_));
	const int fqOffset = static_cast<int>(fastq_write_offset_);
	const size_t qlen = std::min(seqLen, quality_sum_per_base_.size());
	for (size_t i = 0; i < qlen; ++i) {
		int q = static_cast<int>(std::round(static_cast<double>(quality_sum_per_base_[i]) / static_cast<double>(std::max<qual_accum_t>(1, quality_coverage_per_base_[i])))) + fqOffset;
		qual_traf[i] = sanitizeFastqQualityChar(q, fqOffset);
	}

	wr.put('@');
	wr << hdr;
	wr.put('\n');
	wr.write(sequence_.data(), seqLen);
	wr.write("\n+\n", 3);
	wr.write(qual_traf.data(), static_cast<std::streamsize>(qual_traf.size()));
	wr.put('\n');
}

string DNAunique::writeFastQ(bool newHD) {
	string ret = string("");
	this->writeFastQ(ret, newHD);
	return ret;
}

void DNAunique::writeFastQ(string& ret, bool newHD) {
	if (!qual_.empty()) {
		DNA::writeFastQ(ret, newHD);
		return;
	}
	if (sequence_.length() == 0 || length() == 0) {
		return;
	}
	prepSumQuals();
	if (quality_sum_per_base_.empty()) {
		return;
	}
	const size_t seqLen = this->length();
	const string& hdr = newHD ? new_id_ : id_;
	thread_local std::string qual_traf;
	qual_traf.assign(seqLen, static_cast<char>(fastq_write_offset_));
	const int fqOffset = static_cast<int>(fastq_write_offset_);
	const size_t qlen = std::min(seqLen, quality_sum_per_base_.size());
	for (size_t i = 0; i < qlen; ++i) {
		int q = static_cast<int>(std::round(static_cast<double>(quality_sum_per_base_[i]) / static_cast<double>(std::max<qual_accum_t>(1, quality_coverage_per_base_[i])))) + fqOffset;
		qual_traf[i] = sanitizeFastqQualityChar(q, fqOffset);
	}

	ret.reserve(ret.size() + hdr.size() + seqLen + qual_traf.size() + 5);
	ret.push_back('@');
	ret += hdr;
	ret.push_back('\n');
	ret.append(sequence_.data(), seqLen);
	ret += "\n+\n";
	ret.append(qual_traf.data(), qual_traf.size());
	ret.push_back('\n');
}


void DNAunique::writeDerepFastQ(ofstream& wr, bool newHD) {
	if ((sequence_.length() == 0 || length() == 0)) {
		return;
	}
	bool superStrictCanonical = false;
	if (getGlobalOptContainer()->runtime().onlyCanonicalNTs) {
		superStrictCanonical = true;
	}

	const size_t seqLen = length();

	if (quality_sum_per_base_.empty()) {
		prepSumQuals();
	}
	thread_local std::string qualities_avg;
	qualities_avg.assign(seqLen, static_cast<char>(derep_fastq_offset_));
	const int derepOffset = static_cast<int>(derep_fastq_offset_);
	const size_t qlen = std::min(seqLen, quality_sum_per_base_.size());
	for (size_t i = 0; i < qlen; i++) {
		int newbase = static_cast<int>(std::round(static_cast<double>(quality_sum_per_base_[i]) / static_cast<double>(std::max<qual_accum_t>(1, quality_coverage_per_base_[i])))) + derepOffset;
		qualities_avg[i] = sanitizeFastqQualityChar(newbase, derepOffset);
	}
	thread_local std::string seqOut;
	seqOut.resize(seqLen);
	for (size_t i = 0; i < seqLen; ++i) {
		const char base = sequence_[i];
		if (!canonicalDNA(base)) {
			if (superStrictCanonical) {
				seqOut[i] = 'A';//default to A and set minimum legal quality
				qualities_avg[i] = sanitizeFastqQualityChar(derepOffset, derepOffset);
			}
			else {
				seqOut[i] = 'N';
			}

		}
		else {
			seqOut[i] = base;
		}
	}

	wr.put('@');
	if (newHD) {
		wr << new_id_;
	}
	else {
		wr << id_;
	}
	wr.put('\n');
	wr.write(seqOut.data(), static_cast<std::streamsize>(seqLen));
	wr.write("\n+\n", 3);
	wr.write(qualities_avg.data(), static_cast<std::streamsize>(seqLen));
	wr.put('\n');

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





//*******************************************
//*        DNAuniqSet OBJECT
//*******************************************


void DNAuniqSet::addNewDNAuniq(shared_ptr<DNA> dna, shared_ptr<DNA> dna2, shared_ptr<DNA> dnaM, int MrgPos1, int sample_id, bool b_derep_as_fasta_, bool markDereplicated) {
	if (markDereplicated) { dna->setDereplicated(); }
	shared_ptr<DNAunique> new_dna_unique = make_shared<DNAunique>(dna, sample_id);
	new_dna_unique->saveMem();
	if (new_dna_unique == nullptr) { return; }
	if (dna2 != nullptr) {
		auto pair_unique = make_shared<DNAunique>(dna2, sample_id);
		new_dna_unique->attachPair(pair_unique);
	}
	if (dnaM != nullptr) {
		auto merge_unique = make_shared<DNAunique>(dnaM, sample_id);
		new_dna_unique->attachMerge(merge_unique);
	}
	DNUs[MrgPos1] = new_dna_unique;
}

map<int, shared_ptr<DNAunique>>::iterator DNAuniqSet::find(const int x) {
	map<int, shared_ptr<DNAunique>>::iterator r = DNUs.find(x);
	return r;
}
const map<int, shared_ptr<DNAunique>>::iterator DNAuniqSet::end() {
	return DNUs.end();
}
const map<int, shared_ptr<DNAunique>>::iterator DNAuniqSet::begin() {
	return DNUs.begin();
}
uint64_t DNAuniqSet::totalSum() const {
    uint64_t total = 0;
    for (const auto& entry : DNUs) {
        if (!entry.second) { continue; }
        const uint64_t count = entry.second->totalSum();
        if (count > std::numeric_limits<uint64_t>::max() - total) {
            throw std::overflow_error("Dereplication subgroup count overflow");
        }
        total += count;
    }
    return total;
}

void DNAunique::absorbDerepEvidence(const shared_ptr<DNAunique>& other, bool collectQualities) {
    if (collectQualities) {
        sumQualities(other);
        sumMergeQualities(other);
    }
    transferOccurence(other);
}

void DNAuniqSet::addEvidenceTo(const shared_ptr<DNAunique>& target, bool collectQualities) const {
    for (const auto& entry : DNUs) {
        if (entry.second) { target->absorbDerepEvidence(entry.second, collectQualities); }
    }
}

shared_ptr<DNAunique> DNAuniqSet::best(bool addCnts, bool collectQualities) {
    setBest(addCnts, collectQualities);
    return std::move(bestDNU);
}

void DNAuniqSet::setBest(bool addCnts, bool collectQualities) {
    bestDNU.reset();
    auto winner = DNUs.end();
    uint64_t bestCount = 0;
    for (auto it = DNUs.begin(); it != DNUs.end(); ++it) {
        if (!it->second) { continue; }
        if (winner == DNUs.end()) {
            winner = it;
            bestCount = it->second->totalSum();
            continue;
        }
        auto candidateMerge = it->second->getMerge();
        auto currentMerge = winner->second->getMerge();
        float candidateWeight = 1.f, currentWeight = 1.f, lengthRatio = 1.f;
        if (candidateMerge && currentMerge) {
            lengthRatio = static_cast<float>(candidateMerge->length()) / currentMerge->length();
        } else {
            if (!candidateMerge) { candidateWeight = 0.8f; }
            if (!currentMerge) { currentWeight = 0.8f; }
        }
        const float ratio = (static_cast<float>(it->second->totalSum()) * candidateWeight) /
            (static_cast<float>(bestCount) * currentWeight) * lengthRatio;
        if (ratio > 1.f) {
            winner = it;
            bestCount = it->second->totalSum();
        }
    }
    if (winner == DNUs.end()) { return; }
    // Output edits and subgroup aggregation must not mutate the live tracker:
    // producers can resume and subsequent snapshots must see every new count.
    bestDNU = winner->second->cloneForOutput();
    for (const auto& entry : DNUs) {
        if (!entry.second || entry.first == winner->first) { continue; }
        const int64_t distance = static_cast<int64_t>(entry.first) - winner->first;
        if (!addCnts && distance != -1 && distance != 1) { continue; }
        bestDNU->absorbDerepEvidence(entry.second, collectQualities);
    }
}
