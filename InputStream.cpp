/* sdm: simple demultiplexer
Copyright (C) 2013  Falk Hildebrand
email: Falk.Hildebrand@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif

#include <mutex>
#include <thread>
#include "InputStream.h"
#include "Common.h"
#include "FastxReader.h"

// Utilities are implemented in Common.cpp (defined once in Common.cpp)




//compares two DNA entries, decides which one has overall better stats
bool whoIsBetter(shared_ptr<DNA> d1, shared_ptr<DNA> d2, shared_ptr<DNA> dM, 
			shared_ptr<DNA> r1, shared_ptr<DNA> r2, shared_ptr<DNA> rM, 
			float& ever_best, bool forSeed, DNAunique* merge_stats_owner,
			int dSiz , int rSiz) {
	
	if (d1 == nullptr || r1 == nullptr) {
		return false;
	}

	//if (forSeed) { 		return false; 	}
	//check if two primers present
	if (d2 == nullptr) {//hard reason .. only for PacBio etc reads
		if (d1->has2PrimersDetected() && !r1->has2PrimersDetected()) { return true; }
		if (!d1->has2PrimersDetected() && r1->has2PrimersDetected()) { return false; }
    } else {
		//if (r2 != nullptr && d2->getRevPrimDetect() && !r1->getFwdPrimDetect() && !r2->getRevPrimDetect()) {	return true;}
	}
	//check if at least 1 primers present
	if (d1->getFwdPrimDetect() && !r1->getFwdPrimDetect()) { return true; }//hard reason
	if (!d1->getFwdPrimDetect() && r1->getFwdPrimDetect()) { return false; }

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
    if (merge_stats_owner != nullptr) {
		merge_stats_owner->recordWhoIsBetterCompare(d1->getMergeLength(), r1->getMergeLength());
	}
    double mergedFraction = 0.0;
	double avgMergeLen = 0.0;
	if (merge_stats_owner != nullptr) {
		const auto& ms = merge_stats_owner->getWhoIsBetterMergeStats();
		const uint64_t totalCompared = ms.merged_compares + ms.non_merged_compares;
		if (totalCompared > 0) {
			mergedFraction = static_cast<double>(ms.merged_compares) / static_cast<double>(totalCompared);
		}
		avgMergeLen = ms.merged_compares > 0 ? static_cast<double>(ms.accumulated_merge_length) / static_cast<double>(ms.merged_compares) : 0.0;
	}
	bool mergAdv = dMerge && !rMerg;//if d merged but ref did not, go for d, hard filter

	//first check if d1 has merged, but ref did not.. clearly go for d, hard filter
	//if (d1->getMergeLength() != -1 && r1->getMergeLength() == -1) { return true; }

	//hard check on length ratios.. too drastically small , don't use d1
   if (refL <= 0.) { return false; }
	if ((curL) / (refL) < BestLengthRatio && !mergAdv) { return false; }
	//at least 90% length of "good" hit
	//if (r1->getMergeErrors() < 0) {//no merge, can look at read1 only
	//	if (d1->length() / r1->length() < RefLengthRatio) { return false; }
	//}


	float dmergErrSco = 0.f; 	float refMergErrSco = 0.f;
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
		if (maxMerr > 0.f) {
			dmergErrSco /= maxMerr; 	refMergErrSco /= maxMerr;
			dmergErrSco = 1.f - dmergErrSco; refMergErrSco = 1.f - refMergErrSco;
		}
	}

	//choose merged DNA if possible; d2 no longer needed then
	shared_ptr<DNA> dx, rx;
	bool allowMergeGuide = dMerge && rMerg;
	if (allowMergeGuide && dM != nullptr) { dx = dM; 	d2 = nullptr; } 	else { dx = d1; }
	if (allowMergeGuide && rM != nullptr) { rx = rM; r2 = nullptr;}else { rx = r1; }
	float thScore =  dx->getAvgQual(); //*(d1pid / 100)* log((float)curL);
	float rScore =  rx->getAvgQual();// *(refpid / 100)* log((float)refL);//r1->length()
	//if (thScore > rScore) {
		//also check for stable lowest score
		// if (d1->minQual() > r1->minQual() - MinQualDiff) { return true; }
	//}
  float maxScore = max(thScore, rScore);
	if (maxScore > 0.f) {
		thScore /= maxScore; 	rScore /= maxScore;
	}


	/*
	* accumulated error was overall too harsh.. switch to qual ratio instead..
	//checks if the new DNA has a better overall quality
	double dAcSc = dx->getAccumError();	double dLen = dx->mem_length();
	double tAcSc = rx->getAccumError();	double tLen = rx->mem_length();
	if (d2 != nullptr) { dAcSc += d2->getAccumError(); dLen += d2->mem_length(); }
	if (r2 != nullptr) { tAcSc += r2->getAccumError(); tLen += r2->mem_length(); }
	if (dLen <= 0. || tLen <= 0.) {
		return false;
	}
	//dAcSc /= dLen; 	tAcSc /= tLen;
	
	//calculate ratios to compare more easily among metrics
	double ratAccErr = dAcSc / tAcSc;///logRatio(dAcSc, tAcSc);
	*/

	double dAccuE = dx->getAccumError(); double rAccuE = rx->getAccumError();
	double AccErrRatio = pointFiveRatio(dAccuE, rAccuE,0.25);

	double dQ=dx->getAvgQual(); double rQ=rx->getAvgQual();
	if (d2 != nullptr) { dQ += d2->getAccumError();  }
	if (r2 != nullptr) { rQ += r2->getAccumError();  }

	double qualRatio = dQ / rQ; //higher better


	//double ratLength = double(curL) / (double)refL; //higher better
	double ratLength = pointFiveRatio(double(curL), (double)refL); //higher better, but not too much better, as this can be due to chimeras etc
    double ratId = 1.;
	if (refpid > 0.) {
		//ratId = d1pid / refpid * 10. - 9.; //higher better //10-fold weighting
		ratId=pointFiveRatio(double(d1pid), (double)refpid,2.);
	}

	//preference based on merge status, but modulated by how often merged DNAs have been better in the past (mergedFraction)
	double mergePreference = 1.0;
	if (dMerge != rMerg) {
		if (dMerge) {
			mergePreference = 1.0 + mergedFraction;
		}else {
			mergePreference = 1.0 - (0.5 * mergedFraction);
		}
	}

	//also take into account the merge length and if this is the averaged norm:
	double mergeLenDevRatio = 1.;
	if (dMerge && rMerg) {
		double dMLDIFF = (double) abs(avgMergeLen - d1->getMergeLength());
		double rMLDIFF = (double) abs(avgMergeLen - r1->getMergeLength());
		mergeLenDevRatio = (dMLDIFF+ avgMergeLen) / (rMLDIFF+ avgMergeLen);
		//mergeLenDeviation = (double)(rMLDIFF - dMLDIFF) / (double)(avgMergeLen + 1e-12); //positive if d is closer to average merge length, negative if r is closer
	}

	//also capture the size difference.. as higher abudant reads might be more trustworthy, but only if other metrics are not too different
	double logsizeRatio = pointFiveRatio(double(dSiz), double(rSiz),0.5); //higher better for r


	double thresh(1.01f);
	//qualRatio 
	if ( ((AccErrRatio* ratLength * ratId * mergePreference * logsizeRatio * mergeLenDevRatio) ) > thresh) {
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


bool whoIsBetter(shared_ptr<DNAunique> d1,
	shared_ptr<DNAunique> r1,
	float& ever_best, bool forSeed) {
	

	bool ret = whoIsBetter(d1, d1->getPair(), d1->getMerge(),
		r1, r1->getPair(), r1->getMerge(), ever_best, forSeed, d1.get(),
		d1->totalSum(), r1->totalSum());

	return ret;

}

double pointFiveRatio(double v1, double v2, double scale ) {
	
	double ret = 1.0;
	
	if (scale <= 0.) { cerr << "pointFiveRatio::scale not ok: "<<scale<<endl; scale = 0.5; }
	
	ret = (v1 / (v1 + v2)) + 0.5;
	if (scale != 0.5) {
		ret = (-1. * (1. - ret) * scale) + 1;
	}
	return ret;
}

double logRatio(double v1, double v2) {
	double ret = 1.f;
	if (v1 > 0. && v2 > 0.) {
		double logRef = std::log(v2);
		if (std::isfinite(logRef) && std::abs(logRef) > 1e-12) {
			ret = std::log(v1) / logRef; //smaller better, log changes terms around
		}
		else {
			ret = (v1 / (v1+v2)) + 0.5;
		}
	}
	if (ret > 0.) {
		return ret;
	}
	return -1. * ret;
}

double loglogRatio(double v1, double v2) {
	double ret = 1.f;
	if (v1 > 0. && v2 > 0.) {
		double logRef = std::log(std::log(v2));
		if (std::isfinite(logRef) && std::abs(logRef) > 1e-12) {
			ret = std::log(std::log(v1)) / logRef; //smaller better, log changes terms around
		}
		else {
			ret = v1 / v2;
		}
	}
	if (ret > 0.) {
		return ret;
	}
	return -1. * ret;
}





///////////////////////////////////////////////////////////////
//INPUT STREAMER


///////////////////////////////////////////////////////////////
//INPUT STREAMER

InputStreamer::~InputStreamer(){
	allStreamClose();
	for (size_t i = 0; i < pairs_read.size(); i++) {
		cdbg("~InputStreamer: pairsRead" + to_string(i) + ": " + to_string(pairs_read[i]) + " ");
	}
	cdbg("Returned DNA objects: "+to_string(N_DNAreturned) + " pairs: " + to_string(N_DNApairsReturned) + " ");
}

bool InputStreamer::read_fasta_entry(ifbufstream* fasta_is, ifbufstream* quality_is, shared_ptr<DNA> tdn1, shared_ptr<DNA> tdn2, int &cnt){
	if (fasta_is->eof()) return false;
	string line;
	string tseq;
	string tqual;
	string lineQ;
	fasta_is->getlines(line,false);
	if (!qualAbsent) {
		quality_is->getlines(lineQ,false);
	}
	cnt++;

	if (cnt == 1) {
		if (line[0] != '>' && fasta_is) {
			cerr << "ERROR: Line 1 in fasta file does not start with \">\" \n";
			exit(23);
		}
		tdn1->setHeader(line.substr(1));
		if (!fasta_is->getlines(line,true))
			return false;
		if (!qualAbsent) {
			quality_is->getlines(lineQ,true);
		}
		cnt++;
	}

	tdn1->setSequence(tseq);
	size_t lsize = tseq.size();
	vector<qual_score> Iqual(lsize, 11);
	if (!qualAbsent) {
		rtrim(tqual);
		const char* lQ = tqual.c_str();
		uint ii(0);
		qual_score nn(0);

		for (; ii < lsize; ii++) {
			nn = (qual_score)parseInt(&lQ);
			Iqual[ii] = nn;
			if (*lQ == '\0') break;
		}

		if (Iqual.size() != lsize) { cerr << "Unequal fasta/qual length in read_fasta_entry\n"; exit(923); }
	}
	tdn1->setQual(move(Iqual));
	if (!line.empty() && line[0] == '>') {
		tdn2->setHeader(line.substr(1));
	} else {
		tdn2->setHeader("");
	}
	return true;
}

void InputStreamer::IO_Error(string x) {
	fqverMTX.lock();
	cerr << x << endl;
	if (DieOnError) exit(632);
	ErrorLog.push_back(x);
	fqverMTX.unlock();
}

shared_ptr<DNA> InputStreamer::read_fastq_entry_fast(istream & fna, int& lnCnt, bool& corrupt) {
	string line;
	if (!fna) { return NULL; }
	if (!safeGetline(fna, line)) { return NULL; }
	shared_ptr<DNA> tdn = make_shared<DNA>("", "");
	if (line.length() == 0) { return tdn; }
	while (line[0] != '@') {
		IO_Error("ERROR on line " + itos(lnCnt) + ": Could not find '@' when expected (file likely corrupt, trying to recover):\n" + line);
		corrupt = true;
		if (!safeGetline(fna, line)) { corrupt = true;  return NULL; }
	}
	tdn->setHeader(line.substr(1));
	if (!safeGetline(fna, tdn->getSequence())) { corrupt = true;  return NULL; }
	if (!safeGetline(fna, line)) { corrupt = true; return NULL; }
	if (line[0] != '+') { IO_Error("Error input line " + itos(lnCnt + 2) + ": Could not find '+' when expected (file likely corrupt, aborting):\n" + line); corrupt = true; return tdn; }

	vector<qual_score> Iqual(tdn->mem_length(), 0);
	if (!safeGetline(fna, line)) { corrupt = true;  return NULL; }
	uint qcnt(0); uint lline = (uint)line.length();
	for (; qcnt < lline; qcnt++) {
		qual_score q = (qual_score)line[qcnt] - fastQver;
		Iqual[qcnt] = q;
	}

	if (qcnt == tdn->mem_length()) {
		tdn->setQual(move(Iqual));
	} else if (line.length() + qcnt != tdn->length()) {
		corrupt = true;
		IO_Error("Error input line " + itos(lnCnt + 3) + ": More quality positions than nucleotides detected for sequence\n " + tdn->getId());
	}
	lnCnt+=4;
	corrupt = false;
	return tdn;
}

void InputStreamer::minmaxQscore(shared_ptr<DNA> t, bool& print) {
	fqverMTX.lock();
	vector<qual_score> Qs = t->getQual();
	for (size_t i = 0; i < Qs.size(); i++) { minmaxQscore(Qs[i], print); }
	fqverMTX.unlock();
}

qual_score InputStreamer::minmaxQscore(qual_score t,bool& print) {
	if (t < 0) {
		if (fqSolexaFmt ){
            if (t < -5 && print) { cerr << "Unusually low Solexa quality score (" << t << "); setting to 0.\n"; print = false; }
		} else {
          if (fastQver == 64 && t >= -5) {
				fqSolexaFmt = true;
				if (print) { cerr << "Detected negative qualities in +64 range; enabling Solexa interpretation.\n"; print = false; }
			}
			else if (fastQver == 64 && t < -5) {
				fastQver = 33;
				fqSolexaFmt = false;
				if (print) { cerr << "Detected qualities incompatible with +64; switching to Phred+33.\n"; print = false; }
			}
			else { if (print) { cerr << "Unusually low quality score (" << t << "); setting to 0.\n"; print = false; } }
		}
		t = 0;
	}
	if (minQScore > t) { minQScore = t; }
	else if (maxQScore < t) { maxQScore = t; }
	return t;
}

bool InputStreamer::checkInFileStatus() {
	for (uint i = 0; i < 3; i++) {
		if (isFasta) {
			if (fasta_istreams[i] != NULL && !fasta_istreams[i]->eof()) { return true; }
		} else {
			if (fastq_istreams[i] != NULL && !fastq_istreams[i]->eof()) { return true; }
		}
	}
	return false;
}

void InputStreamer::allStreamReset() {
	std::unique_lock<std::shared_mutex> lock(protect);
	resetStats();
	cdbg("Resetting input streams");
	for (uint i = 0; i < 3; i++) {
		if (fasta_istreams[i] != NULL) { fasta_istreams[i]->reset(); }
		if (quality_istreams[i] != NULL ) { quality_istreams[i]->reset(); }
		if (fastq_istreams[i] != NULL ) { fastq_istreams[i]->reset();  }
	}
	cdbg("Done");
}

void InputStreamer::allStreamClose(){
	for (uint i = 0; i < 3; i++){
		if (fasta_istreams[i] != NULL) { delete fasta_istreams[i];  fasta_istreams[i] = NULL;}
		if (quality_istreams[i] != NULL) { delete quality_istreams[i];  quality_istreams[i] = NULL; }
		if (fastq_istreams[i] != NULL) { delete fastq_istreams[i];   fastq_istreams[i] = NULL; }
	}

    if (!isFasta && minQScore < SCHAR_MAX) { /* maxminQualWarns_fq() moved/removed with utilities */ }
}

void InputStreamer::jumpToNextDNA(bool& stillMore, int pos) {
	if (isFasta) {
		stillMore = read_fasta_entry(fasta_istreams[pos], quality_istreams[pos], dnaTemp1[pos], dnaTemp2[pos], lnCnt[pos]);
		dnaTemp1[pos] = dnaTemp2[pos];
		dnaTemp2[pos] = make_shared <DNA>("", "");
	} else {
		fastq_istreams[pos]->jumpLines(4);
		if (fastq_istreams[pos]->eof()) { stillMore = false; }
	}
}

vector<shared_ptr<DNA>> InputStreamer::getDNApairs() {
	vector<shared_ptr<DNA>> ret(3, nullptr);
	ret[0] = getDNA(0);

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

	N_DNApairsReturned++;
	
	return ret;
}


string InputStreamer::current_infiles() {
	string ret(""); int pos(0);
	if (isFasta) {
		fasta_istreams[pos]->getInFile();
	} else {
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
        if (!b1) {
			tmpO->setSize(k);
			if (safe) { protect.unlock(); }
			return false;
		}
		_localRdsRead++;
		if (numPairs > 1) {
			b2 = this->getDNAlines(tmpO->tmp[k][1], 1);
			if (!b2) {
				if (safe) { protect.unlock(); }
				cerr << "Problem: in file " << current_infiles() << " read2 ended before read1 (paired fastq not synchronized).\n";
				tmpO->setSize(k);
				return false;
			}
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
    if (ret.size() < 4) { ret.resize(4); }
	//vector<string>tmpLines(4, "");
	if (isFasta) {//get DNA from fasta + qual_ files
		//cerr << "Can't read fasta in multicore mode!\n";
		//exit(232);
		stillMore = fasta_istreams[pos]->getlines(ret[0],false);//header
		int lnRd = 1;//read 1 header line
		stillMore = fasta_istreams[pos]->getlines(ret[1], true);//fasta
		lnCnt[pos] += lnRd;
		if (quality_istreams[pos] != nullptr) {
			stillMore = quality_istreams[pos]->getlines(ret[2],false);//header
			stillMore = quality_istreams[pos]->getlines(ret[2], true);//qual
		}
	} else {// fastq format
		//fqRead
		stillMore = fastq_istreams[pos]->get4lines(ret);
		/*for (int xx = 0; xx < 4; xx++) {
			stillMore = fastq_istreams[pos]->getline(ret[xx]);
		}*/
        if (stillMore) {
			lnCnt[pos] += 4;
		}
     if (stillMore && fastQver == 0) {
		 fastQver = auto_fq_version(ret);
		 QverSet = true;			
		}
	}
	if (stillMore) {
		pairs_read[pos]++;
	}
	return stillMore;
}

shared_ptr<DNA> InputStreamer::getDNA(int pos){
    // acquire shared lock so multiple readers can proceed concurrently
	std::shared_lock<std::shared_mutex> lock(protect);
	//if (sync) {
	//	while (desync(pos)) {
	//		jumpToNextDNA(stillMore, pos);
	//	}
	//}
	bool stillMore = true;
	if (pos == 1 && numPairs <= 1) {
		return nullptr;
	}
	shared_ptr<DNA> ret(nullptr);
	bool corrupt(false); //corrupt state isn't implemented for fnaread
	bool repairInStream(false);
	if (isFasta) {//get DNA from fasta + qual_ files
		if (fasta_istreams[pos] == nullptr || fasta_istreams[pos]->eof()) {
			stillMore = false;
			return nullptr;
		}
		vector<string>tmpStr(3, "");
		stillMore = fasta_istreams[pos]->getlines(tmpStr[0],false);//header
		int lnRd = 1;//read 1 header line
		stillMore = fasta_istreams[pos]->getlines(tmpStr[1], true);//fasta
		lnCnt[pos] += lnRd;
		if (quality_istreams[pos] != nullptr) {
			stillMore = quality_istreams[pos]->getlines(tmpStr[2],false);//qual
			stillMore = quality_istreams[pos]->getlines(tmpStr[2], true);//qual
		}

		ret = str2DNA(tmpStr, keepPairHD, fastQver,pos);

		if (!stillMore || fasta_istreams[pos]->eof() || (!*(fasta_istreams[pos]))) {
			if (ret != nullptr) { if (!ret->seal() || ret->isEmpty()) { ret = nullptr; } } //delete ret;
			stillMore = false; 
		}
		else if (ret == nullptr || !ret->seal() || ret->isEmpty()) {
			corrupt = true;
		}
	}
	else { //fqRead
		if (fastq_istreams[pos] == nullptr || fastq_istreams[pos]->eof()) {
			stillMore = false;
			return nullptr;
		}
		vector<string>tmpLines(4, "");
		stillMore = getDNAlines(tmpLines,pos);
		if (stillMore) {
			ret = str2DNA(tmpLines, keepPairHD, fastQver,pos);
		}

		if (!stillMore || fastq_istreams[pos]->eof()) {
			if (ret != nullptr) { if (!ret->seal() || ret->isEmpty()) { ret = nullptr; } }
			stillMore = false;
		} else if (ret == nullptr || !ret->seal() || ret->isEmpty()) {
			corrupt = true;
		}

		if (ret != nullptr && !keepPairHD) {//better cut in early
			ret->setHeader(ret->getPositionFreeId());
		}
		if (corrupt) {
			ret = nullptr;
			repairInStream = true;
		}
	}
	if (!corrupt && stillMore) {
		pairs_read[pos]++;
	}
	N_DNAreturned++;

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
	if (pairs == -1) {
		pairs = (int)tfas.size();
		if (pairs > 1) { cerr << "Paired input (\",\" separated) detected\n"; }
	}
	numPairs = pairs;
	string p1(""), p2(""), midp("");
	string xtraMsg = "";
	if (getCurFileN() > 0) {
		xtraMsg = " " + itos(getCurFileN()) + " of " + itos(totalFileNumber);
		if (BCnumber > 1) {
			xtraMsg = ", looking for " + itos(BCnumber) + "BCs.\n";
		}
		else {
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
	}
	else if (tfas.size() == 2) {
		p1 = path + tfas[0];
		p2 = path + tfas[1];
		file1 = tfas[0];
	}
	else if (tfas.size() == 1) {
		p1 = path + fileS;
		file1 = fileS;
		//cerr << "Reading fastq " << p1 << endl;
	}
	if (subsPairs == "1") {
		p2 = "";
	}
	else 	if (subsPairs == "2") {
		p1 = p2; p2 = "";
	}

	if (simu) {

		return fileExists(p1) && fileExists(p2) && fileExists(midp);
	}


	if (verbose && !p1.empty() && !p2.empty()) {
		if (!midp.empty()) {
			//cerr << "Reading paired fastq + MID file" << xtraMsg<<"."<<endl;
		}
		else {
			//cerr << "Reading paired fastq" << xtraMsg << "." << endl;
			//cerr << p1 << " + " << p2 << endl;
		}
	}
	else {
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


bool InputStreamer::setupFastaQual2(string p1, string p2, string file_type) {
	if (p1.empty()) return false;
	int pos = 0;
	if (fasta_istreams[0] != nullptr) pos = 1;

	fasta_istreams[pos] = DBG_NEW ifbufstream(p1, INPUT_BUFFER_SIZE, doTIO);
	if (fasta_istreams[pos]->eof()) {
		cerr << "\nWarning: Could not open or empty " << file_type << " " << p1 << " — skipping this input.\n";
		fasta_istreams[pos] = nullptr;
		return false;
	}

	if (!p2.empty()) {
		quality_istreams[pos] = DBG_NEW ifbufstream(p2, INPUT_BUFFER_SIZE, doTIO);
		if (quality_istreams[pos]->eof()) {
			cerr << "\nWarning: Could not open or empty quality file " << p2 << " — continuing without quality stream.\n";
			quality_istreams[pos] = nullptr;
			qualAbsent = true;
		}
		else {
			qualAbsent = false;
		}
	}
	else {
		quality_istreams[pos] = nullptr;
		qualAbsent = true;
	}

	return true;
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
	string file_type = "";
	size_t bufS = INPUT_BUFFER_SIZE;

	if (!p1.empty()) {
		file_type = "fastq file 1";
		fastq_istreams[0] = DBG_NEW ifbufstream(p1, (size_t)round(bufS * 0.8), doTIO);
		if (fastq_istreams[0]->eof()) {
			cerr << "\nWarning: Could not open or empty " << file_type << " " << p1 << " — skipping this input.\n";
			fastq_istreams[0] = nullptr;
			return false;
		}
	}
	if (!p2.empty()) {
		file_type = "fastq file 2";
		fastq_istreams[1] = DBG_NEW ifbufstream(p2, (size_t)round(bufS * 1.2), doTIO);
		if (fastq_istreams[1]->eof()) {
			cerr << "\nWarning: Could not open or empty " << file_type << " " << p2 << " — skipping this input.\n";
			fastq_istreams[1] = nullptr;
			return false;
		}
	}
	if (!midp.empty()) {
		this->openMIDseqs("", midp);
	}
	return true;
}

string InputStreamer::setupInput(string path, int t, const string& uniqueFastxFile,
	filesStr& files, int& paired, string onlyPair,
	string& mainFilename, bool simulate) {
	string mainFilepath("");
	vector<string> fastqFiles = files.FastqF;
	vector<string> fastaFiles = files.FastaF;
	vector<string> qualityFiles = files.QualF;
	vector<string> midFiles = files.MIDfq;

	if (isFasta) {
		if (fastaFiles[t] != uniqueFastxFile) {
			cerr << "Error in matching FASTA target filenames.\n";
			exit(11);
		}
		this->setupFastaQual(path, fastaFiles[t], qualityFiles[t], paired, onlyPair);
		mainFilepath = path + fastaFiles[t];
		mainFilename = fastaFiles[t];
	}
	else {
		if (fastqFiles[t] != uniqueFastxFile) {
			cerr << "Error in matching target filenames.\n";
			exit(11);
		}
		this->setupFastq(path, fastqFiles[t], paired, onlyPair, simulate, !simulate);
		mainFilepath = path + fastqFiles[t];
		mainFilename = fastqFiles[t];
	}
	if ((size_t)t < midFiles.size()) {
		this->openMIDseqs(path, midFiles[t]);
	}
	return mainFilepath;
}

void InputStreamer::setupFna(string in) {
	int paired = 1;
	if (detectSeqFmt(in) == "-i_fna") {
		setupFastaQual("", in, "", paired, "", false);
	}
	else {
		setupFastq("", in, paired, "", false, false);
	}
}

qual_score InputStreamer::auto_fq_version() {
	if (fastQver != 0) {
		return fastQver;
	}

	return (qual_score)auto_fq_version(minQScore, maxQScore);
}

qual_score InputStreamer::auto_fq_version(qual_score minQScoreIn, qual_score maxQScoreIn) {
  // minQScoreIn/maxQScoreIn are raw ASCII quality characters from FASTQ lines.
	// Prefer Phred+33 unless there is strong evidence for legacy +64.
	if ((int)minQScoreIn < 59) { return 33; }
	if ((int)maxQScoreIn > 74) { return 64; }
	return 33;
}
qual_score InputStreamer::auto_fq_version(const vector<string>& ret) {
	qual_score minAscii = SCHAR_MAX;
	qual_score maxAscii = 0;
	for (unsigned char ch : ret[3]) {
		qual_score q = static_cast<qual_score>(ch);
		if (q < minAscii) { minAscii = q; }
		if (q > maxAscii) { maxAscii = q; }
	}

	if (minAscii < minQScore) { minQScore = minAscii; }
	if (maxAscii > maxQScore) { maxQScore = maxAscii; }

	const bool modernIlluminaHeader =
		(ret[0].find(" 1:N:") != string::npos) ||
		(ret[0].find(" 2:N:") != string::npos) ||
		(ret[0].find(" 1:Y:") != string::npos) ||
		(ret[0].find(" 2:Y:") != string::npos);

	if (modernIlluminaHeader) {
		fastQver = 33;
		fqSolexaFmt = false;
	}
	else {
		fastQver = static_cast<qual_score>(auto_fq_version(minQScore, maxQScore));
	}
	return fastQver;
}


bool InputStreamer::setupFastaQual(string path, string fasta, string qual, int& pairs, string subsPairs, bool simu) {
	allStreamClose();
	resetStats();
	vector<string> fna = splitByCommas(fasta);
	vector<string> qua = qual.empty() ? vector<string>(0) : splitByCommas(qual);

	if (pairs == -1) {
		pairs = (int)fna.size();
	}
	numPairs = pairs;

	fastaFilepathTemp[0] = (fna.empty() ? "" : path + fna[0]);
	fastaFilepathTemp[1] = (fna.size() > 1 ? path + fna[1] : "");
	qualityFilepathTemp[0] = (qua.size() > 0 ? path + qua[0] : "");
	qualityFilepathTemp[1] = (qua.size() > 1 ? path + qua[1] : "");

	if (subsPairs == "1") {
		fastaFilepathTemp[1].clear();
		qualityFilepathTemp[1].clear();
	}
	else if (subsPairs == "2") {
		fastaFilepathTemp[0] = fastaFilepathTemp[1];
		qualityFilepathTemp[0] = qualityFilepathTemp[1];
	}

	if (simu) {
		bool ok = !fastaFilepathTemp[0].empty() && fileExists(fastaFilepathTemp[0]);
		if (!qualityFilepathTemp[0].empty()) ok = ok && fileExists(qualityFilepathTemp[0]);
		if (pairs == 2 && !fastaFilepathTemp[1].empty()) ok = ok && fileExists(fastaFilepathTemp[1]);
		if (pairs == 2 && !qualityFilepathTemp[1].empty()) ok = ok && fileExists(qualityFilepathTemp[1]);
		return ok;
	}

	bool suc = setupFastaQual2(fastaFilepathTemp[0], qualityFilepathTemp[0]);
	if (pairs == 2 && !fastaFilepathTemp[1].empty()) {
		suc = suc && setupFastaQual2(fastaFilepathTemp[1], qualityFilepathTemp[1], "fasta file 2");
	}
	return suc;
}
