#include "ReadMerger.h"


void mergeStats::printHisto(string File) {
	int maxVS = max((int)matchFwd.size(), (int)matchRev.size());
	if (maxVS <= 0) { return; }
	cout << "Print Merge error log of S="<< matchFwd.size()<<endl;
	if (matchFwd.size() == 0 && matchRev.size() == 0) { return; }
	prepStats();
	ofstream temp;
	temp.open(File.c_str(), ios::out);
	if (!temp) { cerr << "Could not open outstream to read merger stat file:\n" << File << endl; exit(478); }

	temp << "Pos\tFracFwd\tmismatchFwd\tmatchFwd\tFracRev\tmismatchRev\tmatchRev\n";
	for (int i = 0; i < maxVS; i++) {
		temp << i << "\t";
		if (i < matchFwd.size()) {
			temp << percFwd[i] << "\t" << mismatchFwd[i] << "\t" << matchFwd[i] << "\t";
		}
		else { temp << "\t\t\t"; }
		if (i < matchFwd.size()) {
			temp << percRev[i] << "\t" << mismatchRev[i] << "\t" << matchRev[i] << "\n";
		}
		else { temp << "\t\t\n"; }
	}
	temp.close();

}
void mergeStats::prepStats() {
	if (matchFwd.size() == 0) {	return;	}
	percFwd.resize(matchFwd.size(), 0.f); percRev.resize(mismatchRev.size(), 0.f);
	maxF = 0.f; maxR = 0.f;
	for (size_t i = 0; i < matchFwd.size(); i++) {
		percFwd[i] = mismatchFwd[i] / (mismatchFwd[i] + matchFwd[i]);
		if (percFwd[i] > maxF) { maxF = percFwd[i]; }
		percFwd[i] = -10 * log10(percFwd[i]);
	}
	//normalize to max val
	//for (size_t i = 0; i < matchFwd.size(); i++) { percFwd[i] /= maxF; }
	//same game for rev scores
	for (size_t i = 0; i < mismatchRev.size(); i++) {
		percRev[i] = mismatchRev[i] / (mismatchRev[i] + matchRev[i]);
		if (percRev[i] > maxR) { maxR = percRev[i]; }
		percRev[i] = -10 * log10(percRev[i]);
	}
	//normalize to max val
	//for (size_t i = 0; i < matchRev.size(); i++) { percRev[i] /= maxR; }

}

void mergeStats::printLogs() {
	if (matchFwd.size() == 0 && matchRev.size() == 0) { return; }
	prepStats();
	cout << "FwdRes::" << maxF << "\n";
	for (size_t i = 0; i < matchFwd.size(); i++) { cout << (int)(percFwd[i] * 100) << " "; }
	cout << "\n\nRevRes::" << maxR << "\n";
	for (size_t i = 0; i < mismatchRev.size(); i++) { cout << (int)(percRev[i] * 100) << " "; }
	cout << "\n\n";
}

void mergeStats::logDistri(int p1, int p2, int overlap, bool same) {
	if (p1>p2) {
		if (mismatchFwd.size() <= p1) { mismatchFwd.resize(p1 + 1, 0); matchFwd.resize(p1 + 1, 0); }
		if (!same) { mismatchFwd[p1]++; }
		matchFwd[p1]++;
	}
	else {
		//int i2 = overlap - p2 - 1;
		if (mismatchRev.size() <= p2) { mismatchRev.resize(p2 + 1, 0); matchRev.resize(p2 + 1, 0); }
		if (!same) { mismatchRev[p2]++; }
		matchRev[p2]++;
	}
}


void qualStats::printHisto(string File) {
	if (r1.size() == 0 && r2.size() == 0) { return; }
	ofstream temp;
	temp.open(File.c_str(), ios::out);
	if (!temp) { cerr << "Could not open outstream to read merger stat file:\n" << File << endl; exit(478); }
	temp << "Pos\tAvgQ_r1\tAvgQ_r2\n";
	int maxVS = max((int)r1.size(), (int)r2.size());
	//int q11 = int(r1[0]);
	for (int i = 0; i < maxVS; i++) {
		temp << i <<"\t";
		if (i < r1.size()) {
			temp <<  (float(int(r1[i])) /((float)N1[i]))  << "\t";
		}
		else { temp << "\t"; }
		if (i < r2.size()) {
			temp << (float(int(r2[i])) / ((float)N2[i])) << "\n";
		}
		else { temp << "\n"; }
	}


	temp.close();
}
void qualStats::logQuals(vector<qual_score> q1, vector<qual_score> q2) {
	if (r1.size() < q1.size()) { r1.resize(q1.size(), 0);N1.resize(q1.size(), 0);}
	if (r2.size() < q2.size()) { r2.resize(q2.size(), 0); N2.resize(q2.size(), 0);}
	for (size_t i = 0; i < q1.size(); i++) {
		r1[i] += (long)q1[i]; N1[i]++;
	}
	for (size_t i = 0; i < q2.size(); i++) {
		r2[i] += (long)q2[i]; N2[i]++;
	}
}




//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////



qual_score ReadMerger::getMergedQual(qual_score q1, qual_score q2,bool same) {
	auto p1 = qualToProb(q1);
	auto p2 = qualToProb(q2);

	//        std::cout << "p1: " << probToQual(p1) << "   p2: " << probToQual(p2) << "   p_new: ";
	double p_new(0);
	if (same) {
		p_new = mergeQProbabilities(p1, p2);
	} else {
		p_new = mismatchQProbability(p1, p2);
	}
	qual_score q_new = int(-10 * log10(p_new) + 0.5);//+0.5: round
//        std::cout << probToQual(p_new) << std::endl;
	return fastMin(q_new, (qual_score)40); //TODO figure out cap
//        return q_new;
}

 inline double ReadMerger::mismatchQProbability(double p1, double p2) {
	 if (p1 > p2) {
		 double tmp = p2;
		 p2 = p1;
		 p1 = tmp;
	 }
	 return (p1 * (1 - (p2 / 3)) / (p1 + p2 - 4 * p1 * p2 / 3));
 }
inline double ReadMerger::mergeQProbabilities(double p1, double p2) {
	//        std::cout << "(p1 * p2 / 3): " << (p1 * p2 / 3) << std::endl;
	//        std::cout << "(1 - p1 - p2 + (4 * p1 * p2 / 3)): " << (1 - p1 - p2 + (4 * p1 * p2 / 3)) << std::endl;

	return ((p1 * p2 / 3) / (1 - p1 - p2 + (4 * p1 * p2 / 3)));
}
string ReadMerger::reverseComplement(std::string str) {
	reverseComplement(const_cast<char *>(str.c_str()), (int)str.length());
	return str;
}

double ReadMerger::percentIdentity(const char* sequence1, const char* sequence2, 
	int length, int mismatches ) {
	int match_count = 0;
	int mismatch_count = 0;
	int base_count = 0;
	if (mismatches == -1) mismatches = length;

	for (int i = 0; i < length; ++i) {

		match_count += sequence1[i] == sequence2[i];
		mismatch_count += sequence1[i] != sequence2[i];
		if (mismatch_count == mismatches) return -1.f;
		++base_count;
	}
	return (double)match_count / base_count;
}

double ReadMerger::percentIdentity(std::string_view sequence1, std::string_view sequence2, 
	int mismatches ) {
	int match_count = 0;
	int mismatch_count = 0;
	int base_count = 0;

	if (sequence1.length() != sequence2.length()) {
		std::cerr << "Sequences must be of the same length." << std::endl;
		return -1.f;
	}

	if (mismatches == -1) { mismatches = (int)sequence1.length(); }

	for (size_t i = 0; i < sequence1.length(); ++i) {

		//match_count+=DNA_IUPAC[256*sequence1[i] + sequence2[i]];
		match_count += sequence1[i] == sequence2[i];
		//mismatch_count += sequence1[i] != sequence2[i];
		if ((i-match_count) >= mismatches) return -1.f;
		++base_count;
	}
	return (double)match_count / base_count;
}
void ReadMerger::reverseComplement(char *str, int len) {
	char *p1 = str;
	char *p2 = str + len - 1;

	while (p1 < p2) {
		char tmp = complement(*p1);
		*p1++ = complement(*p2);
		*p2-- = tmp;
	}
}

void ReadMerger::reverseStringInPlace(char *str, int len) {
	char *p1 = str;
	char *p2 = str + len - 1;

	while (p1 < p2) {
		char tmp = *p1;
		*p1++ = *p2;
		*p2-- = tmp;
	}
}

void ReadMerger::testMergeWithReads(std::istream &is1, std::istream &is2) {
	Benchmark merge_bm("Merge");
	merge_bm.start();
	BufferedFastxReader r1;
	BufferedFastxReader r2;

	FastxRecord record1;
	FastxRecord record2;

	MergeResult result;

	ReadMerger merger;

	const size_t batch_size = 100000;

	bool check_r1 = false;
	bool check_r2 = false;

	size_t counter = 0;
	size_t merge_counter = 0;

	size_t greater95 = 0;
	size_t smaller95 = 0;
	size_t zero = 0;

	while (true) {
		check_r1 = r1.LoadBatch(is1, batch_size);
		check_r2 = r2.LoadBatch(is2, batch_size);
		if (!check_r1 || !check_r2) break;

		while (r1.NextSequence(record1) && r2.NextSequence(record2)) {
			++counter;

			std::cout << record1.to_string() << std::endl;
            std::cout << record2.to_string() << std::endl;

			std::string sequence1 = "TACGTACTATCTTTCTCTCCTGAGTCGTACTGACGTAGTCGCGCGCGGCGCGCAGCTAGTCTACTGTACGTGCATGATCGTGACTGTACGTACTACTCGGTGTCGTACGTACGTAGTCAGTCAGTCACTGATCGTAGCTACAGCGACTTACGTAGTCATCGTACTATCCGATCGTGCAGTACGTACGTATATGCGTACTACTACTGATGCTACGTACGTACTAGTCAGTGTACTACGATGCAGC";
			std::string sequence2 = "GTCGCGCGCGGCGCGCAGCTAGTCTACTGTACGTGCATGATCGTGACTGTACGTACTACTCGGTGTCGTACGTACGTAGTCAGTCAGTCACTGATCGTAGCTACAGCGACTTACGTAGTCATCGTACTATCCGATCGTGCAGTACGTACGTATATGCGTACTACTACTGATGCTACGTACGTACTAGTCAGTGTACTACGATGCAGCGTGCTTGTCGTGTCGTCGTCGTCGTCGTCGCGTCGCTGTCGTGACG";

			std::cout << "test" << std::endl;
			std::cout << record1.sequence << std::endl;
            std::cout << record2.sequence << std::endl;

			auto res = merger.findSeed(record1.sequence, record2.sequence, result);

			std::cout << "res: " << res << std::endl;

			if (res) {
				++merge_counter;
			}
		}
	}
	merge_bm.stop();
	std::cout << counter << " reads" << std::endl;
	std::cout << greater95 << " > 95% id reads" << std::endl;
	std::cout << smaller95 << " < 95% id reads" << std::endl;
	std::cout << zero << " no seed found" << std::endl;
	std::cout << counter << " reads" << std::endl;
	std::cout << (double)greater95 / counter << " merged" << std::endl;

	std::cout << "merged: " << merge_counter << "/" << counter << std::endl;
	std::cout << "merged: " << (double)merge_counter / counter << "%" << std::endl;
	merge_bm.printResults();
}




///// alll important initial routine to find the best place to merge
//////////////////////////////////////////
bool ReadMerger::findSeed(std::string& sequence1, std::string& sequence2, MergeResult &result) {
	seedmap_.clear();
	int hit = 0;

	if (check_reverse_complement_)
		reverse_complement_tmp_ = reverseComplement(sequence2);

	// The max possible overlap
	auto end = std::min(sequence1.length(), sequence2.length()) - seed_size_;

	// add seeds to map
	// Take seed from read 2 (both ends) and if set, also from the reverse complement
	for (auto seed_num = 0; seed_num < seed_count_; ++seed_num) {
		auto offset = seed_margin_ + seed_positions_[seed_num];// + seed_num * seed_dist_; //
		if (sequence2.length() - offset == 0 || sequence2.length() - offset * -1 == 0) {
			//                std::cerr << "seq2 length: " << sequence2.length() << std::endl;
			//                std::cerr << "offset: " << offset << std::endl;
			break;
		}
		if (offset >= end) break;
		seedmap_.insert({ SeedStr(sequence2.c_str() + offset, seed_size_), offset });
		seedmap_.insert({ SeedStr(sequence2.c_str() + sequence2.length() - offset, seed_size_), sequence2.length() - offset });
		if (check_reverse_complement_) {
			seedmap_.insert({ SeedStr(reverse_complement_tmp_.c_str() + offset, seed_size_), offset * -1 });
			seedmap_.insert({ SeedStr(reverse_complement_tmp_.c_str() + reverse_complement_tmp_.size() - offset, seed_size_), (sequence2.length() - offset) * -1 });
			//                auto rev_seed1 = std::string_view(reverse_complement_tmp_.c_str() + offset, seed_size_);
			//                auto rev_seed2 = std::string_view(reverse_complement_tmp_.c_str() + reverse_complement_tmp_.size() - offset, seed_size_);
			//                seedmap_[rev_seed1] = offset * -1;
			//                seedmap_[rev_seed2] = (sequence2.length() - offset) * -1;
		}
	}

	int offset1(-1);
	int offset2(-1);
	int overlap(-1);
	int pos2(-1);

	// Now traverse through read one and look for a match in seedmap_ that indicates a potential seed.
	for (auto pos1 = seed_margin_; pos1 < sequence1.length() - seed_size_; ++pos1) {
		pos2 = -1;
		bool found = false;

		SeedStr seed1(sequence1.c_str() + pos1, seed_size_);
		//            std::string_view seed1(sequence1.c_str() + sequence1.length() - pos1 - seed_size_, seed_size_);

		auto find = seedmap_.find(seed1);
		if (find != seedmap_.end()) {
			found = true;
			pos2 = find->second;
		}
		//            if (!found && seedmap_robin_.find(seed2) != seedmap_robin_.end()) {
		//                found = true;
		//                pos2 = seedmap_robin_.at(seed2);
		//            }

		if (found) {
			////            if (seedmap_robin_.find(seed) != seedmap_robin_.end()) {
			////                pos2 = seedmap_robin_.at(seed);
			bool rev = pos2 < 0;

			pos2 = -1 * pos2 * rev + !rev * pos2;

			int offset = int(pos2 - pos1);

			offset1 = (offset >= 0) * offset;
			offset2 = (offset < 0) * offset * -1;
			//
			overlap =
				(offset >= 0) * int(std::min(sequence2.length() - offset, sequence1.length())) +
				(offset < 0) * int(std::min(sequence1.length() + offset, sequence2.length()));

			if (overlap < min_overlap_) { offset1 = -1; continue; }

			auto pi = percentIdentity(
				sequence1.c_str() + offset2,
				(rev ? reverse_complement_tmp_.c_str() : sequence2.c_str()) + offset1, overlap,
				overlap / 10);
			if (pi > percent_identity_threshold_) {
				result.seed.pos1 = pos1;
				result.seed.pos2 = pos2;
				result.seed.is2reversed = rev;
				result.percent_identity = pi;
				result.overlap = overlap;
				result.offset = offset;
				result.offset1 = offset1;
				result.offset2 = offset2;
				return true;
			}
		}
	}
	return false;
}


/*bool ReadMerger::findSeed(std::string& sequence1, std::string& sequence2) {
	return findSeed(sequence1, sequence2, result);
}
*/

bool ReadMerger::findSeedForMerge(shared_ptr<DNA> dna1, shared_ptr<DNA> dna2) {
	bool didMerge(false);
	if (dna1->length() < 20 || dna2->length() < 20) {//complete garbage sequence.. just ignore
		return false;
	}
	MergeResult res;
	if (findSeed(dna1->getSequence(), dna2->getSequence(), res) ) {
		dna1->merge_seed_pos_ = (int)res.seed.pos1;
		dna1->merge_offset_ = res.offset1;
		dna2->merge_seed_pos_ = (int)res.seed.pos2;
		dna2->merge_offset_ = res.offset2;
		dna2->reversed_merge_ = res.seed.is2reversed;
		return true;
	}
	return false;
}

shared_ptr<DNA> ReadMerger::merge(shared_ptr<DNA> read1, shared_ptr<DNA> read2) {
	if (read1->merge_seed_pos_ == -1 || read2->merge_seed_pos_ == -1) {
		return nullptr;
	}
	//if (read1->QualCtrl.MaxAmb || read2->QualCtrl.MaxAmb) {return nullptr;}
	//merg2MTX.lock();
	if (b_takeStats) {
		qS.logQuals(read1->getQual(), read2->getQual());
	}


	if (read2->reversed_merge_) {
		read2->reverse_transcribe(false);
	}
	// This checks if there are dovetails. When read 2 is reverse transcribed and read1 is offset, 
	//then the only constellation is that there are dovetails
	string Seq1 = read1->getSequence();
	string Seq2 = read2->getSequence();
	vector<qual_score> Qual1 = read1->getQual();
	vector<qual_score> Qual2 = read2->getQual();
	int seq1_length = (int)Seq1.length();
	int seq2_length = (int)Seq2.length();

	int offset1 = read1->merge_offset_;
	int offset2 = read2->merge_offset_;

	bool overlap_only = read2->reversed_merge_ && offset1;

	//offset1 >= 0 & read1->merge_offset_>0
	//offset1 > 0 
	int offset = int(read2->merge_seed_pos_ - read1->merge_seed_pos_);

	size_t overlap =
		(offset >= 0 ) * min(seq2_length - offset1, seq1_length)
		+ (offset < 0) * min(seq1_length - offset2, seq2_length);


	//	overlap =
	//		(offset >= 0) * int(std::min(sequence2.length() - offset, sequence1.length())) +
	//		(offset < 0) * int(std::min(sequence1.length() + offset, sequence2.length()));


	size_t overlap_start = (!overlap_only * ((bool)read1->merge_offset_) * offset1 +
		((bool)offset2) * offset2) + overlap_only * overlap;


	const size_t new_length = !overlap_only * max(read1->merge_offset_ + seq1_length, read2->merge_offset_ + seq2_length) + overlap_only * overlap;

	//        std::cout << "newlength> " << new_length << std::endl;
	//string new_seq = string(new_length + 1,'N');//
	char* new_seq = { new char[new_length + 1] };
	//char new_seq[new_length + 1];
	new_seq[new_length] = '\0';
	std::vector<qual_score> new_qual(new_length);

	int pos1 = 0;
	int pos2 = 0;
	int pos_overlap = 0;

	//        std::cout << "overlap start: " << overlap_start << std::endl;
			// merge first part of read (before overlap
	if (!overlap_only) {
		//can I copy over r1 or r2 into first part?
		//this is for the head
		if (read1->merge_offset_) {
			memcpy(new_seq, Seq2.c_str(), read1->merge_offset_);
			for (size_t i = 0; i < read1->merge_offset_; i++) {
				new_qual[i] = Qual2[i];
			}
			//            overlap_start = read1->merge_offset_;
		}
		else if (read2->merge_offset_) {
			memcpy(new_seq, Seq1.c_str(), read2->merge_offset_);
			for (size_t i = 0; i < read2->merge_offset_; i++) {
				new_qual[i] = Qual1[i];
			}
			//            overlap_start = read2->merge_offset_;
		}
		//this is for the tail
		size_t tail_start = overlap + overlap_start;
		if ((overlap + overlap_start) < read1->merge_offset_ + seq1_length) {
			pos1 = int(read2->merge_offset_ + overlap);
			for (size_t i = tail_start; pos1 < seq1_length && i < new_length; i++, ++pos1) {
				new_seq[i] = Seq1[pos1];
				new_qual[i] = Qual1[pos1];
			}
		}
		else {//read to tail
			pos2 = int(read1->merge_offset_ + overlap);
			for (size_t i = tail_start; i < new_length; i++, ++pos2) {
				new_seq[i] = Seq2[pos2];
				new_qual[i] = Qual2[pos2];
			}
		}
	}
	//        std::cout << "head: " << std::endl;
	//        std::cout << std::string(new_seq) << std::endl;

			// merge overlap
	pos1 = (int)read2->merge_offset_;
	pos2 = (int)read1->merge_offset_;
	pos_overlap = int(overlap_only ? 0 : overlap_start);

	//        std::cout << "overlap: " << overlap<< std::endl;
	//        std::cout << "pos1: " << pos1 << std::endl;
	//        std::cout << "pos2 " << pos2 << std::endl;
	//        std::cout << "pos_overlap " << pos_overlap << std::endl;

	//        std::cout << "s1ov " << std::string_view(Seq1.c_str() + pos1, overlap) << std::endl;
	//        std::cout << "s2ov " << std::string_view(Seq2.c_str() + pos2, overlap) << std::endl;


	//also log quality along read

	int errInOverlap(0);
	qual_score summedMismathcQ(0);

	for (size_t i = 0; i < overlap; i++, pos1++, pos2++, pos_overlap++) {
		//            std::cout << "___" << std::endl;
		char S1 = Seq1[pos1];
		char S2 = Seq2[pos2];


		bool sameNT = S1 == S2;

		//logging.. tmp
		if (b_takeStats) {
			mS.logDistri(pos1, (int)Qual2.size() - pos2, (int)overlap, sameNT);
		}

		new_qual[pos_overlap] = getMergedQual(
			Qual1[pos1], Qual2[pos2], sameNT);
		if (sameNT) {
			// same base at overlap position i
			new_seq[pos_overlap] = S1;

		}
		else { //Oh oh..
			errInOverlap++;
			summedMismathcQ += min(Qual1[pos1],Qual2[pos2]);

			bool S1canonical = canonicalDNA(S1);
			bool S2canonical = !canonicalDNA(S2);

			// Different base -> take the higher qual one
			// here qualities: the bigger the better (as opposed to probabilities)
			if (Qual1[pos1] > Qual2[pos2]
				&& (S1canonical || (!S1canonical && !S2canonical))) {
				// Qual read1 is better at overlap position i
				new_seq[pos_overlap] = S1;
			}
			else if (S2canonical) {
				// Qual read2 is better at overlap position i
				new_seq[pos_overlap] = S2;
			}
			else {
				new_seq[pos_overlap] = 'N';
			}
		}
		//            std::cout << "s1: " << Seq1[pos1] << " q: " << (int)Qual1[pos1] << "  -  s2: " << Seq2[pos2] << " q: " << (int)Qual2[pos2] << "  new: " << new_seq[pos_overlap] << " q: " << (int)new_qual[pos_overlap] << std::endl;
	}

	//        std::cout << "overlap: " << std::endl;
	//        std::cout << std::string(new_seq) << std::endl;

			// merge last part


	//        std::cout << "inmerge seq1:\n" << Seq1 << std::endl;

	//        std::cout << "tail: " << std::endl;
	//        std::cout << "new: " << std::string(new_seq) << std::endl;

	//        std::string merged_sequence(new_seq);
	//        std::string test = "IAMATESTSTRING";

	shared_ptr<DNA> merged_read = make_shared<DNA>();
	merged_read->setSequence(std::string(new_seq));
	delete[] new_seq;
	merged_read->setQual(move(new_qual));
	merged_read->setNewID(read1->getId());// +"_merged");
	merged_read->getEssentialsFts(read1);


	qual_score meanMisMatchQ = (qual_score)round((float)summedMismathcQ / (float)errInOverlap);
	merged_read->setMergeErrors(errInOverlap, meanMisMatchQ);
	read1->setMergeErrors(errInOverlap, meanMisMatchQ);
	read2->setMergeErrors(errInOverlap, meanMisMatchQ);

	merged_read->setMergeLength(merged_read->length());
	read1->setMergeLength(merged_read->length());
	read2->setMergeLength(merged_read->length());

	//backtranslate to make sure correct read again
	if (read2->reversed_merge_) {
		read2->reverse_transcribe(false);
	}

	//filter for Ns
	merged_read->stripLeadEndN();
	if (merged_read->numNonCanonicalDNA(true)) {
		merged_read= nullptr;
	}
	//merg2MTX.unlock();

	return merged_read;
}