#include "ReadMerger.h"

#include <algorithm>


void mergeStats::printHisto(string File) {
	int maxVS = max((int)matchFwd.size(), (int)matchRev.size());
	prepStats();
	ofstream temp;
	temp.open(File.c_str(), ios::out);
	if (!temp) { cerr << "Could not open outstream to read merger stat file:\n" << File << endl; exit(478); }

	temp << "Pos\tFracFwd\tmismatchFwd\tmatchFwd\tFracRev\tmismatchRev\tmatchRev\n";
	for (size_t i = 0; i < (size_t)maxVS; i++) {
		temp << i << "\t";
		if (i < matchFwd.size()) {
			temp << percFwd[i] << "\t" << mismatchFwd[i] << "\t" << matchFwd[i] << "\t";
		}
		else { temp << "\t\t\t"; }
      if (i < matchRev.size()) {
			temp << percRev[i] << "\t" << mismatchRev[i] << "\t" << matchRev[i] << "\n";
		}
		else { temp << "\t\t\n"; }
	}
	temp.close();

}
void mergeStats::prepStats() {
	percFwd.resize(matchFwd.size(), 0.f); percRev.resize(mismatchRev.size(), 0.f);
	maxF = 0.f; maxR = 0.f;
	for (size_t i = 0; i < matchFwd.size(); i++) {
		const float observations = mismatchFwd[i] + matchFwd[i];
		percFwd[i] = observations > 0.f ? mismatchFwd[i] / observations : 0.f;
		if (percFwd[i] > maxF) { maxF = percFwd[i]; }
	}
	//normalize to max val
	//for (size_t i = 0; i < matchFwd.size(); i++) { percFwd[i] /= maxF; }
	//same game for rev scores
	for (size_t i = 0; i < mismatchRev.size(); i++) {
		const float observations = mismatchRev[i] + matchRev[i];
		percRev[i] = observations > 0.f ? mismatchRev[i] / observations : 0.f;
		if (percRev[i] > maxR) { maxR = percRev[i]; }
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
		if ((int)mismatchFwd.size() <= p1) { mismatchFwd.resize(p1 + 1, 0); matchFwd.resize(p1 + 1, 0); }
		if (same) { matchFwd[p1]++; }
		else { mismatchFwd[p1]++; }
	}
	else {
		//int i2 = overlap - p2 - 1;
		if ((int)mismatchRev.size() <= p2) { mismatchRev.resize(p2 + 1, 0); matchRev.resize(p2 + 1, 0); }
		if (same) { matchRev[p2]++; }
		else { mismatchRev[p2]++; }
	}
}


void qualStats::printHisto(string File) {
	ofstream temp;
	temp.open(File.c_str(), ios::out);
	if (!temp) { cerr << "Could not open outstream to read merger stat file:\n" << File << endl; exit(478); }
	temp << "Pos\tAvgQ_r1\tAvgQ_r2\n";
	int maxVS = max((int)r1.size(), (int)r2.size());
	//int q11 = int(r1[0]);
	for (size_t i = 0; i < (size_t)maxVS; i++) {
		temp << i <<"\t";
		if (i < r1.size()) {
			temp <<  (N1[i] ? static_cast<double>(r1[i]) / N1[i] : 0.0)  << "\t";
		}
		else { temp << "\t"; }
		if (i < r2.size()) {
			temp << (N2[i] ? static_cast<double>(r2[i]) / N2[i] : 0.0) << "\n";
		}
		else { temp << "\n"; }
	}


	temp.close();
}
void qualStats::logQuals(const vector<qual_score>& q1, const vector<qual_score>& q2) {
	if (r1.size() < q1.size()) { r1.resize(q1.size(), 0);N1.resize(q1.size(), 0);}
	if (r2.size() < q2.size()) { r2.resize(q2.size(), 0); N2.resize(q2.size(), 0);}
	for (size_t i = 0; i < q1.size(); i++) {
		r1[i] += static_cast<uint64_t>(std::max<int>(0, q1[i])); N1[i]++;
	}
	for (size_t i = 0; i < q2.size(); i++) {
		r2[i] += static_cast<uint64_t>(std::max<int>(0, q2[i])); N2[i]++;
	}
}




//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////



qual_score ReadMerger::getMergedQual(qual_score q1, qual_score q2, bool same) {
    static const auto table = [] {
        vector<qual_score> values((maxSAqualP + 1) * (maxSAqualP + 1) * 2, 0);
        ReadMerger calculator;
        for (int left = 0; left <= maxSAqualP; ++left) {
            for (int right = 0; right <= maxSAqualP; ++right) {
                const double p1 = qualToProb(static_cast<qual_score>(left));
                const double p2 = qualToProb(static_cast<qual_score>(right));
                for (int match = 0; match < 2; ++match) {
                    const double probability = match
                        ? calculator.mergeQProbabilities(p1, p2)
                        : calculator.mismatchQProbability(p1, p2);
                    const double quality = -10 * log10(probability) + 0.5;
                    values[(left * (maxSAqualP + 1) + right) * 2 + match] =
                        static_cast<qual_score>(std::clamp(quality, 0.0, 40.0));
                }
            }
        }
        return values;
    }();
    const int left = std::clamp<int>(q1, 0, maxSAqualP);
    const int right = std::clamp<int>(q2, 0, maxSAqualP);
    return table[(left * (maxSAqualP + 1) + right) * 2 + static_cast<int>(same)];
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
  if (!str.empty()) {
		reverseComplement(&str[0], (int)str.length());
	}
	return str;
}

double ReadMerger::percentIdentity(const char* sequence1, const char* sequence2, 
	int length, int mismatches ) {
	int match_count = 0;
	int mismatch_count = 0;
	if (mismatches < 0) mismatches = length;
	if (length <= 0) return 0.f;

	for (int i = 0; i < length; ++i) {
		if (sequence1[i] == sequence2[i]) {
			++match_count;
		}
		else if (++mismatch_count > mismatches) {
			return -1.f;
		}
	}
    return (double)match_count / length;
}

double ReadMerger::percentIdentity(std::string_view sequence1, std::string_view sequence2, 
	int mismatches ) {
	int match_count = 0;
	int mismatch_count = 0;
 const int length = static_cast<int>(sequence1.length());

	if (sequence1.length() != sequence2.length()) {
		std::cerr << "Sequences must be of the same length." << std::endl;
		return -1.f;
	}

	if (mismatches < 0) { mismatches = length; }
	if (length <= 0) return 0.f;

 for (int i = 0; i < length; ++i) {

      if (sequence1[i] == sequence2[i]) {
			++match_count;
		}
		else if (++mismatch_count > mismatches) {
			return -1.f;
		}
	}
    return (double)match_count / length;
}
void ReadMerger::reverseComplement(char *str, int len) {
 if (str == nullptr || len <= 0) {
		return;
	}
	char *p1 = str;
	char *p2 = str + len - 1;

	while (p1 < p2) {
		char tmp = complement(*p1);
		*p1++ = complement(*p2);
		*p2-- = tmp;
	}
   if (p1 == p2) {
		*p1 = complement(*p1);
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

/*
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

*/


///// alll important initial routine to find the best place to merge
//////////////////////////////////////////
bool ReadMerger::findSeed(std::string& sequence1, std::string& sequence2, MergeResult &result) {
   const size_t seq1_len = sequence1.length();
	const size_t seq2_len = sequence2.length();
	seedmap_.clear();
   if (seq1_len < seed_size_ || seq2_len < seed_size_) {
		return false;
	}
	seedmap_.reserve(seed_count_ * (check_reverse_complement_ ? 4u : 2u));
	const char* sequence1_ptr = sequence1.c_str();
	const char* sequence2_ptr = sequence2.c_str();
	const char* reverse_sequence2_ptr = nullptr;

  if (check_reverse_complement_) {
		reverse_complement_tmp_ = sequence2;
		reverseComplement(&reverse_complement_tmp_[0], static_cast<int>(reverse_complement_tmp_.size()));
		reverse_sequence2_ptr = reverse_complement_tmp_.c_str();
	}

	// The max possible overlap
   auto end = std::min(seq1_len, seq2_len) - seed_size_;

	// add seeds to map
	// Take seed from read 2 (both ends) and if set, also from the reverse complement
   for (size_t seed_num = 0; seed_num < seed_count_; ++seed_num) {
		auto offset = seed_margin_ + seed_positions_[seed_num];// + seed_num * seed_dist_; //
        if (offset >= end) break;
		const size_t forward_pos = offset;
		const size_t reverse_pos = sequence2.length() - offset - seed_size_;
        seedmap_.insert({ SeedStr(sequence2_ptr + forward_pos, seed_size_), static_cast<int>(forward_pos) });
		seedmap_.insert({ SeedStr(sequence2_ptr + reverse_pos, seed_size_), static_cast<int>(reverse_pos) });
        if (check_reverse_complement_) {
			seedmap_.insert({ SeedStr(reverse_sequence2_ptr + forward_pos, seed_size_), -static_cast<int>(forward_pos) - 1 });
			seedmap_.insert({ SeedStr(reverse_sequence2_ptr + reverse_pos, seed_size_), -static_cast<int>(reverse_pos) - 1 });
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
	const int seq1_len_i = static_cast<int>(seq1_len);
	const int seq2_len_i = static_cast<int>(seq2_len);
	for (size_t pos1 = seed_margin_; pos1 <= seq1_len - seed_size_; ++pos1) {
		pos2 = -1;
		bool found = false;

        SeedStr seed1(sequence1_ptr + pos1, seed_size_);
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

			pos2 = rev ? (-pos2 - 1) : pos2;

          int offset = pos2 - static_cast<int>(pos1);

			offset1 = (offset >= 0) * offset;
			offset2 = (offset < 0) * offset * -1;
            if (offset >= 0) {
				overlap = std::min(seq2_len_i - offset, seq1_len_i);
			}
			else {
                overlap = std::min(seq1_len_i + offset, seq2_len_i);
			}

			if (overlap < (int)min_overlap_) { offset1 = -1; continue; }

          auto pi = percentIdentity(
				sequence1_ptr + offset2,
				(rev ? reverse_sequence2_ptr : sequence2_ptr) + offset1, overlap,
				overlap / 10);
			if (pi > percent_identity_threshold_) {
                result.seed.pos1 = static_cast<size_t>(pos1);
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
    if (dna1) { dna1->clearMergeSeed(); }
    if (dna2) { dna2->clearMergeSeed(); }
    if (!dna1 || !dna2 || dna1->mem_length() < 20 || dna2->mem_length() < 20 ||
        dna1->getQual().size() != dna1->mem_length() ||
        dna2->getQual().size() != dna2->mem_length()) {
        return false;
    }
    // Logical mate-quality trims retain the original overlap evidence. Both
    // seed search and construction use this same backing sequence, after any
    // physical barcode/primer cuts. The consensus is filtered separately.
    MergeResult res;
    if (!findSeed(dna1->getSequence(), dna2->getSequence(), res)) {
        return false;
    }
    dna1->merge_seed_pos_ = static_cast<int>(res.seed.pos1);
    dna1->merge_offset_ = res.offset1;
    dna2->merge_seed_pos_ = static_cast<int>(res.seed.pos2);
    dna2->merge_offset_ = res.offset2;
    dna2->reversed_merge_ = res.seed.is2reversed;
    return true;
}

shared_ptr<DNA> ReadMerger::merge(shared_ptr<DNA> read1, shared_ptr<DNA> read2) {
    if (!read1 || !read2 || read1->merge_seed_pos_ < 0 || read2->merge_seed_pos_ < 0 ||
        read1->merge_offset_ < 0 || read2->merge_offset_ < 0) {
        return nullptr;
    }
    const string& seq1 = read1->getSequence();
    const string& storedSeq2 = read2->getSequence();
    const auto& qual1 = read1->getQual();
    const auto& qual2 = read2->getQual();
    const size_t len1 = seq1.size(), len2 = storedSeq2.size();
    const size_t start1 = static_cast<size_t>(read1->merge_offset_);
    const size_t start2 = static_cast<size_t>(read2->merge_offset_);
    if (len1 == 0 || len2 == 0 || qual1.size() != len1 || qual2.size() != len2 ||
        (start1 != 0 && start2 != 0) || start1 >= len2 || start2 >= len1) {
        return nullptr;
    }
    const size_t overlapStart = max(start1, start2);
    const size_t overlapEnd = min(start1 + len1, start2 + len2);
    if (overlapEnd <= overlapStart) { return nullptr; }

    const bool reversed = read2->reversed_merge_;
    // Dereplication shares a merger across workers. Scratch storage must be
    // local to a worker, and neither input's bases nor qualities may be changed.
    thread_local string reverseSeq2;
    std::string_view seq2 = storedSeq2;
    if (reversed) {
        reverseSeq2 = storedSeq2;
        reverseComplement(reverseSeq2.data(), static_cast<int>(reverseSeq2.size()));
        seq2 = reverseSeq2;
    }
    const auto quality2 = [&](size_t pos) {
        return qual2[reversed ? len2 - 1 - pos : pos];
    };

    // For inward-facing mates, R1's 5' end and R2's 5' end delimit the
    // amplicon. Exclude read-through outside those ends, including primers or
    // barcodes already removed from the other mate. Forward overlaps use the union.
    const size_t outputStart = reversed ? start1 : 0;
    const size_t outputEnd = reversed ? start2 + len2 : max(start1 + len1, start2 + len2);
    if (outputEnd <= outputStart) { return nullptr; }
    string sequence(outputEnd - outputStart, 'N');
    vector<qual_score> qualities(sequence.size());
    int errors = 0, mismatchQualSum = 0;

    for (size_t pos = outputStart; pos < outputEnd; ++pos) {
        const bool in1 = pos >= start1 && pos - start1 < len1;
        const bool in2 = pos >= start2 && pos - start2 < len2;
        const size_t out = pos - outputStart;
        if (!in1 && !in2) { return nullptr; }
        if (!in2) {
            sequence[out] = seq1[pos - start1];
            qualities[out] = qual1[pos - start1];
        } else if (!in1) {
            sequence[out] = seq2[pos - start2];
            qualities[out] = quality2(pos - start2);
        } else {
            const size_t p1 = pos - start1, p2 = pos - start2;
            const char base1 = seq1[p1], base2 = seq2[p2];
            const qual_score q1 = qual1[p1], q2 = quality2(p2);
            const bool same = base1 == base2;
            qualities[out] = getMergedQual(q1, q2, same);
            if (same) {
                sequence[out] = base1;
            } else {
                ++errors;
                mismatchQualSum += min(q1, q2);
                const bool canonical1 = canonicalDNA(base1), canonical2 = canonicalDNA(base2);
                sequence[out] = (q1 >= q2)
                    ? ((canonical1 || !canonical2) ? base1 : base2)
                    : ((canonical2 || !canonical1) ? base2 : base1);
            }

        }
    }
    if (b_takeStats.load(std::memory_order_relaxed)) {
        // Consensus construction remains parallel. Commit one read's statistics
        // together, since dereplication shares this merger between workers.
        std::lock_guard<std::mutex> guard(statsMutex);
        qS.logQuals(qual1, qual2);
        for (size_t pos = max(outputStart, overlapStart); pos < min(outputEnd, overlapEnd); ++pos) {
            const size_t p1 = pos - start1, p2 = pos - start2;
            const size_t originalP2 = reversed ? len2 - 1 - p2 : p2;
            mS.logDistri(static_cast<int>(p1), static_cast<int>(originalP2),
                static_cast<int>(overlapEnd - overlapStart), seq1[p1] == seq2[p2]);
        }
    }
    auto merged = make_shared<DNA>();
    merged->setSequence(std::move(sequence));
    merged->setQual(std::move(qualities));
    merged->setHeader(read1->getId());
    merged->getEssentialsFts(read1);
    if (reversed && read2->getRevPrimDetect()) { merged->setRevPrimCut(); }
    const qual_score meanMismatch = errors == 0 ? 0
        : static_cast<qual_score>(round(static_cast<double>(mismatchQualSum) / errors));
    merged->setMergeErrors(errors, meanMismatch);
    read1->setMergeErrors(errors, meanMismatch);
    read2->setMergeErrors(errors, meanMismatch);
    merged->stripLeadEndN();
    merged->setMergeLength(merged->length());
    read1->setMergeLength(merged->length());
    read2->setMergeLength(merged->length());
    return merged;
}

namespace {
template<class T> void addHistogram(vector<T>& target, const vector<T>& source) {
    target.resize(max(target.size(), source.size()), 0);
    for (size_t i = 0; i < source.size(); ++i) { target[i] += source[i]; }
}
}

void ReadMerger::addRMstats(ReadMerger* other) {
    if (!other || other == this) { return; }
    std::scoped_lock guard(statsMutex, other->statsMutex);
    addHistogram(mS.matchFwd, other->mS.matchFwd);
    addHistogram(mS.matchRev, other->mS.matchRev);
    addHistogram(mS.mismatchFwd, other->mS.mismatchFwd);
    addHistogram(mS.mismatchRev, other->mS.mismatchRev);
    addHistogram(qS.r1, other->qS.r1);
    addHistogram(qS.r2, other->qS.r2);
    addHistogram(qS.N1, other->qS.N1);
    addHistogram(qS.N2, other->qS.N2);
}
void ReadMerger::printMergeHisto() {
    std::lock_guard<std::mutex> guard(statsMutex);
    mS.printLogs();
}
void ReadMerger::printMergeHisto(string file, string) {
    std::lock_guard<std::mutex> guard(statsMutex);
    mS.printHisto(file);
}
void ReadMerger::printQualHisto(string file, string) {
    std::lock_guard<std::mutex> guard(statsMutex);
    qS.printHisto(file);
}
void ReadMerger::restStats() {
    std::lock_guard<std::mutex> guard(statsMutex);
    mS = mergeStats();
    qS = qualStats();
}
