#pragma once
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


//input through combined getDNApairs
#define togRe//ad
//input through getDNAlines
#define IOsep
#ifndef _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_DEPRECATE
#endif

//#include "include/iowrap.h"
//#include "FastxReader.h"
#include <functional>
#include <cctype>
#include <locale>
#include <climits>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <string_view>
#include <shared_mutex>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <memory>
#include <atomic>
#include <unordered_map>

#include "ThreadPool.h"
#include "Benchmark.h"
#include "DNAconsts.h"

class DNA {
    friend class DNAunique;
public:
    struct WhoIsBetterMergeStats {
        uint64_t merged_compares = 0;
        uint64_t non_merged_compares = 0;
        uint64_t accumulated_merge_length = 0;
    };

    DNA(string seq, string names) : sequence_(seq), sequence_length_(sequence_.length()),
        id_(names), new_id_(names),
        qual_(0), sample_id_(-1), avg_qual_(-1.f),
        quality_sum_(0), accumulated_error_(0.),
        good_quality_(false), mid_quality_(false), failed_(false), reversed_(false),
        read_position_(-1),
        FtsDetected(),
        fastq_write_offset_(0), id_fixed_(false), tempFloat(0.f), merge_offset_(-1) {}
    DNA() : sequence_(""), sequence_length_(0), id_(""), new_id_(""), qual_(0),
        sample_id_(-1), avg_qual_(-1.f),
        quality_sum_(0), accumulated_error_(0.),
        good_quality_(false), mid_quality_(false), failed_(false), reversed_(false),
        read_position_(-1),
        FtsDetected(),
        fastq_write_offset_(0), id_fixed_(false), tempFloat(0.f), merge_offset_(-1) {
        sequence_.reserve(151);//set to resonable expectation
    }
    //starts with fastx record
    //DNA(FastxRecord*, qual_score & minQScore, qual_score & maxQScore, qual_score fastQver);
    //works directly with 4 lines in fastq format
    DNA(const vector<string>& fq, qual_score fastQver);//fastq input
    DNA(const vector<string>& fas);//this is for fasta only..
    // move-aware constructors to allow zero-copy handoff from temporary buffers
    DNA(vector<string>&& fq, qual_score fastQver);
    DNA(vector<string>&& fas);

    ~DNA() {
        //cout << "destruct" << endl;
    }

    bool operator==(const DNA& i) const {
        return i.getSeqPseudoView() == this->getSeqPseudoView();
    }
    bool operator==(const shared_ptr<DNA>& i) const {
        return i != nullptr && i->getSeqPseudoView() == this->getSeqPseudoView();
    }
    //something wrong with DNA object, just del all info
    void delself() {
        sequence_ = ""; sequence_length_ = 0; id_ = ""; new_id_ = ""; qual_.resize(0);
        good_quality_ = false; mid_quality_ = false; failed_ = true;
    }
    //~DNA(){}
    void appendSequence(const string& s) { sequence_ += s; sequence_length_ = sequence_.length(); }
    void appendQuality(const vector<qual_score>& q);
    void fixQ0(void);
    void setSequence(string& s) {
        sequence_ = s;
        sequence_length_ = sequence_.length();
    }    void setSequence(string&& s) {
        sequence_ = std::move(s);
        sequence_length_ = sequence_.length();
    }


    string& getSequence() { return sequence_; }
    const vector<qual_score>& getQual() const {
        return qual_;
    }
    const vector<qual_score> getQual(int sta, int end) const {
        if (sta == 0 && end == 0) {
            return qual_;
        }
        vector<qual_score>::const_iterator first = qual_.begin() + sta;
        vector<qual_score>::const_iterator last = qual_.begin() + end;
        vector<qual_score> newVec(first, last);
        return newVec;
    }

    string getSeqPseudo() {
        return sequence_.substr(0, sequence_length_);
    }
    std::string_view getSeqPseudoView() const {
        const size_t len = std::min(sequence_length_, sequence_.size());
        return std::string_view(sequence_.data(), len);
    }

    void setQual(vector<qual_score> Q) { qual_ = move(Q); avg_qual_ = -1.f; }
    //void setQual(vector<qual_score>&& Q) { qual_ = Q; avg_qual_ = -1.f; }

    void setAllQual(qual_score q) { for (size_t i = 0; i < qual_.size(); i++) { qual_[i] = q; } avg_qual_ = -1.f; }

    const string& getId() {
        if (id_fixed_) {
            return new_id_;
        }
        return id_;
    }

    shared_ptr<DNA> getDNAsubseq(int start, int end, string id);

    string getPositionFreeId(); // remove /1 /2 #1:0 etc
    void normalizePositionFreeId();
    const string& getOldId() { return id_; }
    string getShortId() { return id_.substr(0, getShorterHeadPos(id_)); }
    string getNewIDshort() { return new_id_.substr(0, getShorterHeadPos(new_id_)); }
    bool seal(bool isFasta = false);


    bool isEmpty() {
        if (id_.empty() && sequence_.empty()) {
            this->setPassed(false);
            return true;
        }
        return false;
    }

    int numACGT();
    // Fused single-pass: computes ambiguity count and checks homonucleotide runs
    // Returns ambiguity count via ambNTs_out, and false if a homo-NT run >= maxHomoNT is found
    bool scanSequenceChecks(int maxHomoNT, int& ambNTs_out);
    void stripLeadEndN();
    int numNonCanonicalDNA(bool);
    float getAvgQual();
    int getMedianQual();
    unsigned int getQsum() { return quality_sum_; }
    float qualWinfloat(uint, float, int&);


    float binomialFilter(int, float);
    //float qualWinfloat_hybr(int,float,int,float,int&);

    bool qualWinPos(unsigned int, float);
    bool qualAccumTrim(double d);
    int qualAccumulate(double d);

    double getAccumError() {
        if (accumulated_error_ == 0.f) {
            for (uint i = 0; i < qual_.size(); i++) {
                if (qual_[i] >= 0) {
                    accumulated_error_ += SAqualP[qual_[i]];
                }
            }
        }
        if (std::isinf((double)accumulated_error_)) {
            accumulated_error_ = 5.f;
        }
        return accumulated_error_;
    }
    int minQual() { int mq = 50; for (uint i = 0; i < qual_.size(); i++) { if (qual_[i] < mq) { mq = qual_[i]; } }return mq; }
    void ntSpecQualScores(vector<long>&, vector<long>&);
    void ntSpecQualScoresMT(vector<long>&, vector<long>&);

    //returns qual filtered length 
    inline uint length() { return (uint)sequence_length_; }
    //returns original length 
    inline uint mem_length() { return (uint)sequence_.length(); }
    bool cutSeq(int start, int stop = -1, bool = false);
    bool cutSeqPseudo(int start) { return cutSeq(start, -1, true); }
    bool HomoNTRuns(int);
    int HomoNTTrim(int);
    int matchSeq(string, int Err, int searchSpace, int startPos, bool exhaustive = false);
    void reverse_compliment(bool reset = true);
    int matchSeqRev(const string&, int Err, int searchSpace, int start = 0, bool = false);
    int matchSeq_tot(const string&, int, int, int&);
    void writeSeq(ostream&, bool singleLine = false);
    void writeSeq_derep(ostream&, bool singleLine = false);
    string writeSeq(bool singleLine = false);
    void writeSeq(string& seq, bool singleLine = false);
    void writeQual(ostream&, bool singleLine = false);
    string writeQual(bool singleLine = false);
    void writeQual(string&, bool singleLine = false);
    void writeFastQ(ostream&, bool = true);
    void writeFastQ_derep(ostream&, bool = true);
    string writeFastQ(bool = true);
    void writeFastQ(string&, bool = true);
    //void writeFastQ(ofbufstream&, bool = true);
    void writeFastQEmpty(ostream&);
    void setNewID(const string& x) { new_id_ = x; }
    void setNewID(string&& x) { new_id_ = std::move(x); }
    void setHeader(const string& x) { id_ = x; new_id_ = x; }
    void setHeader(string&& x) {
        id_ = std::move(x);
        new_id_ = id_;
    }
    void changeHeadPver(int ver);
    void setTA_cut(bool x) { FtsDetected.TA_cut = x; }
    bool getTA_cut() { return FtsDetected.TA_cut; }
    void setBarcodeCut() { FtsDetected.barcode_cut = true; FtsDetected.barcode_detected = true; }
    bool getBarcodeCut() { return FtsDetected.barcode_cut; }
    void setBarcodeDetected(bool x) { FtsDetected.barcode_detected = x; }
    bool getBarcodeDetected() { return FtsDetected.barcode_detected; }
    bool isMIDseq() { if (read_position_ == 3) { return true; } return false; }
    void setMIDseq(bool b) { if (b) { read_position_ = 3; } }
    void setpairFWD() { read_position_ = 0; }
    void setpairREV() { read_position_ = 1; }
    int getReadMatePos() { return (int)read_position_; }

    bool sameHeadPosFree(shared_ptr<DNA>);

    bool sameHead(const string&);
    //inline void reverseTranscribe();
    void setTempFloat(float i) { tempFloat = i; }
    float getTempFloat() { return tempFloat; }
    //void adaptHead(shared_ptr<DNA>,const int,const int);
    void failed() { good_quality_ = false; mid_quality_ = false; failed_ = true; }
    bool control() { if (qual_.size() == 0) { return false; }return true; }

    void setBCnumber(int i, int BCoff) {
        if (i < 0) {
            sample_id_ = i;
            FtsDetected.barcode_detected = false;
        }
        else {
            sample_id_ = i + BCoff;
            FtsDetected.barcode_detected = true;
        }
    }

    //always return BC tag IDX global (no local filter idx accounted for, use getBCoffset() to correct)
    int getBCnumber() const {
        return sample_id_;
    }

    void prepareWrite(int fastQver);
    void reset();
    void resetTruncation() { sequence_length_ = sequence_.length(); }
    void setPassed(bool b);
    void setYellowQual(bool b) {
        mid_quality_ = b; failed_ = false;
    }
    bool isFailed() {
        return failed_;
    }
    bool isGreenQual(void) {
        if (mid_quality_) { return false; }
        return good_quality_;
    }
    bool isYellowQual(void) {
        if (good_quality_) { return false; }
        return mid_quality_;
    }
    bool isGreenYellowQual(void) {
        return (mid_quality_ || good_quality_);
    }
    // Return a non-owning view into the underlying sequence to avoid heap allocations
    std::string_view getSubSeq(int sta, int sto) const {
        if (sta < 0) sta = 0;
        if (sto < 0) sto = 0;
        size_t maxLen = sequence_.size();
        if ((size_t)sta >= maxLen) return std::string_view();
        size_t len = static_cast<size_t>(sto);
        if (sta + (int)len > (int)maxLen) {
            len = maxLen - static_cast<size_t>(sta);
        }
        return std::string_view(sequence_.data() + sta, len);
    }
    // Return a non-owning view of id_ with paired suffix (/1,/2,.1,.2 or trailing space) removed
    std::string_view getPositionFreeIdView() const;
    void resetQualOffset(int off, bool solexaFmt);

    //control & check what happened to any primers (if)
    bool has2PrimersDetected() { return (FtsDetected.reverse && FtsDetected.forward); }
    bool getRevPrimDetect() { return FtsDetected.reverse; }
    bool getFwdPrimDetect() { return FtsDetected.forward; }
    void setRevPrimCut() { FtsDetected.reverse = true; }
    void setFwdPrimCut() { FtsDetected.forward = true; }
    void setDereplicated() { FtsDetected.dereplicated = true; }
    bool isDereplicated() { return FtsDetected.dereplicated; }
    void constellationPairRev(bool b) { FtsDetected.revPairConstellation = b; }
    bool isConstellationPairRev() { return FtsDetected.revPairConstellation; }
    int getMergeLength() { return FtsDetected.mergeLength; }

    //only used in pre best seed step
    //float getSeedScore() { return tempFloat; }
    //void setSeedScore(float i) { tempFloat = (float)i; }

    struct QualStats {
        uint32_t maxL : 1;
        uint32_t PrimerFwdFail : 1;
        uint32_t AvgQual : 1; //sAvgQual
        uint32_t HomoNT : 1;
        uint32_t HomoNTtrimmed : 1;
        uint32_t PrimerRevFail : 1;
        uint32_t minL : 1;
        uint32_t minLqualTrim : 1; //<-sMinQTrim trimmed due to quality
        uint32_t TagFail : 1;
        uint32_t MaxAmb : 1;
        uint32_t QualWin : 1;//sQualWin 
        uint32_t AccErrTrimmed : 1;
        uint32_t QWinTrimmed : 1;  // either of these makes bool Trimmed; 
        uint32_t fail_correct_BC : 1;
        uint32_t suc_correct_BC : 1;
        uint32_t failedDNAread : 1;
        //bool adapterRem; -> setTA_cut
        uint32_t RevPrimFound : 1;
        uint32_t BinomialErr : 1;
        uint32_t dblTagFail : 1;
        QualStats() :
            maxL(false), PrimerFwdFail(false), AvgQual(false), HomoNT(false), HomoNTtrimmed(false),
            PrimerRevFail(false), minL(false), minLqualTrim(false),
            TagFail(false), MaxAmb(false), QualWin(false),
            AccErrTrimmed(false), QWinTrimmed(false),
            fail_correct_BC(false), suc_correct_BC(false),
            failedDNAread(false), RevPrimFound(false),
            BinomialErr(false),
            dblTagFail(false)
        {}
    } QualCtrl;

protected:
    size_t getShorterHeadPos(const string& x, int fastQheadVer = -1);
    //mainly used to mark if rev/Fwd primer was detected
    string xtraHdStr();
    void removeNonCanonicalNTs();
    void removeNegativeQuals();
    size_t getSpaceHeadPos(const string& x);
    //binomial accumulated error calc
    inline float interpolate(int errors1, float prob1, int errors2, float prob2, float alpha);
    float sum_of_binomials(int j, int k, float n, int qual_length, const vector<float>& error_probs, const vector<float>& per_position_accum_probs);
    inline float prob_j_errors(float p, float j, float n);
    inline bool matchDNA(char, char);

    std::string sequence_;
    size_t sequence_length_;
    string id_, new_id_; //original and newly constructed id_
    vector<qual_score> qual_;

public:

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
        mid_quality_ = oD->mid_quality_;
        failed_ = oD->failed_;
        sample_id_ = oD->sample_id_;

        FtsDetected.errInOverlap = oD->FtsDetected.errInOverlap;
        FtsDetected.meanQInOverlapMismatch = oD->FtsDetected.meanQInOverlapMismatch;
        FtsDetected.mergeLength = oD->FtsDetected.mergeLength;

    }

    void setMergeErrors(int eoi, qual_score meanQeoi) {
        FtsDetected.errInOverlap = eoi;
        FtsDetected.meanQInOverlapMismatch = meanQeoi;
    }
    void setMergeLength(int x) { FtsDetected.mergeLength = x; }
    int getMergeErrors() { return FtsDetected.errInOverlap; }
    qual_score getMergeErrorsQual() { return FtsDetected.meanQInOverlapMismatch; }
    bool isReversed() { return reversed_; }

protected:
    int sample_id_;

    //const char* DN;
    float avg_qual_;
    unsigned int quality_sum_;
    double accumulated_error_;
    uint8_t good_quality_ : 1;
    uint8_t mid_quality_ : 1;
    uint8_t failed_ : 1;
    uint8_t reversed_ : 1;


    short read_position_;//-1=unkown; 0=pair1 (fwd primer); 1=pair2 (rev primer); 3=MID seq ;

    qual_score fastq_write_offset_;

    struct ElementsDetection {
        uint8_t forward : 1;
        uint8_t reverse : 1;//primers detected
        uint8_t TA_cut : 1;
        uint8_t barcode_detected : 1;
        uint8_t barcode_cut : 1;
        uint8_t dereplicated : 1;
        uint8_t revPairConstellation : 1;
        //read merging related stats
        int errInOverlap;        qual_score meanQInOverlapMismatch;
        int mergeLength;

        ElementsDetection() : forward(false), reverse(false), TA_cut(false), barcode_detected(false),
            barcode_cut(false), dereplicated(false), revPairConstellation(false),
            errInOverlap(-1), meanQInOverlapMismatch(0), mergeLength(-1) {}
        void transferFrom(ElementsDetection& o) {
            forward = o.forward;
            reverse = o.reverse;
            TA_cut = o.TA_cut;
            barcode_detected = o.barcode_detected;
            barcode_cut = o.barcode_cut;
            dereplicated = o.dereplicated;
            revPairConstellation = o.revPairConstellation;

            errInOverlap = o.errInOverlap;
            mergeLength = o.mergeLength;
            meanQInOverlapMismatch = o.meanQInOverlapMismatch;
        }
        void reset() {
            forward = false;
            reverse = false;
            TA_cut = false;
            barcode_detected = false;
            barcode_cut = false;
            dereplicated = false;
            revPairConstellation = false;
            errInOverlap = -1;
            mergeLength = -1;
            meanQInOverlapMismatch = 0;
        }

    } FtsDetected;



    bool id_fixed_;
    float tempFloat;
};

using sample_count_t = uint32_t;
using read_occ_entry = std::pair<int, sample_count_t>;
typedef std::vector<read_occ_entry> read_occ;

struct DNAHasher
{
    size_t operator()(shared_ptr<DNA> k) const
    {
        return ((hash<std::string_view>()(k->getSeqPseudoView())) >> 1);
    }
};



class DNAunique : public DNA {//used for dereplication
    friend class DNA;
public:
    // Constructors and Destructors
    DNAunique() : DNA(), pair_(0), total_sum_cache_(0), total_sum_cache_valid_(false), derep_fastq_offset_(0) {}
    DNAunique(string s, string x) : DNA(s, x), total_sum_cache_(0), total_sum_cache_valid_(false), derep_fastq_offset_(0) {}

    // Mostly used constructor
    DNAunique(shared_ptr<DNA>d, int BC) : DNA(*d), pair_(0), total_sum_cache_(0), total_sum_cache_valid_(false), derep_fastq_offset_(0) {//best_seed_length_((uint)sequence_.size())
        incrementSampleCounter(BC);
    }
    ~DNAunique() = default;

    void recordWhoIsBetterCompare(int dMergeLength, int rMergeLength);
    const WhoIsBetterMergeStats& getWhoIsBetterMergeStats() const { return who_is_better_merge_stats_; }
    void resetWhoIsBetterMergeStats() { who_is_better_merge_stats_ = WhoIsBetterMergeStats(); }

    //string sequence_;    string id_;
    void Count2Head(bool);
    bool betterPreSeed(shared_ptr<DNA> d1, shared_ptr<DNA> d2);
    void matchedDNA(shared_ptr<DNA>, shared_ptr<DNA>, shared_ptr<DNA>, int, bool);
    void incrementSampleCounter(int sample_id);
    void writeMap(ofstream& os, const string&, vector<int>&, const vector<int>&, vector<uint64_t>&, vector<int>&);
    //inline int getCount() { return count_; }
    uint64_t totalSum();
    //uint getBestSeedLength() { return best_seed_length_; }
    //void setBestSeedLength(uint i) { best_seed_length_ = i; }
    void incrementSampleCounterBy(int sample_id, long count);
    void transferOccurence(shared_ptr<DNAunique>);
    const read_occ& getDerepMap() { return occurence; }
    vector<int> getDerepMapSort(size_t);
    //vector<pair_<int, int>> getDerepMapSort2(size_t wh);
    //void getDerepMapSort(vector<int>&, vector<int>&);

    using qual_accum_t = uint32_t;

    void addToQualSum(qual_accum_t& target, qual_accum_t value);

    void prepSumQuals();
    void sumQualities(const shared_ptr<DNAunique>& dna);
    void sumQualities(const shared_ptr<DNA>& dna);
    void releaseRawQualitiesIfPromoted();

    void prepareDerepQualities(int ofastQver);
    void writeDerepFastQ(ofstream&, bool = true);
    void writeFastQ(ostream&, bool = true);
    string writeFastQ(bool = true);
    void writeFastQ(string&, bool = true);

    void saveMem();

    void attachPair(shared_ptr<DNAunique> dna_unique);
    shared_ptr<DNAunique> getPair(void);

    void attachMerge(shared_ptr<DNAunique> dnamerge);
    shared_ptr<DNAunique> getMerge(void);


    //estimates if one sample occurence covers the unique counts required for sample specific derep min counts
    bool pass_deprep_smplSpc(const vector<int>&);
    //int counts() const { return count_; }

    void takeOver(shared_ptr<DNAunique> dna_unique_old);
    void takeOverDNA(shared_ptr<DNA> dna1, shared_ptr<DNA> dna2, shared_ptr<DNA> dnaMerge, bool keepQualities = true);
    void compactForFastaDerepStorage();
    std::vector<qual_accum_t> transferPerBaseQualitySum();

    std::vector<qual_accum_t> quality_sum_per_base_;

private:
    //int count_;
    //matrixUnit chimeraCnt;
    //int best_seed_length_;

    read_occ occurence;
    mutable uint64_t total_sum_cache_;
    mutable uint8_t total_sum_cache_valid_ : 1;

    inline void invalidateTotalSumCache() {
        total_sum_cache_valid_ = false;
    }

    void promoteQualitiesToSums(bool releaseRawQualities);

    read_occ::iterator findOccurrenceEntry(int sample_id);
    read_occ::const_iterator findOccurrenceEntry(int sample_id) const;

    shared_ptr<DNAunique> pair_;

    shared_ptr<DNAunique> merge_;

    // to calculate mean
    qual_score derep_fastq_offset_;
    WhoIsBetterMergeStats who_is_better_merge_stats_;

    //threadsafe
};




class DNAuniqSet {
public:
    DNAuniqSet() :bestDNU(nullptr), bestSet(false), bestHasMerge(false), totalCnts(0),
        cntsAdded2best(false) {}
    ~DNAuniqSet() {}

    void addNewDNAuniq(shared_ptr<DNA> dna, shared_ptr<DNA> dna2,
        shared_ptr<DNA> dnaM, int MrgPos1, int sample_id, bool b_derep_as_fasta_);
    shared_ptr<DNAunique>& operator[] (int x) {
        return DNUs[x];
    }
    map<int, shared_ptr<DNAunique>>::iterator find(const int x);
    const map<int, shared_ptr<DNAunique>>::iterator end();
    const map<int, shared_ptr<DNAunique>>::iterator begin();
    size_t size() { return DNUs.size(); }
    void setBest(bool addCnts);
    shared_ptr<DNAunique> best(bool addCnts);

    mutex lockMTX;
private:
    map<int, shared_ptr<DNAunique>> DNUs;
    shared_ptr<DNAunique> bestDNU;
    bool bestSet; bool bestHasMerge;
    int totalCnts;
    int cntsAdded2best;
};

// Hash map from sequence string to DNAuniqSet (store by value)
//typedef std::unordered_map<string, unique_ptr<DNAuniqSet>> HashDNA;
typedef robin_hood::unordered_node_map<string, unique_ptr< DNAuniqSet> > HashDNA;



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






shared_ptr<DNA> str2DNA(vector<string>& in, bool keepPairHD, int fastQver, int readpos);
// Move-aware overload to accept temporary buffers without copying
shared_ptr<DNA> str2DNA(vector<string>&& in, bool keepPairHD, int fastQver, int readpos);



