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

#include "diffDNA.h"
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
#include <algorithm>

#include "ThreadPool.h"
#include "Benchmark.h"
#include "DNAconsts.h"

struct CompactDNAReferences;
class CompactDNA;
class CompactDerepSet;

// A rarely set string (HTS BC:Z sample tags): one pointer instead of a 32-byte
// std::string in every read and every stored dereplication object.
class CompactOptionalString {
public:
    CompactOptionalString() = default;
    CompactOptionalString(const CompactOptionalString& o)
        : value_(o.value_ ? std::make_unique<std::string>(*o.value_) : nullptr) {}
    CompactOptionalString(CompactOptionalString&&) noexcept = default;
    CompactOptionalString& operator=(const CompactOptionalString& o) {
        if (this != &o) { value_ = o.value_ ? std::make_unique<std::string>(*o.value_) : nullptr; }
        return *this;
    }
    CompactOptionalString& operator=(CompactOptionalString&&) noexcept = default;
    CompactOptionalString& operator=(std::string s) {
        if (s.empty()) { value_.reset(); } else { value_ = std::make_unique<std::string>(std::move(s)); }
        return *this;
    }
    const std::string& str() const { return value_ ? *value_ : emptyString(); }
    operator const std::string&() const { return str(); }
    bool empty() const { return !value_; }
private:
    static const std::string& emptyString() { static const std::string empty; return empty; }
    std::unique_ptr<std::string> value_;
};

class DNA {
    friend class CompactDNA;
    friend class DNAunique;
public:
    using WhoIsBetterMergeStats = SeedComparisonHistory;

    DNA(string seq, string names) : sequence_(std::move(seq)),
        id_(std::move(names)), sequence_length_(sequence_.size()) {}
    DNA() { sequence_.reserve(151); }
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
        sequence_ = ""; sequence_length_ = 0; id_ = ""; resetNewId(); qual_.resize(0);
        good_quality_ = false; mid_quality_ = false; failed_ = true;
    }
    //~DNA(){}
    void appendSequence(const string& s) {
        sequence_ += s;
        sequence_length_ = sequence_.length();
        invalidateQualityCaches();
    }
    void appendQuality(const vector<qual_score>& q);
    void fixQ0(void);
    void setSequence(string& s) {
        sequence_ = s;
        sequence_length_ = sequence_.length();
        invalidateQualityCaches();
    }    void setSequence(string&& s) {
        sequence_ = std::move(s);
        sequence_length_ = sequence_.length();
        invalidateQualityCaches();
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

    const string& sampleBarcode() const { return sample_barcode_; }
    void setSampleBarcode(string barcode) { sample_barcode_ = std::move(barcode); }

    void setQual(vector<qual_score> Q) {
        qual_ = move(Q);
        invalidateQualityCaches();
    }
    //void setQual(vector<qual_score>&& Q) { qual_ = Q; avg_qual_ = -1.f; }

    void setAllQual(qual_score q) {
        for (size_t i = 0; i < qual_.size(); i++) { qual_[i] = q; }
        invalidateQualityCaches();
    }

    const string& getId() {
        if (id_fixed_) {
            return newId();
        }
        return id_;
    }

    shared_ptr<DNA> getDNAsubseq(int start, int end, string id);

    string getPositionFreeId(); // remove /1 /2 #1:0 etc
    void normalizePositionFreeId();
    const string& getOldId() { return id_; }
    string getShortId() { return id_.substr(0, getShorterHeadPos(id_)); }
    string getNewIDshort() { return newId().substr(0, getShorterHeadPos(newId())); }
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
    uint64_t getQsum() { return quality_sum_; }
    float qualWinfloat(uint, float, int&);


    float binomialFilter(int, float);
    //float qualWinfloat_hybr(int,float,int,float,int&);

    bool qualWinPos(unsigned int, float);
    bool qualAccumTrim(double d);
    int qualAccumulate(double d);

    double getAccumError();
    int minQual() { int mq = 50; for (uint i = 0; i < qual_.size(); i++) { if (qual_[i] < mq) { mq = qual_[i]; } }return mq; }
    void ntSpecQualScores(vector<uint64_t>&, vector<uint64_t>&);
    void ntSpecQualScoresMT(vector<uint64_t>&, vector<uint64_t>&);

    //returns qual filtered length 
    inline uint length() { return (uint)sequence_length_; }
    //returns original length 
    inline uint mem_length() { return (uint)sequence_.length(); }
    bool cutSeq(int start, int stop = -1, bool = false);
    bool cutSeqPseudo(int start) { return cutSeq(start, -1, true); }
    bool HomoNTRuns(int);
    int HomoNTTrim(int);
    struct PrimerMatchDetails {
        int start = -1;
        int mismatches = 0;
        bool fullLength = false;
        string observed;
    };
    PrimerMatchDetails getPrimerMatchDetails(const string& primer, int start) const;
    int matchSeq(string, int Err, int searchSpace, int startPos, bool exhaustive = false);
    // Read-only view overload: reuse primer matching while considering orientation.
    static int matchSeq(std::string_view sequence, std::string_view primer,
        int Err, int searchSpace, int startPos, bool exhaustive = false);
    void reverse_compliment(bool reset = true);
    int matchSeqRev(const string&, int Err, int searchSpace, int start = 0, bool = false);
    int matchSeq_tot(const string&, int, int, int&, bool endScan = false);
    struct BarcodeMatch {
        int start = -1, stop = -1, errors = 0;
        bool ambiguous = false; // Equally good matches imply different trim coordinates.
    };
    // ONT overload: best full-barcode edit match, with bounded start positions.
    static BarcodeMatch matchSeq_tot(std::string_view sequence, std::string_view barcode,
        int error, int maxPos);
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
    void setNewID(const string& x) { new_id_ = x; new_id_is_original_ = false; }
    void setNewID(string&& x) { new_id_ = std::move(x); new_id_is_original_ = false; }
    void setHeader(const string& x) { id_ = x; resetNewId(); }
    void setHeader(string&& x) {
        id_ = std::move(x);
        resetNewId();
    }
    void changeHeadPver(int ver);
    // Rewrites an existing mate tag ("/1" "/2" or Illumina " 1:"/" 2:") to this mate number,
    // e.g. after R1/R2 were swapped; headers without a tag are unchanged.
    void setMateNumber(int mate);
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
    void resetTruncation() {
        sequence_length_ = sequence_.length();
        invalidateQualityCaches();
    }
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
    };

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


public:


    void clearMergeSeed() {
        merge_seed_pos_ = -1;
        seed_length_ = -1;
        merge_offset_ = -1;
        reversed_merge_ = false;
    }

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
    // Marks this mate with an accepted consensus's merge length and overlap
    // errors. ReadMerger::merge leaves its inputs untouched, so a rejected
    // consensus never makes a read look merged.
    void adoptMergeStats(const DNA& merged) {
        FtsDetected.errInOverlap = merged.FtsDetected.errInOverlap;
        FtsDetected.meanQInOverlapMismatch = merged.FtsDetected.meanQInOverlapMismatch;
        FtsDetected.mergeLength = merged.FtsDetected.mergeLength;
    }
    int getMergeErrors() { return FtsDetected.errInOverlap; }
    qual_score getMergeErrorsQual() { return FtsDetected.meanQInOverlapMismatch; }
    bool isReversed() { return reversed_; }

protected:

    //const char* DN;

    void invalidateQualityCaches() {
        avg_qual_ = -1.f;
        quality_sum_ = 0;
        accumulated_error_ = 0.;
    }




    struct ElementsDetection {
        uint8_t forward : 1;
        uint8_t reverse : 1;//primers detected
        uint8_t TA_cut : 1;
        uint8_t barcode_detected : 1;
        uint8_t barcode_cut : 1;
        uint8_t dereplicated : 1;
        uint8_t revPairConstellation : 1;
        //read merging related stats
        qual_score meanQInOverlapMismatch;
        int errInOverlap;
        int mergeLength;

        ElementsDetection() : forward(false), reverse(false), TA_cut(false), barcode_detected(false),
            barcode_cut(false), dereplicated(false), revPairConstellation(false),
            meanQInOverlapMismatch(0), errInOverlap(-1), mergeLength(-1) {}
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

    };



    // Large aligned fields first; byte-sized state shares the final word.
    CompactOptionalString sample_barcode_;
    string sequence_;
    string id_, new_id_;
    vector<qual_score> qual_;
    size_t sequence_length_ = 0;
    uint64_t quality_sum_ = 0;
    double accumulated_error_ = 0.;
    ElementsDetection FtsDetected;
public:
    QualStats QualCtrl;
    int merge_seed_pos_ = -1;
    int seed_length_ = -1;
    int merge_offset_ = -1;
protected:
    int sample_id_ = -1;
    float avg_qual_ = -1.f;
    float tempFloat = 0.f;
    short read_position_ = -1; // -1 unknown; 0 forward; 1 reverse; 3 MID
    qual_score fastq_write_offset_ = 0;
    uint8_t good_quality_ : 1 = false;
    uint8_t mid_quality_ : 1 = false;
    uint8_t failed_ : 1 = false;
    uint8_t reversed_ : 1 = false;
    bool id_fixed_ = false;
    bool new_id_is_original_ = true;
public:
    bool reversed_merge_ = false;
protected:
    // Parsing constructors adopt a buffer and must not allocate a reserve first.
    struct AdoptBuffer {};
    explicit DNA(AdoptBuffer) {}
    const string& newId() const { return new_id_is_original_ ? id_ : new_id_; }
    string& mutableNewId() {
        if (new_id_is_original_) { new_id_ = id_; new_id_is_original_ = false; }
        return new_id_;
    }
    void resetNewId() { new_id_.clear(); new_id_is_original_ = true; }

};

using sample_count_t = uint64_t;
using read_occ_entry = diffDNA::SampleCount;
using read_occ = diffDNA::SampleCounts;

struct DNAHasher
{
    size_t operator()(shared_ptr<DNA> k) const
    {
        return ((hash<std::string_view>()(k->getSeqPseudoView())) >> 1);
    }
};



// Shared arithmetic for selected representatives and detached merge history.
// Uniform coverage avoids an allocation and a per-base increment stream.
class DNAQualityEvidence {
    friend class CompactDNA;
public:
    using qual_accum_t = uint64_t;
    std::vector<qual_accum_t> quality_sum_per_base_;
    static void addToQualSum(qual_accum_t& target, qual_accum_t value);
    void initializeQualityEvidence(const vector<qual_score>& qualities, size_t length);
    void addQualityEvidence(const vector<qual_score>& qualities, size_t length);
    void addQualityEvidence(const DNAQualityEvidence& other);
    void clearQualityEvidence();
    qual_accum_t qualityCoverage(size_t position) const {
        return quality_coverage_per_base_.empty() ? uniform_quality_coverage_ : quality_coverage_per_base_[position];
    }
protected:
    void buildAverageQualityString(string& output, size_t length, int offset) const;
private:
    void extendQualityEvidence(size_t length, bool uniformSource);
    std::vector<qual_accum_t> quality_coverage_per_base_;
    qual_accum_t uniform_quality_coverage_ = 0;
};

// A stored mate (R2) or merged consensus of a dereplicated unique: the read, its
// quality evidence and the derep FASTQ writers. Sample counts, mate/merge links and
// selection history exist only on the representative (DNAunique), saving memory per pair.
class DerepRead : public DNA, public DNAQualityEvidence {
    friend class DNA;
    friend class CompactDNA;
    friend class DNAunique;
public:
    using qual_accum_t = uint64_t;
    DerepRead() : DNA(), derep_fastq_offset_(0) {}
    DerepRead(string s, string x) : DNA(std::move(s), std::move(x)), derep_fastq_offset_(0) {}
    explicit DerepRead(const shared_ptr<DNA>& d) : DNA(*d), derep_fastq_offset_(0) {}

    void prepSumQuals();
    void sumQualities(const shared_ptr<DerepRead>& dna);
    void sumQualities(const shared_ptr<DNA>& dna);
    void releaseRawQualitiesIfPromoted();

    void prepareDerepQualities(int ofastQver);
    void writeDerepFastQ(ostream&, bool = true);
    void writeFastQ(ostream&, bool = true);
    string writeFastQ(bool = true);
    void writeFastQ(string&, bool = true);

    void saveMem();
    shared_ptr<DerepRead> cloneForOutput() const { return std::make_shared<DerepRead>(*this); }
    void compactForFastaDerepStorage();
    std::vector<qual_accum_t> transferPerBaseQualitySum();

    // Stored R2/merged objects omit headers that can be rebuilt from R1's header
    // (R2: trailing "/1" -> "/2", or identical; merged: identical). Only exact
    // matches are dropped, so restoring is lossless.
    static string mateHeaderFrom(const string& r1Header);
    const string& storedHeader() const { return newId(); }
    void restoreMateHeader(const string& r1Header);

protected:
    void promoteQualitiesToSums(bool releaseRawQualities);
    void clearStoredHeader() {
        string().swap(id_);
        string().swap(new_id_);
        new_id_is_original_ = true;
    }
    // to calculate mean
    qual_score derep_fastq_offset_;
};

class DNAunique : public DerepRead {//used for dereplication
    friend class DNA;
    friend class CompactDNA;
public:
    // Constructors and Destructors
    DNAunique() : DerepRead(), pair_(0) {}
    DNAunique(string s, string x) : DerepRead(std::move(s), std::move(x)) {}

    // Mostly used constructor
    DNAunique(shared_ptr<DNA>d, int BC) : DerepRead(d), pair_(0) {//best_seed_length_((uint)sequence_.size())
        incrementSampleCounter(BC);
    }
    ~DNAunique() = default;

    void recordWhoIsBetterCompare(int dMergeLength, int rMergeLength);
    const WhoIsBetterMergeStats& getWhoIsBetterMergeStats() const { return who_is_better_merge_stats_; }
    void resetWhoIsBetterMergeStats() { who_is_better_merge_stats_ = WhoIsBetterMergeStats(); }

    //string sequence_;    string id_;
    void Count2Head(bool);
    bool betterPreSeed(shared_ptr<DNA> d1, shared_ptr<DNA> d2);
    void matchedDNA(shared_ptr<DNA>, shared_ptr<DNA>, shared_ptr<DNA>, int, bool, bool collectMergedQualities = true, bool markDereplicated = true);
    // Coarse clusters retain only the chosen read's qualities and attachments.
    void matchedCoarseDNA(std::shared_ptr<DNA> dna, std::shared_ptr<DNA> mate,
        std::shared_ptr<DNA> merged, int sample) {
        matchedDNA(std::move(dna), std::move(mate), std::move(merged), sample, true, false);
        merge_quality_evidence_.reset();
    }
    void incrementSampleCounter(int sample_id);
    void writeMap(ostream& os, const string&, vector<uint64_t>&, const vector<int>&, vector<uint64_t>&, vector<int>&);
    //inline int getCount() { return count_; }
    uint64_t totalSum();
    //uint getBestSeedLength() { return best_seed_length_; }
    //void setBestSeedLength(uint i) { best_seed_length_ = i; }
    void incrementSampleCounterBy(int sample_id, uint64_t count);
    void transferOccurence(shared_ptr<DNAunique>);
    const read_occ& getDerepMap() const { return differences_.sampleCounts(); }
    vector<uint64_t> getDerepMapSort(size_t);
    //vector<pair_<int, int>> getDerepMapSort2(size_t wh);
    //void getDerepMapSort(vector<int>&, vector<int>&);

    using DerepRead::sumQualities;
    // Exact overload so DNAunique arguments do not make the base overloads ambiguous.
    void sumQualities(const shared_ptr<DNAunique>& dna) { DerepRead::sumQualities(std::static_pointer_cast<DerepRead>(dna)); }

    // Replaces the complete count/variant state, as used to create an output snapshot.
    void attachDiffDNA(std::shared_ptr<diffDNA> differences) { differences_ = differences ? *differences : diffDNA{}; }
    const diffDNA& getDerepCounts() const { return differences_; }
    // Preserve the optional, owning snapshot API for nucleotide observations.
    std::shared_ptr<const diffDNA> getDiffDNA() const {
        return differences_.hasSequenceData() ? std::make_shared<const diffDNA>(differences_) : nullptr;
    }
    void recordDiffDNA(std::string_view sequence, int sample = -1, uint64_t copies = 1);
    void recordDiffDNA(std::string_view sequence, const diffDNA::SampleCounts& samples);
    void recordDiffQuality(DNA& dna, int sample) { differences_.addQuality(dna, 1, sample); }
    void recordDiffPairQuality(DNA& r1, DNA& r2, int sample) { differences_.addPairQuality(r1, r2, 1, sample); }
    void recordDiffPair(std::string_view r1, std::string_view r2, const diffDNA::SampleCounts& samples, const diffDNA::QualityRecord& quality) {
        differences_.addPair(r1, r2, samples, quality);
    }
    void recordDiffDNA(std::string_view sequence, const diffDNA::SampleCounts& samples, const diffDNA::QualityRecord& quality) {
        differences_.add(sequence, samples, quality);
    }
    void clearDerepCounts() { differences_.clear(); }
    // Recovered deferred reads: counts only, no sequence or quality evidence.
    void addSampleCount(int sample, uint64_t copies) { differences_.addSample(sample, copies); }
    shared_ptr<DNAunique> cloneForOutput() const;
    void sumMergeQualities(const shared_ptr<DNAunique>& other);
    void absorbDerepEvidence(const shared_ptr<DNAunique>& other, bool collectQualities);

    void attachPair(shared_ptr<DerepRead> dna_unique);
    shared_ptr<DerepRead> getPair(void);
    void dropDerivableMateHeaders();

    void attachMerge(shared_ptr<DerepRead> dnamerge);
    shared_ptr<DerepRead> getMerge(void);


    //estimates if one sample occurence covers the unique counts required for sample specific derep min counts
    bool pass_deprep_smplSpc(const vector<int>&);
    //int counts() const { return count_; }

    void takeOver(shared_ptr<DNAunique> dna_unique_old);
    void takeOverDNA(shared_ptr<DNA> dna1, shared_ptr<DNA> dna2, shared_ptr<DNA> dnaMerge, bool keepQualities = true);


private:
    //int count_;
    //matrixUnit chimeraCnt;
    //int best_seed_length_;

    diffDNA differences_;
    struct MergeQualityEvidence : DNAQualityEvidence {
        // Preserve the unpromoted case without allocating eight-byte sums yet.
        vector<qual_score> raw;
        MergeQualityEvidence() = default;
        explicit MergeQualityEvidence(DNA& dna) {
            const auto& qualities = dna.getQual();
            const size_t n = std::min(static_cast<size_t>(dna.length()), qualities.size());
            raw.assign(qualities.begin(), qualities.begin() + n);
        }
        void add(DNA& dna) {
            if (dna.getQual().empty() || !dna.length()) { return; }
            if (quality_sum_per_base_.empty()) { initializeQualityEvidence(raw, raw.size()); }
            addQualityEvidence(dna.getQual(), dna.length());
            vector<qual_score>().swap(raw);
        }
    };

    shared_ptr<DerepRead> pair_;

    shared_ptr<DerepRead> merge_;
    // Retain accepted consensus evidence if the chosen representative loses its merge.
    shared_ptr<MergeQualityEvidence> merge_quality_evidence_;

    WhoIsBetterMergeStats who_is_better_merge_stats_;

    //threadsafe
};




// Merge-offset subgroups of one dereplication key. Nearly every key has a single
// subgroup, so a sorted vector (same ascending order as std::map) replaces a tree
// node allocation per key. Callers must not keep iterators across operator[].
class DerepSubgroupMap {
public:
    using value_type = std::pair<int, shared_ptr<DNAunique>>;
    using iterator = std::vector<value_type>::iterator;
    using const_iterator = std::vector<value_type>::const_iterator;
    iterator begin() { return entries_.begin(); }
    iterator end() { return entries_.end(); }
    const_iterator begin() const { return entries_.begin(); }
    const_iterator end() const { return entries_.end(); }
    size_t size() const { return entries_.size(); }
    iterator find(int key) {
        auto it = lowerBound(key);
        return it != entries_.end() && it->first == key ? it : entries_.end();
    }
    shared_ptr<DNAunique>& operator[](int key) {
        auto it = lowerBound(key);
        if (it == entries_.end() || it->first != key) { it = entries_.emplace(it, key, nullptr); }
        return it->second;
    }
private:
    iterator lowerBound(int key) {
        return std::lower_bound(entries_.begin(), entries_.end(), key,
            [](const value_type& entry, int k) { return entry.first < k; });
    }
    std::vector<value_type> entries_;
};

class DNAuniqSet {
public:
    DNAuniqSet();
    explicit DNAuniqSet(std::shared_ptr<const CompactDNAReferences> references);
    ~DNAuniqSet();
    bool isCompact() const { return bool(compact_); }
    std::vector<diffDNA> compactVariants() const;

    // One ordinary observation, retaining merge-offset selection and quality evidence.
    void addObservation(shared_ptr<DNA> dna, shared_ptr<DNA> mate, shared_ptr<DNA> merged,
        int mergePosition, int sample, bool asFasta, bool collectMergedQualities, bool markDereplicated = true);
    void addNewDNAuniq(shared_ptr<DNA> dna, shared_ptr<DNA> dna2,
        shared_ptr<DNA> dnaM, int MrgPos1, int sample_id, bool b_derep_as_fasta_, bool markDereplicated = true);
    shared_ptr<DNAunique>& operator[] (int x) {
        return DNUs[x];
    }
    DerepSubgroupMap::iterator find(const int x);
    DerepSubgroupMap::iterator end();
    DerepSubgroupMap::iterator begin();
    size_t size() const;
    uint64_t totalSum() const;
    void addEvidenceTo(const shared_ptr<DNAunique>& target, bool collectQualities) const;
    void setBest(bool addCnts, bool collectQualities = true);
    shared_ptr<DNAunique> best(bool addCnts, bool collectQualities = true);

private:
    // (An unused per-key mutex was removed: shards already serialize updates.)
    std::unique_ptr<CompactDerepSet> compact_;
    DerepSubgroupMap DNUs;
    shared_ptr<DNAunique> bestDNU;

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



