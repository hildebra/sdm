#pragma once
// Lossless internal storage for ordinary exact dereplicates.
// SPDX-License-Identifier: GPL-3.0-or-later
#include "DNA.h"
#include <optional>

struct CompactDNAReferences {
    std::shared_ptr<const std::string> first, mate, merged;
    CompactDNAReferences(DNA& r1, DNA* r2, DNA* consensus);
};

// No DNA object survives an update. Nucleotides share coarse references; the
// selected qualities and ordinary quality evidence use lossless packed arrays.
class CompactDNA {
public:
    CompactDNA(const CompactDNAReferences& references);
    void observe(const std::shared_ptr<DNA>& dna, const std::shared_ptr<DNA>& mate,
        const std::shared_ptr<DNA>& merged, int sample, bool asFasta,
        bool collectMergedQualities, bool markDereplicated, const CompactDNAReferences& references);
    std::shared_ptr<DNAunique> materialize(const CompactDNAReferences& references) const;
    uint64_t totalSum() const { return variants_.readCount(); }
    const diffDNA& variants() const { return variants_; }
private:
    struct Writer;
    struct Reader;
    static void packEvidence(Writer&, const DNAQualityEvidence&);
    static void unpackEvidence(Reader&, DNAQualityEvidence&);
    static void packRead(Writer&, const DerepRead&, const CompactDNAReferences&, int level);
    static void packUnique(Writer&, const DNAunique&, const CompactDNAReferences&);
    static void unpackRead(Reader&, DerepRead&, const CompactDNAReferences&, int level);
    static std::shared_ptr<DNAunique> unpackUnique(Reader&, const CompactDNAReferences&);
    void store(DNAunique&, const CompactDNAReferences&);
    std::string selected_;
    diffDNA variants_; // Full unmerged nucleotide variants and occurrences, no per-variant qualities.
    std::optional<SeedReadStats> firstStats_, mateStats_, mergedStats_;
    SeedComparisonHistory history_;
};

class CompactDerepSet {
public:
    explicit CompactDerepSet(std::shared_ptr<const CompactDNAReferences> references) : references_(std::move(references)) {}
    void addObservation(std::shared_ptr<DNA> dna, std::shared_ptr<DNA> mate, std::shared_ptr<DNA> merged,
        int mergePosition, int sample, bool asFasta, bool collectMergedQualities, bool markDereplicated);
    size_t size() const { return subgroups_.size(); }
    uint64_t totalSum() const;
    std::shared_ptr<DNAunique> best(bool addCounts, bool collectQualities) const;
    void addEvidenceTo(const std::shared_ptr<DNAunique>& target, bool collectQualities) const;
    std::vector<diffDNA> variants() const;
private:
    std::shared_ptr<const CompactDNAReferences> references_;
    std::map<int, CompactDNA> subgroups_;
};
