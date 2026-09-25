#pragma once
// Unified sparse sample counts and optional nucleotide variants.
// SPDX-License-Identifier: GPL-3.0-or-later
#include <cstdint>
#include <functional>
#include <iosfwd>
#include <memory>
#include <span>
#include <string>
#include <string_view>
#include <vector>
#include "include/robin_hood.h"
#include "SeedStats.h"

class diffDNA {
public:
    using SampleCount = std::pair<int, uint64_t>;
    using SampleCounts = std::vector<SampleCount>;
    // Counts-only mode: no reference, variant table, or auxiliary allocation.
    diffDNA() = default;
    explicit diffDNA(std::shared_ptr<const std::string> reference, bool trackSamples = false, bool storeQualities = false, std::shared_ptr<const std::string> mateReference = nullptr, bool allowPartialPairs = false);
    explicit diffDNA(std::string reference, bool trackSamples = false, bool storeQualities = false);
    struct QualityRead {
        SeedReadStats stats;
        std::vector<qual_score> qualities; // Full backing sequence; stats.logical retains the filtering length.
    };
    struct QualitySeed : QualityRead {
        std::shared_ptr<const QualityRead> mate; // Both reads belong to one selected fragment.
    };
    struct QualityRecord {
        std::shared_ptr<const QualitySeed> seed; // Immutable payload shared by output snapshots.
        SeedComparisonHistory history;          // Observation history survives seed replacement.
    };
    // The same lossless nucleotide codec is used by compact ordinary seeds.
    static std::string encodeDifference(std::string_view reference, std::string_view sequence);
    static std::string decodeDifference(const std::string& reference, std::string_view patch);
    bool storesPairs() const;
    const std::string& mateReference() const;
    void addPair(std::string_view r1, std::string_view r2, uint64_t copies = 1, int sample = -1);
    void addPair(std::string_view r1, std::string_view r2, const SampleCounts& samples, const QualityRecord& quality);
    void addPairQuality(DNA& r1, DNA& r2, uint64_t copies = 1, int sample = -1);
    // Views/references are borrowed for this callback only. Quality is null after binary reload.
    void forEachPair(const std::function<void(std::string_view, std::string_view, uint64_t, const SampleCounts*, const QualityRecord*)>& visitor, bool ordered = false) const;
    void writeFastq(std::ostream& r1, std::ostream& r2, std::string_view cluster, int qualityOffset = 33) const;
    // Every exact variant, decoded, in the order of the ordered visits and
    // exports. Its pointers borrow from this object and stay valid while it
    // is unchanged, so ranges of one snapshot can be written concurrently.
    struct Fragment {
        std::string first, second; // second: mate, empty for single-end
        uint64_t copies;
        const SampleCounts* samples;
        const QualityRecord* quality;
    };
    struct SortedFragments {
        std::vector<Fragment> fragments;     // storage order
        std::vector<const Fragment*> order;  // sorted
        size_t size() const { return order.size(); }
        const Fragment& operator[](size_t i) const { return *order[i]; }
    };
    // Runs task(0 .. count-1), possibly concurrently, and returns once all have finished.
    using ParallelFor = std::function<void(size_t count, const std::function<void(size_t)>& task)>;
    // With `parallel`, decoding and sorting are split into tasks; the order is the same.
    SortedFragments sortedFragments(const ParallelFor* parallel = nullptr) const;
    // Records [begin, end) of writeFastq's output, from a sortedFragments() snapshot.
    void writeFastq(std::ostream& r1, std::ostream& r2, std::string_view cluster, const SortedFragments& fragments,
        size_t begin, size_t end, int qualityOffset = 33) const;
    void writeFastq(std::ostream& out, std::string_view cluster, const SortedFragments& fragments,
        size_t begin, size_t end, int qualityOffset = 33) const;
    void writeFasta(std::ostream& r1, std::ostream& r2, std::string_view cluster) const;
    void addSample(int sample, uint64_t copies = 1);
    // Retain one real selected quality vector per full processed single-end variant.
    void addQuality(DNA& dna, uint64_t copies = 1, int sample = -1);
    void add(std::string_view sequence, const SampleCounts& samples, const QualityRecord& quality);
    bool storesQualities() const;
    void forEachWithQualities(const std::function<void(std::string_view, uint64_t, const SampleCounts*, const QualityRecord&)>& visitor, bool ordered = false) const;
    void writeFastq(std::ostream& out, std::string_view cluster, int qualityOffset = 33) const;
    void add(std::string_view sequence, uint64_t copies = 1, int sample = -1);
    void add(std::string_view sequence, const SampleCounts& samples);
    // A counts-only destination keeps only counts. A variant destination also
    // rebases and merges the source variants; missing provenance is rejected.
    void merge(const diffDNA& other);
    void clear();
    const SampleCounts& sampleCounts() const;
    bool hasSequenceData() const { return sequences_ != nullptr; }
    bool tracksSamples() const;
    void forEachWithSamples(const std::function<void(std::string_view, uint64_t, const SampleCounts&)>& visitor) const;
    const std::string& reference() const;
    uint64_t readCount() const { return reads_; }
    uint64_t referenceCount() const;
    size_t uniqueCount() const;
    size_t payloadBytes() const;
    // Visit each exact nucleotide subcluster once, with its full abundance.
    void forEach(const std::function<void(std::string_view, uint64_t)>& visitor, bool ordered = false) const;
    void writeFasta(std::ostream& out, std::string_view cluster) const;
    void write(std::ostream& out, std::string_view cluster) const;
    // Reads single-end versions 1/2 or paired version 3. No qualities/sample IDs are serialized.
    // Returns nullptr only at clean EOF between blocks.
    static std::shared_ptr<diffDNA> read(std::istream& in, std::string& cluster);
private:
    struct SequenceData;
    void detachSequences();
    void record(std::string_view sequence, uint64_t copies, std::span<const SampleCount> samples, DNA* read = nullptr, const QualityRecord* quality = nullptr, std::string_view mateSequence = {}, DNA* mate = nullptr);
    void visitPatches(const std::function<void(const std::string&, uint64_t, const SampleCounts*, const QualityRecord*)>& visitor) const;
    void visit(const std::function<void(std::string_view, uint64_t, const SampleCounts*, const QualityRecord*)>& visitor, bool ordered = false) const;
    void visitFragment(const std::function<void(std::string_view, std::string_view, uint64_t, const SampleCounts*, const QualityRecord*)>& visitor, bool ordered = false) const;
    // In reassignment mode, occurrences_ is a lazily rebuilt aggregate view;
    // the per-variant sparse counts are the underlying source of truth.
    mutable SampleCounts occurrences_;
    uint64_t reads_ = 0;
    std::shared_ptr<SequenceData> sequences_; // Optional, shared copy on write.
};
