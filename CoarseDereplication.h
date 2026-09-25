#pragma once
// SPDX-License-Identifier: GPL-3.0-or-later
#include "DNA.h"
#include "diffDNA.h"
#include <map>

// Streaming, same-start substitution clustering. The caller serializes add();
// probe() is read-only and may run concurrently with other probes.
class CoarseDereplication {
public:
    static constexpr uint64_t identityScale = 100000000; // 100 percent
    static constexpr size_t seedLength = 20;
    static uint64_t parseIdentity(const std::string& percent);
    static unsigned parseHomopolymer(const std::string& value);
    static bool seedWord(std::string_view sequence, size_t position, unsigned homopolymer, uint64_t& word);
    struct Group {
        std::shared_ptr<const std::string> reference, mateReference;
        size_t searchLength = 0, mateSearchLength = 0;
        std::shared_ptr<const CompactDNAReferences> compactReferences;
        DNAuniqSet representatives;
        std::string_view key() const { return std::string_view(*reference).substr(0, searchLength); }
        std::string_view mateKey() const { return mateReference ? std::string_view(*mateReference).substr(0, mateSearchLength) : std::string_view{}; }
    };
    struct Stats {
        uint64_t queries = 0, comparisons = 0, fallbackQueries = 0, indexedSeeds = 0, exactHits = 0, rejectedExactMisses = 0;
        void add(const Stats& other) {
            queries += other.queries; comparisons += other.comparisons; fallbackQueries += other.fallbackQueries;
            indexedSeeds += other.indexedSeeds; exactHits += other.exactHits; rejectedExactMisses += other.rejectedExactMisses;
        }
    };
    struct SeedKey {
        uint64_t word; size_t position;
        bool operator==(const SeedKey& other) const { return word == other.word && position == other.position; }
    };
    struct Cursor { const std::vector<size_t>* ids; size_t offset; };
    // Per-caller search buffers, so that probes can search concurrently.
    struct SearchScratch {
        std::vector<SeedKey> seeds;
        std::vector<const std::vector<size_t>*> lists;
        std::vector<Cursor> heap;
    };
    // Result of a read-only candidate search against a snapshot of the index.
    // Anchor IDs grow with creation order and search returns the smallest
    // matching ID, so a snapshot match stays the first-created match; without
    // one only anchors created after the snapshot (IDs >= snapshot) can match.
    struct Probe {
        size_t group = std::numeric_limits<size_t>::max();
        size_t snapshot = 0;
        bool searched = false;
        Stats stats;
    };
    // Representative/count/variant update for one admitted observation. It
    // only touches its target set, so updates of different targets may run in
    // parallel (the caller serializes updates of one target).
    struct Update {
        enum Kind : uint8_t { Observation, Matched, MatchedCoarse, RecordOnly };
        DNAuniqSet* target = nullptr;
        Kind kind = RecordOnly;
        bool paired = false, derepAsFasta = true, collectMergedQualities = false;
        int sample = 0, mergePosition = 0;
        std::shared_ptr<DNA> dna, mate, merged;
        std::string_view fullSequence;
    };
    CoarseDereplication(uint64_t identity, unsigned homopolymer, bool storeDiffs, bool reassign = false, bool storeQualities = false, bool exactOutput = false);
    struct ExactAssignment {
        size_t group;
        // Only the requested output level owns counts, variants, and seeds.
        std::unique_ptr<DNAuniqSet> representatives;
    };
    using ExactAssignments = robin_hood::unordered_map<std::string, ExactAssignment, robin_hood::hash<std::string_view>, std::equal_to<>>;
    size_t storagePartitionHash(std::string_view sequence) const;
    bool outputsExact() const { return exactOutput_; }
    bool hasExactKey(std::string_view key) const {
        return exactAssignments_.find(key, robin_hood::is_transparent_tag{}) != exactAssignments_.end();
    }
    bool storageOnly() const { return exactOutput_ && !storeQualities_ && !storeDiffs_ && !reassign_; }
    const ExactAssignments& exactGroups() const { return exactAssignments_; }
    struct Refinement {
        std::vector<std::shared_ptr<DNAunique>> representatives;
        uint64_t movedReads = 0, movedVariants = 0;
        Stats searches;
    };
    // Independent output snapshot: live admission anchors and counts stay valid.
    Refinement refine(bool searchWithMerge, int searchLimit) const;
    bool reassigns() const { return reassign_; }
    // Returns the cluster ID, or npos for a rejected read without an existing match.
    // With a probe, the candidate search it performed is reused (see Probe).
    // With `deferred`, representative updates are returned instead of applied;
    // new anchors and exact keys are still created here.
    size_t add(std::string_view search, std::string_view fullSequence,
        std::shared_ptr<DNA> dna, std::shared_ptr<DNA> mate, std::shared_ptr<DNA> merged,
        int sample, bool passed, int searchLimit = -1,
        bool derepAsFasta = true, bool collectMergedQualities = false, int mergePosition = 0,
        const Probe* probe = nullptr, std::vector<Update>* deferred = nullptr);
    // The read-only part of add(): exact-key lookup and candidate search.
    Probe probe(std::string_view search, const std::shared_ptr<DNA>& dna, const std::shared_ptr<DNA>& mate,
        bool passed, int searchLimit, SearchScratch& scratch) const;
    void apply(const Update& update) const;
    size_t find(std::string_view sequence, std::string_view mate = {});
    size_t mismatchBudget(size_t length) const;
    void clear();
    const std::vector<std::unique_ptr<Group>>& groups() const { return groups_; }
    const Stats& stats() const { return stats_; }
    bool storesDiffs() const { return storeDiffs_; }
    bool storesQualities() const { return storeQualities_; }
    static constexpr size_t npos = std::numeric_limits<size_t>::max();
private:
    struct SeedHash {
        size_t operator()(const SeedKey& key) const {
            return robin_hood::hash<uint64_t>{}(key.word) ^ robin_hood::hash<size_t>{}(key.position);
        }
    };
    // Length bins bound the comparison length without needing one table for
    // every ONT read length. Bin zero contains lengths 1..127.
    static constexpr size_t lengthBin = 128;
    struct Bucket {
        size_t minLength = std::numeric_limits<size_t>::max(), maxLength = 0;
        std::vector<size_t> ids;
        robin_hood::unordered_map<SeedKey, std::vector<size_t>, SeedHash> postings;
    };
    uint64_t identity_;
    unsigned homopolymer_;
    bool storeDiffs_, reassign_, storeQualities_, exactOutput_;
    size_t search(std::string_view sequence, bool nearest, size_t original, std::string_view mate = {}, bool paired = false);
    // Const search core; only anchors with IDs >= minId are considered.
    size_t search(std::string_view sequence, bool nearest, size_t original, std::string_view mate, bool paired,
        SearchScratch& scratch, Stats& stats, size_t minId, bool countQuery) const;
    void collectSeeds(std::string_view sequence, std::vector<SeedKey>& seeds) const;
    void collectSeeds(std::string_view sequence) { collectSeeds(sequence, scratch_.seeds); }
    void indexGroup(size_t id);
    // Mate key shared by add() and probe().
    std::string_view mateKeyFor(const std::shared_ptr<DNA>& mate, bool paired, int searchLimit) const;
    ExactAssignment* exactLookup(std::string_view sequence, bool paired);
    // Authoritative ordinary dereplication by effective R1/search sequence.
    // Entries survive until clear(): an R2 mismatch must never split an R1
    // dereplicate merely because its lookup entry was evicted.
    size_t find(std::string_view sequence, std::string_view mate, ExactAssignment*& exactMatch, bool paired, bool searchCandidates = true);
    ExactAssignments exactAssignments_;
    Stats stats_;
    std::vector<std::unique_ptr<Group>> groups_;
    std::map<size_t, Bucket> buckets_;
    SearchScratch scratch_;
};
