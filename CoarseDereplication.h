#pragma once
// SPDX-License-Identifier: GPL-3.0-or-later
#include "DNA.h"
#include "diffDNA.h"
#include <map>

// Streaming, same-start substitution clustering. The caller serializes access.
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
        DNAuniqSet representatives;
        std::string_view key() const { return std::string_view(*reference).substr(0, searchLength); }
        std::string_view mateKey() const { return mateReference ? std::string_view(*mateReference).substr(0, mateSearchLength) : std::string_view{}; }
    };
    struct Stats { uint64_t queries = 0, comparisons = 0, fallbackQueries = 0, indexedSeeds = 0; };
    CoarseDereplication(uint64_t identity, unsigned homopolymer, bool storeDiffs, bool reassign = false, bool storeQualities = false);
    struct Refinement {
        std::vector<std::shared_ptr<DNAunique>> representatives;
        uint64_t movedReads = 0, movedVariants = 0;
        Stats searches;
    };
    // Independent output snapshot: live admission anchors and counts stay valid.
    Refinement refine(bool searchWithMerge, int searchLimit) const;
    bool reassigns() const { return reassign_; }
    // Returns the cluster ID, or npos for a rejected read without an existing match.
    size_t add(std::string_view search, std::string_view fullSequence,
        std::shared_ptr<DNA> dna, std::shared_ptr<DNA> mate, std::shared_ptr<DNA> merged,
        int sample, bool passed, int searchLimit = -1);
    size_t find(std::string_view sequence, std::string_view mate = {});
    size_t mismatchBudget(size_t length) const;
    void clear();
    const std::vector<std::unique_ptr<Group>>& groups() const { return groups_; }
    const Stats& stats() const { return stats_; }
    bool storesDiffs() const { return storeDiffs_; }
    bool storesQualities() const { return storeQualities_; }
    static constexpr size_t npos = std::numeric_limits<size_t>::max();
private:
    struct SeedKey {
        uint64_t word; size_t position;
        bool operator==(const SeedKey& other) const { return word == other.word && position == other.position; }
    };
    struct SeedHash {
        size_t operator()(const SeedKey& key) const {
            return robin_hood::hash<uint64_t>{}(key.word) ^ robin_hood::hash<size_t>{}(key.position);
        }
    };
    // Length bins bound the comparison length without needing one table for
    // every ONT read length. Bin zero contains lengths 1..127.
    static constexpr size_t lengthBin = 128;
    struct Bucket {
        std::vector<size_t> ids;
        robin_hood::unordered_map<SeedKey, std::vector<size_t>, SeedHash> postings;
    };
    uint64_t identity_;
    unsigned homopolymer_;
    bool storeDiffs_, reassign_, storeQualities_;
    size_t search(std::string_view sequence, bool nearest, size_t original, std::string_view mate = {});
    void collectSeeds(std::string_view sequence);
    void indexGroup(size_t id);
    Stats stats_;
    std::vector<std::unique_ptr<Group>> groups_;
    std::map<size_t, Bucket> buckets_;
    std::vector<SeedKey> querySeeds_;
    std::vector<const std::vector<size_t>*> lists_;
    struct Cursor { const std::vector<size_t>* ids; size_t offset; };
    std::vector<Cursor> heap_;
};
