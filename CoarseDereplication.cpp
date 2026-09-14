// SPDX-License-Identifier: GPL-3.0-or-later
#include "CoarseDereplication.h"
#include "ReadMerger.h"
#include <algorithm>
#include <climits>
#include <stdexcept>

uint64_t CoarseDereplication::parseIdentity(const std::string& value) {
    uint64_t integer = 0, fraction = 0, scale = 1;
    bool decimal = false, digit = false, fractionDigit = false;
    for (char c : value) {
        if (c == '.' && !decimal && digit) { decimal = true; continue; }
        if (c < '0' || c > '9') { throw std::invalid_argument("-derepIdentity must be a percentage in (0,100], with at most six decimal places"); }
        digit = true;
        if (decimal) {
            if (scale == 1000000) { throw std::invalid_argument("-derepIdentity accepts at most six decimal places"); }
            scale *= 10; fraction = fraction * 10 + c - '0'; fractionDigit = true;
        } else {
            integer = integer * 10 + c - '0';
            if (integer > 100) { throw std::invalid_argument("-derepIdentity must be in (0,100]"); }
        }
    }
    const uint64_t result = integer * 1000000 + fraction * (1000000 / scale);
    if (!digit || (decimal && !fractionDigit) || result == 0 || result > identityScale) {
        throw std::invalid_argument("-derepIdentity must be a percentage in (0,100]");
    }
    return result;
}
unsigned CoarseDereplication::parseHomopolymer(const std::string& value) {
    if (value.empty() || value.size() > 2 || value.find_first_not_of("0123456789") != std::string::npos) {
        throw std::invalid_argument("-derepSeedHomopolymer must be an integer from 2 to 20");
    }
    const unsigned result = static_cast<unsigned>(std::stoul(value));
    if (result < 2 || result > seedLength) { throw std::invalid_argument("-derepSeedHomopolymer must be from 2 to 20"); }
    return result;
}
CoarseDereplication::CoarseDereplication(uint64_t identity, unsigned homopolymer, bool storeDiffs, bool reassign, bool storeQualities)
    : identity_(identity), homopolymer_(homopolymer), storeDiffs_(storeDiffs), reassign_(reassign), storeQualities_(storeQualities) {
    if (!identity || identity > identityScale || homopolymer < 2 || homopolymer > seedLength) {
        throw std::invalid_argument("Invalid coarse dereplication configuration");
    }
}
size_t CoarseDereplication::mismatchBudget(size_t length) const {
    // Split the product to keep the inclusive decimal boundary exact without overflow.
    const uint64_t errors = identityScale - identity_;
    return (length / identityScale) * errors + ((length % identityScale) * errors) / identityScale;
}
bool CoarseDereplication::seedWord(std::string_view sequence, size_t position, unsigned homopolymer, uint64_t& word) {
    if (position > sequence.size() || sequence.size() - position < seedLength) { return false; }
    word = 0;
    unsigned run = 0; char previous = 0;
    for (size_t i = position; i < position + seedLength; ++i) {
        const char base = sequence[i];
        run = base == previous ? run + 1 : 1; previous = base;
        if (run >= homopolymer) { return false; }
        unsigned code;
        switch (base) { case 'A': code = 0; break; case 'C': code = 1; break;
            case 'G': code = 2; break; case 'T': code = 3; break; default: return false; }
        word = (word << 2) | code;
    }
    return true;
}
void CoarseDereplication::collectSeeds(std::string_view sequence) {
    querySeeds_.clear();
    for (size_t position = 0; position + seedLength <= sequence.size(); position += seedLength) {
        uint64_t word;
        if (seedWord(sequence, position, homopolymer_, word)) { querySeeds_.push_back({word, position}); }
    }
}
void CoarseDereplication::indexGroup(size_t id) {
    auto& bucket = buckets_[groups_[id]->searchLength / lengthBin];
    bucket.ids.push_back(id);
    for (const auto& seed : querySeeds_) { bucket.postings[seed].push_back(id); ++stats_.indexedSeeds; }
}
size_t CoarseDereplication::find(std::string_view sequence, std::string_view mate) { return search(sequence, false, npos, mate); }
size_t CoarseDereplication::search(std::string_view sequence, bool nearest, size_t original, std::string_view mate) {
    if (sequence.empty()) { return nearest ? original : npos; }
    if (sequence.size() > INT_MAX || mate.size() > INT_MAX) { throw std::length_error("Coarse dereplication sequence exceeds supported length"); }
    ++stats_.queries;
    size_t winner = nearest ? original : npos;
    bool fallback = false;
    ReadMerger verifier;
    auto identityTo = [&](size_t id, bool bounded) {
        ++stats_.comparisons;
        const auto mateAnchor = groups_[id]->mateKey();
        if (mate.empty() != mateAnchor.empty()) { throw std::invalid_argument("Cannot mix single and paired cluster comparisons"); }
        uint64_t errors = 0, overlap = 0;
        auto compare = [&](std::string_view query, std::string_view anchor) {
            const size_t length = std::min(query.size(), anchor.size());
            const double identity = verifier.percentIdentity(query.data(), anchor.data(), static_cast<int>(length),
                bounded ? static_cast<int>(mismatchBudget(length)) : -1);
            if (identity < 0.) { return false; }
            // Integral match counts avoid floating-point tie/boundary drift.
            errors += length - static_cast<uint64_t>(identity * length + 0.5); overlap += length;
            return true;
        };
        if (!compare(sequence, groups_[id]->key()) || (!mate.empty() && !compare(mate, mateAnchor))) {
            return std::pair<uint64_t, uint64_t>{1, 0};
        }
        // Pair proximity is combined errors / combined overlap, with each mate
        // independently required to pass its identity threshold when bounded.
        return std::pair<uint64_t, uint64_t>(errors, overlap);
    };
    auto bestIdentity = nearest ? identityTo(original, false) : std::pair<uint64_t, uint64_t>{0, 1};
    if (nearest && bestIdentity.first == 0) { return winner; }
    collectSeeds(sequence);
    for (const auto& entry : buckets_) {
        const auto& bucket = entry.second;
        const size_t lowerLength = entry.first * lengthBin;
        const size_t upperLength = lowerLength + lengthBin - 1;
        const size_t prefix = std::min(sequence.size(), lowerLength);
        const size_t required = mismatchBudget(std::min(sequence.size(), upperLength)) + 1;
        lists_.clear();
        for (const auto& seed : querySeeds_) {
            if (seed.position + seedLength > prefix) { break; }
            auto it = bucket.postings.find(seed);
            lists_.push_back(it == bucket.postings.end() ? nullptr : &it->second);
        }
        auto verify = [&](size_t id) {
            if (!nearest && id >= winner) { return true; }
            if (nearest && id == original) { return false; }
            const auto identity = identityTo(id, true);
            if (identity.first > identity.second) { return false; }
            if (nearest) {
                // Keep the original cluster on ties. Otherwise use the earliest
                // cluster ID for equal best targets, independent of length bins.
                const uint64_t candidateErrors = identity.first * bestIdentity.second;
                const uint64_t bestErrors = bestIdentity.first * identity.second;
                if (candidateErrors < bestErrors || (candidateErrors == bestErrors && winner != original && id < winner)) {
                    bestIdentity = identity; winner = id;
                }
                return false; // Nearest search must consider every candidate.
            }
            winner = id; return true;
        };
        if (lists_.size() < required) {
            // Low complexity, ambiguous bases, short overlap, or a low identity
            // threshold must affect speed only, never candidate sensitivity.
            fallback = true;
            for (size_t id : bucket.ids) { if (verify(id)) { break; } }
        } else {
            std::partial_sort(lists_.begin(), lists_.begin() + required, lists_.end(),
                [](const auto* a, const auto* b) { return (a ? a->size() : 0) < (b ? b->size() : 0); });
            // Merge sorted posting lists lazily. Do not copy E+1 potentially
            // huge lists into a temporary vector just to find the first match.
            heap_.clear();
            for (size_t i = 0; i < required; ++i) {
                if (lists_[i] && !lists_[i]->empty()) { heap_.push_back({lists_[i], 0}); }
            }
            auto later = [](const Cursor& a, const Cursor& b) {
                return (*a.ids)[a.offset] > (*b.ids)[b.offset];
            };
            std::make_heap(heap_.begin(), heap_.end(), later);
            size_t previous = npos;
            while (!heap_.empty()) {
                std::pop_heap(heap_.begin(), heap_.end(), later);
                auto cursor = heap_.back(); heap_.pop_back();
                const size_t id = (*cursor.ids)[cursor.offset++];
                if (id != previous && verify(id)) { break; }
                previous = id;
                if (cursor.offset != cursor.ids->size()) {
                    heap_.push_back(cursor); std::push_heap(heap_.begin(), heap_.end(), later);
                }
            }
        }
    }
    stats_.fallbackQueries += fallback;
    return winner;
}
size_t CoarseDereplication::add(std::string_view search, std::string_view fullSequence,
    std::shared_ptr<DNA> dna, std::shared_ptr<DNA> mate, std::shared_ptr<DNA> merged, int sample, bool passed, int searchLimit) {
    if (!dna || search.empty() || sample < 0) { return npos; }
    if (!fullSequence.starts_with(search)) { throw std::invalid_argument("Dereplication key must be a prefix of its sequence"); }
    const bool paired = storeQualities_ && bool(mate);
    if (paired && mate->length() == 0) { return npos; } // A mate removed by trimming cannot form a complete pair.
    if (storeQualities_ && (merged || dna->getMergeLength() >= 0 || dna->getQual().size() < fullSequence.size() ||
        dna->getSeqPseudoView() != fullSequence || (mate && (mate->getMergeLength() >= 0 || mate->getQual().size() < mate->length())))) {
        throw std::invalid_argument("-derepStoreQuals requires complete unmerged reads and qualities for each mate");
    }
    if (storeQualities_ && !groups_.empty() && paired != bool(groups_[0]->mateReference)) {
        throw std::invalid_argument("-derepStoreQuals cannot mix single and paired fragments");
    }
    if (paired) { passed = passed && mate->isGreenQual(); }
    const auto mateKey = paired ? mate->getSeqPseudoView().substr(0, searchLimit > 0 ? std::min(size_t(mate->length()), size_t(searchLimit)) : mate->length()) : std::string_view{};
    size_t id = find(search, mateKey);
    if (id == npos) {
        if (!passed) { return npos; }
        id = groups_.size();
        auto group = std::make_unique<Group>();
        group->reference = std::make_shared<const std::string>((storeDiffs_ || reassign_ || storeQualities_) ? fullSequence : search);
        group->searchLength = search.size();
        if (paired) { group->mateReference = std::make_shared<const std::string>(mate->getSeqPseudoView()); group->mateSearchLength = mateKey.size(); }
        // One representative across all merge offsets; no per-base quality sums.
        group->representatives.addNewDNAuniq(dna, mate, merged, 0, -1, true);
        if (storeDiffs_ || reassign_ || storeQualities_) { group->representatives[0]->attachDiffDNA(std::make_shared<diffDNA>(group->reference, reassign_, storeQualities_, group->mateReference)); }
        groups_.push_back(std::move(group));
        indexGroup(id);
    } else if (!paired || passed) {
        // Reuse whoIsBetter without averaging qualities across alleles. Sample
        // and variant counts are recorded together below, after selection.
        groups_[id]->representatives[0]->matchedCoarseDNA(dna, mate, merged, -1);
    }
    auto& representative = groups_[id]->representatives[0];
    if (paired) { representative->recordDiffPairQuality(*dna, *mate, sample); }
    else if (storeQualities_) { representative->recordDiffQuality(*dna, sample); }
    else { representative->recordDiffDNA(fullSequence, sample); }
    return id;
}
CoarseDereplication::Refinement CoarseDereplication::refine(bool searchWithMerge, int searchLimit) const {
    if (!reassign_) { throw std::logic_error("Final reassignment requires sample-tracked variants"); }
    Refinement result;
    CoarseDereplication representatives(identity_, homopolymer_, false);
    auto key = [searchLimit](std::string_view sequence) {
        return sequence.substr(0, searchLimit > 0 ? std::min(sequence.size(), size_t(searchLimit)) : sequence.size());
    };
    for (const auto& group : groups_) {
        auto dna = group->representatives[0]->cloneForOutput();
        auto source = searchWithMerge && dna->getMerge() ? dna->getMerge() : dna;
        auto entry = std::make_unique<Group>();
        entry->reference = std::make_shared<const std::string>(source->getSeqPseudoView());
        entry->searchLength = key(*entry->reference).size();
        if (group->mateReference) {
            entry->mateReference = std::make_shared<const std::string>(dna->getPair()->getSeqPseudoView());
            entry->mateSearchLength = key(*entry->mateReference).size();
        }
        if (!entry->searchLength || entry->searchLength > INT_MAX) {
            throw std::length_error("Invalid final representative search length");
        }
        // Insert all selected representatives, even if they now match each other.
        representatives.collectSeeds(entry->key());
        representatives.groups_.push_back(std::move(entry));
        representatives.indexGroup(representatives.groups_.size() - 1);
        dna->attachDiffDNA((storeDiffs_ || storeQualities_) ? std::make_shared<diffDNA>(group->reference, false, storeQualities_, group->mateReference) : nullptr);
        result.representatives.push_back(std::move(dna));
    }
    for (size_t from = 0; from < groups_.size(); ++from) {
        auto differences = groups_[from]->representatives[0]->getDiffDNA();
        auto moveVariant = [&](std::string_view sequence, std::string_view mate, uint64_t count, const diffDNA::SampleCounts& samples, const diffDNA::QualityRecord* quality) {
            const size_t to = representatives.search(key(sequence), true, from, key(mate));
            auto& destination = result.representatives[to];
            // Output differences are rebuilt once against the destination's
            // immutable initial reference, without retaining per-read metadata.
            if (!mate.empty()) { destination->recordDiffPair(sequence, mate, samples, *quality); }
            else if (quality) { destination->recordDiffDNA(sequence, samples, *quality); }
            else { destination->recordDiffDNA(sequence, samples); }
            if (to != from) { result.movedReads += count; ++result.movedVariants; }
        };
        if (differences->storesPairs()) {
            differences->forEachPair([&](auto sequence, auto mate, auto count, const auto* samples, const auto* quality) {
                moveVariant(sequence, mate, count, *samples, quality);
            });
        } else if (storeQualities_) {
            differences->forEachWithQualities([&](auto sequence, auto count, const auto* samples, const auto& quality) {
                moveVariant(sequence, {}, count, *samples, &quality);
            });
        } else {
            differences->forEachWithSamples([&](auto sequence, auto count, const auto& samples) {
                moveVariant(sequence, {}, count, samples, nullptr);
            });
        }
    }
    // Empty groups (possible with truncated/prefix representatives) emit nothing.
    for (auto& dna : result.representatives) {
        if (dna->getDerepMap().empty()) { dna.reset(); }
    }
    result.searches = representatives.stats();
    return result;
}
void CoarseDereplication::clear() {
    buckets_.clear(); groups_.clear(); stats_ = {};
    querySeeds_.clear(); lists_.clear(); heap_.clear();
}
