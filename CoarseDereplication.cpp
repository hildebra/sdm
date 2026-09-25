// SPDX-License-Identifier: GPL-3.0-or-later
#include "CoarseDereplication.h"
#include "CompactDNA.h"
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
CoarseDereplication::CoarseDereplication(uint64_t identity, unsigned homopolymer, bool storeDiffs, bool reassign, bool storeQualities, bool exactOutput)
    : identity_(identity), homopolymer_(homopolymer), storeDiffs_(storeDiffs), reassign_(reassign), storeQualities_(storeQualities), exactOutput_(exactOutput) {
    if (exactOutput && reassign) { throw std::invalid_argument("Reassignment requires coarse cluster output"); }
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
size_t CoarseDereplication::storagePartitionHash(std::string_view sequence) const {
    // Reuse the candidate index's eligible 20-mers. Nearby keys can share a
    // compression reference instead of being scattered by whole-key hashing.
    // This only routes storage; exact key equality/admission remains unchanged.
    size_t selected = std::numeric_limits<size_t>::max();
    bool found = false;
    for (size_t position = 0; position + seedLength <= sequence.size(); position += seedLength) {
        uint64_t word;
        if (!seedWord(sequence, position, homopolymer_, word)) { continue; }
        selected = std::min(selected, SeedHash{}(SeedKey{word, position})); found = true;
    }
    return found ? selected : robin_hood::hash<std::string_view>{}(sequence);
}
void CoarseDereplication::collectSeeds(std::string_view sequence, std::vector<SeedKey>& seeds) const {
    seeds.clear();
    for (size_t position = 0; position + seedLength <= sequence.size(); position += seedLength) {
        uint64_t word;
        if (seedWord(sequence, position, homopolymer_, word)) { seeds.push_back({word, position}); }
    }
}
void CoarseDereplication::indexGroup(size_t id) {
    const size_t length = groups_[id]->searchLength;
    auto& bucket = buckets_[length / lengthBin];
    bucket.minLength = std::min(bucket.minLength, length);
    bucket.maxLength = std::max(bucket.maxLength, length);
    bucket.ids.push_back(id);
    for (const auto& seed : scratch_.seeds) { bucket.postings[seed].push_back(id); ++stats_.indexedSeeds; }
}
size_t CoarseDereplication::find(std::string_view sequence, std::string_view mate) {
    ExactAssignment* exactMatch;
    return find(sequence, mate, exactMatch, !mate.empty());
}
CoarseDereplication::ExactAssignment* CoarseDereplication::exactLookup(std::string_view sequence, bool paired) {
    const auto found = exactAssignments_.find(sequence, robin_hood::is_transparent_tag{});
    if (found == exactAssignments_.end()) { return nullptr; }
    if (paired != bool(groups_[found->second.group]->mateReference)) {
        throw std::invalid_argument("Cannot mix single and paired cluster comparisons");
    }
    return &found->second;
}
size_t CoarseDereplication::find(std::string_view sequence, std::string_view mate, ExactAssignment*& exactMatch, bool paired, bool searchCandidates) {
    exactMatch = exactLookup(sequence, paired);
    if (exactMatch) {
        ++stats_.queries; ++stats_.exactHits;
        return exactMatch->group;
    }
    // Default output only admits failed observations at an already established
    // exact key. Searching approximate anchors cannot change that decision.
    if (!searchCandidates) {
        ++stats_.queries; ++stats_.rejectedExactMisses;
        return npos;
    }
    // A new R1 still uses the coarse matcher. Mere lookups never create an
    // exact assignment: only an admitted observation owns a dereplicate.
    return search(sequence, false, npos, mate, paired);
}
size_t CoarseDereplication::search(std::string_view sequence, bool nearest, size_t original, std::string_view mate, bool paired) {
    return search(sequence, nearest, original, mate, paired, scratch_, stats_, 0, true);
}
size_t CoarseDereplication::search(std::string_view sequence, bool nearest, size_t original, std::string_view mate, bool paired,
    SearchScratch& scratch, Stats& stats, size_t minId, bool countQuery) const {
    if (sequence.empty()) { return nearest ? original : npos; }
    if (sequence.size() > INT_MAX || mate.size() > INT_MAX) { throw std::length_error("Coarse dereplication sequence exceeds supported length"); }
    stats.queries += countQuery;
    size_t winner = nearest ? original : npos;
    // Exact-parent output only needs any anchor within budget (a compression
    // reference): first fit, cheapest. Coarse-cluster output uses the nearest
    // anchor at admission, so membership does not depend on arrival order
    // beyond ties (which go to the earliest anchor).
    const bool bestFit = !nearest && !exactOutput_;
    bool fallback = false;
    auto identityTo = [&](size_t id, bool bounded) {
        ++stats.comparisons;
        const auto mateAnchor = groups_[id]->mateKey();
        if (paired != bool(groups_[id]->mateReference)) { throw std::invalid_argument("Cannot mix single and paired cluster comparisons"); }
        uint64_t errors = 0, overlap = 0;
        auto compare = [&](std::string_view query, std::string_view anchor) {
            const size_t length = std::min(query.size(), anchor.size());
            if (!length) { return true; } // A logical trim can hide the entire retained mate.
            const double identity = ReadMerger::percentIdentity(query.data(), anchor.data(), static_cast<int>(length),
                bounded ? static_cast<int>(mismatchBudget(length)) : -1);
            if (identity < 0.) { return false; }
            // Integral match counts avoid floating-point tie/boundary drift.
            errors += length - static_cast<uint64_t>(identity * length + 0.5); overlap += length;
            return true;
        };
        if (!compare(sequence, groups_[id]->key()) || (paired && !compare(mate, mateAnchor))) {
            return std::pair<uint64_t, uint64_t>{1, 0};
        }
        // Pair proximity is combined errors / combined overlap, with each mate
        // independently required to pass its identity threshold when bounded.
        return std::pair<uint64_t, uint64_t>(errors, overlap);
    };
    auto bestIdentity = nearest ? identityTo(original, false) : std::pair<uint64_t, uint64_t>{0, 1};
    if (nearest && bestIdentity.first == 0) { return winner; }
    auto& querySeeds = scratch.seeds;
    auto& lists = scratch.lists;
    auto& heap = scratch.heap;
    collectSeeds(sequence, querySeeds);
    for (const auto& entry : buckets_) {
        const auto& bucket = entry.second;
        // Use populated search-key bounds, not the nominal 128-base interval:
        // a bin of 200 bp anchors guarantees 200 bp of overlap, not just 128.
        // Every seed must fit the shortest overlap, while E bounds errors in
        // the longest overlap to preserve sensitivity for mixed-length keys.
        const size_t prefix = std::min(sequence.size(), bucket.minLength);
        const size_t required = mismatchBudget(std::min(sequence.size(), bucket.maxLength)) + 1;
        lists.clear();
        for (const auto& seed : querySeeds) {
            if (seed.position + seedLength > prefix) { break; }
            auto it = bucket.postings.find(seed);
            lists.push_back(it == bucket.postings.end() ? nullptr : &it->second);
        }
        auto verify = [&](size_t id) {
            if (!nearest && !bestFit && id >= winner) { return true; }
            if (nearest && id == original) { return false; }
            // After an exact hit only a smaller ID (a tie) can still win; IDs
            // ascend within a bucket, so stop this bucket at the first larger one.
            if (bestFit && winner != npos && bestIdentity.first == 0 && id >= winner) { return true; }
            const auto identity = identityTo(id, true);
            if (identity.first > identity.second) { return false; }
            if (nearest || bestFit) {
                // Keep the original cluster on ties. Otherwise use the earliest
                // cluster ID for equal best targets, independent of length bins.
                const uint64_t candidateErrors = identity.first * bestIdentity.second;
                const uint64_t bestErrors = bestIdentity.first * identity.second;
                if (winner == npos || candidateErrors < bestErrors ||
                    (candidateErrors == bestErrors && winner != original && id < winner)) {
                    bestIdentity = identity; winner = id;
                }
                // Nearest/best-fit searches consider every candidate; an exact
                // hit ends this bucket (later IDs in it cannot tie-break).
                return bestFit && identity.first == 0;
            }
            winner = id; return true;
        };
        // ID lists are ascending: skip anchors older than minId.
        auto firstFrom = [minId](const std::vector<size_t>& ids) {
            return minId ? static_cast<size_t>(std::lower_bound(ids.begin(), ids.end(), minId) - ids.begin()) : size_t(0);
        };
        if (lists.size() < required) {
            // Low complexity, ambiguous bases, short overlap, or a low identity
            // threshold must affect speed only, never candidate sensitivity.
            fallback = true;
            for (size_t i = firstFrom(bucket.ids); i < bucket.ids.size(); ++i) { if (verify(bucket.ids[i])) { break; } }
        } else {
            std::partial_sort(lists.begin(), lists.begin() + required, lists.end(),
                [](const auto* a, const auto* b) { return (a ? a->size() : 0) < (b ? b->size() : 0); });
            // Merge sorted posting lists lazily. Do not copy E+1 potentially
            // huge lists into a temporary vector just to find the first match.
            heap.clear();
            for (size_t i = 0; i < required; ++i) {
                if (!lists[i]) { continue; }
                const size_t start = firstFrom(*lists[i]);
                if (start < lists[i]->size()) { heap.push_back({lists[i], start}); }
            }
            auto later = [](const Cursor& a, const Cursor& b) {
                return (*a.ids)[a.offset] > (*b.ids)[b.offset];
            };
            std::make_heap(heap.begin(), heap.end(), later);
            size_t previous = npos;
            while (!heap.empty()) {
                std::pop_heap(heap.begin(), heap.end(), later);
                auto cursor = heap.back(); heap.pop_back();
                const size_t id = (*cursor.ids)[cursor.offset++];
                if (id != previous && verify(id)) { break; }
                previous = id;
                if (cursor.offset != cursor.ids->size()) {
                    heap.push_back(cursor); std::push_heap(heap.begin(), heap.end(), later);
                }
            }
        }
    }
    stats.fallbackQueries += fallback;
    return winner;
}
std::string_view CoarseDereplication::mateKeyFor(const std::shared_ptr<DNA>& mate, bool paired, int searchLimit) const {
    return paired ? mate->getSeqPseudoView().substr(0, searchLimit > 0 ? std::min(size_t(mate->length()), size_t(searchLimit)) : mate->length()) : std::string_view{};
}
CoarseDereplication::Probe CoarseDereplication::probe(std::string_view search, const std::shared_ptr<DNA>& dna,
    const std::shared_ptr<DNA>& mate, bool passed, int searchLimit, SearchScratch& scratch) const {
    // Mirrors the lookup that add() performs; validation stays in add().
    Probe result;
    result.snapshot = groups_.size();
    if (!dna || search.empty()) { return result; }
    const bool paired = storeQualities_ && bool(mate);
    if (paired && mate->mem_length() == 0) { return result; }
    if (paired && !exactOutput_) { passed = passed && mate->isGreenQual(); }
    if (!(passed || !exactOutput_)) { return result; }
    if (exactAssignments_.find(search, robin_hood::is_transparent_tag{}) != exactAssignments_.end()) { return result; }
    result.group = this->search(search, false, npos, mateKeyFor(mate, paired, searchLimit), paired,
        scratch, result.stats, 0, true);
    result.searched = true;
    return result;
}
void CoarseDereplication::apply(const Update& update) const {
    auto& representatives = *update.target;
    if (update.kind == Update::Observation) {
        representatives.addObservation(update.dna, update.mate, update.merged, update.mergePosition, update.sample,
            update.derepAsFasta, update.collectMergedQualities);
        return;
    }
    if (update.kind == Update::Matched) {
        representatives[0]->matchedDNA(update.dna, update.mate, update.merged, -1, update.derepAsFasta, update.collectMergedQualities);
    } else if (update.kind == Update::MatchedCoarse) {
        representatives[0]->matchedCoarseDNA(update.dna, update.mate, update.merged, -1);
    }
    auto& representative = representatives[0];
    if (update.paired) { representative->recordDiffPairQuality(*update.dna, *update.mate, update.sample); }
    else if (storeQualities_) { representative->recordDiffQuality(*update.dna, update.sample); }
    else { representative->recordDiffDNA(update.fullSequence, update.sample); }
}
size_t CoarseDereplication::add(std::string_view search, std::string_view fullSequence,
    std::shared_ptr<DNA> dna, std::shared_ptr<DNA> mate, std::shared_ptr<DNA> merged, int sample, bool passed, int searchLimit,
    bool derepAsFasta, bool collectMergedQualities, int mergePosition, const Probe* probe, std::vector<Update>* deferred) {
    if (!dna || search.empty() || sample < 0) { return npos; }
    if (!fullSequence.starts_with(search)) { throw std::invalid_argument("Dereplication key must be a prefix of its sequence"); }
    const bool paired = storeQualities_ && bool(mate);
    if (paired && mate->mem_length() == 0) { return npos; } // A physically removed mate cannot form a complete pair.
    if (storeQualities_ && (merged || dna->getMergeLength() >= 0 || dna->getQual().size() < dna->mem_length() ||
        dna->getSeqPseudoView() != fullSequence || (mate && (mate->getMergeLength() >= 0 || mate->getQual().size() < mate->mem_length())))) {
        throw std::invalid_argument("-derepStoreQuals requires complete unmerged reads and qualities for each mate");
    }
    if (storeQualities_ && !groups_.empty() && paired != bool(groups_[0]->mateReference)) {
        throw std::invalid_argument("-derepStoreQuals cannot mix single and paired fragments");
    }
    // Exact output follows ordinary R1 admission. A failed R2 must not
    // discard a passing R1 or prevent it from competing as representative.
    if (paired && !exactOutput_) { passed = passed && mate->isGreenQual(); }
    const auto mateKey = mateKeyFor(mate, paired, searchLimit);
    ExactAssignment* exactMatch;
    size_t id;
    if (!probe || !probe->searched) {
        id = find(search, mateKey, exactMatch, paired, passed || !exactOutput_);
    } else if ((exactMatch = exactLookup(search, paired))) {
        // An earlier commit created this exact key after the probe; the probe's
        // search is discarded and accounting matches a serial exact hit.
        ++stats_.queries; ++stats_.exactHits;
        id = exactMatch->group;
    } else {
        stats_.add(probe->stats);
        id = probe->group;
        if (probe->snapshot < groups_.size()) {
            // Anchors created after the snapshot: first fit only needs them when
            // the snapshot had no match; best fit (cluster output) may find a
            // nearer one among them, so it searches everything again.
            if (!exactOutput_) { id = this->search(search, false, npos, mateKey, paired, scratch_, stats_, 0, false); }
            else if (id == npos) { id = this->search(search, false, npos, mateKey, paired, scratch_, stats_, probe->snapshot, false); }
        }
    }
    // Admission of an exact R1 dereplicate must not depend on whether a
    // different, nearby R1 happened to create a coarse anchor first.
    if (!passed && (id == npos || (exactOutput_ && !exactMatch))) { return npos; }
    const bool newGroup = id == npos;
    if (newGroup) {
        id = groups_.size();
        auto group = std::make_unique<Group>();
        group->reference = std::make_shared<const std::string>(storeQualities_ ? std::string_view(dna->getSequence()) : (storeDiffs_ || reassign_) ? fullSequence : search);
        group->searchLength = search.size();
        if (paired) { group->mateReference = std::make_shared<const std::string>(mate->getSequence()); group->mateSearchLength = mateKey.size(); }
        if (storageOnly()) {
            group->compactReferences = std::make_shared<const CompactDNAReferences>(*dna, mate.get(), merged.get());
            // The immutable search key is a view of an already shared backing
            // reference; do not allocate a third copy for unmerged groups.
            group->reference = fullSequence == dna->getSeqPseudoView() ?
                group->compactReferences->first : group->compactReferences->merged;
        }
        groups_.push_back(std::move(group));
        // Index the seeds of this key (a probe searched with its own scratch).
        collectSeeds(search);
        indexGroup(id);
    }
    const bool newExact = !exactMatch;
    if (newExact) {
        exactMatch = &exactAssignments_.emplace(std::string(search), ExactAssignment{id,
            storageOnly() ? std::make_unique<DNAuniqSet>(groups_[id]->compactReferences) :
                exactOutput_ ? std::make_unique<DNAuniqSet>() : nullptr}).first->second;
    }
    auto& group = *groups_[id];
    auto& representatives = exactOutput_ ? *exactMatch->representatives : group.representatives;
    Update update;
    update.target = &representatives;
    update.paired = paired;
    update.derepAsFasta = derepAsFasta;
    update.collectMergedQualities = collectMergedQualities;
    update.sample = sample;
    update.mergePosition = mergePosition;
    update.fullSequence = fullSequence;
    if (storageOnly()) {
        // Internal search groups must not replace ordinary merge-offset subgroups
        // or their counts, selected pairs, and per-base quality accumulation.
        update.kind = Update::Observation;
    } else if (exactOutput_ ? newExact : newGroup) {
        // Reuse ordinary paired seed selection, now within the exact R1 key.
        // Variants share immutable coarse references without duplicating counts
        // or maintaining a second set of coarse representative qualities.
        // A new target is initialized here, before anyone else can update it.
        representatives.addNewDNAuniq(dna, mate, merged, 0, -1, true);
        if (storeDiffs_ || reassign_ || storeQualities_) {
            representatives[0]->attachDiffDNA(std::make_shared<diffDNA>(group.reference, reassign_, storeQualities_, group.mateReference));
        }
        update.kind = Update::RecordOnly;
    } else if (!paired || passed) {
        // Exact output must collect the same quality evidence as ordinary
        // dereplication, including averaged main FASTQ output. Variant export
        // remains a separate mode with selected per-variant qualities.
        update.kind = exactOutput_ && !storeQualities_ ? Update::Matched : Update::MatchedCoarse;
    } else {
        update.kind = Update::RecordOnly;
    }
    update.dna = std::move(dna);
    update.mate = std::move(mate);
    update.merged = std::move(merged);
    if (deferred) { deferred->push_back(std::move(update)); }
    else { apply(update); }
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
            // Full HQ variants can extend past their logical filtering lengths.
            // Reassignment must use the same logical/search view as admission.
            const auto searchSequence = quality ? sequence.substr(0, quality->seed->stats.logical) : sequence;
            const auto searchMate = quality && quality->seed->mate ? mate.substr(0, quality->seed->mate->stats.logical) : mate;
            const size_t to = representatives.search(key(searchSequence), true, from, key(searchMate), bool(groups_[from]->mateReference));
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
    exactAssignments_.clear();
    buckets_.clear(); groups_.clear(); stats_ = {};
    scratch_ = {};
}
