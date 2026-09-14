// SPDX-License-Identifier: GPL-3.0-or-later
#include "diffDNA.h"
#include "DNA.h"
#include "SeedSelection.h"
#include <type_traits>
#include <algorithm>
#include <istream>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <variant>
#include <optional>

namespace {
template<class Put> void putVar(uint64_t value, Put put) {
    while (value >= 128) { put(static_cast<uint8_t>(value) | 128); value >>= 7; }
    put(static_cast<uint8_t>(value));
}
template<class Get> uint64_t getVar(Get get) {
    uint64_t value = 0;
    for (unsigned shift = 0; shift < 64; shift += 7) {
        const uint8_t byte = get();
        if (shift == 63 && byte > 1) { throw std::runtime_error("diffDNA integer overflow"); }
        value |= uint64_t(byte & 127) << shift;
        if (!(byte & 128)) { return value; }
    }
    throw std::runtime_error("Invalid diffDNA integer");
}
size_t varBytes(uint64_t value) {
    size_t bytes = 1; while (value >= 128) { value >>= 7; ++bytes; } return bytes;
}
size_t asSize(uint64_t value) {
    if (value > std::numeric_limits<size_t>::max()) { throw std::runtime_error("diffDNA length overflow"); }
    return static_cast<size_t>(value);
}
template<class Get> std::string decode(const std::string& reference, Get get) {
    const size_t length = asSize(getVar(get));
    const size_t count = asSize(getVar(get));
    if (count > length || (length > reference.size() && length - reference.size() > count)) {
        throw std::runtime_error("Invalid diffDNA length/change count");
    }
    std::string sequence = reference.substr(0, length);
    // Grow suffixes only as their explicitly stored bases arrive. An invalid
    // length in a damaged input must not cause an unbounded allocation.
    size_t next = 0;
    for (size_t c = 0; c < count; ++c) {
        const size_t gap = asSize(getVar(get));
        if (next >= length || gap >= length - next) { throw std::runtime_error("Invalid diffDNA position"); }
        next += gap;
        const char base = static_cast<char>(get());
        if (next < sequence.size()) { sequence[next] = base; }
        else if (next == sequence.size()) { sequence.push_back(base); }
        else { throw std::runtime_error("Missing diffDNA suffix bases"); }
        ++next;
    }
    if (sequence.size() != length) { throw std::runtime_error("Missing diffDNA suffix bases"); }
    return sequence;
}
}

namespace {
std::string encodePatch(std::string_view ref, std::string_view sequence) {
    if (sequence == ref) { return {}; }
    std::string patch;
    size_t count = 0;
    for (size_t i = 0; i < sequence.size(); ++i) { count += i >= ref.size() || sequence[i] != ref[i]; }
    auto put = [&](uint8_t byte) { patch.push_back(static_cast<char>(byte)); };
    putVar(sequence.size(), put); putVar(count, put);
    size_t next = 0;
    for (size_t i = 0; i < sequence.size(); ++i) {
        if (i < ref.size() && sequence[i] == ref[i]) { continue; }
        putVar(i - next, put); put(static_cast<uint8_t>(sequence[i])); next = i + 1;
    }
    return patch;
}
std::string decodePatch(const std::string& reference, std::string_view patch) {
    if (patch.empty()) { return reference; }
    size_t cursor = 0;
    auto get = [&]() -> uint8_t {
        if (cursor == patch.size()) { throw std::runtime_error("Truncated diffDNA patch"); }
        return static_cast<uint8_t>(patch[cursor++]);
    };
    auto sequence = decode(reference, get);
    if (cursor != patch.size()) { throw std::runtime_error("Trailing diffDNA patch data"); }
    return sequence;
}
std::pair<std::string_view, std::string_view> splitPairPatch(std::string_view patch) {
    if (patch.empty()) { return {}; }
    size_t cursor = 0;
    auto get = [&]() -> uint8_t {
        if (cursor == patch.size()) { throw std::runtime_error("Truncated diffDNA pair key"); }
        return static_cast<uint8_t>(patch[cursor++]);
    };
    const size_t first = asSize(getVar(get));
    if (first > patch.size() - cursor) { throw std::runtime_error("Invalid diffDNA pair key length"); }
    return {patch.substr(cursor, first), patch.substr(cursor + first)};
}
uint64_t totalCounts(std::span<const diffDNA::SampleCount> counts) {
    uint64_t sum = 0;
    for (const auto& entry : counts) {
        if (entry.first < 0) { throw std::invalid_argument("diffDNA requires nonnegative sample IDs"); }
        if (entry.second > std::numeric_limits<uint64_t>::max() - sum) { throw std::overflow_error("diffDNA sample total overflow"); }
        sum += entry.second;
    }
    return sum;
}
void addSparse(diffDNA::SampleCounts& counts, int sample, uint64_t copies) {
    if (sample < 0 || !copies) { return; }
    auto found = std::find_if(counts.begin(), counts.end(), [sample](const auto& entry) { return entry.first == sample; });
    if (found == counts.end()) { counts.emplace_back(sample, copies); }
    else {
        if (copies > std::numeric_limits<uint64_t>::max() - found->second) { throw std::overflow_error("diffDNA sample count overflow"); }
        found->second += copies;
    }
}
}
struct diffDNA::SequenceData {
    template<class Counts> struct Qualified {
        Counts counts{};
        QualityRecord quality;
    };
    using Counts = robin_hood::unordered_map<std::string, uint64_t>;
    using Samples = robin_hood::unordered_map<std::string, SampleCounts>;
    using QualityCounts = robin_hood::unordered_map<std::string, Qualified<uint64_t>>;
    using QualitySamples = robin_hood::unordered_map<std::string, Qualified<SampleCounts>>;
    SequenceData(std::shared_ptr<const std::string> ref, bool trackSamples, bool storeQualities, std::shared_ptr<const std::string> mateRef = nullptr) : reference(std::move(ref)), mateReference(std::move(mateRef)) {
        if (storeQualities) {
            if (trackSamples) { variants.emplace<QualitySamples>(); }
            else { variants.emplace<QualityCounts>(); }
        } else if (trackSamples) { variants.emplace<Samples>(); }
    }
    template<class T> static const SampleCounts* samples(const T& value) {
        if constexpr (std::is_same_v<T, SampleCounts>) { return &value; }
        else if constexpr (requires { value.counts; }) { return samples(value.counts); }
        else { return nullptr; }
    }
    template<class T> static uint64_t count(const T& value) {
        if constexpr (std::is_same_v<T, uint64_t>) { return value; }
        else if constexpr (std::is_same_v<T, SampleCounts>) { return totalCounts(value); }
        else { return count(value.counts); }
    }
    template<class T> static const QualityRecord* quality(const T& value) {
        if constexpr (requires { value.quality; }) { return &value.quality; }
        else { return nullptr; }
    }
    template<class T> static void addCounts(T& value, uint64_t copies, std::span<const SampleCount> samples) {
        if constexpr (std::is_same_v<T, uint64_t>) { value += copies; }
        else if constexpr (std::is_same_v<T, SampleCounts>) {
            for (const auto& entry : samples) { addSparse(value, entry.first, entry.second); }
        } else { addCounts(value.counts, copies, samples); }
    }
    static void selectQuality(QualityRecord& target, DNA* read, DNA* mate, const QualityRecord* retained) {
        if (retained && !target.seed) { target = *retained; return; }
        if (read && target.seed && (!read->isGreenQual() || (mate && !mate->isGreenQual()))) { return; }
        const SeedReadStats incoming = read ? SeedReadStats(*read) : retained->seed->stats;
        std::optional<SeedReadStats> mateStats;
        if (mate) { mateStats.emplace(*mate); }
        const SeedReadStats* incomingMate = mate ? &*mateStats :
            (retained && retained->seed->mate ? &retained->seed->mate->stats : nullptr);
        const bool passing = incoming.passing && (!incomingMate || incomingMate->passing);
        SeedComparisonHistory history = target.history;
        if (retained) { history.merge(retained->history); }
        bool replace = !target.seed;
        if (!replace && passing) {
            if (!target.seed->stats.passing || (target.seed->mate && !target.seed->mate->stats.passing)) { replace = true; }
            else {
                float best = -1.f; // Derep representative rules, not external-reference seed selection.
                replace = seed_selection::better<const SeedReadStats>(&incoming, incomingMate, nullptr,
                    &target.seed->stats, target.seed->mate ? &target.seed->mate->stats : nullptr, nullptr, best, false, &history);
            }
        }
        if (replace) {
            if (retained) { target.seed = retained->seed; }
            else {
                const auto& qualities = read->getQual();
                std::shared_ptr<const QualityRead> selectedMate;
                if (mate) {
                    const auto& q2 = mate->getQual();
                    selectedMate = std::make_shared<const QualityRead>(QualityRead{*incomingMate,
                        std::vector<qual_score>(q2.begin(), q2.begin() + incomingMate->logical)});
                }
                target.seed = std::make_shared<const QualitySeed>(QualitySeed{{incoming,
                    std::vector<qual_score>(qualities.begin(), qualities.begin() + incoming.logical)}, std::move(selectedMate)});
            }
        }
        target.history = history;
    }
    std::shared_ptr<const std::string> reference, mateReference;
    // Only nucleotide-only, untracked mode uses a separate reference counter.
    // Other modes store the reference under an empty patch key in the same map.
    uint64_t identical = 0;
    std::variant<Counts, Samples, QualityCounts, QualitySamples> variants;
};

diffDNA::diffDNA(std::shared_ptr<const std::string> reference, bool trackSamples, bool storeQualities, std::shared_ptr<const std::string> mateReference) {
    if (!reference || (mateReference && (reference->empty() || mateReference->empty()))) { throw std::invalid_argument("diffDNA requires complete references"); }
    sequences_ = std::make_shared<SequenceData>(std::move(reference), trackSamples, storeQualities, std::move(mateReference));
}
diffDNA::diffDNA(std::string reference, bool trackSamples, bool storeQualities)
    : diffDNA(std::make_shared<const std::string>(std::move(reference)), trackSamples, storeQualities) {}
bool diffDNA::storesPairs() const { return sequences_ && sequences_->mateReference; }
const std::string& diffDNA::mateReference() const {
    if (!storesPairs()) { throw std::logic_error("diffDNA has no mate reference"); }
    return *sequences_->mateReference;
}
bool diffDNA::tracksSamples() const {
    return sequences_ && (std::holds_alternative<SequenceData::Samples>(sequences_->variants) ||
        std::holds_alternative<SequenceData::QualitySamples>(sequences_->variants));
}
bool diffDNA::storesQualities() const {
    return sequences_ && (std::holds_alternative<SequenceData::QualityCounts>(sequences_->variants) ||
        std::holds_alternative<SequenceData::QualitySamples>(sequences_->variants));
}
const std::string& diffDNA::reference() const {
    if (!sequences_) { throw std::logic_error("Counts-only diffDNA has no sequence reference"); }
    return *sequences_->reference;
}
void diffDNA::detachSequences() {
    if (sequences_ && !sequences_.unique()) { sequences_ = std::make_shared<SequenceData>(*sequences_); }
}
void diffDNA::clear() {
    if (sequences_) {
        // Starting empty avoids copying a retained snapshot's entire variant map.
        if (!sequences_.unique()) { sequences_ = std::make_shared<SequenceData>(sequences_->reference, tracksSamples(), storesQualities(), sequences_->mateReference); }
        else {
            std::visit([](auto& variants) { variants.clear(); }, sequences_->variants);
            sequences_->identical = 0;
        }
    }
    occurrences_.clear(); reads_ = 0;
}
const diffDNA::SampleCounts& diffDNA::sampleCounts() const {
    if (tracksSamples() && occurrences_.empty() && reads_) {
        SampleCounts aggregate;
        std::visit([&](const auto& variants) {
            for (const auto& variant : variants) {
                if (const auto* samples = SequenceData::samples(variant.second)) {
                    for (const auto& entry : *samples) { addSparse(aggregate, entry.first, entry.second); }
                }
            }
        }, sequences_->variants);
        std::sort(aggregate.begin(), aggregate.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        occurrences_ = std::move(aggregate);
    }
    return occurrences_;
}
void diffDNA::addSample(int sample, uint64_t copies) {
    if (sample < 0 || !copies) { return; }
    if (sequences_) { throw std::logic_error("Nucleotide-enabled diffDNA requires a sequence with each sample count"); }
    // Keep ordinary exact dereplication out of the nucleotide/span machinery.
    if (copies > std::numeric_limits<uint64_t>::max() - reads_) { throw std::overflow_error("diffDNA read count overflow"); }
    addSparse(occurrences_, sample, copies);
    reads_ += copies;
}
void diffDNA::add(std::string_view sequence, uint64_t copies, int sample) {
    const SampleCount entry(sample, copies);
    record(sequence, copies, sample >= 0 ? std::span<const SampleCount>(&entry, 1) : std::span<const SampleCount>());
}
void diffDNA::add(std::string_view sequence, const SampleCounts& samples) {
    record(sequence, totalCounts(samples), samples);
}
void diffDNA::addPair(std::string_view r1, std::string_view r2, uint64_t copies, int sample) {
    if (!storesPairs()) { throw std::logic_error("diffDNA pair storage was not enabled"); }
    const SampleCount entry(sample, copies);
    record(r1, copies, sample >= 0 ? std::span<const SampleCount>(&entry, 1) : std::span<const SampleCount>(), nullptr, nullptr, r2);
}
void diffDNA::addPair(std::string_view r1, std::string_view r2, const SampleCounts& samples, const QualityRecord& quality) {
    if (!storesPairs()) { throw std::logic_error("diffDNA pair storage was not enabled"); }
    record(r1, totalCounts(samples), samples, nullptr, &quality, r2);
}
void diffDNA::addPairQuality(DNA& r1, DNA& r2, uint64_t copies, int sample) {
    if (!storesPairs() || !storesQualities()) { throw std::logic_error("diffDNA paired quality retention was not enabled"); }
    const SampleCount entry(sample, copies);
    record(r1.getSeqPseudoView(), copies, sample >= 0 ? std::span<const SampleCount>(&entry, 1) : std::span<const SampleCount>(), &r1, nullptr, r2.getSeqPseudoView(), &r2);
}
void diffDNA::addQuality(DNA& dna, uint64_t copies, int sample) {
    if (!storesQualities()) { throw std::logic_error("diffDNA quality retention was not enabled"); }
    const SampleCount entry(sample, copies);
    record(dna.getSeqPseudoView(), copies, sample >= 0 ? std::span<const SampleCount>(&entry, 1) : std::span<const SampleCount>(), &dna);
}
void diffDNA::add(std::string_view sequence, const SampleCounts& samples, const QualityRecord& quality) {
    record(sequence, totalCounts(samples), samples, nullptr, &quality);
}
void diffDNA::record(std::string_view sequence, uint64_t copies, std::span<const SampleCount> samples,
    DNA* read, const QualityRecord* quality, std::string_view mateSequence, DNA* mate) {
    if (copies > std::numeric_limits<uint64_t>::max() - reads_) { throw std::overflow_error("diffDNA read count overflow"); }
    if (!copies) { return; }
    if (tracksSamples() && samples.empty()) { throw std::invalid_argument("Sample-tracked diffDNA requires sample IDs"); }
    if (storesPairs() && sequence.empty()) { throw std::invalid_argument("diffDNA requires both mates to be nonempty"); }
    if (sequences_ && storesPairs() != !mateSequence.empty()) { throw std::invalid_argument("diffDNA requires consistent complete single or paired observations"); }
    if (storesQualities()) {
        auto validRead = [](DNA& d, std::string_view sequence) {
            return d.getMergeLength() < 0 && d.getSeqPseudoView() == sequence && !sequence.empty() && d.getQual().size() >= sequence.size();
        };
        auto validQuality = [](const QualityRead& q, std::string_view sequence) {
            return q.stats.logical == sequence.size() && q.qualities.size() == sequence.size() && q.stats.mergeLength < 0;
        };
        if (read) {
            if (!validRead(*read, sequence) || storesPairs() != bool(mate) || (mate && !validRead(*mate, mateSequence))) {
                throw std::invalid_argument("Quality retention requires complete unmerged reads and qualities for each mate");
            }
        } else if (!quality || !quality->seed || !validQuality(*quality->seed, sequence) ||
            storesPairs() != bool(quality->seed->mate) || (quality->seed->mate && !validQuality(*quality->seed->mate, mateSequence))) {
            throw std::invalid_argument("Quality-retaining diffDNA requires matching selected qualities");
        }
    }
    if (sequences_) {
        const auto& ref = reference();
        const bool identical = sequence == ref && (!storesPairs() || mateSequence == mateReference());
        std::string patch = encodePatch(ref, sequence);
        if (storesPairs() && !identical) {
            const auto first = std::move(patch); patch.clear();
            putVar(first.size(), [&](uint8_t byte) { patch.push_back(static_cast<char>(byte)); });
            patch += first; patch += encodePatch(mateReference(), mateSequence);
        }
        detachSequences();
        if (!tracksSamples() && !storesQualities() && identical) { sequences_->identical += copies; }
        else {
            std::visit([&](auto& variants) {
                auto [it, inserted] = variants.try_emplace(patch);
                try {
                    if constexpr (requires { it->second.quality; }) { SequenceData::selectQuality(it->second.quality, read, mate, quality); }
                    SequenceData::addCounts(it->second, copies, samples);
                } catch (...) {
                    if (inserted) { variants.erase(it); }
                    throw;
                }
            }, sequences_->variants);
        }
        if (tracksSamples()) { occurrences_.clear(); }

    }
    if (!tracksSamples()) {
        for (const auto& entry : samples) { addSparse(occurrences_, entry.first, entry.second); }
    }
    reads_ += copies;
}
void diffDNA::merge(const diffDNA& other) {
    if (other.reads_ > std::numeric_limits<uint64_t>::max() - reads_) { throw std::overflow_error("diffDNA merge count overflow"); }
    if (!other.reads_) { return; }
    if (!sequences_) {
        // Exact dereplication merges counts frequently; avoid copying its sparse
        // vector. The total-overflow check also covers every doubled entry.
        if (this == &other) {
            for (auto& entry : occurrences_) { entry.second *= 2; }
            reads_ *= 2;
        } else {
            for (const auto& entry : other.sampleCounts()) { addSparse(occurrences_, entry.first, entry.second); }
            reads_ += other.reads_;
        }
    } else {
        // A snapshot makes self-merging and shared variant data well-defined.
        const diffDNA source(other);
        if (!source.sequences_) { throw std::invalid_argument("Cannot merge unknown nucleotide observations into variant storage"); }
        if (storesPairs() != source.storesPairs()) { throw std::invalid_argument("Cannot merge single and paired variants"); }
        if (storesQualities() && !source.storesQualities()) { throw std::invalid_argument("Cannot merge variants without retained qualities"); }
        if (tracksSamples()) {
            if (!source.tracksSamples()) { throw std::invalid_argument("Cannot merge variants without per-variant sample provenance"); }
            source.visitFragment([&](auto sequence, auto mate, auto copies, const SampleCounts* samples, const QualityRecord* quality) {
                record(sequence, copies, *samples, nullptr, quality, mate);
            });
        } else {
            source.visitFragment([&](auto sequence, auto mate, auto copies, const SampleCounts*, const QualityRecord* quality) {
                record(sequence, copies, {}, nullptr, quality, mate);
            });
            for (const auto& entry : source.sampleCounts()) { addSparse(occurrences_, entry.first, entry.second); }
        }
    }
}
uint64_t diffDNA::referenceCount() const {
    if (!sequences_) { return 0; }
    if (!tracksSamples() && !storesQualities()) { return sequences_->identical; }
    return std::visit([](const auto& variants) -> uint64_t {
        const auto found = variants.find(std::string{});
        return found == variants.end() ? 0 : SequenceData::count(found->second);
    }, sequences_->variants);
}
size_t diffDNA::uniqueCount() const {
    if (!sequences_) { return 0; }
    return std::visit([](const auto& variants) { return variants.size(); }, sequences_->variants) + (sequences_->identical != 0);
}
void diffDNA::visitPatches(const std::function<void(const std::string&, uint64_t, const SampleCounts*, const QualityRecord*)>& visitor) const {
    if (!sequences_) { return; }
    if (sequences_->identical) { visitor({}, sequences_->identical, nullptr, nullptr); }
    std::visit([&](const auto& variants) {
        for (const auto& entry : variants) {
            visitor(entry.first, SequenceData::count(entry.second), SequenceData::samples(entry.second), SequenceData::quality(entry.second));
        }
    }, sequences_->variants);
}
size_t diffDNA::payloadBytes() const {
    size_t bytes = 0;
    visitPatches([&](const auto& patch, uint64_t count, const auto*, const auto*) {
        if (patch.empty()) { return; } // Reference observations live in the header.
        const size_t countBytes = varBytes(count) + (storesPairs() ? varBytes(patch.size()) : 0);
        if (patch.size() > std::numeric_limits<size_t>::max() - countBytes ||
            bytes > std::numeric_limits<size_t>::max() - countBytes - patch.size()) { throw std::overflow_error("diffDNA payload size overflow"); }
        bytes += countBytes + patch.size();
    });
    return bytes;
}
void diffDNA::visit(const std::function<void(std::string_view, uint64_t, const SampleCounts*, const QualityRecord*)>& visitor) const {
    if (storesPairs()) { throw std::logic_error("Paired variants require a paired visitor or export"); }
    visitFragment([&](auto sequence, auto, auto copies, const auto* samples, const auto* quality) { visitor(sequence, copies, samples, quality); });
}
void diffDNA::forEachPair(const std::function<void(std::string_view, std::string_view, uint64_t, const SampleCounts*, const QualityRecord*)>& visitor) const {
    if (!storesPairs()) { throw std::logic_error("diffDNA pair storage was not enabled"); }
    visitFragment(visitor);
}
void diffDNA::visitFragment(const std::function<void(std::string_view, std::string_view, uint64_t, const SampleCounts*, const QualityRecord*)>& visitor) const {
    visitPatches([&](const std::string& patch, uint64_t copies, const SampleCounts* samples, const QualityRecord* quality) {
        if (patch.empty()) { visitor(reference(), storesPairs() ? mateReference() : std::string_view{}, copies, samples, quality); return; }
        if (storesPairs()) {
            const auto parts = splitPairPatch(patch);
            const auto r1 = decodePatch(reference(), parts.first), r2 = decodePatch(mateReference(), parts.second);
            visitor(r1, r2, copies, samples, quality);
        } else {
            const auto sequence = decodePatch(reference(), patch);
            visitor(sequence, {}, copies, samples, quality);
        }
    });
}
void diffDNA::forEach(const std::function<void(std::string_view, uint64_t)>& visitor) const {
    visit([&](std::string_view sequence, uint64_t count, const SampleCounts*, const QualityRecord*) { visitor(sequence, count); });
}
void diffDNA::forEachWithSamples(const std::function<void(std::string_view, uint64_t, const SampleCounts&)>& visitor) const {
    if (!tracksSamples()) { throw std::logic_error("diffDNA sample provenance was not retained"); }
    visit([&](std::string_view sequence, uint64_t count, const SampleCounts* samples, const QualityRecord*) { visitor(sequence, count, *samples); });
}
namespace {
std::string subclusterLabel(std::string_view cluster) {
    std::string label(cluster.substr(0, cluster.find_first_of(" \t\r\n")));
    const auto size = label.find(";size=");
    if (size != std::string::npos) {
        const auto end = label.find(';', size + 6);
        label.erase(size, end == std::string::npos ? std::string::npos : end - size + 1);
    }
    if (label.empty()) { label = "cluster"; }
    return label;
}
void writeQualityRecord(std::ostream& out, std::string_view sequence, const diffDNA::QualityRead& quality, const std::string& label, int offset) {
    if (quality.qualities.size() != sequence.size()) { throw std::logic_error("Incomplete diffDNA qualities"); }
    DNA expanded(std::string(sequence), label);
    expanded.setQual(quality.qualities); expanded.prepareWrite(offset); expanded.writeFastQ(out, false);
}
}
void diffDNA::forEachWithQualities(const std::function<void(std::string_view, uint64_t, const SampleCounts*, const QualityRecord&)>& visitor) const {
    if (!storesQualities()) { throw std::logic_error("diffDNA qualities were not retained"); }
    visit([&](auto sequence, auto count, const SampleCounts* samples, const QualityRecord* quality) {
        visitor(sequence, count, samples, *quality);
    });
}
void diffDNA::writeFastq(std::ostream& out, std::string_view cluster, int qualityOffset) const {
    if (qualityOffset < 33 || qualityOffset > 126) { throw std::invalid_argument("Invalid FASTQ quality offset"); }
    const auto label = subclusterLabel(cluster);
    size_t subcluster = 0;
    forEachWithQualities([&](auto sequence, auto count, const SampleCounts*, const QualityRecord& quality) {
        if (!quality.seed) { throw std::logic_error("Missing diffDNA qualities"); }
        writeQualityRecord(out, sequence, *quality.seed, label + ".sub" + std::to_string(++subcluster) + ";size=" + std::to_string(count) + ";", qualityOffset);
    });
    if (!out) { throw std::runtime_error("Could not write diffDNA subcluster FASTQ"); }
}
void diffDNA::writeFastq(std::ostream& r1, std::ostream& r2, std::string_view cluster, int qualityOffset) const {
    if (!storesQualities()) { throw std::logic_error("diffDNA qualities were not retained"); }
    if (qualityOffset < 33 || qualityOffset > 126) { throw std::invalid_argument("Invalid FASTQ quality offset"); }
    const auto label = subclusterLabel(cluster); size_t subcluster = 0;
    forEachPair([&](auto first, auto second, auto count, const auto*, const QualityRecord* quality) {
        if (!quality || !quality->seed || !quality->seed->mate) { throw std::logic_error("Missing paired diffDNA qualities"); }
        const auto name = label + ".sub" + std::to_string(++subcluster) + ";size=" + std::to_string(count) + ";";
        writeQualityRecord(r1, first, *quality->seed, name, qualityOffset);
        writeQualityRecord(r2, second, *quality->seed->mate, name, qualityOffset);
    });
    if (!r1 || !r2) { throw std::runtime_error("Could not write paired diffDNA FASTQ"); }
}
void diffDNA::writeFasta(std::ostream& r1, std::ostream& r2, std::string_view cluster) const {
    const auto label = subclusterLabel(cluster); size_t subcluster = 0;
    forEachPair([&](auto first, auto second, auto count, const auto*, const auto*) {
        const auto name = label + ".sub" + std::to_string(++subcluster) + ";size=" + std::to_string(count) + ";";
        r1 << '>' << name << '\n' << first << '\n'; r2 << '>' << name << '\n' << second << '\n';
    });
    if (!r1 || !r2) { throw std::runtime_error("Could not write paired diffDNA FASTA"); }
}
void diffDNA::writeFasta(std::ostream& out, std::string_view cluster) const {
    if (!sequences_) { throw std::logic_error("Counts-only diffDNA has no nucleotide sequences to export"); }
    const auto label = subclusterLabel(cluster);
    size_t subcluster = 0;
    forEach([&](std::string_view sequence, uint64_t count) {
        out << '>' << label << ".sub" << ++subcluster << ";size=" << count << ";\n";
        out.write(sequence.data(), sequence.size()); out.put('\n');
    });
    if (!out) { throw std::runtime_error("Could not write diffDNA subcluster FASTA"); }
}
void diffDNA::write(std::ostream& out, std::string_view cluster) const {
    const auto& ref = reference(); // Reject counts-only objects before writing a header.
    out.write(storesPairs() ? "SDMDIFF3" : "SDMDIFF2", 8);
    auto put = [&](uint8_t byte) { out.put(static_cast<char>(byte)); };
    putVar(cluster.size(), put); out.write(cluster.data(), cluster.size());
    putVar(ref.size(), put); out.write(ref.data(), ref.size());
    if (storesPairs()) { const auto& mate = mateReference(); putVar(mate.size(), put); out.write(mate.data(), mate.size()); }
    putVar(reads_, put); putVar(referenceCount(), put); putVar(payloadBytes(), put);
    visitPatches([&](const auto& patch, uint64_t count, const auto*, const auto*) {
        if (!patch.empty()) { putVar(count, put); if (storesPairs()) { putVar(patch.size(), put); } out.write(patch.data(), patch.size()); }
    });
    if (!out) { throw std::runtime_error("Could not write diffDNA output"); }
}
std::shared_ptr<diffDNA> diffDNA::read(std::istream& in, std::string& cluster) {
    if (in.peek() == std::char_traits<char>::eof()) {
        if (in.bad() || !in.eof()) { throw std::runtime_error("Could not read diffDNA input"); }
        return nullptr;
    }
    char magic[8]; in.read(magic, 8);
    const bool version1 = in && std::string_view(magic, 8) == "SDMDIFF1";
    const bool version3 = in && std::string_view(magic, 8) == "SDMDIFF3";
    if (!in || (!version1 && !version3 && std::string_view(magic, 8) != "SDMDIFF2")) { throw std::runtime_error("Invalid diffDNA header"); }
    auto get = [&]() -> uint8_t {
        const int byte = in.get();
        if (byte == std::char_traits<char>::eof()) { throw std::runtime_error("Truncated diffDNA input"); }
        return static_cast<uint8_t>(byte);
    };
    auto readString = [&]() {
        const size_t length = asSize(getVar(get));
        std::string result;
        for (size_t remaining = length; remaining;) {
            char chunk[4096]; const size_t n = std::min(remaining, sizeof(chunk));
            in.read(chunk, n);
            if (!in) { throw std::runtime_error("Truncated diffDNA string"); }
            result.append(chunk, n); remaining -= n;
        }
        return result;
    };
    cluster = readString();
    const auto ref = std::make_shared<const std::string>(readString());
    auto result = std::make_shared<diffDNA>(ref, false, false, version3 ? std::make_shared<const std::string>(readString()) : nullptr);
    const uint64_t reads = getVar(get), identical = getVar(get);
    if (identical > reads) { throw std::runtime_error("Invalid diffDNA counts"); }
    if (version3) { result->addPair(result->reference(), result->mateReference(), identical); }
    else { result->add(result->reference(), identical); }
    size_t remaining = asSize(getVar(get));
    auto payloadGet = [&]() -> uint8_t {
        if (!remaining) { throw std::runtime_error("Truncated diffDNA payload"); }
        --remaining; return get();
    };
    while (remaining) {
        const uint64_t copies = version1 ? 1 : getVar(payloadGet);
        if (!copies || copies > reads - result->readCount()) { throw std::runtime_error("Invalid diffDNA subcluster count"); }
        if (version3) {
            const size_t length = asSize(getVar(payloadGet));
            if (!length || length > remaining) { throw std::runtime_error("Invalid diffDNA pair payload length"); }
            std::string patch;
            for (size_t i = 0; i < length; ++i) { patch.push_back(static_cast<char>(payloadGet())); }
            const auto parts = splitPairPatch(patch);
            const auto r1 = decodePatch(result->reference(), parts.first), r2 = decodePatch(result->mateReference(), parts.second);
            result->addPair(r1, r2, copies);
        } else { result->add(decode(result->reference(), payloadGet), copies); }
    }
    if (result->readCount() != reads) { throw std::runtime_error("Inconsistent diffDNA observation count"); }
    return result;
}
