// SPDX-License-Identifier: GPL-3.0-or-later
#include "CompactDNA.h"
#include "SeedSelection.h"
#include <cstring>
#include <bit>
#include <stdexcept>
#include <type_traits>

CompactDNAReferences::CompactDNAReferences(DNA& r1, DNA* r2, DNA* consensus)
    : first(std::make_shared<const std::string>(r1.getSequence())),
      mate(std::make_shared<const std::string>(r2 ? r2->getSequence() : std::string{})),
      merged(consensus ? std::make_shared<const std::string>(consensus->getSequence()) : first) {}

// Private process-local serialization, never a file format. Preserve the exact
// scalar types, including lazy caches/flags; arrays use raw, run-length or base-offset
// bit packing, choosing the smallest representation. Integers in runs are varints.
struct CompactDNA::Writer {
    std::string bytes;
    template<class T> void pod(const T& value) {
        static_assert(std::is_trivially_copyable_v<T>);
        bytes.append(reinterpret_cast<const char*>(&value), sizeof(T));
    }
    void var(uint64_t value) {
        while (value >= 128) { bytes.push_back(static_cast<char>(value | 128)); value >>= 7; }
        bytes.push_back(static_cast<char>(value));
    }
    void signedVar(int64_t value) { var((uint64_t(value) << 1) ^ uint64_t(-(value < 0))); }
    void string(std::string_view value) { var(value.size()); bytes.append(value); }
    template<class T> void array(const std::vector<T>& values) {
        Writer runs;
        for (size_t i = 0; i < values.size();) {
            size_t j = i + 1; while (j < values.size() && values[j] == values[i]) { ++j; }
            runs.var(j - i); runs.var(static_cast<std::make_unsigned_t<T>>(values[i])); i = j;
        }
        Writer packed;
        using Unsigned = std::make_unsigned_t<T>;
        uint64_t minimum = 0, maximum = 0;
        if (!values.empty()) {
            minimum = maximum = static_cast<Unsigned>(values.front());
            for (const auto value : values) {
                minimum = std::min(minimum, uint64_t(static_cast<Unsigned>(value)));
                maximum = std::max(maximum, uint64_t(static_cast<Unsigned>(value)));
            }
        }
        const unsigned width = std::bit_width(maximum - minimum);
        packed.var(minimum); packed.pod(static_cast<uint8_t>(width));
        uint8_t byte = 0; unsigned used = 0;
        for (const auto value : values) {
            uint64_t difference = static_cast<Unsigned>(value) - minimum;
            for (unsigned remaining = width; remaining;) {
                const unsigned take = std::min(8 - used, remaining);
                byte |= uint8_t(difference & ((1u << take) - 1)) << used;
                difference >>= take; used += take; remaining -= take;
                if (used == 8) { packed.pod(byte); byte = 0; used = 0; }
            }
        }
        if (used) { packed.pod(byte); }
        var(values.size());
        const size_t rawSize = values.size() * sizeof(T);
        const uint8_t encoding = packed.bytes.size() < rawSize && packed.bytes.size() < runs.bytes.size() ? 2 :
            runs.bytes.size() < rawSize ? 1 : 0;
        pod(encoding);
        if (encoding == 2) { bytes += packed.bytes; }
        else if (encoding == 1) { bytes += runs.bytes; }
        else if (!values.empty()) { bytes.append(reinterpret_cast<const char*>(values.data()), rawSize); }
    }
};
struct CompactDNA::Reader {
    std::string_view bytes;
    template<class T> T pod() {
        static_assert(std::is_trivially_copyable_v<T>);
        if (bytes.size() < sizeof(T)) { throw std::logic_error("Truncated compact DNA state"); }
        T result; std::memcpy(&result, bytes.data(), sizeof(T)); bytes.remove_prefix(sizeof(T)); return result;
    }
    uint64_t var() {
        uint64_t result = 0;
        for (unsigned shift = 0; shift < 64; shift += 7) {
            const auto byte = pod<uint8_t>();
            if (shift == 63 && byte > 1) { break; }
            result |= uint64_t(byte & 127) << shift;
            if (!(byte & 128)) { return result; }
        }
        throw std::logic_error("Invalid compact DNA integer");
    }
    int64_t signedVar() { const auto value = var(); return int64_t(value >> 1) ^ -int64_t(value & 1); }
    std::string_view string() {
        const auto size = var();
        if (size > bytes.size()) { throw std::logic_error("Invalid compact DNA string"); }
        const auto value = bytes.substr(0, size); bytes.remove_prefix(size); return value;
    }
    template<class T> void array(std::vector<T>& values) {
        const auto size = var(); const auto encoding = pod<uint8_t>();
        values.clear();
        if (encoding == 2) {
            const auto minimum = var(); const unsigned width = pod<uint8_t>();
            if (width > sizeof(T) * 8 || (width && size > uint64_t(bytes.size()) * 8 / width)) {
                throw std::logic_error("Invalid compact DNA packed array");
            }
            uint8_t byte = 0; unsigned available = 0;
            values.reserve(size);
            for (uint64_t i = 0; i < size; ++i) {
                uint64_t difference = 0;
                for (unsigned shift = 0; shift < width;) {
                    if (!available) { byte = pod<uint8_t>(); available = 8; }
                    const unsigned take = std::min(available, width - shift);
                    difference |= uint64_t(byte & ((1u << take) - 1)) << shift;
                    byte >>= take; available -= take; shift += take;
                }
                values.push_back(static_cast<T>(minimum + difference));
            }
        } else if (encoding == 1) {
            while (values.size() < size) {
                const auto count = var(); const auto value = var();
                if (!count || count > size - values.size()) { throw std::logic_error("Invalid compact DNA run"); }
                values.insert(values.end(), count, static_cast<T>(value));
            }
        } else {
            if (encoding != 0) { throw std::logic_error("Invalid compact DNA array encoding"); }
            if (size > bytes.size() / sizeof(T)) { throw std::logic_error("Invalid compact DNA array"); }
            values.resize(size);
            if (size) { std::memcpy(values.data(), bytes.data(), size * sizeof(T)); }
            bytes.remove_prefix(size * sizeof(T));
        }
    }
};
void CompactDNA::packEvidence(Writer& out, const DNAQualityEvidence& d) {
    out.array(d.quality_sum_per_base_); out.array(d.quality_coverage_per_base_); out.var(d.uniform_quality_coverage_);
}
void CompactDNA::unpackEvidence(Reader& in, DNAQualityEvidence& d) {
    in.array(d.quality_sum_per_base_); in.array(d.quality_coverage_per_base_); d.uniform_quality_coverage_ = in.var();
}
void CompactDNA::packRead(Writer& out, const DerepRead& d, const CompactDNAReferences& refs, int level) {
    const auto& ref = level == 1 ? refs.mate : level == 2 ? refs.merged : refs.first;
    out.string(diffDNA::encodeDifference(*ref, d.sequence_));
    out.string(d.sample_barcode_.str()); out.string(d.id_); out.string(d.new_id_); out.array(d.qual_);
    out.var(d.sequence_length_); out.var(d.quality_sum_); out.pod(d.accumulated_error_);
    const auto& f = d.FtsDetected;
    const uint8_t detected = f.forward | (f.reverse << 1) | (f.TA_cut << 2) | (f.barcode_detected << 3) |
        (f.barcode_cut << 4) | (f.dereplicated << 5) | (f.revPairConstellation << 6);
    out.pod(detected); out.pod(f.meanQInOverlapMismatch); out.signedVar(f.errInOverlap); out.signedVar(f.mergeLength);
    out.pod(d.QualCtrl);
    out.signedVar(d.merge_seed_pos_); out.signedVar(d.seed_length_); out.signedVar(d.merge_offset_); out.signedVar(d.sample_id_);
    out.pod(d.avg_qual_); out.pod(d.tempFloat); out.pod(d.read_position_); out.pod(d.fastq_write_offset_);
    const uint8_t flags = d.good_quality_ | (d.mid_quality_ << 1) | (d.failed_ << 2) | (d.reversed_ << 3) |
        (d.id_fixed_ << 4) | (d.new_id_is_original_ << 5) | (d.reversed_merge_ << 6);
    out.pod(flags); out.pod(d.derep_fastq_offset_);
    packEvidence(out, d);
}
void CompactDNA::packUnique(Writer& out, const DNAunique& d, const CompactDNAReferences& refs) {
    // The top-level diff object owns the exact parent's authoritative counts;
    // the selection history is stored separately in history_.
    packRead(out, d, refs, 0);
    out.pod(bool(d.pair_)); if (d.pair_) { packRead(out, *d.pair_, refs, 1); }
    out.pod(bool(d.merge_)); if (d.merge_) { packRead(out, *d.merge_, refs, 2); }
    out.pod(bool(d.merge_quality_evidence_));
    if (d.merge_quality_evidence_) { out.array(d.merge_quality_evidence_->raw); packEvidence(out, *d.merge_quality_evidence_); }
}
void CompactDNA::unpackRead(Reader& in, DerepRead& d, const CompactDNAReferences& refs, int level) {
    const auto& ref = level == 1 ? refs.mate : level == 2 ? refs.merged : refs.first;
    d.sequence_ = diffDNA::decodeDifference(*ref, in.string());
    d.sample_barcode_ = std::string(in.string()); d.id_ = in.string(); d.new_id_ = in.string(); in.array(d.qual_);
    d.sequence_length_ = in.var(); d.quality_sum_ = in.var(); d.accumulated_error_ = in.pod<double>();
    const auto detected = in.pod<uint8_t>(); auto& f = d.FtsDetected;
    f.forward = detected & 1; f.reverse = (detected >> 1) & 1; f.TA_cut = (detected >> 2) & 1;
    f.barcode_detected = (detected >> 3) & 1; f.barcode_cut = (detected >> 4) & 1;
    f.dereplicated = (detected >> 5) & 1; f.revPairConstellation = (detected >> 6) & 1;
    f.meanQInOverlapMismatch = in.pod<qual_score>(); f.errInOverlap = in.signedVar(); f.mergeLength = in.signedVar();
    d.QualCtrl = in.pod<DNA::QualStats>();
    d.merge_seed_pos_ = in.signedVar(); d.seed_length_ = in.signedVar(); d.merge_offset_ = in.signedVar(); d.sample_id_ = in.signedVar();
    d.avg_qual_ = in.pod<float>(); d.tempFloat = in.pod<float>(); d.read_position_ = in.pod<short>(); d.fastq_write_offset_ = in.pod<qual_score>();
    const auto flags = in.pod<uint8_t>();
    d.good_quality_ = flags & 1; d.mid_quality_ = (flags >> 1) & 1; d.failed_ = (flags >> 2) & 1; d.reversed_ = (flags >> 3) & 1;
    d.id_fixed_ = flags & 16; d.new_id_is_original_ = flags & 32; d.reversed_merge_ = flags & 64;
    d.derep_fastq_offset_ = in.pod<qual_score>();
    unpackEvidence(in, d);
}
std::shared_ptr<DNAunique> CompactDNA::unpackUnique(Reader& in, const CompactDNAReferences& refs) {
    auto result = std::make_shared<DNAunique>(); auto& d = *result;
    unpackRead(in, d, refs, 0);
    if (in.pod<bool>()) { d.pair_ = std::make_shared<DerepRead>(); unpackRead(in, *d.pair_, refs, 1); }
    if (in.pod<bool>()) { d.merge_ = std::make_shared<DerepRead>(); unpackRead(in, *d.merge_, refs, 2); }
    if (in.pod<bool>()) {
        d.merge_quality_evidence_ = std::make_shared<DNAunique::MergeQualityEvidence>();
        in.array(d.merge_quality_evidence_->raw); unpackEvidence(in, *d.merge_quality_evidence_);
    }
    return result;
}
CompactDNA::CompactDNA(const CompactDNAReferences& refs) : variants_(refs.first, false, false, refs.mate, true) {}
void CompactDNA::store(DNAunique& dna, const CompactDNAReferences& refs) {
    firstStats_.emplace(dna);
    mateStats_.reset(); if (dna.pair_) { mateStats_.emplace(*dna.pair_); }
    mergedStats_.reset(); if (dna.merge_) { mergedStats_.emplace(*dna.merge_); }
    history_ = dna.who_is_better_merge_stats_;
    Writer out; packUnique(out, dna, refs); selected_ = std::move(out.bytes); selected_.shrink_to_fit();
}
std::shared_ptr<DNAunique> CompactDNA::materialize(const CompactDNAReferences& refs) const {
    Reader in{selected_}; auto result = unpackUnique(in, refs);
    if (!in.bytes.empty()) { throw std::logic_error("Trailing compact DNA state"); }
    result->differences_.merge(variants_); // Counts-only snapshot for ordinary exporters.
    result->who_is_better_merge_stats_ = history_;
    return result;
}
void CompactDNA::observe(const std::shared_ptr<DNA>& dna, const std::shared_ptr<DNA>& mate,
    const std::shared_ptr<DNA>& merged, int sample, bool asFasta, bool collectMergedQualities,
    bool markDereplicated, const CompactDNAReferences& refs) {
    if (selected_.empty()) {
        DNAuniqSet ordinary;
        ordinary.addObservation(dna, mate, merged, 0, sample, asFasta, collectMergedQualities, markDereplicated);
        store(*ordinary[0], refs);
    } else if (!asFasta) {
        auto selected = materialize(refs);
        selected->matchedDNA(dna, mate, merged, -1, asFasta, collectMergedQualities, markDereplicated);
        store(*selected, refs);
    } else {
        // Exact repeats normally need only cached selection statistics. Use the
        // ordinary selector/history, materializing the pair only on replacement.
        if (markDereplicated) { dna->setDereplicated(); }
        if (dna->isGreenQual()) {
            const SeedReadStats incoming(*dna);
            std::optional<SeedReadStats> mateStats, mergedStats;
            if (mate) { mateStats.emplace(*mate); }
            if (merged) { mergedStats.emplace(*merged); }
            float best = -1.f;
            if (seed_selection::better<const SeedReadStats>(&incoming, mateStats ? &*mateStats : nullptr,
                mergedStats ? &*mergedStats : nullptr, &*firstStats_, mateStats_ ? &*mateStats_ : nullptr,
                mergedStats_ ? &*mergedStats_ : nullptr, best, false, &history_)) {
                auto selected = materialize(refs);
                selected->takeOverDNA(dna, mate, merged, true);
                store(*selected, refs);
            }
        }
    }
    variants_.addPair(dna->getSequence(), mate ? std::string_view(mate->getSequence()) : std::string_view{}, 1, sample);
}
void CompactDerepSet::addObservation(std::shared_ptr<DNA> dna, std::shared_ptr<DNA> mate, std::shared_ptr<DNA> merged,
    int mergePosition, int sample, bool asFasta, bool collectMergedQualities, bool markDereplicated) {
    auto [it, inserted] = subgroups_.try_emplace(mergePosition, *references_);
    it->second.observe(dna, mate, merged, sample, asFasta, collectMergedQualities, markDereplicated, *references_);
}
uint64_t CompactDerepSet::totalSum() const {
    uint64_t result = 0;
    for (const auto& [position, subgroup] : subgroups_) { DNAQualityEvidence::addToQualSum(result, subgroup.totalSum()); }
    return result;
}
std::shared_ptr<DNAunique> CompactDerepSet::best(bool addCounts, bool collectQualities) const {
    // Reuse ordinary subgroup weighting, tie handling and evidence aggregation.
    // Only this exact key is expanded; output snapshots never expand the index.
    DNAuniqSet ordinary;
    for (const auto& [position, subgroup] : subgroups_) { ordinary[position] = subgroup.materialize(*references_); }
    return ordinary.best(addCounts, collectQualities);
}
void CompactDerepSet::addEvidenceTo(const std::shared_ptr<DNAunique>& target, bool collectQualities) const {
    for (const auto& [position, subgroup] : subgroups_) { target->absorbDerepEvidence(subgroup.materialize(*references_), collectQualities); }
}
std::vector<diffDNA> CompactDerepSet::variants() const {
    std::vector<diffDNA> result;
    for (const auto& [position, subgroup] : subgroups_) { result.push_back(subgroup.variants()); }
    return result;
}
