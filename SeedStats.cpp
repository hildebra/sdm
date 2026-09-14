// SPDX-License-Identifier: GPL-3.0-or-later
#include "SeedStats.h"
#include "DNA.h"
#include <limits>
#include <stdexcept>
SeedReadStats::SeedReadStats(DNA& d)
    : errors(d.getAccumError()), logical(d.length()), backing(d.mem_length()),
      average(d.getAvgQual()), identity(d.getTempFloat()), ambiguous(d.numNonCanonicalDNA(true)),
      mergeLength(d.getMergeLength()), mergeErrors(d.getMergeErrors()), overlapQuality(d.getMergeErrorsQual()),
      forward(d.getFwdPrimDetect()), reverse(d.getRevPrimDetect()), passing(d.isGreenQual()) {}
void SeedComparisonHistory::recordWhoIsBetterCompare(int candidateLength, int referenceLength) {
    for (int length : {candidateLength, referenceLength}) {
        if (length >= 0) { ++merged_compares; accumulated_merge_length += static_cast<uint64_t>(length); }
        else { ++non_merged_compares; }
    }
}
void SeedComparisonHistory::merge(const SeedComparisonHistory& other) {
    const auto limit = std::numeric_limits<uint64_t>::max();
    if (other.merged_compares > limit - merged_compares || other.non_merged_compares > limit - non_merged_compares ||
        other.accumulated_merge_length > limit - accumulated_merge_length) { throw std::overflow_error("Seed comparison history overflow"); }
    merged_compares += other.merged_compares; non_merged_compares += other.non_merged_compares;
    accumulated_merge_length += other.accumulated_merge_length;
}
