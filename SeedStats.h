#pragma once
// Compact inputs to the existing representative-selection rules.
// SPDX-License-Identifier: GPL-3.0-or-later
#include <cstdint>
#include "DNAconsts.h"
class DNA;
struct SeedComparisonHistory {
    uint64_t merged_compares = 0, non_merged_compares = 0, accumulated_merge_length = 0;
    void recordWhoIsBetterCompare(int candidateLength, int referenceLength);
    const SeedComparisonHistory& getWhoIsBetterMergeStats() const { return *this; }
    void merge(const SeedComparisonHistory& other);
};
struct SeedReadStats {
    explicit SeedReadStats(DNA& dna);
    double errors;
    uint32_t logical, backing;
    float average, identity;
    int ambiguous, mergeLength, mergeErrors;
    qual_score overlapQuality;
    bool forward, reverse, passing;
    bool has2PrimersDetected() const { return forward && reverse; }
    bool getFwdPrimDetect() const { return forward; }
    bool getRevPrimDetect() const { return reverse; }
    double getAccumError() const { return errors; }
    float getAvgQual() const { return average; }
    float getTempFloat() const { return identity; }
    uint32_t length() const { return logical; }
    uint32_t mem_length() const { return backing; }
    int numNonCanonicalDNA(bool) const { return ambiguous; }
    int getMergeLength() const { return mergeLength; }
    int getMergeErrors() const { return mergeErrors; }
    qual_score getMergeErrorsQual() const { return overlapQuality; }
};
