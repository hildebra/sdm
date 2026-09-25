#pragma once
// Shared modern selector. Borrowed read views stay alive throughout the call.
// SPDX-License-Identifier: GPL-3.0-or-later
#include "DNAconsts.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
double pointFiveRatio(double, double, double);
namespace seed_selection {
template<class Read, class History>
bool better(Read* d1, Read* d2, Read* dM, Read* r1, Read* r2, Read* rM,
    float& ever_best, bool forSeed, History* merge_stats_owner, uint64_t dSiz = 2, uint64_t rSiz = 2) {
    using namespace std;
if (d1 == nullptr || r1 == nullptr) {
		return false;
	}
	if (d2 == nullptr) {
		if (d1->has2PrimersDetected() && !r1->has2PrimersDetected()) { return true; }
		if (!d1->has2PrimersDetected() && r1->has2PrimersDetected()) { return false; }
    } else {
	}
	if (d1->getFwdPrimDetect() && !r1->getFwdPrimDetect()) { return true; }
	if (!d1->getFwdPrimDetect() && r1->getFwdPrimDetect()) { return false; }
	double d1pid(1.), refpid(1.);
	if (ever_best >=0) {
		d1pid = (double) d1->getTempFloat(); refpid = (double)r1->getTempFloat();
		if (d1pid > ever_best) {
			ever_best = (float) d1pid;
		}
		if (d1pid < (refpid - 0.3f) || d1pid < (ever_best - 0.5f)) { return false; }
	}
	double ratId = 1.;
	if (refpid > 0.) {
		if (forSeed) {
			ratId = (double(d1pid) / double(refpid)) * 50 - 49.0;
		}
		else {
			ratId = pointFiveRatio(double(d1pid), (double)refpid, 2.0);
		}
	}
	bool dMerge(true), rMerg(true);
	double curL = (double)d1->getMergeLength();
	if (curL < 0) {
		curL = (double) d1->length();
		dMerge = false;
	}
	double refL = (double)r1->getMergeLength();
	if (refL < 0) {
		refL = (double) r1->length();
		rMerg = false;
	}
    if (merge_stats_owner != nullptr) {
		merge_stats_owner->recordWhoIsBetterCompare(d1->getMergeLength(), r1->getMergeLength());
	}
    double mergedFraction = 0.0;
	double avgMergeLen = 0.0;
	if (merge_stats_owner != nullptr) {
		const auto& ms = merge_stats_owner->getWhoIsBetterMergeStats();
		const uint64_t totalCompared = ms.merged_compares + ms.non_merged_compares;
		if (totalCompared > 0) {
			mergedFraction = static_cast<double>(ms.merged_compares) / static_cast<double>(totalCompared);
		}
		avgMergeLen = ms.merged_compares > 0 ? static_cast<double>(ms.accumulated_merge_length) / static_cast<double>(ms.merged_compares) : 0.0;
	}
	bool mergAdv = dMerge && !rMerg;
   if (refL <= 0.) { return false; }
	if ((curL) / (refL) < BestLengthRatio && !mergAdv) { return false; }
	int d1Ncnt = d1->numNonCanonicalDNA(true); int r1Ncnt = r1->numNonCanonicalDNA(true);
	float CanoScore = 1.f;
	if (d1Ncnt > 0 || r1Ncnt > 0) {
		CanoScore = 1.f - (float)(d1Ncnt - r1Ncnt) / (float)(d1Ncnt + r1Ncnt + 10);
	}
	// Overlap errors of the two consensus reads (mismatches in the overlap plus
	// their mean quality), normalized so that the cleaner merge scores 1.
	// overlapRatio > 1 favors the candidate; 1 unless both reads merged.
	float dmergErrSco = 0.f; 	float refMergErrSco = 0.f;
	double overlapRatio = 1.;
	if (r1->getMergeLength() >=0 && d1->getMergeLength() >=0) {
		if (d1->getMergeErrors() > 0) {
			dmergErrSco = (float)d1->getMergeErrors();
			dmergErrSco += d1->getMergeErrorsQual()/30;
		}
		if (r1->getMergeErrors() > 0) {
			refMergErrSco = (float)r1->getMergeErrors() + r1->getMergeErrorsQual() / 30;
		}
        float maxMerr = max(refMergErrSco, dmergErrSco);
		if (maxMerr > 0.f) {
			dmergErrSco /= maxMerr; 	refMergErrSco /= maxMerr;
			dmergErrSco = 1.f - dmergErrSco; refMergErrSco = 1.f - refMergErrSco;
			overlapRatio = pointFiveRatio(dmergErrSco, refMergErrSco, 0.5);
		}
	}
	Read* dx = nullptr; Read* rx = nullptr;
	bool allowMergeGuide = dMerge && rMerg;
	if (allowMergeGuide && dM != nullptr) { dx = dM; 	d2 = nullptr; } 	else { dx = d1; }
	if (allowMergeGuide && rM != nullptr) { rx = rM; r2 = nullptr;}else { rx = r1; }
	// Expected errors per base: the sum alone penalized longer reads.
	double dAccuE = dx->getAccumError() / std::max(1.0, (double)dx->length());
	double rAccuE = rx->getAccumError() / std::max(1.0, (double)rx->length());
	double AccErrRatio = pointFiveRatio(rAccuE, dAccuE,0.25);
	double dQ=dx->getAvgQual(); double rQ=rx->getAvgQual();
	if (d2 != nullptr) { dQ += d2->getAvgQual();  }
	if (r2 != nullptr) { rQ += r2->getAvgQual();  }
	double qualRatio = pointFiveRatio(dQ, rQ, 0.5); 
	double ratLength = pointFiveRatio(double(curL), (double)refL, 0.5); 
	double mergePreference = 1.0;
	if (dMerge != rMerg) {
		if (forSeed) {
			mergePreference = dMerge ? 1.35 : 0.75;
		}
		else if (dMerge) {
			mergePreference = 1.0 + mergedFraction;
		}
		else {
			mergePreference = 1.0 - (0.5 * mergedFraction);
		}
	}
	double mergeLenDevRatio = 1.;
	if (!forSeed && dMerge && rMerg) {
		double dMLDIFF = (double) abs(avgMergeLen - d1->getMergeLength());
		double rMLDIFF = (double) abs(avgMergeLen - r1->getMergeLength());
		mergeLenDevRatio = (rMLDIFF + std::max(1.0, avgMergeLen)) / (dMLDIFF + std::max(1.0, avgMergeLen));
	}
	double logsizeRatio = pointFiveRatio(double(dSiz), double(rSiz),0.5); 
	double thresh(1.01f);
	auto clampScore = [](double v) {
		if (!std::isfinite(v)) {
			return 1.0;
		}
		if (v < 0.2) {
			return 0.2;
		}
		if (v > 2.5) {
			return 2.5;
		}
		return v;
	};
	// Weights sum to 1.28 (seed) and 1.40 (dereplication), as before the
	// overlap-error term was added; its weight came off wAccErr/wMergeLenDev.
	double wCano, wAccErr, wQual, wLen, wId, wMergePref, wSize, wMergeLenDev, wOverlap;
	if (forSeed) {
		wCano = 0.08;
		wAccErr = 0.34;
		wQual = 0.14;
		wLen = 0.08;
		wId = 0.42;
		wMergePref = 0.10;
		wSize = 0.02;
		wMergeLenDev = 0.00;
		wOverlap = 0.10;
	}
	else {
		wCano = 0.08;
		wAccErr = 0.28;
		wQual = 0.16;
		wLen = 0.08;
		wId = 0.0;
		wMergePref = 0.34;
		wSize = 0.04;
		wMergeLenDev = 0.30;
		wOverlap = 0.12;
	}
	double cumulativeQualityScore =
		std::pow(clampScore(CanoScore), wCano) *
		std::pow(clampScore(AccErrRatio), wAccErr) *
		std::pow(clampScore(qualRatio), wQual) *
		std::pow(clampScore(ratLength), wLen) *
		std::pow(clampScore(ratId), wId) *
		std::pow(clampScore(mergePreference), wMergePref) *
		std::pow(clampScore(logsizeRatio), wSize) *
		std::pow(clampScore(mergeLenDevRatio), wMergeLenDev) *
		std::pow(clampScore(overlapRatio), wOverlap);
	if (cumulativeQualityScore > thresh) {
		return true;
	}
	return false;
}
}
