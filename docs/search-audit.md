# Primer and barcode search audit

This audit accompanied SDM 3.51's ONT extension and used the current working checkout, including its existing uncommitted changes. It focused on search bounds, orientation, barcode ambiguity, trim coordinates and sequence/quality consistency. It was not a review of every filtering or clustering algorithm.

## Reuse and ONT changes

The existing `DNA::matchSeq` implementation is shared through a read-only sequence-view overload. ONT orientation selection uses this same substitution/IUPAC primer matcher on bounded end views. Final trimming and reversal still use `cutPrimer`, `cutPrimerRev`, `cutSeq`, `reverse_compliment` and `BCintoHead`. Existing statistics aggregation, output routing, dereplication and matrix code remain in use.

The original barcode matcher counts substitutions only. No suitable existing bounded barcode indel matcher was found. A `DNA::matchSeq_tot` overload adds full-barcode edit matching, an exact-match fast path, an error-banded dynamic program and ambiguous trim-coordinate reporting. ONT assignment decisions are localized in `Filters::detectCutBC`; they do not invoke the legacy run-wide barcode-orientation guessing/restart logic.

## Findings corrected

| Finding | Consequence before the fix | Correction and validation |
| --- | --- | --- |
| Legacy `scanBC` received the barcode error count by value. | A corrected barcode could reach the caller with an incorrect `bc_diffs` value; the final candidate's mismatch count was also not saved independently. | Pass by reference and retain the selected candidate's count. Test the observed barcode and one-error return value. |
| Rear-barcode approximate fallback called a forward-only search. | A mismatching rear barcode could be missed, or a front barcode could incorrectly satisfy a rear search. | Extend the existing substitution matcher with rear start positions and use it in the fallback. Test both rear detection and rejection of a front-only hit. |
| Mixed-length barcode maps required the read to fit the longest barcode. | A short read containing a valid shorter barcode could be rejected before searching. | Guard against the applicable minimum length and bound the longest-barcode fast path. Test a 6-base read in a 4-/8-base barcode map. |
| The preliminary read-length check always used mate 1's barcode length. | A shorter, valid mate 2 barcode could be skipped. | Check the requested mate's minimum barcode length. |
| Barcode scans used backing sequence length after a logical trim. | Bases already removed from the logical read could still contribute a barcode hit. | Use logical read length consistently, including the substitution helper and rear search. Test a barcode found only in a trimmed suffix. |
| Reverse primer matching did not validate empty patterns or explicit rear offsets before indexing. | Negative offsets could address beyond the read; small offsets could admit a partial primer. | Require a nonempty complete primer, clamp the starting offset and lower bound, and retain the existing exclusive lower-bound convention. Exercise empty, short, negative-offset, past-end and position-zero cases. |
| Several primer helpers indexed per-sample primer arrays without validating the index. | Absent or invalid primer indices could access invalid memory; detection could report a nonexistent primer. | Validate indices and empty sequences in the existing search/cut helpers. Tests cover missing indices and empty primers. |
| The first empty primer entry disabled primer processing for every sample. | Later samples with configured primers could retain primers. | Check whether any configured primer is present; test a map with an empty first primer row and a populated second row. |
| Retained reverse primers left sequence beyond the primer and were not marked detected. | Adapter/flank bases could remain despite locating the primer, and repeated processing could search it again. | Trim outside the retained primer and mark detection in the existing helpers. Test retained sequence and detection state; ONT also runs the helpers with `keepPrimerSeq=1`. |
| Retaining a barcode skipped the existing header assignment call. | `keepBarcodeSeq=1` could lose the sample prefix despite successful matching. | Separate physical cutting from `BCintoHead`. Test that both the sample prefix and original bases survive. |
| Barcode-map rebuilding retained stale cache entries. | Reused filter configurations could search stale barcode strings/indices. | Clear forward/reverse lookup caches before rebuilding. Test a changed barcode panel. |

These are targeted corrections to the shared code. Ordinary non-ONT substitution budgets, legacy search-window sizes and first-hit selection remain in place. ONT's stronger ambiguity policy is opt-in.

## Limits retained deliberately

- Primers still use substitutions and IUPAC matching, with the existing primary/alternate controls. There is no primer-indel alignment, concatemer splitting or primer alignment score model.
- Legacy exact barcode matching remains first-hit oriented. It does not gain ONT's exhaustive rejection of multiple exact sample hits in an end window.
- Legacy primer searches retain their first/last occurrence behavior and historical window conventions. Repeated primer motifs within a window are not a general structural-layout parser.
- A symmetric barcode pair does not establish biological strand. ONT uses available primer evidence and leaves primer-free reads in their input orientation.
- `keepBarcodeSeq=1` preserves barcode bases at the barcode-processing step; subsequent primer trimming can still remove sequence exterior to a primer, as in the existing pipeline.
- Existing primary/alternate filtering and primer diagnostic conventions were not redesigned. ONT assignment counts precede subsequent primer/quality rejection and should not be interpreted as final output counts.

## Regression coverage

`tests/TestSearch.cpp` includes an independent, full-matrix Levenshtein reference that checks 2,000 deterministic barcode searches for hit presence, minimum edit distance, valid coordinates and tied trim endpoints. It also directly tests the shared-code defects above, including IUPAC primer matching.

`tests/test_ont_pipeline.py` verifies assignment, physical trimming and the actual position-dependent quality strings in both orientations. It also checks full per-sample FASTQ, dereplicated abundance, HQ sequences, sample/run grouping and the existing matrix path. Both are included in `make test`. The existing regression and paired-read pipeline suites provide coverage outside ONT mode. These are synthetic correctness checks, not evidence of real ONT kit accuracy.

Validation of this change passed with the static HTSlib 1.24 build: `make test`, 19 thorough-audit pipeline cases, 34 output/audit pipeline cases, the merged paired-gzip suite and the prefix-dereplication suite. The updated search unit test and full ONT pipeline also passed with AddressSanitizer and UndefinedBehaviorSanitizer, including leak detection. The documented six-read example was executed and its FASTQ, dereplication and HQ outputs checked. A four-thread run of the 4,000-read/24-barcode fixture took about 0.2 seconds in this workspace; that timing is only a local sanity check.
