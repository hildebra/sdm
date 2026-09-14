# Unified dereplication sample counts

Measured 2026-09-13, comparing the saved implementation immediately before the sample-count refactor with the final implementation. Every `DNAunique` now embeds `diffDNA` as its count store; optional nucleotide variants live behind a shared pointer with copy on write. Ordinary exact/prefix dereplication retains no nucleotide reference or variant table in `diffDNA`. For final reassignment, one map of variant-to-sparse-sample-counts replaces the previous separate variant-count and variant-sample maps. Cluster sample counts are derived lazily from those entries.

## Method

The existing production benchmark generated one million 200 bp Q30 reads, one at a time, from 50 artificial templates in five lineages, using one sample and one processing thread. The observed substitution-error fraction was 0.001001. It is a rough substitution simulation, not an ONT indel or biological-accuracy evaluation. Both builds used GCC 13.3.0, `-O3`, C++20, `WITH_HTS=0`, `STATIC=0`, and `derep_reserve=1`.

Three runs per case and implementation started fresh processes. Before/after runs alternated within each case, reversing order in the second repetition. No builds or other tests from this task ran during measurement; unrelated host load was uncontrolled. Timing includes generation, DNA construction, production grouping, final reassignment when enabled, and normal FASTA/HQ/map output plus enabled binary differences. It excludes FASTQ parsing, upstream filtering, exact-subcluster FASTA export, and cleanup. Memory is total-process Linux `VmHWM`, not estimated live-object storage.

## Results

Values below are medians of three runs.

| Mode | Before time (s) | After time (s) | Before peak RSS (MiB) | After peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| Exact, equal lengths | 2.220 | 2.027 | 68.63 | 68.59 |
| 97%, counts only | 2.278 | 2.354 | 5.46 | 5.39 |
| 97%, differences | 2.806 | 2.833 | 8.98 | 8.96 |
| 97%, reassignment | 3.051 | 2.893 | 17.15 | 13.52 |
| 97%, reassignment + differences | 3.062 | 3.170 | 20.88 | 17.54 |

Reassignment peak RSS decreased by about **21% without difference export** and **16% with export**. Count-only and difference-only process peaks changed little. On this build, `sizeof(DNAunique)` decreased from 424 to 416 bytes and the base `sizeof(diffDNA)` from 104 to 48 bytes. These sizes exclude sparse entries, reference strings, optional variant tables, and allocator overhead; a smaller object does not necessarily change its allocator size class.

Runtime changes were small relative to the observed run-to-run variation. The final exact-mode median was lower, but an earlier comparison had the opposite sign, so these runs do not establish an exact-mode speedup. The dedicated exact counter path updates the sparse vector directly. A separate hot-counter check with seven alternating runs measured medians of 0.0447 seconds before and 0.0517 seconds after for ten million increments; this small absolute cost cannot explain the larger whole-run timing variation.

| Mode | Before time range (s) | After time range (s) |
| --- | ---: | ---: |
| Exact, equal lengths | 2.017–2.568 | 1.991–2.227 |
| 97%, counts only | 2.260–2.599 | 2.152–3.260 |
| 97%, differences | 2.620–2.925 | 2.759–3.033 |
| 97%, reassignment | 3.006–3.178 | 2.850–3.104 |
| 97%, reassignment + differences | 3.017–3.099 | 2.854–3.217 |

Every run verified that output abundances sum to one million. Before and after produced the same 44,728 exact groups or 47 coarse groups. Reassignment moved the same 333,957 reads; enabled sidecar byte sizes and simulated error fractions also matched in every paired run. The single-sample fixture does not measure scaling with many samples.

## Validation and reproduction

The final private portable build passed `TestCoarseDerep`, `TestRegression`, `TestInputHelpers`, `TestSearch`, command-line difference/subcluster exports with and without reassignment, and the ONT pipeline regressions. The dedicated coarse/diffDNA suite also passed AddressSanitizer, UndefinedBehaviorSanitizer, and leak detection. New tests compare variants and sparse samples against an independent observation tally, exercise rebasing and self-merges, reject missing provenance and count overflow, and verify exact-mode storage and retained-snapshot independence. The `SDMDIFF2` nucleotide-only format is unchanged.

Raw measurements: [unified-derep-counts.csv](benchmarks/unified-derep-counts.csv). Build `TestCoarseBenchmark` separately for the two source versions, then alternate these cases; the arguments are identity, reads, differences, length jitter, length, quality, and reassignment:

```sh
make WITH_HTS=0 STATIC=0 .build/portable-0/TestCoarseBenchmark
.build/portable-0/TestCoarseBenchmark exact 1000000 0 0 200 30 0
.build/portable-0/TestCoarseBenchmark 97 1000000 0 0 200 30 0
.build/portable-0/TestCoarseBenchmark 97 1000000 1 0 200 30 0
.build/portable-0/TestCoarseBenchmark 97 1000000 0 0 200 30 1
.build/portable-0/TestCoarseBenchmark 97 1000000 1 0 200 30 1
```
