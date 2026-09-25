# Full-length 16S-like Q20 simulation

Measured 2026-09-12 as a follow-up to the [99% identity evaluation](99-identity-evaluation.md). This reuses the same candidate indexes, verifier, centroid ordering, and 99% threshold. It adds a rough synthetic workload; production SDM behavior is unchanged.

**Result:** packed disjoint 20-mers remain the fastest tested strategy. However, strict 99% clustering merges few independent Q20 reads, because their combined sequencing errors put most same-template pairs below 99% identity.

Follow-up: [pseudo-dereplication for lower SDM memory and reduced Savont input](pseudo-dereplication-evaluation.md) evaluates looser streaming grouping and capped retention of real reads.

## Simulation

- **5,000 reads**, generated from **50 templates**, with 100 reads per template.
- Every read is **1,500 bp**, so the average is exactly 1,500 bp.
- Each base independently receives a substitution with probability **0.01**. This continues the earlier substitutions-only comparison; insertions and deletions are outside this experiment.
- Templates have a shared backbone and nine artificial 60 bp variable blocks. Five lineages each contain ten related variants. Mean template divergence is **1.35% within a lineage** and **13.34% between lineages**.
- This conserved/variable structure is a simple performance fixture, not a biological 16S sequence or sequencer-specific error model. No reference templates are supplied to the clustering algorithm.

The generated reads contained **74,966 substitutions**, averaging **14.993 errors/read**, for an observed error rate of **0.99955%**. Individual reads had 3–31 substitutions. All 5,000 reads were exact unique sequences.

## Runtime and candidate counts

Four interleaved repetitions per method, using one thread, GCC 13.3, C++20, `-O3`, and the portable build on the same Intel Core Ultra 7 258V / WSL2 host as before. Values are medians. Timing covers the candidate index, seed extraction, lookup, and verification; it excludes simulation, exact duplicate collection, packing, and output aggregation. Exact collection plus packing separately took about **20 ms**.

| Strategy | Median time | Candidate checks | Centroid/seed postings | Clusters |
| --- | ---: | ---: | ---: | ---: |
| Disjoint 20-mers, byte verification | 0.181 s | 475,695 | 364,125 | 4,855 |
| **Disjoint 20-mers, packed verification** | **0.082 s** | **475,695** | 364,125 | 4,855 |
| Sparse 15-mers, byte verification | 0.739 s | 2,394,977 | 218,475 | 4,855 |
| Sparse 15-mers, packed verification | 0.265 s | 2,394,977 | 218,475 | 4,855 |
| Sparse 15-mers, seed rejection then byte verification | 0.811 s | 2,394,977 | 218,475 | 4,855 |
| 15-mer minimizers, `w=50`, unpositioned, packed | 2.638 s | 11,955,576 | 295,783 | 4,855 |
| 15-mer minimizers, `w=50`, positional, packed | 2.731 s | 11,955,576 | 295,783 | 4,855 |
| 15-mer minimizers, `w=66`, positional, packed | 2.354 s | 11,955,576 | 224,740 | 4,855 |
| Require first four 15-mers, then byte verification **(lossy)** | 0.267 s | 2,437,234 | 220,230 | 4,894 |

At 1,500 bp the mismatch budget is **15**, so the fixed-seed methods query the **16 rarest disjoint seeds**. The 20-mer index stores 75 words per centroid; the sparse 15-mer grid stores 45. Both enumerate the union of the selected posting lists. The minimizer methods enumerate all selected-word lists to preserve sensitivity.

Packed 20-mers were about **3.2 times faster** than packed sparse 15-mers and required **five times fewer candidate checks**, at the cost of about **67% more postings**. Postings are not total memory measurements.

The conserved backbone made every tested minimizer strategy visit the same candidate pairs as exhaustive clustering. Position information did not help because the shared sequence was already at corresponding coordinates. Thus minimizer extraction and posting-list traversal added cost without reducing verification work in this fixture. This result applies to these index strategies, not all possible uses of minimizers.

Unlike at 200 bp, the `w=66` variant has a common-window guarantee here: `1500 >= (15+1)*(15+66-1) = 1280`. Its speed disadvantage is therefore not a sensitivity problem at this length.

Individual timings were noisy, including one 20-mer packed run at 0.149 s versus 0.078–0.083 s in the other three. These numbers compare algorithms on a small synthetic workload and do not measure whole-pipeline slowdown or establish linear scaling.

## Clustering behavior at Q20

All eight methods intended to retain valid matches produced **exactly the same assignments as exhaustive greedy clustering**: **4,855 clusters from 5,000 reads**. Only 145 reads joined an earlier centroid. Individual templates were spread across **88–100 clusters**. No cross-template assignments occurred in this fixture.

The first-four mandatory-match variant changed **95 assignments** and produced 39 additional clusters. Assignment differences include downstream changes caused by altered centroid selection; they are not a direct pairwise false-negative count.

The fragmentation is expected under this error model. A Q20 read has 1% expected divergence from its error-free template, but two independently noisy copies have expected divergence:

```text
2*p*(1-p) + (2/3)*p*p = 1.9867%, for p = 0.01
```

That is about **29.8 mismatches per 1,500 bp pair**, versus the **15 mismatches** allowed at 99% identity. In 10,000 sampled same-template read pairs, the observed mean was **29.743 mismatches**, and only **27 pairs (0.27%)** passed 99% identity. Those pairs are sampled from the generated reads, not supplied to clustering.

Even against their known error-free templates, only **2,799 of 5,000 reads (55.98%)** passed 99%, because the per-read error count varies around 15. Therefore 99% clustering of raw Q20 reads should not be interpreted as recovering the original template groups in this model. A different threshold or a denoising/consensus step would require a separate evaluation of grouping and separation accuracy.

## Validation and reproduction

The complete 5,000-read result was checked against exhaustive clustering, which took 1.48 s outside the timed repetitions. All lossless index variants matched it. An additional 10,000 long-read comparisons checked the packed and existing byte verifiers against direct mismatch counting.

The [simulation program](../.claude/archive/benchmarks/Test99OntSimulation.cpp) includes and reuses [Test99Candidates.cpp](../.claude/archive/benchmarks/Test99Candidates.cpp). It generates its own templates and reads from fixed RNG seeds, so no input files or external database are needed:

```sh
make -j4 WITH_HTS=0 STATIC=0 BUILD_DIR=/tmp/sdm-99-ont-build \
  /tmp/sdm-99-ont-build/Test99OntSimulation
/tmp/sdm-99-ont-build/Test99OntSimulation > candidates.csv 2> diagnostics.txt
```

The optional first argument changes reads per template; the default is 100. Raw results are saved as [candidate timings](benchmarks/99-identity-ont-q20-candidates.csv) and [simulation/validation diagnostics](benchmarks/99-identity-ont-q20-diagnostics.txt). The private build did not replace the shared SDM executable or modify production IO code.
