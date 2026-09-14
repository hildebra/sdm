# Selected single-end qualities in diffDNA

**Paired support is now implemented:** see [exact paired variants](paired-diff-quality.md). The measurements below describe the preceding single-end implementation.

Implemented and evaluated on 2026-09-13. Enable with `-derepStoreQuals 1`; the default is `0`. The flag is also available as a tab-separated options-file entry, with command-line precedence.

```sh
sdm -i_fastq ont.fastq -o_dereplicate clusters.fna -options sdm_options.txt \
    -derepIdentity 97 -derepStoreQuals 1 -threads 1
```

The output includes `clusters.subclusters.fq`, with one reconstructed full processed sequence and one selected observation's qualities per exact nucleotide variant, counted with `;size=N;`. The usual parent output still contains one representative per coarse cluster. Original read names and the qualities of losing observations are discarded; this is not a raw-read archive or averaged-quality consensus.

The flag enables exact variants in memory without requiring `-derepStoreDiffs 1`. That separate option writes the existing nucleotide-only `SDMDIFF2` sidecar. Reloading `.diff` cannot reconstruct the retained qualities. `-derepSubclusterFasta 1` continues to require `-derepStoreDiffs 1` and exports the same nucleotide variants as FASTA.

## Implementation and scope

- **One reference, independent quality replacements.** Each coarse cluster retains its immutable initial sequence. Variant keys encode processed length and differing bases against that reference. Selecting a new parent representative or a better quality witness does not change the reference or recalculate its nucleotide patches.
- **Shared representative selector.** `SeedSelection.h` contains the modern scoring routine used by both the existing `DNA` wrapper and compact `SeedReadStats`. The classic external-reference seed selector remains unchanged. The compact record preserves field precision, primer flags, expected errors, mean quality, logical and backing lengths, and the ambiguity count of the original backing sequence. This matters when trimming hides a tail that the existing selector still examines.
- **One quality vector per exact variant.** A winning observation replaces its quality vector and scoring summaries together. Comparison history and copy/sample counts survive replacement. The vector stores only qualities corresponding to the full processed sequence. Duplicate observations compare against the compact summary without expanding the old nucleotide sequence. Parent representative selection still runs separately; skipping that comparison could change established behavior.
- **Eligibility is preserved.** A passing observation takes priority over a failed witness. A failed observation admitted to an existing passing cluster can supply the first witness for its exact variant, but subsequent failed observations only add counts. Failed reads cannot establish a new coarse cluster.
- **Snapshots and reassignment work.** Immutable quality payloads are shared across snapshots; later replacement leaves earlier output snapshots intact. Optional final reassignment moves the complete variant, its selected qualities, and all sample occurrences together. Rebased patches use the destination's initial reference. Merging already compressed records compares their surviving witnesses; it cannot replay discarded observations under a different selection policy.
- **Streaming export.** FASTQ export expands one sequence at a time and reuses the existing DNA writer and configured quality offset. Generated subcluster labels replace the parent's abundance with the exact variant's abundance. Exports include groups below the parent minimum-copy cutoff, and follow the existing per-run naming scheme.

The initial implementation described here accepts unmerged single-end input with complete qualities. The subsequent [paired extension](paired-diff-quality.md) accepts complete unmerged pairs as well. A previously merged FASTQ supplied as ordinary single-end input can use this mode, but its original mate linkage is then unavailable. At identity 100 the flag selects the streaming same-start exact-prefix backend. Default 100% dereplication remains on its existing backend when the optional variant modes are disabled.

Matching remains same-start, substitutions-only, across the shorter search key. Homopolymer-containing search seeds are excluded by the existing policy; homopolymer bases still participate in final identity checks. Internal indels can substantially reduce grouping/compression. The positional delta encoding reconstructs such observations exactly, but a shifted suffix can require many stored changes. Neither the storage change nor the tests add an indel-aware clustering algorithm.

## Validation

The private portable build used GCC 13.3, C++20 and `-O3`. Tests passed with the following coverage:

- The retained-quality C++ test compares selected payloads, scoring summaries, comparison history, and counts against uncompressed `DNAunique` selection over 2,400 observations. It covers 15–1,500 bp variants, ambiguity, primer flags, hidden trimmed tails, failed observations, reference variants, and varying quality vectors.
- An independent copy of the previous selector agrees with the refactored production wrapper on 24,000 single/paired/merged comparison cases, including external-reference and classic selection. This was an additional local refactoring check; the maintained tests exercise the shared production path.
- The full FASTQ pipeline test uses 600 distinct noisy sequences, with three observations each, lengths approximately 1,400–1,600 bp, 1% simulated substitutions, explicit homopolymer insertion/deletion cases, and Q10/Q20/Q30 profiles. All 1,800 observations and selected per-base qualities match the source oracle at 97%, 98%, 99%, and 100%; also with final reassignment and with four threads plus a 700 bp search cap. It validates full sequences beyond that cap and sample-map totals. Quality/error generators are independent, and this is a rough fixture rather than a calibrated ONT model.
- Additional C++ cases cover FASTQ offset 64, binary compatibility, absent-quality rejection, paired/merged rejection, overflow, reference rebasing, self-merges, retained snapshots, resumed input, per-run filenames, below-cutoff exports, and CLI/options-file precedence.
- Existing input-helper, search, targeted regression, ONT pipeline, and coarse CLI checks pass. AddressSanitizer, UndefinedBehaviorSanitizer, and leak detection pass for the coarse/diffDNA suite.

Run the maintained pipeline test directly with:

```sh
python3 tests/test_diff_quality_pipeline.py /path/to/sdm
```

It is included in `make test`, alongside `TestCoarseDerep` and the existing regression suites. Builds used private directories and did not replace the workspace's root executable.

## Benchmarks

The [raw CSV](benchmarks/diff-quality-production.csv) contains 42 runs: three independent processes per configuration, interleaved in a deterministic shuffled order. Host: Intel Core Ultra 7 258V, x86-64 Linux/WSL, GCC 13.3, portable/shared-library build, one processing thread. The simulation is the existing 50-template conserved/variable-block model in `tests/CoarseSimulation.h`. Nucleotide errors are independent substitutions at nominal Q20 (1%) or Q30 (0.1%). Quality profiles alternate at nominal Q minus six, nominal Q, and nominal Q plus six in **both** compared modes, to exercise replacements. These profiles are not calibrated to each simulated error draw.

Times include simulation, DNA construction, production batching/clustering, parent FASTA/HQ/map output, and every enabled sidecar. They exclude FASTQ input parsing, upstream filtering, fixture cleanup, and the subsequent output validation scan. Quality-enabled times include reconstructing and writing the full `.subclusters.fq` file. Peak RSS is the process's Linux `VmHWM`, including allocator/application overhead, sampled after output validation. All runs use `derep_reserve=1`; every run verifies parent abundance, and quality-enabled runs additionally verify FASTQ sequence/quality lengths and total exact-variant abundance.

### 20,000 reads, nominal Q20 errors, 1,500 bp

These inputs contain 20,000 distinct processed nucleotide sequences, so every observation needs its own selected quality vector. Both columns below enable binary nucleotide differences; the quality column additionally retains qualities and writes exact-subcluster FASTQ. Values are medians of three runs.

| Identity | Coarse clusters | Deltas only: seconds | Deltas + qualities: seconds | Deltas only: peak MiB | Deltas + qualities: peak MiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| 97% | 156 | 0.889 | 1.069 | 10.99 | 46.05 |
| 98% | 5,572 | 7.811 | 8.128 | 42.73 | 77.71 |
| 99% | 19,400 | 15.559 | 15.942 | 122.04 | 161.50 |

The median increases are approximately 20%, 4%, and 2.5%, respectively. The smaller percentages at 98–99% reflect the existing candidate-search cost dominating the run; they do not establish negligible quality overhead. Timing varies on this desktop host: the quality-enabled 97% runs span 0.921–1.501 s, and the 99% runs span 15.640–19.914 s. Treat small differences cautiously.

The raw quality vectors alone require 28.61 MiB here. Actual extra peak memory is about 35–39 MiB because payload objects, reference-variant entries, comparison history, and allocation overhead also occupy space. Each FASTQ sidecar is approximately 60.6 MB; nucleotide compression in memory does not make the expanded FASTQ smaller.

The intended minimal quality configuration, **97% with quality retention but no binary delta output**, takes **1.006 s and 46.06 MiB**. Ordinary length-matched exact dereplication takes **0.357 s and 110.51 MiB** on the same reads; it retains 20,000 representatives. Thus the new mode uses about **58% less peak memory**, at roughly **2.8 times the runtime** in this fixture. Coarse counts-only mode is smaller still: 0.817 s and 7.47 MiB, but does not preserve all exact variants/qualities.

At 99%, quality retention is larger than ordinary exact dereplication in this fixture. Strict identity admits too few noisy reads together to amortize cluster/reference/index overhead. The mode is therefore useful for stronger grouping, not a universal memory reduction at every threshold.

### Other configurations

| Input/configuration | Deltas only: seconds / MiB | Deltas + qualities: seconds / MiB | Result |
| --- | ---: | ---: | --- |
| Q20, 20,000 reads, terminal lengths 1,400–1,600 bp, 97% | 1.382 / 12.27 | 1.680 / 47.92 | 147 coarse clusters; 20,000 exact variants |
| Q30, 1,000,000 reads, 200 bp, 97% | 2.944 / 9.13 | 3.272 / 24.40 | 47 coarse clusters; 44,728 exact variants |

The Q30 case adds approximately 11% to the median runtime and retains 8.53 MiB of raw quality bytes. Exact repeats increment counts and can improve the selected witness without allocating another permanent quality vector. Its FASTQ sidecar has 44,728 records totaling about 19.2 MB.

On the Q20 1,500 bp case, enabling **final reassignment plus quality retention**, with binary delta export disabled, takes **5.589 s and 57.64 MiB**, moving 7,375 observations. This is substantially more work than the 1.006 s quality-retention-only run because reassignment searches the frozen final representatives and rebases variants. It remains independently optional and off by default.

Internal-indel reconstruction is covered by the full pipeline validation fixture, not the performance simulation above. Real ONT length/error/abundance distributions and indel rates may change both grouping and memory use considerably.

### Reproduce

Build the benchmark without replacing the root executable:

```sh
make WITH_HTS=0 STATIC=0 BUILD_DIR=.build/quality-bench .build/quality-bench/TestCoarseBenchmark

# identity reads diffs length-jitter length nominal-Q reassign retain-qualities
.build/quality-bench/TestCoarseBenchmark 97 20000 1 0 1500 20 0 0
.build/quality-bench/TestCoarseBenchmark 97 20000 1 0 1500 20 0 1
.build/quality-bench/TestCoarseBenchmark 98 20000 1 0 1500 20 0 1
.build/quality-bench/TestCoarseBenchmark 99 20000 1 0 1500 20 0 1
.build/quality-bench/TestCoarseBenchmark 97 1000000 1 0 200 30 0 1
```

Supplying the final argument enables the variable-quality fixture in both the `0` and `1` cases; omit it only when reproducing the older uniform-quality benchmark. Repeat each configuration three times in interleaved order. The CSV records all arguments and the observed simulated error fraction.

## Paired design choices from the initial evaluation

The existing shared selector already accepts R1, R2, and merged scoring views. The substantial design choice is what constitutes an exact variant and which fragment data must survive; the scoring formula does not need another implementation.

| Option | Representation and output | Tradeoff |
| --- | --- | --- |
| Premerge, then use the single-end mode | Treat each merged sequence and its consensus qualities as one observation; exact merged variants get FASTQ output. | Already possible for separately prepared merged FASTQ. Original mates, failed merges, and fragment linkage are absent. Supporting native merges in SDM would additionally need the original scoring summaries and explicit handling of unmerged input. |
| R1-keyed subclusters with one selected complete pair | Keep the current R1 variant key and attach the selected fragment's R2 delta, both quality vectors, and both scoring summaries. Export synchronized R1/R2 files. | Smallest native paired extension. Different R2 sequences sharing identical R1 collapse into one R1 subcluster, so the output does not preserve exact pair-level variants. |
| Exact `(R1, R2)` subclusters | Use a length-delimited pair of delta keys, immutable references for both mates, one quality vector/summary per mate, and one shared sparse fragment-count vector. | Recommended when retaining unique paired sequences matters. Preserves variants differing only in R2, at the cost of more unique keys and a second quality vector. Candidate grouping can still start from R1, with an explicit policy for verifying R2. |
| Native paired and merged observations in one mode | Add fragment state, mate/merged presence, orientation and overlap metadata, plus separate score views and output rules for successful/failed merges. | Most complete, but requires broader changes to admission, reassignment, and synchronized output. Whether exact identity means the two raw mates or the merged consensus must be decided explicitly. |

For a native paired implementation, select and replace the **entire fragment's payload atomically**. Independently choosing the best R1 and best R2 could synthesize a pair that was never observed. Keep comparison history and sample abundance at fragment level, count each fragment once, and use matching generated IDs/counts in both output files. Per-mate deltas can reuse the current encoder/decoder; avoid putting two complete sample-counting `diffDNA` objects underneath a fragment.

The exact pair-keyed option has since been implemented; see [paired behavior and validation](paired-diff-quality.md). The initial recommendation was to use exact pair-keyed variants as the next native paired mode. Start with unmerged, complete pairs, reuse two compact scoring views and the current selector, and add native merged/mixed input afterward. This keeps pair fidelity explicit and makes synchronization, reassignment, missing-mate handling, and count conservation testable before adding merge-state complexity. Runtime/memory for paired mode have not been measured; raw retained quality bytes will scale with `len(R1) + len(R2)` per selected unique pair.
