# Coarse storage with ordinary output semantics

Updated 2026-09-15. This supersedes the earlier recommendation to enable full-variant retention automatically in coarse LotuS workflows. With only `derepIdentity` changed from 100 to 97, the contract is ordinary exact dereplicates and ordinary representative seed candidates. `derepStoreQuals=1` remains an explicit full-variant policy, rather than a prerequisite for keeping representative qualities.

## Corrections

The supplied quality-accumulation patch is applied. It passes the dereplicated-output format and merged-quality setting into the coarse path and restores ordinary `matchedDNA` accumulation for exact output. Deferred prefix failures use ordinary quality accumulation when no variant export is selected. The regression reproduces the earlier main FASTQ quality difference: before that patch, main/rest FASTQ differed even though the representative HQ pairs agreed.

An additional bug remained: every observation of an exact key went into coarse subgroup zero. Ordinary storage instead keeps merge-offset subgroups and selects among them using their abundance and merged length, then gathers their evidence. A rare, high-quality merged observation could therefore replace the representative of a more abundant unmerged subgroup only in coarse storage, changing parent IDs, HQ mates and downstream seeds without changing total abundance.

Storage-only updates now share `DNAuniqSet::addObservation` with ordinary updates and retain the original merge position. This preserves subgroup selection, sample counts, per-base quality sums/coverage, merged quality evidence, and the selected full representative pair. Export uses the existing ordinary snapshots, prefix recovery, copy cutoffs, map and main/merged/rest/HQ writers. Internal search-group IDs are not exposed to the seed reader. Empty fixed search keys follow the ordinary tracker path; explicit variant validation remains separate.

These changes apply to storage-only exact output (`derepStoreQuals=0`, `derepStoreDiffs=0`, no coarse-cluster output or reassignment). Explicit full-variant and coarse-cluster policies retain their separate behavior. In particular, native merged search with explicitly retained full variants remains unsupported; that limitation does not apply to storage-only mode.

## Compact persistent representation

Default coarse storage now uses `CompactDNA` within `DNAuniqSet`. The index still owns the ordinary exact search keys, but **does not keep a full `DNA`/`DNAunique` object per key**. Each coarse group supplies immutable full R1/R2/merged references; its search key borrows the appropriate reference. Exact parents retain the ordinary merge-offset subgroups.

With multiple workers, storage-only routing uses the smallest hash among the candidate index's eligible positioned 20-mers, falling back to the full-key hash when no seed is eligible. This lets nearby keys share more references than hashing the whole sequence, while identical keys always reach the same owner. It does not require all 97% matches to land in one partition. Explicit variant/coarse-cluster modes keep their established routing.

Each subgroup stores:

- A nucleotide-only `diffDNA` containing full available unmerged fragment variants, their multiplicities, and aggregate sample counts. Missing/physically emptied mates are representable internally; they do not change ordinary admission. Variant qualities are not allocated.
- One losslessly packed selected pair/consensus, expressed as deltas from the shared references, with its original IDs, full observed quality vectors, logical lengths, filtering/merging metadata and comparison history. Quality arrays choose raw, run-length or base-offset bit packing, whichever is smallest, without rounding or quantization.
- Ordinary main-FASTQ sums and per-base coverage, including detached merged-quality evidence, when main FASTQ is requested. Packing retains their full integer precision.

FASTA-mode repeats use cached `SeedReadStats` and the existing `seed_selection::better` rule. Only a winning replacement needs the selected pair reconstructed for the ordinary `takeOverDNA` routine. Main-FASTQ updates reconstruct one subgroup and call `matchedDNA`, preserving ordinary accumulation instead of implementing separate arithmetic. This path can cost more CPU; timings are reported alongside memory.

At export, only the current exact parent's merge-offset subgroups are reconstructed. Ordinary subgroup weighting, prefix recovery, copy cutoffs, main/merged/rest/map writers and HQ seed files remain in use. No collection-wide expansion or coarse-parent merge is needed. Internal nucleotide variants do not add `.subN` HQ candidates. `derepStoreQuals=1` still explicitly requests the separate full-variant quality policy.

Deferred failed-prefix observations still use the ordinary evidence-only buffer until recovery, and empty fixed keys follow the ordinary tracker. The compact representation described here applies to admitted coarse exact-key storage; the fixed-length memory benchmark does not measure a large deferred-prefix buffer.

Physical technical cuts are already part of the supplied read objects. HQ output restores their full available backing sequence and observed qualities; it does not restore physically removed barcodes/primers. The existing seed-extension parity checks cover this input contract.

## Acceptance coverage

`tests/test_storage_only_parity.pl` reuses the 1,000-pair, two-sample mixed-quality fixture and compares **ordinary 100% with storage-only 97%**. It supplies variable observed quality strings, R2-only variants, a dominant unmerged subgroup with a rare better merged observation, failed/mid-quality repeats, failures preceding passing anchors, short prefixes, below-cutoff parents and separate sequencing runs.

There are **48 comparisons**: eight processing profiles, main FASTA/main FASTQ, and 1/4/12 workers. Profiles exercise capped R1 search, capped merged search, merged output with R1-only search, uncapped prefix recovery with/without merging, truncation/quality trimming with/without merging, and cumulative sequencing-run output. Both modes use the same binary, inputs, filters and minimum-copy cutoff.

The checks compare complete IDs, sequences, quality strings and abundances, permitting only record/map-entry ordering differences. They validate that:

- Main + merged + rest records partition the exact parent map with matching abundance headers and per-sample counts.
- HQ contains one linked representative pair per parent, with full available sequences and observed qualities from the same input pair, including below-cutoff parents. No `.subN` candidates appear.
- Search-policy logs confirm merged-key use in merged-search cases. `noSearchWithMerge` is omitted there: its presence selects R1 search even when its value is `0`.
- Normal seed commands agree after replacing the output directory, without `seedSubclusters`. Seed statistics, final FASTQ seeds and matrices agree as well.
- Results agree across the tested worker counts when representatives are unambiguous.

The seed comparison runs SDM with deterministic parent assignments; it does not execute external LotuS clustering/mapping tools. The separate LotuS integration regression described by the caller remains the check of the complete pipeline.

The older executable reproduces main FASTQ accumulation differences; the executable containing only the supplied accumulation patch still reproduces merge-subgroup differences. The corrected path passes the 48 comparisons. Native tests separately cover subgroup selection, accumulation, exact ties and fixed-key admission ordering. A further 1,000-pair compact-storage test verifies empty persistent DNA maps, shared coarse references, reconstruction of every admitted nucleotide variant/count, changing representatives, variable logical lengths and quality coverage, failed/empty mates and mixed merged/unmerged subgroups, for both main formats.

## Ties and ordering

IDs are not normalized away to make comparisons pass. During testing, equal-quality observations `S0_fzero_empty_r2_4/1` and `S0_fzero_empty_r2_14/1` supplied different IDs under different schedules while their sequence and quality strings were identical. The primary parity fixture now has an unambiguous best observation in that subgroup; the observed tie is retained as evidence of the ordering limitation.

The native tie test supplies two equal observations in both orders: both ordinary and coarse storage keep the first ID, and reversing the input reverses that ID. The failed-read test similarly verifies one admitted observation when failures precede the first passing anchor versus three when that anchor comes first. Prefix recovery is independently tested through the ordinary deferred path.

The existing representative selector also uses comparison history and a replacement threshold; it is not a total, order-independent ranking. This update preserves that policy. Arbitrary parallel runs cannot be promised byte-identical IDs or admission outcomes in order-sensitive cases merely because output granularity is the same. Deterministic total ordering or deterministic admission scheduling would be a separate change to ordinary behavior; one-worker runs with the same input order provide the controlled comparison here.

## Peak memory, measured separately

`.claude/archive/benchmarks/benchmark_storage_only.py` runs ordinary and 97% storage as fresh processes, alternates mode order and measures whole-process peak RSS using `/usr/bin/time`. Each measured run also passes a strict canonical checksum of parent IDs, main/merged/rest sequences and quality strings, sample maps and both observed HQ quality strings. Input generation and validation are outside the measured process. Sorting changes record order only.

Local GCC 13.3, `-O3`, C++20, static portable build, **250,000 paired 200 bp fragments and 50,000 exact keys**, three repetitions per configuration. Merging is disabled in both modes; search is explicitly capped R1, minimum copies is 1, `derepStoreQuals=0` and `derepStoreDiffs=0`. `derep_reserve=1` prevents a large default reservation from inflating the ordinary comparison. The first observation per key has basewise random Q40–45, later observations Q25–30, so representatives are unambiguous and the quality strings resist simple run-length encoding. Closely related families contain 100 exact keys differing in the first four bases; the control uses unrelated random pairs.

Median results (negative RSS change means less memory):

| Input diversity | Main format | Workers | Ordinary 100%, MiB | Compact 97%, MiB | RSS change | Ordinary, s | Compact, s | 97% groups |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Families of 100 nearby keys | FASTA | 1 | 157.68 | 136.25 | -13.6% | 1.007 | 1.489 | 500 |
| Families of 100 nearby keys | FASTA | 12 | 168.80 | 154.30 | -8.6% | 0.614 | 0.806 | 4,506 |
| Families of 100 nearby keys | FASTQ | 1 | 233.30 | 159.75 | -31.5% | 1.211 | 3.099 | 500 |
| Families of 100 nearby keys | FASTQ | 12 | 244.55 | 175.80 | -28.1% | 0.734 | 1.479 | 4,506 |
| Unrelated keys | FASTA | 1 | 157.80 | 180.85 | +14.6% | 0.955 | 10.914 | 50,000 |
| Unrelated keys | FASTA | 12 | 168.93 | 197.39 | +16.8% | 0.533 | 0.830 | 50,000 |
| Unrelated keys | FASTQ | 1 | 233.43 | 204.16 | -12.5% | 1.179 | 12.837 | 50,000 |
| Unrelated keys | FASTQ | 12 | 244.68 | 217.95 | -10.9% | 0.683 | 1.431 | 50,000 |

All 48 measured runs retain **50,000 external parents and 50,000 HQ pairs** with matching IDs, sequences, qualities, per-sample counts and mate relationships. Internal group count does not affect output granularity. Raw observations and executable/platform metadata are in [the CSV](benchmarks/compact-storage-memory.csv) and [metadata](benchmarks/compact-storage-memory.json).

For the related-read fixture, compact storage reduces peak RSS by about **14%/9% for main FASTA** and **32%/28% for main FASTQ**, at 1/12 workers. Seed-based partitioning retains 4,506 internal references with 12 workers; whole-key partitioning retained 25,413 on this fixture before that change. The admitted exact-key store contains packed records and nucleotide patches, rather than full DNA objects plus an extra index.

These are not universal savings. The unrelated FASTA control consumes about **15–17% more memory**: almost every exact key needs a separate coarse reference and postings. FASTQ saves roughly 11–13% even in that control because sums/coverage can be packed. Many distinct R2/tail variants, broad quality distributions, deferred-prefix buffers, merged input and other read lengths can change the tradeoff. This benchmark does not measure the user's HPC dataset.

Compression also costs CPU. On the related fixture, wall time is roughly 1.3–1.5 times ordinary for FASTA and 2.0–2.6 times for FASTQ in these short runs. Single-worker unrelated coarse search is much slower; its candidate-search overhead dominates. The existing worker pool and exact-key reuse remain in use. Main FASTQ currently reconstructs a subgroup for each quality update to reuse ordinary arithmetic. Lower peak memory is reported separately from runtime and output parity.

The earlier [storage-only CSV](benchmarks/storage-only-memory.csv) and [metadata](benchmarks/storage-only-memory.json) describe the previous full-DNA representation on a smaller, constant-quality fixture. They remain historical measurements, not measurements of this compact implementation.

## Simulated Q20 long reads

A separate [100,000-read Q20 benchmark](benchmarks/compact-storage-ont-q20.md) uses 1,500 bp single reads from 50 closely related templates, variable quality scores and approximately 1% substitution errors. Compact storage reduced peak RSS by about 30% with one worker and 27% with 12 workers while retaining all 100,000 exact output records. It also exposed limited parallel scaling: one partition received 70% of the reads. See that report for memory/runtime measurements, parity checks and simulation limits.

## Ordinary 100% regression check

A separate alternating before/after benchmark exercises the default ordinary paired path with no identity/retention overrides, 250,000 pairs, 50,000 exact keys and the normal reservation setting. Across three repetitions, median wall time was 1.499 → 1.513 s with one worker and 0.588 → 0.572 s with 12 workers. These short runs show no consistent slowdown; changes are within about 3%. Peak RSS rose by 0.82 MiB in both cases (164.36 → 165.18 MiB and 175.73 → 176.55 MiB), consistent with the extra optional compact-state pointer/allocation rounding in the exact-key wrapper. This comparison validates ordinary sequence/count output, independently of the stricter same-binary 100%/97% output checks above. See [raw measurements](benchmarks/compact-storage-ordinary-check.csv) and [metadata](benchmarks/compact-storage-ordinary-check.json).

## Run the checks

```sh
make WITH_HTS=0 test
perl tests/test_storage_only_parity.pl ./sdm --output /tmp/sdm-storage-parity
python3 .claude/archive/benchmarks/benchmark_storage_only.py ./sdm --output /tmp/sdm-storage-memory \
    --reads 250000 --templates 50000 --threads 1 12 --repeats 3 \
    --formats fa fq --quality-profile variable
```

Use fresh output directories. The full portable `make WITH_HTS=0 test` suite passed with these changes. It includes ordinary-storage parity and the existing explicit-variant and seed-extension regressions. See [LotuS3 IO instructions](lotus3-dereplication-io.md) for the corrected flag and file contract.
