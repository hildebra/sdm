# Quality-aware exact variants in diffDNA: feasibility

**Update:** the optional single-end mode is now implemented as `-derepStoreQuals 1`. See [production behavior, validation, and measurements](diff-quality-retention.md). The results below document the preceding feasibility investigation.

Evaluation date: 2026-09-13. This is an investigation and isolated prototype; production storage, options, and output behavior have not changed.

**Feasible for one selected FASTQ per exact nucleotide variant.** The current nucleotide deltas already reconstruct its FASTA. Adding one selected quality vector and a compact scoring record per variant permits the current representative-selection logic to operate without reconstructing the stored sequence. This needs a localized refactor of the scoring interface and variant entries. Single-end ONT is the simplest initial scope; paired/merged output and persistent quality-bearing sidecars add work.

## What the current selector actually reads

The active dereplication path calls `whoIsBetter(..., forSeed=false)` from `DNAunique::matchedDNA`; the old `betterPreSeed` member is an unused aborting stub. The active selector is in [InputStream.cpp](../InputStream.cpp); accessors and replacement behavior are in [DNA.h](../DNA.h) and [DNA.cpp](../DNA.cpp).

| Information | Retention needed |
| --- | --- |
| Sequence identity and processed length | Existing immutable reference and exact delta key, including length and suffix bases |
| Output qualities | One `vector<qual_score>` belonging to the selected observation, matching the full processed variant length |
| Quality scores | Cached mean quality as `float` and expected-error sum as `double`, calculated with the existing DNA accessors |
| Eligibility and primer information | Passing/failed status and forward/reverse-primer flags; preserve admission/replacement gates |
| Sequence summary | Logical length and full-backing ambiguity count; backing length for the classic seed selector |
| Paired/merged comparisons | Presence plus separate R1, R2, and merged summaries, merge length and overlap error metadata |
| Comparison context/history | Existing merge-comparison counters, caller-provided abundance, and identity/`ever_best` context where applicable |

No active scoring branch needs the actual nucleotide string after these summaries are available. The sequence scan is `numNonCanonicalDNA(true)`: it counts across the **backing sequence**, while mean quality and expected errors use the logical trimmed prefix. Capture that ambiguity count before discarding the backing sequence; deriving it from the shorter delta reconstruction can change a decision. Average quality alone is insufficient: reads with the same average can have different expected-error sums.

Preserve field precision and the order of the existing checks. In particular, primer/identity rejection precedes the merge-history update, while some later rejections still update that history. Replacing a selected observation updates its qualities and associated summaries together. Sample counts and comparison history describe accumulated observations and survive replacement. Identity against a later external reference is comparison context and must be refreshed for that reference.

## Proposed integration

1. Extract the statistical comparison into a shared helper accepting compact read summaries and explicit comparison history/context. Keep the DNA-based function as a wrapper. The new diff-backed path uses the same helper rather than maintaining another scoring formula.
2. Add an optional quality-bearing variant-record mode inside `diffDNA::SequenceData`. Keep a single map key per exact variant, carrying its existing count or sparse samples and its selected quality/statistics payload. Include a payload for the reference sequence itself; the current reference-only counter is insufficient for quality retention. Keep count-only/default modes unallocated for these additions.
3. After clustering and exact-delta lookup, initialize the first observation or compare a later observation against the stored scores. Allocate/copy the quality vector on insertion or a winning replacement. Rejected duplicates only update counts/history. Reuse incoming summaries across cluster and variant comparisons where possible. Moving qualities out of an incoming shared DNA requires proven ownership; copying on selection is the safer initial path.
4. Reuse `diffDNA::forEach`/the existing decoder at export, expanding one sequence at a time. Attach its selected qualities to a temporary DNA and reuse the current FASTQ writer and offset handling. Generate a subcluster header such as `parent.sub7;size=42;`. Validate equal sequence/quality lengths and require available qualities. This emits the selected observation's qualities; it is distinct from the averaged-quality output used by ordinary exact FASTQ dereplication.
5. Carry the quality payload through snapshots and final reassignment alongside the exact variant's sample counts. Rebasing a nucleotide patch leaves its quality vector intact because the represented sequence is unchanged. Shared immutable quality payloads avoid copying every quality vector when a retained snapshot forces the variant map to detach. Coalescing variants compares the retained selected observations; the discarded qualities cannot support replaying every original read under a different later selection policy.

Keep the initial cluster reference immutable. A new high-quality read of an existing exact variant changes only its selected quality/statistics record. A new cluster representative can identify another variant while every patch stays relative to the original reference. Stable IDs or owned payload handles are preferable to retaining addresses into a mutable hash table.

Retaining the existing full representative per coarse cluster is a useful first step: it limits changes to normal parent output and existing callers. Per-variant selection then adds a comparison for duplicates in addition to the current cluster-level selection. The selector depends on context/history, so skipping every cluster comparison merely because the read lost its local variant comparison would require a separate equivalence argument.

For paired reads, the emitted merged vector alone cannot reproduce selection: when one or both reads are unmerged, the selector may use R1/R2 summaries instead. Full paired FASTQ reconstruction also needs the corresponding mate nucleotide observations and qualities, which current deltas do not retain. That extension is materially larger than single-end quality retention.

The existing `SDMDIFF2` format has no qualities. Same-run FASTQ export can use the optional in-memory payload directly. Reconstruction in a later process requires a versioned format extension or an additional explicitly associated sidecar; old nucleotide-only files cannot recover discarded qualities. Persistent restart of seed comparison also needs its metadata/history and defined comparison semantics.

## Prototype checks

An isolated probe mechanically instantiated the current selector body for either DNA objects or a compact accessor-compatible statistics record. The production selector was also compiled from the same frozen source and used as the decision oracle. **24,000 comparisons matched in decisions, `ever_best` updates, and all three merge-history counters.** Cases included logical trimming with hidden backing tails, ambiguous bases, primer flags, paired and mixed merged/unmerged candidates, varied abundances, and modern dereplication, modern seed, and classic seed modes.

Thirty exact variants, including shortened reads, suffixes, and ambiguous symbols, were reconstructed through the existing `diffDNA` decoder. Their FASTQ output using retained winning qualities matched the selected DNA records exactly. These are checks of the comparison adapter and reconstruction primitives, not a complete implementation of the new mode: admission gates, variant metadata merging, quality-bearing serialization, and concurrent snapshot updates still need integration tests.

The scalar scoring record occupied **40 bytes** in this build. A quality vector header occupies 24 bytes and the existing merge-history record 24 bytes. Together that is 88 bytes before the quality bytes, eligibility flags, optional mate/merge summaries, containers, ownership handles, and allocator overhead. History may be unnecessary per variant in a restricted unmerged-only mode, but omitting it should be an explicit scope decision.

## Comparison microbenchmark

Medians of three fresh runs, GCC `-O3`, C++20, one million comparisons per row. Each run used 256 prepared single-end reads; per-read mean quality and per-base qualities varied, so both accepted and rejected decisions occurred. DNA quality caches were warmed, as they commonly are after filtering. The same extracted selector body was instantiated for DNA and statistics. The incoming-statistics column includes the incoming read's backing-sequence ambiguity scan; the cached column assumes both summaries already exist. All three strategies produced identical decision counts.

| Length | Full DNA objects (s) | Incoming stats + stored reference stats (s) | Both stats cached (s) |
| --- | ---: | ---: | ---: |
| 200 bp | 0.286 | 0.177 | 0.082 |
| 1500 bp | 1.329 | 0.766 | 0.078 |

This supports inexpensive scoring with stored summaries and no repeated decompression. It does **not** predict an end-to-end speedup: the experiment excludes delta lookup, allocating/copying selected quality vectors, sample updates, clustering, and FASTQ output. Both cached and incoming-statistics strategies were run after the DNA strategy within each process; cache warming, specialization, and uncontrolled host load can affect the exact ratios. The cached case is a lower bound unless incoming summaries are already available. Full mode overhead also depends on how many variants are new, how often winners change, and whether cluster-level comparisons are retained in addition to variant comparisons.

## Memory and output tradeoff

`qual_score` is already one byte. Quality values vary at matching nucleotide positions, so nucleotide differences do not provide a lossless encoding of quality differences. The straightforward design retains **U × L quality bytes**, where U is the number of exact variants and L their average processed length.

| Illustrative workload | Retained exact variants | Quality bytes alone | Additional 88-byte records |
| --- | ---: | ---: | ---: |
| Earlier one-million-read, Q30, 200 bp simulation | 44,728 | 8.53 MiB | 3.75 MiB |
| Earlier 20,000-read, Q20, 1,500 bp simulation | 20,000 | 28.61 MiB | 1.68 MiB |

These are arithmetic storage estimates using earlier unique counts, not measured process peaks for the proposed mode. Add existing deltas, sample counts, cluster references, map/allocator overhead, and any additional mate state. They explain why the mode can remain much smaller than retaining a full DNA object and nucleotide string per variant, while using appreciably more memory than nucleotide-only differences. In Q20 full-length data, most error patterns were unique, so almost every read needs a retained quality vector. FASTQ export also writes both bases and qualities for every variant; the resulting output volume can dominate reconstruction cost.

Raw comparison measurements: [diff-quality-selector.csv](benchmarks/diff-quality-selector.csv). The workspace probe and frozen selector are under `.build/diff-quality-eval/`. Frozen InputStream source SHA-256: `a3d912f663508a06b358322ca9dcabde49987c74ab6a1aff7966c649313c7499`.

Recommendation: proceed with an **optional single-end quality-retention mode**, preserving one real selected quality vector per exact variant, shared scoring logic, and an immutable delta reference. Use same-run subcluster FASTQ export initially. Benchmark the integrated insertion/replacement path before committing to a percentage runtime target; paired/mate output and persistent quality-bearing files can then build on the same record structure.
