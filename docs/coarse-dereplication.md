# Streaming coarse dereplication

SDM can group noisy reads while retaining one high-quality representative and all assigned counts. `-derepIdentity 97.5` is an example setting for stronger compression of Q20 reads, not a universal ONT recommendation. The default remains 100% exact/prefix dereplication with no nucleotide deltas. There is no dependency on Savont.

```sh
sdm -i_fastq ont.fastq -o_dereplicate clusters.fna -options sdm_options.txt \
    -derepIdentity 97 -derepSeedHomopolymer 4 -derepStoreDiffs 1 \
    -derepSubclusterFasta 1 -threads 1
```

These options also accept tab-separated entries in the options file; command-line values win. `derepIdentity` accepts a decimal percentage greater than zero and at most 100, with up to six fractional digits. `derepSeedHomopolymer` accepts integers 2 through 20. `derepStoreDiffs` accepts 0 or 1 and defaults to 0. `derepSubclusterFasta` also defaults to 0; setting it to 1 requires `derepStoreDiffs=1` and exports exact subclusters as FASTA.

## Matching and retained state

With only `derepIdentity=97` enabled, ordinary outputs are the contract: retain the configured merged/R1 search policy, exact keys, prefix recovery, admission, copy cutoffs and sample counts. Main FASTQ uses ordinary per-base quality accumulation; standard HQ files contain one full representative pair with its observed qualities per ordinary parent. Representative qualities do not require `derepStoreQuals=1`, and the normal seed command needs no `seedSubclusters`. Explicit full-variant retention is a separate policy. See the [storage-only parity and memory assessment](storage-only-dereplication.md).

By default, coarse groups are internal candidate indexes and immutable compression references. The main FASTA/FASTQ, merged/rest files, `.map` and representative HQ pairs preserve ordinary exact search-key dereplicates, search/truncation caps, prefix recovery and copy cutoffs. Storage-only mode preserves the configured merged/R1 search policy and ordinary merge-offset subgroup selection. With explicit full-variant quality retention, the search source is currently R1 and R2 differences do not split that parent; the existing paired better-seed selection can replace both mates. Select `-derepCoarseClusters 1` only when the output should describe coarse clusters instead. `-derepReassign 1` also explicitly selects coarse-cluster output and performs final reassignment. `-derepGlobal 1` controls search scope, not output granularity.

Each newly admitted exact R1 key is assigned to an immutable coarse reference through the candidate index. Later exact-key hits bypass that search. In default exact-output mode, a non-passing exact-key miss is rejected immediately: approximate matching cannot admit a new failed exact parent. The report counts these shortcuts separately from candidate searches; prefix-mode deferred recovery is unchanged. Each ordinary merge-offset subgroup of an exact key owns one packed selected pair and its sample counts; nucleotide-only fragment variants share coarse references and store differences and occurrence counts. Default storage does not retain a full DNA object per exact key. `DNAunique::matchedDNA` and `whoIsBetter` retain the existing seed selection. Default output does not also retain a coarse representative or duplicate its variant/quality counts. Selecting a new seed does not move a search anchor or rebase stored differences.

Memory scales with packed exact representatives, their selected quality strings/accumulation evidence, shared coarse references and distinct nucleotide patches. Main FASTQ retains ordinary sums and coverage losslessly. One exact parent is reconstructed at a time for export; temporary read batches remain bounded. Compression depends on sequence and quality diversity, and coarse indexing adds overhead when keys are unrelated. Measured default-output peak RSS is reported in the [compact-storage benchmark](storage-only-dereplication.md#peak-memory-measured-separately). Explicit variant-output snapshots share immutable quality payloads. Text export reconstructs and sorts one parent's variants at a time to assign `.subN` independently of the compression reference and hash insertion order; it does not merge coarse groups.

Comparisons start at zero and cover `L = min(query key length, anchor key length)`. A match allows at most `floor(L * (100 - identity) / 100)` substitutions, inclusively. The decimal threshold is represented as an integer, avoiding floating-point boundary errors. Verification reuses the existing bounded `ReadMerger::percentIdentity` comparison and exits when the mismatch budget is exceeded. Bases, including ambiguity symbols, are compared literally.

There is no indel alignment, offset search, reverse-complement search, or homopolymer compression. Consequently, ONT homopolymer indels can still split clusters. Excluding homopolymers from seeds addresses candidate selection, not alignment. Reads should already have consistent orientation and starts.

Explicit positive `derepSrchLen` and `TruncateSequenceLength` cap search keys, including merged keys. The smaller applicable cap wins. Without such a cap, the full logical sequence is compared; `minSeqLength` alone does not shorten this mode's keys. `derepPrefix` also controls exact-parent output with streaming storage: automatic variable-length prefix mode and forced prefix mode reuse the ordinary final consolidation and deferred-recovery pass. Forced prefix mode overrides the search-only cap, while respecting `TruncateSequenceLength`. Explicit coarse-parent output supersedes this exact-prefix policy. Setting `derepStoreDiffs=1` or `derepStoreQuals=1` at identity 100 also selects streaming storage with exact effective R1 output keys. `derepReassign=1` requests coarse output instead.

Greedy grouping depends on input/batch order. Use one processing thread for repeatable runs. Only distance to the initial anchor is bounded: two members may differ by up to approximately twice the mismatch budget. If a short read establishes a cluster, unmatched longer suffixes do not affect admission. Choosing a different output representative does not change those membership semantics. These boundary effects concern internal grouping, or coarse-cluster output when explicitly requested; default output preserves distinct effective R1 keys.

Failed/alternate-tier reads contribute if an existing passing exact R1 dereplicate matches them (an approximate cluster in explicit coarse-output mode). They cannot create parents. With exact output and prefix mode, failures are buffered for recovery against passing prefixes at output; otherwise unmatched failures are discarded. Existing minimum-copy rules and sample-count map output apply to the requested output level. No matched reads are downsampled. As in the existing dereplication mechanism, a paired fragment contributes one parent observation.

Retention currently changes merged-search behavior and excludes physically empty mates. HQ export preserves full available backing mates and their observed qualities, including hidden tails. See the [four-mode behavioral audit](derep-behavior-parity.md) for corrected admission rules, remaining limits and count reconciliation.

## Search seeds and fallback

SDM indexes complete, non-overlapping 20-mers at positions 0, 20, 40, and so on. A seed is excluded if it contains an ambiguous base or a run of at least `derepSeedHomopolymer` identical bases. The default is four. Filtering is based on the contents of the seed, so an identical word has identical eligibility in query and reference. Homopolymer bases still participate in verification.

Postings are keyed by both the exact 40-bit word and its sequence position; hash collisions cannot imply a sequence match. Anchors are binned into 128-base length intervals. Each bin tracks the actual minimum and maximum search-key lengths of its stored anchors, including keys capped for clustering. For each bin, the query chooses the `E+1` rarest eligible, disjoint seeds lying wholly inside the overlap guaranteed by the shortest stored key, where `E` bounds the mismatch budget for the longest stored key. Using populated bounds avoids forced exhaustive scans for ordinary 200–255 bp amplicons merely because their bin starts at 128 bp. An absent posting list has frequency zero and still counts as a selected seed.

At most `E` substitutions can disrupt `E` disjoint seeds. Thus every admissible candidate must appear in at least one selected posting list. SDM verifies candidate IDs in creation order and accepts the first passing anchor. If there are fewer than `E+1` eligible seeds, it scans that length bin directly. This preserves sensitivity for short, low-complexity, ambiguous, and variable-length reads; a low threshold can therefore be substantially slower. There is no four-seed gate or arbitrary candidate cap that could silently reject valid assignments.

With identity below 100 and multiple workers, SDM now maintains 64 independent coarse indexes using the existing hash shards and a shared work pool. Processing workers submit their partition tasks and help execute any ready partition while waiting for their batch. Only one worker updates a given index at a time; queued tasks for a busy index do not block work on other ready indexes. No additional worker threads are created, and a submitting call finishes only after its entire batch has completed. Each index owns its search scratch, anchors, variants, and counts. Storage-only mode routes by the smallest hash of an eligible positioned 20-mer, reusing the candidate-index seed rules, with whole-key hashing as the fallback when there is no eligible seed. This improves reference sharing across nearby keys. Explicit variant modes keep whole-key hashing. Explicit coarse-cluster output (`-derepCoarseClusters 1`) always uses the global index, because approximate matches routed to different partitions would never be grouped; partitioning is a memory-saving arrangement for exact-parent output. Identical effective search keys reach one owner, independently of R2. This preserves normal R1 dereplication and paired better-seed selection. An abundant R1 can concentrate work in one index, although its repeats bypass candidate search. Approximate matches may live in different partitions, which are gathered for output without merging them. More internal anchors are possible, but exact R1 output counts and minimum-copy decisions do not depend on those boundaries. With explicit coarse-cluster output, representative counts and minimum-copy decisions can differ. Exact retained variants and selected quality witnesses remain available; fixed-key failed repeats contribute within their partition, while default prefix output reconciles compatible exact prefixes across partitions.

`-derepGlobal 1` selects the previous single global index when global grouping is required. One worker, streaming 100% identity, `-derepReassign 1` or `-derepCoarseClusters 1` also retains that path. Global coarse updates remain serialized; ordinary read processing remains parallel. Reassignment continues to consider all global representatives. The legacy exact-key mode retains its existing sharded updates. See [options and partition semantics](options.md#streaming-coarse-dereplication) and the [work-pool measurements](benchmarks/coarse-work-pool.md). The report counts exact-key reuse, candidate searches, partition tasks processed for other workers, and peak simultaneously active partitions; that peak is not sustained CPU utilization.

## Optional final reassignment

`-derepReassign 1` selects coarse-cluster output and enables one final pass at the configured `derepIdentity`; the default is `0`. It works with or without `-derepStoreDiffs 1`. For example:

```sh
sdm -i_fastq reads.fastq -o_dereplicate clusters.fna -options sdm_options.txt \
    -derepIdentity 97 -derepReassign 1 -threads 1
```

At output, SDM freezes the selected high-quality representatives and builds a temporary index using the same seed and length-bin implementation as streaming admission. Each exact nucleotide variant is reconstructed once and searched against these representatives. The candidate search uses bounded substitution verification and the existing homopolymer exclusions, with exhaustive fallback when too few seeds survive. It considers every qualifying candidate, rather than stopping at the first acceptable representative. An exact match to the current representative needs no candidate search. There is no fixed neighbor cap or heuristic similarity gate that could exclude a better qualifying target.

A variant moves, with all its copies, only when the target meets the configured threshold and has strictly higher identity than its current representative. Comparisons use the same start, shorter key, merged/R1 source selection, and search/truncation limits as admission. Identities with different overlap lengths are compared as exact integer ratios. Equal identities keep the original group. If several other targets share the best improved identity, the earliest-created target wins. For example, at 97% on 200 bp reads, a variant five substitutions from its current representative and three from another moves; a four-versus-four tie stays.

Representatives stay fixed for the whole pass. It does not iterate, create new clusters, or merge whole clusters. A member whose current representative is outside the configured radius remains assigned if there is no strictly better qualifying target. Therefore this step improves assignments but does not guarantee that every final member lies within the requested radius of its representative. Existing input-order dependence and possible loss of close biological variants still apply.

Reassignment requires exact variants and their sample abundances during streaming. SDM automatically enables the nucleotide portion of `diffDNA`, using one map from encoded variants to sparse sample-count vectors. This replaces the variant-to-count map used without reassignment. Exact repeats increment the corresponding sample entry; variant multiplicity is its sample-count sum. The aggregate cluster sample vector is rebuilt lazily when needed and invalidated on additions, rather than updated as a second set of sample counters on every read. No read names are retained; selected qualities require the separate `derepStoreQuals` option. Variant data uses copy on write for retained snapshots. With reassignment off, per-variant sample vectors are absent. The temporary sample details are not serialized: the `SDMDIFF2` file format stays nucleotide/count-only and cannot by itself restore the sample information needed for another reassignment run.

The pass constructs an independent output snapshot while producers are paused. It redistributes sample counts, and rebuilds nucleotide differences against the destination's immutable initial reference only when delta export is enabled. It then applies minimum-copy filtering and writes matching parent abundances, sample maps, and optional variant exports. Empty groups are omitted. Live input anchors, variant membership, representatives, and counts are unchanged, so a later snapshot or resumed input does not compound previous moves. Per-run output uses that run's current representatives.

This option adds memory for exact variants, per-variant sample counts, a temporary representative index, and the output snapshot; enabled delta or quality export additionally needs the rebased output variants. It can cost considerably more than coarse clustering with no variants retained. The report gives the number of moved reads, moved exact source subclusters, and representative comparisons. See the [Q30 200 bp measurements](coarse-dereplication-q30.md).

## `diffDNA` observations and sidecar

Every `DNAunique` embeds a base `diffDNA` by value. It owns the sparse sample occurrences and one 64-bit total, replacing the previous separate occurrence vector and total cache in `DNAunique`. The existing sample-map, count, and minimum-copy APIs delegate to this object.

| Mode | Base counts | Optional nucleotide data |
| --- | --- | --- |
| Ordinary 100% exact/prefix, or coarse without differences/qualities/reassignment | Sparse sample vector and total | Absent: no reference or variant table allocation |
| Nucleotide differences enabled | Sparse sample vector and total | Immutable reference plus one counter per exact variant |
| Final reassignment enabled | Total and lazy aggregate sample view | Immutable reference plus one sparse sample vector per exact variant; multiplicity derived from its entries |
| Selected qualities enabled | Same counts as the corresponding mode above | Each variant additionally retains selected qualities, scoring summaries, and comparison history |

The [unified-count benchmark](unified-derep-counts.md) compares memory and runtime before and after this refactor on one million Q30, 200 bp reads.

The base object requires no separate heap allocation. Sparse entries allocate only as samples occur. Nucleotide data, when enabled, is allocated separately and shared with copy on write, so a retained output snapshot stays independent if streaming resumes. Count-only merges add sample occurrences; nucleotide-aware merges also rebase and combine variants, rejecting sources whose sequence or required sample provenance is unavailable. Clearing the object clears counts and variants together.

`DNAunique::getDerepCounts()` exposes the always-present base; `getDerepMap()` refers to its aggregate sparse counts. `getDiffDNA()` retains the optional owning snapshot API and returns null when nucleotide storage is disabled. `attachDiffDNA()` replaces the complete count/variant state; it is used when constructing empty output snapshots, not to add metadata to an already-counted group. `diffDNA::readCount()` is zero for an empty object. The legacy `DNAunique::totalSum()` fallback of one for uncounted raw/mate intermediates is preserved.

With nucleotide storage enabled, a full reference string is shared with the cluster and remains valid when the representative changes. The reference need not be the sequence written to the representative FASTA/FASTQ. Per-read names are discarded. Qualities are stored only with `derepStoreQuals=1`, once per exact variant.

Without reassignment or quality retention, an observation equal to the complete reference increments its 64-bit counter along with the base sample counts. For every other observation, its processed length and changed positions/bases form an exact binary key. `diffDNA` stores that key once and increments its 64-bit counter whenever the same full nucleotide sequence occurs again. With reassignment enabled, sparse sample entries replace these variant counters, including for the reference. Hash collisions compare the complete keys. Original read order is discarded, but each distinct combination of changes and its multiplicity are preserved. Differently sized sequences remain distinct even when one is an exact prefix of another. Shorter observations retain their lengths; longer observations retain their additional suffix bases. Ambiguous symbols are retained literally. This is sufficient to reconstruct every stored processed nucleotide observation and its multiplicity.

Without quality retention, the stored sequence is the logical sequence supplied to the search: the merged read when merge-based search is used and available, otherwise R1. A search-only cap does not discard the rest of that logical sequence. With quality retention, variant identity and HQ export instead include full backing mates and their observed qualities, while scoring retains logical lengths separately. Physically cut sequence is outside the record, as are R2 observations when paired quality retention is disabled. Sample abundances are exported in the ordinary `.map`; the nucleotide-only sidecar does not carry sample IDs, even when per-variant sample counts are retained in memory for reassignment.

When enabled, `clusters.fna` or `clusters.fq` gets a binary `clusters.diff` sidecar. It includes all output groups, including those below the minimum-copy threshold. Its cluster labels match the count-bearing dereplicate/map labels. The main output cannot itself use the reserved `.diff` sidecar filename. Per-sequencing-run files use the same run suffix as the main output. A fresh snapshot rewrites the matching sidecar; representative changes do not require rewriting the in-memory observations. With final reassignment enabled, exported differences are rebuilt against each destination cluster's initial reference; live observations remain unchanged.

A file consists of consecutive self-contained blocks:

1. Eight ASCII bytes `SDMDIFF2`.
2. Unsigned varint cluster-label byte length, then label bytes.
3. Unsigned varint reference byte length, then reference nucleotide bytes.
4. Unsigned varint total observation count, then count identical to the complete reference.
5. Unsigned varint payload byte length, then payload bytes.

All integers use unsigned LEB128, at most ten bytes for 64 bits. Each distinct non-reference payload record contains its multiplicity, processed length, and number of changed bases as varints, then one `(position gap varint, nucleotide byte)` pair per change. The first position is zero-based; later gaps are relative to one past the previous changed position. Added suffix positions must be represented explicitly. Records with shortened length and no changed bases are valid.

`diffDNA::read` reads and validates one block, returning null only at clean EOF. `diffDNA::forEach` reconstructs observations one at a time and passes `(sequence_view, multiplicity)` to a callback; callers must copy a view if they want to retain it. Each exact subcluster is delivered once with its full count, including non-reference variants. Iteration order is unspecified. This API does not reconstruct FASTQ metadata.

The reader also accepts the earlier `SDMDIFF1` encoding, whose payload records had implicit multiplicity one, and combines repeated patterns on import. New single-end writes use version 2. Paired quality retention writes linked nucleotide pairs using version 3; the same reader accepts both formats.

### Selected qualities for single-end subclusters

`-derepStoreQuals 1` retains nucleotide variants and selected qualities and writes exact variants into the normal `.1.hq.fq` output (and `.2.hq.fq` for pairs). These replace the parent-only HQ records; there is no separate subcluster FASTQ copy. HQ files accumulate across sequencing runs just as in ordinary dereplication. Main representative output and `.map` describe exact R1 parents by default; coarse-parent IDs/counts require `derepCoarseClusters=1` or `derepReassign=1`. Use `-seedSubclusters 1` when reading variant HQ files for seed extension; see the [LotuS3 integration instructions](lotus3-dereplication-io.md). It is off by default and requires unmerged single-end input or complete pairs with qualities. Paired variants and synchronized mate output are described in [paired quality retention](paired-diff-quality.md). Binary deltas remain independently controlled by `-derepStoreDiffs`; `SDMDIFF2` does not contain qualities. See [quality-retention implementation and benchmarks](diff-quality-retention.md).

A quality-enabled variant map attaches a `QualityRecord` to each existing patch/count entry, including the reference under the empty patch. Its immutable `QualitySeed` holds one quality vector and compact `SeedReadStats`; the record also holds comparison history. Sample provenance uses the same sparse counts when reassignment is enabled. Snapshots share immutable quality payloads, and a replacement allocates a new payload only when the incoming observation wins. Reassignment carries the quality payload and sample counts together while rebasing the nucleotide patch. No additional full `DNA` object is retained per variant.

`diffDNA::forEachWithQualities` exposes reconstructed sequences, multiplicities, optional sample provenance, and the selected payload. `diffDNA::writeFastq` reconstructs one sequence at a time and reuses the existing DNA FASTQ writer. The normal `forEach` and FASTA export also work on these objects. Reloading a `.diff` file only restores nucleotide variants; it cannot recover the quality payload.

### Exact subcluster FASTA

With `-derepStoreDiffs 1 -derepSubclusterFasta 1`, SDM writes `clusters.subclusters.fna` in addition to the normal dereplicate representatives and binary differences. Each exact variant within each output parent gets a record such as:

```text
>read_17.sub1;size=12;
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
>read_17.sub2;size=5;
TCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
```

These counts are specific to the exact variants, and sum to the output parent's count. Export uses the fixed reference plus the stored changes, so a replaced representative does not affect reconstruction. Subclusters are local to their output parent, and their numeric suffixes are not persistent identifiers. Existing minimum-copy filtering applies to output parents, not independently to these exact variants; the export includes groups written to `.rest` too. The file uses the same run suffix as other per-run outputs.

`diffDNA::writeFasta(ostream, clusterLabel)` provides the same export directly from an object. It removes the parent's `;size=` annotation and adds each subcluster's own abundance. An existing binary sidecar can be exported later without rerunning clustering:

```cpp
std::ifstream input("clusters.diff", std::ios::binary);
std::ofstream output("clusters.subclusters.fna");
std::string label;
while (auto differences = diffDNA::read(input, label)) {
    differences->writeFasta(output, label);
}
```

FASTA export is disabled by default because it expands the distinct variants into full sequences on disk. In-memory storage retains the compact differences and counters during export.


## Validation

`tests/test_paired_diff_quality_pipeline.py` compares complete 100%/97% FASTA, map, and FASTQ records on 1,000 paired fragments (2,000 reads), ignoring record order. It covers 200–255 bp R1s, quality winners, divergent R2s, capped tails, below-cutoff dereplicates, retained qualities on/off, and 1/4/12 workers. Independent input tallies validate seeds and abundances.

`tests/TestCoarseDerep.cpp` checks the production index against an independent exhaustive first-anchor scan, including length-bin boundaries, short sequences, ambiguous bases, homopolymers, and identity thresholds from 0.1 to 100%. It also checks exact decimal boundaries, representative replacement without chaining, exact-variant counting, subcluster FASTA abundance, delta reconstruction and truncated-file rejection, concurrent sample counts, repeated snapshots, reset/finalization, option-file precedence, and search caps. Unified-count tests additionally compare per-variant and aggregate sample abundances against an independent tally, exercise rebasing merges and self-merges, validate count-only exact mode and overflow handling, and check snapshot independence and sidecar compatibility.

```sh
make WITH_HTS=0 STATIC=0 .build/portable-0/TestCoarseDerep
.build/portable-0/TestCoarseDerep
```

The standard `make test` target includes these tests alongside the existing dereplication and pipeline checks. The earlier algorithm/simulation reports describe experiments, not a guarantee of production speed or compression on biological ONT data.

## Production Q20 benchmark

For current default exact-output behavior, see the [100,000-read compact-storage Q20 benchmark](benchmarks/compact-storage-ont-q20.md). The historical coarse-parent measurements below use a different output contract.

The streaming benchmark generates 20,000 reads from 50 artificial templates in five lineages, with conserved and variable blocks and independent substitution errors at approximately 1%. The primary cases use 1,500 bp reads. Reads are generated one at a time; the complete input is never retained. This follows the earlier rough 16S-like model, not a validated biological or ONT homopolymer/indel model.

Results below are medians of three separate process runs on the same simulated input. The build used GCC 13.3, `-O3`, C++20, portable/shared-library mode on x86-64 Linux, and one processing thread. Both modes used `derep_reserve=1` to avoid attributing a large default hash reservation to exact-mode storage. Times include simulation, DNA construction, batching/clustering, and normal FASTA/HQ/map output plus the binary sidecar when enabled. They exclude FASTQ input parsing, upstream filtering, optional exact-subcluster FASTA export, and fixture cleanup. Memory is total process peak RSS, including allocator and application overhead, not a sum of live sequence bytes.

These historical measurements used coarse-parent output (now selected with `-derepCoarseClusters 1`). They do not measure the current default exact-output memory or runtime.

| Identity | Counted differences | Coarse representatives | Time (s) | Peak RSS (MiB) | Binary differences (MiB) |
| --- | --- | ---: | ---: | ---: | ---: |
| Exact, equal lengths | Off | 20,000 | 0.309 | 109.97 | — |
| 99% | Off | 19,400 | 14.375 | 119.03 | — |
| 98% | Off | 5,572 | 4.627 | 39.32 | — |
| 97.5% | Off | 983 | 0.852 | 12.76 | — |
| 97.5% | On | 983 | 1.025 | 16.70 | 2.66 |
| 97% | Off | 156 | 0.652 | 12.24 | — |
| 97% | On | 156 | 0.813 | 12.24 | 1.73 |

A separately generated case with uniformly varied terminal lengths of 1,400–1,600 bp, at 97.5% with differences, retained 950 representatives and took 2.116 s at 18.61 MiB peak RSS. Its sidecar occupied 3.97 MiB. This varies terminal length only; it does not simulate internal indels.

For this fixture, 97% with differences reduced peak RSS about ninefold while taking about 2.6 times the exact-mode runtime. Equal peak RSS with differences on and off at 97% does not mean deltas have zero memory cost: process peak RSS includes other allocations and allocator reuse. Exact repeated variants save additional space by sharing one patch/count, but Q20 full-length error patterns are usually distinct.

The 99% setting provided little compression and was much slower, with slightly higher memory than exact dereplication. A strict threshold leaves many noisy anchors to compare and index. This implementation should therefore not be presented as a cheap replacement for exact dereplication at 99% on Q20 data. Stronger compression at 97–97.5% also permits more biological variants to share a coarse representative; this benchmark does not measure taxonomic accuracy. With differences enabled, the processed nucleotide variants and their exact counts remain recoverable.

Raw measurements: [coarse-production-q20.csv](benchmarks/coarse-production-q20.csv). Reproduce an individual case with:

```sh
make WITH_HTS=0 STATIC=0 .build/portable-0/TestCoarseBenchmark
.build/portable-0/TestCoarseBenchmark exact 20000 0 0
.build/portable-0/TestCoarseBenchmark 97 20000 1 0
.build/portable-0/TestCoarseBenchmark 97.5 20000 1 100
```

The arguments are identity (`exact` or a percentage), read count, binary-difference storage (`0|1`), and terminal-length jitter (`0..200`). The input seed is fixed; each invocation starts a fresh process. The benchmark uses the production `Dereplicate` path and verifies that emitted parent abundances sum to the input count.
