# Dereplication partitioning and parallelization investigation

Investigated on 2026-09-13. The recommendation is to reuse SDM's existing full-key hash shards for fixed-length 100% matching, including the newer `diffDNA` mode. A rare-20-mer partition is unnecessary for exact equality. Approximate clustering needs a separate reconciliation design: merging local representatives alone can violate the configured identity threshold for their members.

This work adds an experimental benchmark and a shared test simulator. It does not add a production partitioning option.

## What SDM already does

The ordinary exact-key backend already implements the proposed general arrangement:

1. Build the effective dereplication key, respecting logical sequence length, merged/R1 search selection, and search/truncation caps.
2. Assign it with `robin_hood::hash<string>{}(srchSeq) % 64`.
3. Collect up to 256 pending reads per producer and group that batch by shard.
4. Acquire only the relevant shard mutex while updating its keys, sample counts, qualities, and representatives.

Identical keys always enter the same shard. Hash collisions merely put unrelated keys in the same shard; the hash table still checks full key equality. Sequencing errors do not undermine this property: a changed sequence is a different exact key. No sequence-clustering merge between these shards is necessary; output gathering and sorting remain.

The 64 shards are independent of the worker count. Retaining more logical shards than workers is generally more flexible than changing the partition count whenever `threads` changes. A very abundant exact key still occupies one shard, so adding shards cannot by itself parallelize all updates to that one key.

Code locations: `Dereplicate::addDNA`, `flush_batch`, and the shard allocation in [containers.cpp](../containers.cpp); batch/shard declarations in [containers.h](../containers.h).

The coarse backend takes a different path. It holds `coarseMutex_` for an entire producer batch, including seed generation, candidate verification, representative selection, sample counting, and optional differences. Its shared search scratch vectors also make simultaneous calls to the current index unsafe. Importantly, `derepIdentity=100` with `derepStoreDiffs=1` or `derepReassign=1` currently selects this coarse backend too.

## Actual 100% mode: measured opportunity

These cases call the production `Dereplicate` API, using 200,000 pre-generated Q30 reads of 200 bp. Both paths retained exactly 22,245 nucleotide sequences and counted every input observation. Values are medians of three processes.

| Current backend | 1 worker | 4 workers | 8 workers |
| --- | ---: | ---: | ---: |
| Ordinary exact keys, existing 64 shards | 0.376 s | 0.200 s | 0.194 s |
| 100% with `diffDNA`, current global coarse lock | 3.811 s | 4.865 s | 5.899 s |

Times include DNA construction, admission and output writing, but exclude input generation and FASTQ parsing. The delta cases additionally write `.diff`. The slow part is admission: the delta path performs approximately **20.6 million bounded candidate comparisons** for these 200,000 reads. Its output stage takes only about 0.08–0.12 seconds. This identifies an avoidable candidate-search cost as well as serialization; it does not isolate mutex overhead from cache, scheduling, allocation, and representative-order effects.

The practical first change should therefore bypass the coarse k-mer candidate search when equal clustering-key lengths and 100% identity are guaranteed, and attach `diffDNA` to records owned by the existing hash shards. Hash the effective key, not the raw full read: bases beyond a truncation/search cap must not split identical clustering outputs.

This should reuse the existing batch routing and mutexes while preserving the newer mode's behavior: one representative across merge offsets, selected-read qualities rather than allele-averaged qualities, immutable difference references, exact full processed variant counts, and consistent sample/output snapshots. Simply falling through to the legacy update helper would not preserve all of those details. Startup already has `Filters::fixedDereplicationLength` and search-cap checks to reuse.

The startup guarantee must cover the entries that actually reach this backend. In particular, shorter failed reads can currently join a passing coarse prefix cluster. Such entries need a correct slow path or an explicit equality-only mode; a minimum length guaranteed only for passing reads is not sufficient to silently discard their existing prefix behavior.

An additional optimization worth evaluating independently is an exact lookup of immutable coarse anchors. Repeating an anchor's full key can return its cluster immediately: that anchor was created only after it failed to match every older anchor, and later anchors cannot supersede it under first-created-match semantics. This could help repeated anchors even in approximate mode, though it adds an index and does not remove the need to synchronize updates.

## 100% prefix matching is a separate routing case

Two reads can be 100% identical across the shorter length without having identical full keys. Hashing the complete short and long strings independently does not ensure that they enter the same block. Their minimum/rarest 20-mer can also differ because the longer sequence contributes additional words, even with no sequencing errors.

A fixed first-k-base routing key is safe if every compared key has at least k bases: a compatible prefix pair shares those bases. Shorter entries require cross-block handling, and conserved amplicon starts can put most reads in one block. Homopolymer/ambiguity exclusions would need a common fallback for a shared routing prefix; they must not make compatible reads choose different later words.

The existing default 100% prefix backend already collects exact keys using local shard locks and performs the prefix reconciliation at output. That arrangement remains useful. It should not be replaced by independent whole-key prefix clusters without retaining the final reconciliation step and its handling of shorter/failed reads.

## Approximate partition prototype

The isolated [TestCoarsePartition.cpp](../tests/TestCoarsePartition.cpp) reuses the real `CoarseDereplication` index and representative/counting methods. It compares a shared index locked in batches of 256 with independently owned bucket indexes. Bucket tasks are assigned to workers; a bucket has one owner at a time and needs no internal mutex in this prototype. The approximate cases retain no `diffDNA`.

After local clustering, the prototype sorts local initial anchors by their first input position and clusters those anchors globally. It then checks every original member against its resulting global anchor. This deliberately measures the proposed inexpensive representative-only merge and exposes when it is insufficient. It does not implement the repair pass that would be needed in production.

Routing methods:

| Name | Experimental rule |
| --- | --- |
| `roundrobin` | Read index modulo bucket count; ideal read-count balance, no sequence locality. |
| `fullhash` | Hash the complete sequence into a fixed number of buckets. Exact duplicates stay together. |
| `min20` | Select the minimum-hash eligible disjoint 20-mer plus its position, then hash that signature into a bucket. This is one whole-read signature, not a sliding-window minimizer implementation. |
| `min8` | Select using the first eight bases of eligible 20-mers, testing a coarser signature. Verification still uses the ordinary 20-mer index and full identity check. |
| `rare20` | Count eligible position-specific 20-mers in the first 4,096 reads, freeze the table, and choose the lowest-frequency word, breaking ties by hash. Unseen words have frequency zero. |
| `supported20` | The same frozen rarity table, but only words seen at least three times can supply the signature. |

All search 20-mers exclude runs of four identical bases and ambiguous bases. Reads without an eligible routing seed use one fallback bucket. The rarity table is frozen so the same exact sequence does not change destination as the run progresses. The current coarse index's changing posting-list frequencies are useful for candidate selection but should not directly define a permanent owner bucket.

### Q20, 1,500 bp, 97% clustering

20,000 reads from the earlier rough 50-template model; medians of three processes. Partitioned rows use eight workers and eight buckets.

| Method | Route + cluster + representative merge | Largest bucket | Local clusters | Final clusters | Members outside final 97% anchor radius |
| --- | ---: | ---: | ---: | ---: | ---: |
| Shared index, one worker | 0.707 s | 100% | 156 | 156 | 0% |
| Shared index, eight workers | 0.719 s | 100% | 125 | 125 | 0% |
| Round-robin buckets | 0.275 s | 12.5% | 492 | 106 | 14.39% |
| Full-key hash | 0.222 s | 13.0% | 480 | 86 | 14.47% |
| Minimum-hash 20-mer | 0.766 s | 77.2% | 463 | 127 | 3.81% |
| Coarser eight-base signature | 0.589 s | 72.5% | 329 | 108 | 5.04% |
| Rarest 20-mer | 0.507 s | 12.8% | 496 | 95 | 13.31% |
| Supported rare 20-mer | 0.431 s | 13.7% | 489 | 90 | 15.72% |

For the shared eight-worker run, summed mutex waiting was approximately 4.71 thread-seconds during 0.72 seconds of admission. This is cumulative waiting across workers, not an extra 4.71 seconds of wall time. Different worker schedules change greedy admission order, which explains the different shared-index cluster counts and prevents interpreting timings as a strictly identical sequence of comparisons.

Rare routing spread reads well but did not keep noisy copies together: the fraction of same-template pairs in the same bucket was about 12.5%, essentially chance for eight buckets. Requiring three observations was insufficient here. A recurring sequencing-error word in a conserved region can exceed that support threshold. The minimum-hash 20-mer retained more locality, but one conserved signature dominated the load. Making the signature shorter did not solve that imbalance.

Increasing to 64 logical buckets did not remove the hot minimum-hash signature: the largest bucket still held 74.7% of reads. Supported-rare routing with 64 buckets balanced the work, but retained 2,883 local clusters and increased representative-merge time to about 0.13 seconds. More buckets trade locality and temporary cluster memory for scheduling flexibility; they cannot split one hot signature without an additional rule.

### Q30, 200 bp

With 100,000 reads at 97%, shared single-worker admission took 0.213 seconds and retained 23 clusters. Eight full-hash buckets plus representative-only merge took 0.057 seconds, temporarily retained 121 clusters, and finished with 22; 0.20% of members exceeded the final anchor radius. Minimum-hash 20-mer routing took 0.119 seconds, with 29.9% of reads in its largest bucket. The coarser signature put 50.9% in its largest bucket.

At 99%, the shared baseline retained 814 clusters in 0.192 seconds. Eight full-hash buckets plus representative-only merge took 0.095 seconds but put 4.16% of reads outside the final 99% radius; minimum-hash and rare routing produced 7.89% and 8.06%, respectively. Thus even the short, high-quality case needs a proper final membership check if the configured radius is to remain meaningful.

These partition timings are **before correcting invalid memberships**. Their smaller final cluster counts are not evidence of equivalent compression accuracy. All observations are counted once, but conserving counts alone does not preserve the identity contract.

## Why representative-only merging is insufficient

At 200 bp, 97% allows six substitutions. Suppose B differs from A at six positions, and C differs from B at another six positions. A and B match at 97%; B and C also match at 97%; A and C match at 94%. A local cluster containing B and C cannot be transferred wholesale to A merely because B matches A.

For equal lengths, keeping a local maximum mismatch radius R gives a sufficient condition for a safe wholesale merge: `distance(local_anchor, target_anchor) + R <= allowed_mismatches`. If this condition fails, the merge may still be valid, but the summary cannot establish it. Variable lengths require additional care because minimum-length comparisons leave suffixes unconstrained.

With counted exact variants and sample provenance retained, a final pass could reconstruct and check variants individually, then decide which local anchors must remain. Without that information, options include conservative radius-certified merges, retaining unresolved local groups, or rereading the input. Existing `derepReassign` is useful code to reuse, but currently fixes the representative set and does not by itself consolidate buckets or resolve removal of a representative. Recomputing a single global clustering from all raw reads would recover the original order-dependent result at the cost of repeating the serialized work.

There is also a general limitation: for a fixed length admitting at least one substitution, the graph connecting sequences by one-base changes is connected. A deterministic rule that always places every allowed near-match in the same single bucket therefore collapses to one bucket. Approximate partitioning needs cross-bucket information, a repair phase, or an explicitly different clustering contract. This limitation does not apply to exact sequence equality.

## Literature and integration direction

The general binning idea is well established for exact counting. [KMC3](https://academic.oup.com/bioinformatics/article/33/17/2759/3796399) partitions by signatures into many bins and uses sampling to improve bin balance. [Gerbil](https://pmc.ncbi.nlm.nih.gov/articles/PMC5374613/) uses minimizer-based partitioning and balances work between hashers. These support using fixed logical blocks and managing skew; their exact k-mer counting guarantees do not imply that approximate whole reads stay together.

[Mash's MinHash formulation](https://mash.readthedocs.io/en/latest/distances.html) relates shared signatures to k-mer-set Jaccard similarity. Under a simplified independent-substitution model with mostly unique 20-mers, `J ≈ p^20 / (2 - p^20)`, where p is pairwise nucleotide identity. An ideal one-MinHash routing signature agrees with probability J. With an independent hash into eight buckets, the approximate co-location probability is `J + (1-J)/8`: about 73%, 56%, and 45% for 99%, 98%, and 97% sequence identity. These are explanatory approximations, not guarantees or measured sensitivities for this disjoint-seed prototype.

Recommended sequence of work:

1. **Extend the existing exact-key shards to fixed-length 100% `diffDNA` mode.** Preserve the current representative, delta, failed-read, and output semantics; avoid the coarse seed search where equality suffices.
2. **Keep variable-length prefix reconciliation explicit.** Whole-key collection can stay sharded, but compatible unequal lengths still need reconciliation or a safe shared-prefix routing rule.
3. **For approximate mode, first test parallel searches against a frozen index with per-worker scratch.** Apply results in a controlled commit phase and resolve reads not matched by that snapshot against newly created anchors. This can parallelize expensive verification without multiplying retained cluster sets. New-anchor publication and representative/count updates need separate treatment; an all-new-cluster workload will still have a serial component.
4. **Treat independent approximate buckets as an optional alternative requiring measured repair.** Use bounded queues and more logical blocks than workers, monitor skew, and budget for local-cluster duplication and member-level reconciliation. Rare words alone did not provide the desired locality in the Q20 experiment.

## Reproduction and measurement limits

The benchmark uses the same [CoarseSimulation.h](../tests/CoarseSimulation.h) as the earlier production benchmark; its Q30 and Q20 cluster counts and error rates were checked after extracting the shared helper. Machine: Intel Core Ultra 7 258V, eight logical CPUs, Linux under WSL/Microsoft virtualization, GCC 13.3, `-O3 -std=c++20`, `WITH_HTS=0 STATIC=0`. CPU placement, frequency, and background load were not controlled.

Reads are preloaded solely to apply identical inputs to the partitions and to inspect every membership after merging. Prototype peak RSS therefore includes the full input and validation structures and must not be compared with the earlier streaming memory measurements. Routing is performed serially in this experiment, including the optional 4,096-read frequency sample; a production frontend could parallelize signature calculation after freezing the routing rule. Approximate timings exclude final FASTA/map/delta output and the exhaustive validation, while the production exact rows include output. Compare within these categories, not across them as equivalent end-to-end pipelines.

```sh
make WITH_HTS=0 STATIC=0 BUILD_DIR=.build/coarse-dev-20260912 -j4 \
    .build/coarse-dev-20260912/TestCoarsePartition
python3 tests/benchmark_coarse_partition.py \
    --binary .build/coarse-dev-20260912/TestCoarsePartition
```

The runner covers 45 configurations, each repeated three times. All 135 runs checked observation accounting; shared-index baselines also checked that all members meet their own immutable anchor's threshold. Raw timings, bucket locality, cluster inflation, and membership violations are in [coarse-partition.csv](benchmarks/coarse-partition.csv).
