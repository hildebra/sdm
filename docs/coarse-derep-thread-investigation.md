# Paired coarse dereplication: single-core CPU usage

**Parallel implementation:** the default coarse path below 100% now uses 64
independent indexes with multiple workers, without final consolidation. Retained
pairs route using both mates. Explicit global mode, reassignment, and 100%
matching retain a global index. See [implementation and measurements](benchmarks/coarse-parallel.md).
The original investigation below describes the earlier serial implementation.

**Follow-up:** the length-bin pruning issue described below is now fixed using
actual stored search-key bounds. See the [fix and measurements](benchmarks/coarse-length-bounds.md).
The original investigation and its measurements below describe the pre-fix code;
at that stage the shared coarse mutex and R1-only candidate index were unchanged.

Investigated 2026-09-14 against the current workspace's SDM 3.52 beta source.
The motivating command uses `derepIdentity=97`, `derepStoreQuals=1`,
`paired=2`, all three native merge flags disabled, `derepPerSR=0`,
`suppressOutput=1`, and `threads=12`.

**The current coarse backend serializes its expensive work.** The thread flag
does configure twelve workers, but all workers share one coarse-index mutex.
A local paired FASTQ reproduction reached approximately 105% process CPU with
twelve threads. The HPC executable, input files, and options file were not
available, so these measurements establish the behavior of the local code,
not the precise runtime breakdown of that remote run.

## Execution path

- `IO.cpp`, `separateByFile`: configures `ThreadPool` from `runtimeOpts.threads`.
  Input file groups are processed in a serial outer loop; `read_sequences`
  submits blocks from the current input to compute workers.
- `containers.cpp`, `Dereplicate` constructor: identity below 100 or enabled
  quality retention selects `CoarseDereplication`. The ordinary exact backend's
  64 hash shards are bypassed by the coarse update path.
- `Dereplicate::addDNA` collects up to 256 pending entries per producer.
  `flush_batch` holds `coarseMutex_` across **every** `coarse_->add` in that batch.
- `CoarseDereplication::add` searches candidates, verifies the two mates,
  creates anchors, selects representatives, records paired nucleotide variants,
  selects retained qualities, and updates sample counts while that mutex is held.
- `writeDereplDNA` sorts and exports clusters serially, including the two exact
  subcluster FASTQ files. `suppressOutput=1` suppresses primary ordinary output;
  it does not disable the dereplication outputs. Minimum-copy filtering occurs
  after admission and does not reduce the cost of searching incoming reads.

The current index cannot safely accept simultaneous searches simply by removing
the mutex. Its seed, posting-list, and heap scratch vectors and statistics are
shared, and insertion mutates the index and cluster storage.

## Additional search bottlenecks

### Broad length bins disable seed pruning for common amplicon lengths

In `CoarseDereplication::search`, the length-bin width is 128. The code uses the
bin's nominal lower bound as the usable seed prefix and its nominal upper bound
to bound the permitted mismatch count, even if every anchor in the bin has the
same actual length.

For a 250 bp query and anchors in the 128–255 bp bin at 97%:

- Only seeds entirely inside the first 128 bases are eligible: at most six
  disjoint 20-mers, before homopolymer exclusions.
- The mismatch budget is `floor(250 * 0.03) = 7`; the completeness guarantee
  therefore requires eight eligible seeds.
- Six cannot satisfy eight, so the index scans the bin in anchor order until
  a match is found or all anchors are rejected.

The same forced fallback occurs for equal-length queries and anchors anywhere
from 200 through 255 bp at 97%. This does not imply every possible dataset or
every length bin falls back. In both 250 bp fixtures below, 199,999 of 200,000
queries used fallback; the first query encounters no populated bucket.

Using actual minimum/maximum anchor lengths within each bucket, or refining
populated bins, is a candidate optimization that can keep exhaustive matching
semantics. It needs the existing variable-length/oracle tests: shorter anchors
must still be searchable, and insufficient eligible seeds must still fall back.

### Paired matching indexes only R1

`indexGroup` builds seed postings from R1. `search` checks R2 only when verifying
each candidate. With many similar R1 anchors distinguished by R2, even effective
R1 seed pruning can leave long candidate lists. A mate-aware candidate index
could avoid many of those comparisons while still requiring both mates to pass.

This matters separately from the length-bin fallback. Fixing bin bounds alone
will not resolve a workload where almost every R1 passes against every anchor.

## Measurements

The new [CLI benchmark](../.claude/archive/benchmarks/benchmark_paired_coarse_threads.py) reuses the
existing quality pipeline's FASTQ parser. Each process reads 200,000 uncompressed
pairs, each mate 250 bp, with Q30 qualities. It uses the motivating command's
dereplication, output suppression, minimum-copy, merging, and threading settings,
with a synthetic mapping file and permissive filtering that accepts all pairs.
Input generation and output validation are outside the timed process.

There are two deliberate synthetic workloads. The first has 50 unrelated parent
pairs. The second has 1,000 parents sharing R1 and differing in R2, stressing
candidate verification. Each is sampled from 10,000 generated variants with
zero to three substitutions per mate. These are diagnostic fixtures, not a
model of the unavailable Apong community or its abundance distribution.

Medians of three separate processes; CPU is `(user + system time) / wall time`,
where 100% means one fully used logical CPU:

| Workload | Threads | Wall time | Average process CPU |
| --- | ---: | ---: | ---: |
| 50 parent pairs | 1 | 1.490 s | 100.0% |
| 50 parent pairs | 4 | 1.021 s | 173.2% |
| 50 parent pairs | 12 | 1.008 s | 177.4% |
| 1,000 parents, shared R1 | 1 | 15.769 s | 99.9% |
| 1,000 parents, shared R1 | 4 | 14.058 s | 105.2% |
| 1,000 parents, shared R1 | 12 | 16.273 s | 105.1% |

Every run preserved all exact pair abundances, synchronized mate names/counts,
the Q30 quality witnesses, and the sample-map total of 200,000 fragments.
The two fixtures retained 50/1,000 coarse clusters and 7,524/7,750 exact pairs.
Representative names and greedy admission order can vary across worker schedules;
byte-identical representative outputs are not claimed.

For direct attribution, a separate temporary instrumented build measured the
mutex wait/hold duration per batch and the duration of `find(search, mateKey)`.
It retained the production algorithm and used timers only in temporary source
copies. A twelve-thread shared-R1 run measured:

| Measurement | Value |
| --- | ---: |
| Process wall time | 15.159 s |
| Average process CPU | 104.85% |
| Time holding the coarse mutex | 15.024 s |
| Candidate search time inside that lock | 14.529 s |
| Summed mutex waiting across workers | 160.260 thread-seconds |
| Candidate comparisons | 99,561,619 |

Search consumed approximately 96% of wall time. The other work inside the lock,
including quality/variant/count updates, used approximately 0.495 seconds in
this fixture. Summed waiting is simultaneous waiting by different workers; it
must not be added to elapsed time. Timer overhead and scheduling mean these
instrumented values should not replace the uninstrumented medians above.

Raw data: [ordinary fixture](benchmarks/paired-coarse-threads.csv),
[shared-R1 fixture](benchmarks/paired-coarse-threads-shared-r1.csv), and
[instrumented processes](benchmarks/paired-coarse-threads-instrumented.csv).
The instrumented [twelve-thread log](benchmarks/paired-coarse-threads-instrumented-1-12.log)
contains both coarse diagnostics and timer totals. The
[instrumentation patch](benchmarks/paired-coarse-thread-instrumentation.patch)
records the temporary modifications; it is not applied to the production source.

## Implementation priorities

1. Improve candidate pruning: tighter actual length bounds and, for retained
   pairs, candidate selection that can exploit R2. This attacks wasted work
   whether one or twelve threads are requested.
2. Parallelize searches against a frozen index with per-worker scratch. Existing
   immutable anchors and verification logic are useful building blocks. Commit
   assignments and new anchors in controlled order; a read unmatched in the
   snapshot must be checked against newly admitted anchors before it seeds a
   cluster or is rejected. Existing snapshot matches cannot be superseded by
   later-created anchors under first-created-match semantics.
3. Separate index publication from cluster updates. Different clusters can own
   independent representative/variant/count updates, reusing the existing
   selection and `diffDNA` helpers. An abundant single cluster can still be a
   bottleneck and needs measurement rather than assumed linear scaling.

Independent per-thread coarse indexes followed by representative-only merging
are not an equivalent shortcut. The earlier
[partitioning investigation](derep-parallel-partitioning.md) demonstrated members
falling outside their final anchor's threshold. Exact-key sharding alone also
does not solve 97% matching because compatible nonidentical reads may hash to
different shards.

Changing `derepIdentity` to 100 while retaining `derepStoreQuals=1` still selects
the current coarse backend. Disabling quality retention changes paired matching
and removes the requested exact-pair quality output; it is not a semantics-neutral
threading fix. `derepPerSR` and search caps likewise change the operation's scope
or matching behavior rather than removing the global mutex.

## Proposed 20-mer bins and optional cluster merging

**Updated requirement:** in the dereplicated-read workflow, the 97% groups are
internal storage, and assigning a variant to a different coarse group is
acceptable. A future binned backend should therefore keep its bins independent
by default and spend no time consolidating them at the end. Stable routing must
keep each exact processed sequence/pair in one owner bin, with its fragment
counts and selected qualities preserved. Final merging is only for an explicit
request for consolidated clusters; the verification discussion below applies
only to that optional operation. No binned backend or new merge option is
implemented by the length-bound fix.

The proposed reuse of the existing dereplication shards is feasible as an
alternative coarse backend. Each logical shard can own a `CoarseDereplication`
instance; producer batches can reuse the current grouping and per-shard locking.
The independent index objects also own independent search scratch. Keep more
logical shards than workers, such as the existing 64, rather than tying the
number of shards to `threads`.

`collectSeeds` currently generates multiple position-specific disjoint 20-mers,
not one guaranteed routing key. A stable choice from those words can route a
read to one owner. If word frequencies influence that choice, use a frozen
sampling table; the growing index's current posting lengths would change routing
as the run proceeds. No eligible word needs an explicit fallback. For the paired
case, routing based only on conserved R1 can leave one very large shard; using
R2 or a paired signature deserves comparison with full-pair hashing.

A fixed routing rule keeps exact duplicates together, but does not guarantee
that approximate matches share an owner. This is acceptable for the default
independent-bin workflow; an explicitly requested consolidation pass can handle
cross-bin groups afterward. The existing `TestCoarsePartition`
prototype already measures minimum-hash and rarity-based 20-mer routing, local
coarse indexes, and representative-only merging. Its earlier results show both
hot bins and identity violations from the unverified merge, so it is the natural
place to extend the experiment rather than write a separate routing system.

If consolidated clusters are explicitly requested with retained pairs, a
conservative merger could:

1. Find candidate destinations for a local cluster's immutable anchor.
2. Reuse `diffDNA::forEachPair` to verify every retained exact pair against a
   destination's immutable anchors, with the configured cap/identity on each mate.
3. If every pair qualifies, reuse `diffDNA::merge` to rebase variants and combine
   paired quality records and aggregate sample counts, plus the existing
   representative-selection helpers to choose output witnesses.
4. If any pair does not qualify, keep the local cluster as a separate cluster
   (or try another destination). Do not transfer it solely on anchor similarity.

This conservative whole-cluster merge does **not** require per-variant sample
provenance: all members and their aggregate sample vector move together.
`derepStoreQuals=1` already retains the exact paired sequences needed for checking.
Splitting a local cluster between multiple final destinations is a more complete
reconciliation strategy, but then each exact variant needs its sparse sample
vector. `diffDNA` already supports this through `trackSamples`; the current
command does not enable it by default. The traversal/rebasing operations in
`CoarseDereplication::refine` can be reused, though its present algorithm only
reassigns against a fixed representative set and is not a bucket merger.

The final verification should process distinct retained variants rather than
replaying all raw observations. In this investigation's shared-R1 fixture, that
means approximately 7,750 exact pairs instead of 200,000 input pairs if routing
keeps duplicates together. This is an opportunity, not a measured merged-backend
speedup; target search, failed merge attempts, load skew, and extra local anchors
still cost time and memory. Quality-selector history also makes a new processing
order potentially observable. The backend can preserve fragment/sample totals,
paired quality witnesses, and the anchor-radius rule without promising identical
clusters or representatives to the existing greedy global order.

## Reproduction and limits

Measurements used GCC 13.3, `-O3 -std=c++20`, a static portable build without
HTSlib, and local temporary files under `/tmp` on an eight-logical-CPU WSL host.
The twelve-thread case therefore oversubscribes this host; the single global
lock and approximately one-core utilization remain the relevant observations.
CPU frequency and background activity were not controlled. Ordinary gzip input
can add a further sequential decompression bottleneck; indexed BGZF follows a
different path in `InputFileStream.cpp`. Neither compression nor remote storage
was involved in these reproductions.

```sh
make WITH_HTS=0 STATIC=1 -j4
python3 .claude/archive/benchmarks/benchmark_paired_coarse_threads.py --binary ./sdm \
    --output docs/benchmarks/paired-coarse-threads.csv
python3 .claude/archive/benchmarks/benchmark_paired_coarse_threads.py --binary ./sdm \
    --templates 1000 --shared-r1 \
    --output docs/benchmarks/paired-coarse-threads-shared-r1.csv
```

Only the benchmark and investigation artifacts were added. Production source
and the existing `sdm` executable were not changed by this investigation.
