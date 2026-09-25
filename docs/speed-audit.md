# Speed audit: where sdm still spends its time

Audited 2026-09-24 on the working tree after the multithreading work in
[multithreading-fixes.md](multithreading-fixes.md). No production source was
changed. Proposed changes were measured in scratch copies of the build, or as
stand-alone microbenchmarks against the build's own object files.

Items 1–5 and 8 are implemented and re-benchmarked in
[speed-fixes.md](speed-fixes.md).

**Summary.** The largest remaining costs are per-read computation and compression,
not locking. In order of expected gain:

| # | Change | Helps | Measured | Expected end to end | Output |
| --- | --- | --- | --- | --- | --- |
| 1 | Linear-time binomial error filter | every run with the binomial filter on (LotuS default) | function 2.2× (two-term) to 5.4× (column-wise) | 1.15× measured with a 2× prototype; about 1.35× at 5× | identical |
| 2 | libdeflate for BGZF output (and BGZF input) | `.gz` output, BGZF input | deflate + CRC 1.8×, 17% smaller files; inflate 3.0× | 1.5× (1 thread) to 1.8× (8 threads) on BGZF-bound runs; about 1.2× for LotuS with `.gz` output | same decompressed content |
| 3 | SIMD candidate verification, fewer candidates | 97%, global and 100% coarse modes | verification loop 5.8× | about 1.2–1.35× (the search is about a third of CPU) | identical |
| 4 | Free reads on the thread that allocated them (or use a thread-caching allocator) | partitioned derep at 8 threads | malloc lock waits are 11–19% of worker samples | 5–15% at 8 threads | identical |
| 5 | Parallel, cheaper sort of the derep output | derep output tail | the sort is 15–30% of the main thread's samples at 8 threads | 0.1–0.4 s per run | identical |
| 6 | Cheaper per-read passes | all modes | quality characters 1.8×, average quality 1.5× | up to 1.1× on plain rewrites, a few % elsewhere | identical |
| 7 | Reuse the merged read (demultiplexing and derep both merge) | LotuS merged modes | merging is 6–10% of CPU | 3–5% | identical |
| 8 | Reader thread: batch wake-ups; longer term, workers parse raw chunks | plain rewrites at 8 threads (reader 73% busy) | – | up to about 1.3× for those | identical |
| 9 | Profile-guided optimization | everything | `-march=x86-64-v3` and LTO: none; PGO 1.01–1.12× | up to about 1.1× (in-sample training) | identical |

"Identical" was checked byte for byte against the current build wherever a
prototype exists (items 1, 3's comparison kernel, 6, 9). For the other items it
follows from the design (same arithmetic, same order).

## Method

- **Profiles.** All 11 scenarios from the multithreading audit, at 1 and 8
  threads, with the `LD_PRELOAD` sampling profiler (2 ms samples, no lock-hold
  tracking) and the final build.
- **Microbenchmarks** (`speedbench.cpp`):
  - The current code and a proposed variant run on the same 100,000 real apong
    reads.
  - Batches alternate between the two within one process, so clock changes hit both
    equally.
  - Every result of the variant is checked against the current code.
- **End-to-end A/B** (`concab.sh`):
  - Both binaries run at the same time, each pinned to its own vCPU, with the pins
    swapped every repeat.
  - The host was throttled for most of the audit (single runs up to 4× slower than
    earlier in the day), and sequential A/B runs were unusable.
  - Running both at once exposes them to the same host state; the reported value is
    the median of the per-pair ratios.
- Hosts, datasets and scenario definitions are as in
  [multithreading-audit.md](multithreading-audit.md).

## Where the time goes

Share of single-thread CPU (profile self time) for the functions that matter:

| Function | rewrite | rewrite `.gz` | pe_lotus | lotus `.gz` out | se_demux | apong_fa | apong_fq | 97% coarse | global coarse | 100% coarse |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| zlib (deflate/inflate, unsymbolized) + `crc32` | 1 | **76** | – | **36** | 1 | – | – | – | – | – |
| binomial filter (`sum_of_binomials`, `binomialFilter`, `prob_j_errors`) | – | – | – | – | – | **32** | **33** | 10 | 9 | 15 |
| coarse candidate verification (`percentIdentity`, `std::min`, `mateKey`) | – | – | – | – | – | – | – | **33** | **35** | 14 |
| accumulated error (`getAccumError`) | 14 | 2 | 6 | 3 | 6 | 2 | 2 | 2 | 1 | 3 |
| FASTQ quality characters (`buildFastqQualityString`) | 14 | 3 | 2 | 1 | – | – | – | – | – | – |
| average quality, `seal`, record construction | 21 | 4 | 6 | 3 | 8 | 4 | 4 | 2 | 2 | 5 |
| sequence checks (`scanSequenceChecks`, `numNonCanonicalDNA`) | – | – | 11 | 7 | 12 | 3 | 3 | 3 | 3 | 5 |
| merging (`percentIdentity`, `merge`, `getMergedQual`) | – | – | 11 | 6 | – | 9 | 8 | – | – | – |
| malloc (internals, arena, `malloc`/`free`) | 2 | 1 | 12 | 8 | 16 | 9 | 8 | 11 | 12 | 16 |

At 8 threads three limits remain:

- **Worker CPU** in the filtering, merging and derep modes: workers are 73–83%
  busy, so the per-read items above translate directly.
- **Compression** for BGZF output: the compressor pool is busy, and `crc32` alone is
  30% of its time.
- **The reader** for plain rewrites: it is 73% busy, and the output writer thread
  60%.
  - 35–50% of the reader's running time is line splitting and copying (`memchr`,
    `memcpy`, string assign).
  - 5–28% is waking a worker for every block (`pthread_cond_signal`).
  - 10–17% is zero-filling new buffers.

Glibc malloc lock waits are 11–19% of worker samples at 8 threads in `pe_lotus`,
`se_demux` and `apong_coarse` (1% in the global mode, whose reads already go back to
their producer).

## Findings and proposals

### 1. Binomial error filter: quadratic where linear suffices

`DNA::binomialFilter` computes the distribution of the number of errors in a read
(a Poisson-binomial distribution) row by row, one row per error count.

- **The waste:** each cell calls `sum_of_binomials`, which loops over all j+1
  smaller error counts and calls `prob_j_errors(p, i, n = 1)` for each. For a
  single Bernoulli trial that probability is exactly 0 for i ≥ 2, so only two
  terms ever contribute. The filter is O(L·E²) instead of O(L·E).
- **The cost:** with the LotuS options (`BinErrorModelMaxExpError 2.7`, so up to 5
  rows) this is 32% of single-thread CPU in LotuS-style runs, and 10–15% in the
  coarse modes.

Two exact rewrites:

- **Two-term sum.** Stop the loop at `min(j, 1)` when `n == 1`. The dropped terms add
  exactly `0.0f`, so every float result is unchanged. Measured: function 2.2×;
  end to end (prototype in `DNA.cpp`) 1.05× (`apong_fa`) and 1.03× (`apong_fq`).
- **Column-wise.** Row e at position k needs only rows e and e−1 at k−1, so all
  rows can advance together along the read.
  - Two small columns replace the per-read, zero-filled rows × length table.
  - Every value is still computed by the same float operations, and the
    early-stopping row is chosen afterwards from the same sums.
  - Stand-alone it measured 5.4× faster; the prototype compiled inside `DNA.cpp`
    was only about 2× faster, and gave 1.15× (`apong_fa`) and 1.14× (`apong_fq`)
    end to end.
  - Inside `DNA.cpp` the compiler turned the per-position `memcpy` of the column
    into a library call. Swapping two column pointers fixed half of the gap.
  - A careful version (fixed-size columns, no calls in the loop) should reach the
    stand-alone speed: about 1.35× end to end on LotuS-style runs.

The prototypes' outputs are byte-identical to the current build on `apong_fa`,
`apong_fq` and `apong_coarse`. That holds for builds without FMA. A build for an
FMA-capable target needs `-ffp-contract=off` (see [Build options](#build-options)),
or the compiler may fuse multiply-adds and change the last bits of these
probabilities.

### 2. Compression: libdeflate for BGZF

BGZF blocks are compressed with zlib (`deflateInit2`, level 1, raw deflate) plus
zlib's `crc32`. On a 64 MB FASTQ sample in 65,280-byte BGZF blocks:

| Operation | zlib 1.3 | libdeflate | Speedup |
| --- | ---: | ---: | ---: |
| Level-1 deflate + `crc32` | 1,950 ms | 1,081 ms | 1.8× |
| Compressed size | 25.3% | 21.0% | 17% smaller |
| `crc32` alone | 97 ms | 26 ms | 3.7× |
| Inflate (BGZF input) | 1,334 ms | 451 ms | 3.0× |

- **Where it applies:** BGZF blocks are small and independent, which suits
  libdeflate's whole-buffer API exactly. `ParallelBgzfOstream` and the BGZF input
  path can switch to it without changing the container format.
- **Expected:** about 1.5× for BGZF-bound runs at 1 thread (zlib is 76% there),
  and up to 1.8× at 8 threads where compression is the limit. About 1.2× for LotuS
  runs with `.gz` output.
- **Output bytes** differ, but the files stay valid BGZF with identical
  decompressed content.
- **Dependency:** libdeflate is MIT-licensed and widely packaged; only the runtime
  library is installed on this host.
- **Plain (non-BGZF) gzip input** is one stream that must be inflated serially;
  libdeflate cannot stream. zlib-ng in compatibility mode (drop-in for zlib) or
  ISA-L's igzip would speed that up. Not measured, as neither is installed.

### 3. Coarse derep: candidate verification

In the coarse modes about a third of CPU is candidate verification in
`CoarseDereplication::search`: 8.8 million bounded comparisons for 26,955 candidate
searches on apong, about 330 per search. Amplicon sequences share conserved seeds,
so the seed filter admits many anchors.

- **SIMD comparison.** `ReadMerger::percentIdentity` compares byte by byte with an
  early exit, which the compiler cannot vectorize.
  - An AVX2 version counts mismatches 32 bytes at a time and exits after the same
    block.
  - It returned the same value for every pair tested (620k pairs, 3% budget) and
    ran 5.8× faster.
  - It needs runtime dispatch or an AVX2 build (item 9).
- **Fewer comparisons.** An anchor within E mismatches must share at least k−E of
  the query's k disjoint seeds.
  - The search currently merges only the E+1 shortest posting lists, then verifies
    every id they contain.
  - Counting seed hits per id over more lists and skipping ids below k−E would
    reject most false candidates before any byte comparison.
  - Ids are still visited in increasing order, so the first-created match stays the
    winner.
- **Small items:**
  - `search` constructs a full `ReadMerger` (seed map and buffers) per call only to
    use `percentIdentity`, which could be static.
  - Each check dereferences group → shared reference string → data; storing the
    anchor keys contiguously would avoid those cache misses.

### 4. Memory allocation across threads

At 8 threads, glibc arena lock waits are 11–19% of worker samples in the partitioned
derep modes (`pe_lotus`, `se_demux`, `apong_coarse`).

- **Cause:** `PartitionWorkPool` runs each shard's task on whichever worker is
  waiting. Reads freed there, such as a representative replaced by a better one or
  entries released after a block is applied, were allocated by another worker.
  Every such free takes the owning arena's lock.
- **Fix used for the global mode:** it hands committed entries and applied updates
  back to their producer, and its lock waits are down to 1%. The same scheme would
  apply to `apply_blocks_ordered`, and to the entries a representative replacement
  drops.
- **Alternative:** a thread-caching allocator (mimalloc or jemalloc) handles
  cross-thread frees without a global lock, would help every mode, and needs no
  bookkeeping in sdm. Neither is installed here, so this is not measured.
- **Tunable only:** `GLIBC_TUNABLES=glibc.malloc.tcache_count=1024` left peak RSS
  unchanged. Its speed effect was lost in the host noise (0.8–1.1×).

### 5. Serial tail: sorting the derep output

`writeDereplDNA` sorts all output groups on the main thread by count, then key.

- **Why it's slow:** many groups tie on count, and their keys are similar amplicons
  with long common prefixes, so most comparisons are long string compares.
- **Cost:** at 8 threads the sort is 15% (`pe_lotus`) to 30% (`se_demux`) of the
  main thread's samples.
- **Fix:**
  - Sort handles or indices in parallel runs and merge them. `diffDNA::sortedFragments`
    already does this for variants.
  - Or bucket by count first and sort the buckets in parallel.
  - Comparing the first 8 key bytes as a big-endian integer settles most ties
    before a full compare.
- **Order:** keys are unique, so the order stays exactly the same.

### 6. Per-read passes

Each read is scanned 6–10 times: parsing, `seal`, average quality, accumulated
error, sequence checks, window quality, binomial filter.

| Change | Function speedup | Share | End to end |
| --- | ---: | --- | --- |
| Quality characters from a 256-entry table instead of a clamp per character | 1.8× | 14% of plain rewrites | about 6% there |
| `getAvgQual` sums into a local (the member accumulator is reloaded every iteration, because `signed char` reads may alias it) | 1.5× | 8% (rewrites), 1–4% elsewhere | 1–3% |
| `getAccumError`, `scanSequenceChecks` with local accumulators | none | – | – |
| Remove the per-value mutex in `ReportStats::addMedian2Histo` | 5 ns per lock | – | negligible (the profiler overstated it) |

Also:

- **`getAccumError`** is bound by the latency of one sequential `double` add per
  base. Faster summation orders change the last bits, which could flip a decision
  exactly at a threshold, so keep it exact, or compute it while decoding the
  qualities.
- **`scanSequenceChecks`** (ambiguous bases and homopolymer runs, 8–12% in
  `pe_lotus`/`se_demux`) can be done with SIMD compares of `s[i]` vs `s[i-1]` and a
  table-free ACGT test, exactly.
- **`shared_ptr<DNA>` passed by value** through the per-read path (for example
  `collectstats::addPreFilt`, `ReportStats::addDNAStats`, `Dereplicate::addDNA`)
  costs two atomic updates per call: 5–8% of samples in `pe_rewrite`. Passing
  `const shared_ptr&` or `DNA&` where ownership is not taken removes it.

### 7. Merging

- **Built twice:** with `-merge_pairs_demulti 1` and `-merge_pairs_derep 1` (the
  LotuS defaults), `ReadMerger::merge` builds the consensus read once for
  demultiplexed output and again in `Dereplicate::addDNA`. Caching the merged read
  on the pair would save one of them. `merge` plus `getMergedQual` are 6–10% of CPU
  in those runs.
- **Seed search:** `findSeed` hashes every read-1 position, then verifies hits with
  the same byte-by-byte `percentIdentity` as item 3, which would benefit equally.

### 8. Reader thread

For plain rewrites the reader is the limit at 8 threads.

- **Zero-filling:** `OutputStreamer` is rebuilt for every input file, and its
  demultiplexing stream buffers are zero-filled each time. That is 9% of the reader
  thread in `apong_fa`. Allocating them uninitialised, or retaining them across
  files as the main outputs already are, removes this.
- **Wake-ups:** `pthread_cond_signal` per submitted block (and per read-ahead
  chunk) costs 5–28% of the reader. Waking a worker only when none is spinning or
  already awake, or submitting several blocks per wake-up, cuts it.
- **Longer term:** hand workers raw input chunks cut at record boundaries in both
  mates, and let them split lines themselves. That removes one full copy of every
  read, and the whole `memchr` split, from the serial path.

### 9. Build options

Full rebuilds of the final sources, timed with the concurrent A/B at one thread
(4 pin-swapped pairs each, median ratio against the current `-O3` build):

| Build | apong_fa | pe_lotus | 97% coarse | se_demux |
| --- | ---: | ---: | ---: | ---: |
| `-march=x86-64-v3 -ffp-contract=off` (AVX2, BMI2, FMA allowed but not contracted) | 1.00× | 0.99× | 1.00× | 0.99× |
| LTO (`-flto=auto`) | 0.99× | 1.00× | 0.98× | 0.98× |
| PGO (trained on apong_fa, 97% coarse, pe_lotus, se_demux, pe_rewrite) | 1.07× | 1.12× | 1.01× | 1.03× |

All outputs were identical to the current build. The run logs differ only in the
output directory they name.

- **`-march` and LTO** give nothing. The hot loops are latency chains, or have early
  exits the compiler cannot vectorize, so the gains need source changes (items 1
  and 3).
  - Once such kernels exist, use per-function targets with runtime dispatch
    (`__attribute__((target))` or `target_clones`), so the static binary stays
    portable.
  - Keep `-ffp-contract=off` on any FMA-capable target, so filter probabilities stay
    bit-identical across builds.
- **PGO** is worth up to 12%, but these numbers are in-sample: the training runs
  used the same workloads. A `make pgo` target that trains on bundled small inputs
  would give somewhat less on other data, at the cost of a two-stage build.

## Reproduction

Tools are in
[.claude/archive/speed-audit/tools](../.claude/archive/speed-audit/tools/); raw data
is in [benchmarks/speed-audit](benchmarks/speed-audit/).

- `profall.sh BIN "SCENARIOS" "THREADS" OUT`: light profiles of every scenario.
  `matrix.pl` tabulates function shares across the reports.
- `speedbench.cpp`: the microbenchmarks. Build it inside a build tree against its
  object files and `libdeflate.so.0`; the build command is in the file header.
- `bincmp.cpp`: the linked `binomialFilter` against the same column-wise code in
  another translation unit. This isolates the code-generation effect in item 1.
- `mkexp.sh`, `mkbin.sh`: scratch copies of a build for patch experiments.
  `cmpall.sh` byte-compares every output of two binaries.
- `concab.sh "SCENARIOS" THREADS REPEATS A B`: concurrent, pin-swapped A/B timing.
  `abmin.pl` summarizes sequential `run_ab.sh` runs by their minimum.
- `flagbuilds.sh`: the `-march`, LTO and PGO builds, including PGO training.

Raw data in [benchmarks/speed-audit](benchmarks/speed-audit/):

- `reports/`: the 22 profile reports.
- `profiled-runs.txt`: wall, user and sys time and peak RSS of the profiled runs.
- `results.txt`: every microbenchmark and end-to-end number quoted here.
