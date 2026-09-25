# Speed fixes: binomial filter, libdeflate, SIMD verification, reader, derep frees and sort

Implemented 2026-09-24 from items 1–5 and 8 of [speed-audit.md](speed-audit.md),
on top of [multithreading-fixes.md](multithreading-fixes.md). Nothing is committed.
Items 4 and 5 came second and are measured against the build with items 1, 2, 3 and
8 (see [Items 4 and 5](#items-4-and-5-against-the-first-round)).

**Summary.** Items 1, 2, 3 and 8 against the previous build:

- **LotuS-style runs** (binomial filter on) are 1.37–1.41× faster at 1, 4 and 8
  threads.
- **Coarse derep modes** are 1.14–1.35× faster.
- **BGZF output** is 1.55–1.78× faster and 11% smaller; LotuS runs with `.gz`
  output are 1.20–1.28× faster.
- **The reader thread** uses 19–49% less CPU. That shows in wall time only for
  plain rewrites, at 1.07× (by minimum; the medians are within noise).
- **Other modes** are unchanged, within ±3%.
- **Outputs:** all are identical. The exceptions are compressed `.gz` bytes and
  `.gzi` offsets; the decompressed content is identical.

Items 4 and 5 against that build:

- **Arena lock waits** at 8 threads fall by 73% (`pe_lotus`), 76% (`se_demux`),
  63% (`apong_fa`) and 25% (`apong_coarse`).
- **The derep output sort** runs on the thread pool. In a microbenchmark it is
  2.2–2.7× faster than `std::sort` from 66,000 groups up, with no extra memory per
  group.
- **Wall time:**
  - The three `pe_lotus` scenarios and `se_demux` are 1.08–1.11× faster at
    4 threads, and 1.10–1.40× at 8 threads (minimum of 5 runs).
  - `apong_fa`, `apong_fq` and 97% coarse gain 1.04–1.06× at 4 threads; their
    8-thread gains are within noise.
  - Nothing changes at 1 thread, or in the global and 100% modes.
- **Outputs** are identical, and peak memory is unchanged at 1×, 2× and 4× input.

## What changed

| Item | Change | Files |
| --- | --- | --- |
| 1 | Column-wise binomial error filter; two-term row-wise fallback | `DNA.cpp` |
| 2 | libdeflate for BGZF blocks, detected at build time, zlib fallback | `OutputFileStream.cpp`, `InputFileStream.cpp`, `Makefile`, `README.md`, `licenses/` |
| 3 | AVX2 bounded mismatch count; stateless `percentIdentity` | `ReadMerger.{h,cpp}`, `CoarseDereplication.cpp` |
| 4 | Reads of blocks applied on another producer's thread go back to their producer to be freed | `containers.{h,cpp}` |
| 5 | In-place parallel sort of the derep output groups; shards collected in parallel | `containers.{h,cpp}`, `tests/TestCoarseDerep.cpp` |
| 8 | Spinning idle workers and fewer wake-ups; read-ahead hysteresis; output buffers without zero-fill; no per-record line reset | `ThreadPool.h`, `InputFileStream.cpp`, `OutputStreamer.{h,cpp}`, `InputStream.cpp` |

### 1. Binomial error filter

`DNA::binomialFilter` now computes all error-count rows in one pass along the read.

- **Recurrence:** row e at position k is `(1 - p) * row_e + p * row_(e-1)`, with
  both inputs taken at position k−1. Descending e updates one small column in
  place.
- **Codegen:** the kernel is a template on the row count, dispatched through a
  table for 1–16 rows (`BinErrorModelMaxExpError` up to 13). Loop bounds are then
  compile-time constants and the loop calls no functions. A column copy with a
  runtime-size `memcpy` had made an earlier prototype 2× slower inside `DNA.cpp`.
- **Stopping:** the error count is chosen afterwards, with the same row-by-row
  accumulation and stopping rule as before.
- **Larger limits:** the existing row-wise code stays. Its `sum_of_binomials` now
  stops after the two terms that can be non-zero for one Bernoulli trial; the
  dropped terms added exactly `0.0f`.
- **Exactness:** every value is computed by the same float operations as before, so
  results are bit-identical (on builds without FMA contraction; see
  [speed-audit.md](speed-audit.md#9-build-options)).
- **Speed:** as linked, 500,000 calls on real apong reads take 441 ms instead of
  4,682 ms (10.6×; measured while other jobs ran). That is the stand-alone speed of
  the same kernel.

### 2. libdeflate for BGZF

`compressBgzfBlock`, the CRC of stored blocks, and `inflateBgzfBlock` use libdeflate
when `SDM_WITH_LIBDEFLATE` is defined.

- **Compression:** level 1, as with zlib's `Z_BEST_SPEED`.
- **Decompression:** `libdeflate_gzip_decompress_ex`, which also checks the gzip
  header, CRC-32 and ISIZE, and must consume the whole block.
- **State:** compressor and decompressor are allocated once per thread. zlib's
  state was set up and torn down for every 64 KB block.
- **Container:** unchanged. Files stay valid BGZF with identical decompressed
  content and a consistent `.gzi` index (every entry checked against the blocks).
  They are about 11% smaller.
- **Build:**
  - `make` defines `SDM_WITH_LIBDEFLATE` when a test program links
    (`WITH_LIBDEFLATE=auto`); for `STATIC=1` that needs `libdeflate.a`.
    `WITH_LIBDEFLATE=0` keeps zlib, and `1` requires libdeflate.
  - A stamp file rebuilds the two stream objects when the setting changes.
  - A `#` inside `$(shell)` needs escaping before GNU make 4.3 and must not be
    escaped after, so the test program takes it from a variable.
- **Documentation:** the README lists the dependency and options, and its MIT
  license is in `licenses/libdeflate-LICENSE`.
- **Scope:** plain (non-BGZF) gzip input still streams through zlib; libdeflate
  cannot stream.

### 3. SIMD candidate verification

`ReadMerger::percentIdentity` counts mismatches with AVX2, 32 bytes per step.

- **Early exit:** it stops after the first block whose running count exceeds the
  budget. It returns −1 exactly when the byte-by-byte loop would, and otherwise the
  same match fraction.
- **Portability:** runtime dispatch (`__builtin_cpu_supports`) keeps the scalar loop
  for CPUs without AVX2 and for non-x86 builds, so the static binary stays portable.
- **Stateless:** the function is now `static`, so `CoarseDereplication::search` no
  longer constructs a `ReadMerger` (seed map and buffers) per call.
- **Merging:** read merging's seed verification uses the same function, and gains
  too.
- **Not done:** the seed-count filter proposed to reduce the number of comparisons
  (about 330 per search).

### 4. Freeing reads on the thread that allocated them

- **Cause:** in the partitioned derep modes, a producer that finds sealed blocks
  ready applies all of them in admission order (`apply_blocks_ordered`), including
  blocks other producers sealed. Their entries were then freed on the applying
  thread, and every free locked the arena of the producer that had allocated the
  read.
- **Fix:** after applying, every block sealed by another producer goes back to it
  through `return_reads`, the scheme the global mode already used. The owner frees
  it at its next `endBlock`, or `flush_all_batches` does at the end. The drain path
  hands back every block. Blocks of the applying producer itself are freed right
  away, in its own arena.
- **Memory:** a producer holds at most the returned entries of blocks applied since
  its last `endBlock`, which the ordered-admission window bounds. On apong input
  scaled 1×, 2× and 4×, peak RSS at 8 threads stays within −3 to +17 MiB of the
  previous build, with no trend (`apong_fa`, `apong_coarse`, `apong_id100`).
- **Lock waits:** 2 ms stack samples in glibc's internal lock wait and wake at
  8 threads, median of 3 profiles per build:

  | Scenario | Before | After | Change |
  | --- | ---: | ---: | ---: |
  | pe_lotus | 1,549 | 414 | −73% |
  | se_demux | 784 | 192 | −76% |
  | apong_fa | 429 | 157 | −63% |
  | apong_coarse | 756 | 570 | −25% |
  | apong_global | 62 | 86 | noise (runs 49–114) |
  | apong_id100 | 80 | 83 | noise (runs 70–96) |

  The global and 100% scenarios do not use `apply_blocks_ordered`.
- **What is left:** in `pe_lotus`, `se_demux` and `apong_fa`, 69–89% of the
  remaining waits are the end-of-run teardown. `clear_shards` destroys the shards
  in parallel, and the `DNAuniqSet` and `DerepRead` objects it frees were allocated
  by the producers.
  - In `apong_coarse`, 88% are the same teardown, mostly the quality seeds of
    `diffDNA::SequenceData`.
  - About 10% there are vector frees in `CoarseDereplication::add` and `indexGroup`
    while blocks are applied.
- **Tried and reverted:** overwriting a representative's quality seed in place
  when a better read replaces it (`SequenceData::selectQuality`). Outputs stayed
  identical, but the seed frees did not go down: a deeper profile showed they come
  from the teardown. The change was dropped.

### 5. Parallel derep output sort

- **Sort:** `sortOnPool` in `containers.cpp` sorts in place, level by level:
  - Each level splits every segment at its median with `std::nth_element`, all
    segments of the level in parallel.
  - That stops at one segment per pool worker, and the segments are then sorted in
    parallel with `std::sort`.
  - Below 32,768 items, at 1 thread, or when called from a pool worker, it is plain
    `std::sort`.
- **Where:** the final count-then-key order in `writeDereplDNA`, and the key sort
  of prefix mode (`outputGroups`).
- **Collection:** `outputGroups` also collects the shards in parallel, each into a
  precomputed slice, in the same order as the serial pass. In prefix mode the
  second sort, of the anchors, is now a reversal: the backward walk over the
  key-sorted groups yields them in descending key order.
- **Exactness:** keys are unique, so both orders are strict total orders, and any
  correct sort gives exactly `std::sort`'s result.
- **Memory:** nothing per item. A first version sorted index arrays in parallel
  runs and merged them. It cost 8 bytes per group, which grows with the data, and
  it was slower from 250,000 groups up, so it was replaced.
- **Speed:** microbenchmark at 8 threads on items shaped like the output groups
  (250 bp keys with a shared prefix, 70% with count 1), minimum of 5 runs:

  | Groups | `std::sort` | Parallel | Speedup |
  | ---: | ---: | ---: | ---: |
  | 33,000 | 4.6 ms | 3.2 ms | 1.4× |
  | 66,000 | 11.4 ms | 5.3 ms | 2.2× |
  | 250,000 | 70 ms | 30 ms | 2.3× |
  | 1,000,000 | 368 ms | 159 ms | 2.3× |
  | 4,000,000 | 2.06 s | 0.75 s | 2.7× |

  It is memory-bound: each comparison follows pointers into scattered sequences.
  The index version reached 1.4–2.1×.
- **Not done:** the audit's other ideas, bucketing by count first and comparing the
  first 8 key bytes as an integer.

### 8. Reader thread

- **Worker wake-ups:**
  - Idle `ThreadPool` workers spin for up to 50 µs before they sleep.
  - `submit` wakes a sleeper only when more tasks are queued than workers are
    spinning.
  - A spinning worker leaves the spin under the queue lock and re-checks the queue,
    so a submitter that relied on it cannot lose a task.
  - Wake-ups had cost the reader 5–28% of its time.
- **Read-ahead:** the gzip read-ahead thread is woken once half of its queue (now
  8 × 256 KB) is free, instead of once per chunk. The reader is woken only when it
  can be waiting on an empty queue.
- **Output buffers:** output-stream buffers (`OutputBuffer`) use an allocator that
  does not zero-fill. This removes megabytes of zeroing per stream whenever
  `OutputStreamer` is rebuilt for an input file.
- **Line reset:** `InputStreamer::getDNAlines` clears the four line strings of a
  record instead of assigning an empty string to each.
- **Not done:** handing raw input chunks to the workers, so they split lines
  themselves. See [Remaining](#remaining).

## Checks

- **Outputs at 1 thread:** every output of all 11 scenarios is byte-identical to the
  previous build. The exceptions are `.gz` files (compared by decompressed content
  and `gzip -t`) and `.gzi` indexes (validated against the blocks).
- **Derep determinism:** outputs are byte-identical at 1, 4 and 8 threads (8 twice)
  in `apong_fa`, `apong_fq`, `apong_coarse`, `apong_global`, `apong_id100`,
  `pe_lotus` and `se_demux`.
- **Test suites:**
  - Release builds with libdeflate, and zlib-only (`WITH_LIBDEFLATE=0`), each pass
    23/23.
  - The ThreadSanitizer build passes 23/23 suites and the 9 scenario runs with no
    warnings.
  - The default static build links against `libdeflate.a`.
- **Items 4 and 5,** against the build with items 1, 2, 3 and 8:
  - Every output of all 11 scenarios is byte-identical at 1 thread (`.gz` by
    content).
  - Derep outputs are byte-identical at 1, 4 and 8 threads (8 twice) in the same 7
    scenarios. `pe_lotus` (66,521 groups) and `se_demux` (58,054) go through the
    parallel sort at 4 and 8 threads; the apong sets (about 26,000) stay below its
    threshold.
  - On the 2× and 4× apong inputs at 8 threads, every derep output equals the
    previous build's.
  - New test `testParallelOutputSort` (`tests/TestCoarseDerep.cpp`): 90,000 exact
    groups go through the parallel sort, since the pool defaults to the hardware
    thread count, and the output must be in `std::sort`'s order.
  - Release (libdeflate), zlib-only and ThreadSanitizer builds pass 23/23 suites;
    TSan also passes the 9 scenario runs; no warnings.

## Benchmarks

Same host and scenarios as the audit, previous build (A) against this one (B).

- **1 and 4 threads:** both binaries ran at the same time on their own vCPUs, with
  pins swapped between pairs. The value is the median of 4 per-pair ratios. This
  cancels the host's throttling (see [speed-audit.md](speed-audit.md#method)).
  - Correction: at 4 threads each binary had only 2 vCPUs (a bug in `concab.sh`,
    fixed since), so the 4 threads were oversubscribed 2×. The table was earlier
    described as 4 vCPUs each.
- **8 threads:** that needs the whole machine, so A and B alternated in 3 repeats.
  Shown are minimum wall times, with user CPU as medians.

| Scenario | 1 thread | 4 threads | 8 threads (A → B, s) | 8 threads | User CPU at 8 threads (s) |
| --- | ---: | ---: | --- | ---: | --- |
| apong_fa | **1.38×** | **1.39×** | 1.56 → 1.14 | **1.37×** | 11.0 → 8.3 |
| apong_fq | **1.39×** | **1.38×** | 1.81 → 1.28 | **1.41×** | 12.3 → 8.3 |
| 97% coarse | 1.35× | 1.17× | 1.43 → 1.25 | 1.14× | 6.3 → 5.3 |
| global coarse | 1.34× | 1.26× | 2.40 → 1.88 | 1.28× | 10.9 → 7.8 |
| 100% coarse | 1.27× | 1.21× | 1.66 → 1.40 | 1.19× | 6.8 → 5.4 |
| pe_rewrite_gz (BGZF output) | **1.78×** | **1.59×** | 3.29 → 2.12 | **1.55×** | 21.4 → 11.9 |
| pe_lotus_gzout | 1.28× | 1.21× | 5.93 → 4.96 | 1.20× | 38.6 → 30.9 |
| pe_rewrite | 1.00× | 1.01× | 0.75 → 0.70 | 1.07× | 3.9 → 3.7 |
| pe_lotus | 0.99× | 1.01× | 3.67 → 3.62 | 1.01× | 23.3 → 22.0 |
| pe_lotus_gzin | 1.01× | 1.01× | 4.39 → 4.28 | 1.03× | 26.6 → 26.5 |
| se_demux | 1.00× | 1.01× | 1.19 → 1.16 | 1.03× | 6.1 → 6.7 |

`pe_lotus`, `pe_lotus_gzin` and `se_demux` use neither the binomial filter, BGZF
output nor the coarse search, so the table shows the others are not slowed down.

Reader (main thread) CPU at 8 threads, light profiles of both builds:

| Scenario | A → B | Main change |
| --- | --- | --- |
| pe_lotus_gzin | 1.38 → 0.71 s | read-ahead wake-ups (`pthread_cond_signal` 0.38 → 0.11 s) and copying |
| apong_fa | 0.50 → 0.32 s | zero-filling of output buffers (0.25 → 0.11 s in libc string/fill routines) |
| pe_rewrite | 0.58 → 0.47 s | line handling and wake-ups |

Only plain rewrites are still limited by the reader at 8 threads.

## Items 4 and 5 against the first round

The build with items 1, 2, 3 and 8 (A) against the one that adds items 4 and 5 (B),
measured the same way as above, with two differences:

- At 4 threads each binary now has 4 vCPUs of its own.
- At 8 threads each binary ran 5 times; the table shows minimum wall times.

Only scenarios with dereplication are included; `pe_rewrite` and `pe_rewrite_gz`
do not run that code.

| Scenario | 1 thread | 4 threads | 8 threads (A → B, s) | 8 threads | User CPU at 8 threads (s) |
| --- | ---: | ---: | --- | ---: | --- |
| pe_lotus | 1.01× | **1.10×** | 3.33 → 2.88 | **1.16×** | 22.4 → 19.7 |
| pe_lotus_gzin | 0.99× | **1.11×** | 3.98 → 3.62 | **1.10×** | 26.0 → 24.3 |
| pe_lotus_gzout | 0.99× | 1.08× | 4.69 → 4.17 | **1.12×** | 29.9 → 26.4 |
| se_demux | 1.01× | **1.11×** | 1.18 → 0.84 | **1.40×** | 6.4 → 5.5 |
| apong_fa | 1.00× | 1.06× | 1.06 → 0.95 | 1.12× | 7.8 → 8.4 |
| apong_fq | 0.99× | 1.05× | 1.14 → 1.04 | 1.10× | 8.0 → 8.3 |
| 97% coarse | 1.00× | 1.04× | 1.13 → 1.10 | 1.03× | 5.3 → 5.0 |
| global coarse | 0.99× | 1.01× | 1.77 → 1.71 | 1.04× | 7.8 → 7.8 |
| 100% coarse | 0.99× | 1.00× | 1.30 → 1.36 | 0.96× | 5.6 → 5.2 |

- **Expected nulls:** at 1 thread neither change is active, and the global and 100%
  modes use neither `apply_blocks_ordered` nor, at about 26,000 groups, the parallel
  sort. Their 0.96–1.04× shows the noise level at 8 threads.
- **Clear gains:** the three `pe_lotus` variants and `se_demux` gain in minimum and
  median alike (medians 1.13–1.26×). `se_demux`'s 1.40× is its minimum; its median
  is 1.24×.
- **Within noise:** the apong gains at 8 threads, 1.03–1.12× by minimum, are not
  confirmed by the medians (0.95–1.10×). At 4 threads, where both binaries run at
  once, they are 1.04–1.06×.
- **Conditions:** another job shared the machine during part of the first run.
  - The 4- and 8-thread runs were repeated once it was idle; the earlier runs are
    kept in the raw data, marked as overlapping.
  - The 1-thread pairs were kept: each pair runs both binaries at once, and all were
    within 0.99–1.01×.

## Remaining

- **Reader line copying:** in `pe_rewrite` at 8 threads, 0.29 of the reader's
  0.47 s is splitting and copying lines into the block's strings.
  - Handing workers raw chunks cut at record boundaries in both mates would move
    that copy into the workers, while the reader still scans for newlines.
  - That is worth at most about 1.15× on plain rewrites; the output writer thread
    (60% busy) is the next limit.
  - It touches every input format, the mate-synchronization checks and
    malformed-record recovery, so it was not done here.
- **Wake-ups at long block intervals:** in plain rewrites a block arrives about
  every 300 µs, longer than the 50 µs spin, so the reader still wakes sleeping
  workers (0.06 s).
  - Keeping one worker spinning, which wakes the next when it takes a task, would
    move that cost off the reader.
  - Longer spins would compete with the reader and writer for the 8 vCPUs.
- **CPU cost of spinning:** the 50 µs spin adds CPU wherever workers go idle often.
  `se_demux` used 10% more user CPU at 8 threads (6.1 → 6.7 s) at unchanged wall
  time. With 3 runs on this host that is within noise, but it matches the expected
  cost. If it matters on shared machines, shorten the spin or make it adaptive.
- **Teardown frees:** 69–89% of the arena lock waits left after item 4 are in the
  parallel `clear_shards` at the end of a run, which frees reads the producers
  allocated.
  - Reads do not record which producer allocated them, so the teardown cannot hand
    them back.
  - A thread-caching allocator (mimalloc, jemalloc) would remove them, and the rest,
    without bookkeeping in sdm. It is still not measured here.
  - Skipping the teardown at program exit would too, but it hides real leaks from
    leak checkers.
- **Audit items not done here:**
  - 6 and 7: per-read passes, merging the pair once;
  - 9: PGO;
  - the seed-count filter of item 3.

## Reproduction

Tools are in [.claude/archive/speed-audit/tools](../.claude/archive/speed-audit/tools/),
and raw data in [benchmarks/speed-audit/after](benchmarks/speed-audit/after/).

- **Tools:**
  - `verify_sp.sh`: determinism, release and zlib-only test suites, TSan
    (`B=<build dir>`).
  - `bench_sp.sh`: the three benchmark runs above (`A=`, `B=` build dirs, `SCEN=`).
  - `cmpall.sh`: all outputs, `.gz` by content.
  - `gzicheck.pl FILE.gz`: validates a BGZF file's blocks and its `.gzi` index.
  - `lockprof.sh` and `lockcallers.pl`: glibc lock-wait samples and their callers.
  - `memscale45.sh`: peak RSS of two builds on scaled inputs.
  - `sortbench.cpp`: the sort microbenchmark.
  - `quiet_bench.sh`: waits until WSL is idle, then runs the 4- and 8-thread
    comparisons and logs any other load seen meanwhile.
- **Raw data:**
  - `concab_t1.log`, `concab_t4.log`: every pair.
  - `ab_sp.csv`: the 8-thread runs.
  - `profile_*_t8.txt`: profiles of the new build.
- **Raw data for items 4 and 5,** in
  [benchmarks/speed-audit/items45](benchmarks/speed-audit/items45/):
  - `concab_t1.log`, `concab_t4.log`, `ab_sp45.csv`: the benchmark runs in the
    table. `quiet_bench.sh` reran the last two once the machine was idle.
  - `*_overlapping_load*`: the earlier runs, made while another job was running.
  - `lockprof_rep.log`, `locks/`: lock-wait profiles, with the top callers per run.
    `locks/apong_coarse_sp45_t8_deep.txt` has 16 frames per chain.
  - `memscale45.csv`: peak RSS on 1×, 2× and 4× input.
