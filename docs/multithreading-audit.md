# Multithreading audit: where the workers wait

Audited 2026-09-24 against the working tree on `master` (c91f5f0 plus uncommitted
changes, sdm 3.60 beta). No production source was changed. This report measures
the complete pipeline, including input, filtering, merging, dereplication and output,
not only the dereplication locks already covered in
[coarse-derep-thread-investigation.md](coarse-derep-thread-investigation.md) and
[benchmarks/coarse-work-pool.md](benchmarks/coarse-work-pool.md).

Mitigations 1–8 have since been implemented, and the results re-benchmarked, in
[multithreading-fixes.md](multithreading-fixes.md).

**Summary.** Dereplication is no longer the main multithreading bottleneck in the
default modes. On these datasets its shard and work-pool locks cost under 2% of
worker time. The largest serialization is now in **output**, and several other
serial points follow:

1. **Output.** Every read pair takes up to five locks. `dualMtx` in
   `dualOfBufStream::write2` is held while buffers are flushed to disk or
   compressed. Workers spend 23–61% of their time waiting on it whenever paired
   FASTQ is written. Demultiplexed output has the same pattern through one global
   `dmltMTX`.
2. **glibc malloc.** The reader allocates the header and sequence strings of every
   read. Workers free them, so the main arena lock bounces between threads. This
   costs about 23% of the reader's samples in one run. A tunable-only experiment
   cuts wall time by 4–15%.
3. **The single reader thread.** It does all file I/O, gzip inflate, FASTQ line
   splitting and a header check that allocates two strings per pair. With gzip
   input the reader is saturated at 8 threads.
4. **BGZF output.** Compression uses a pool of `threads/4` threads (2 at 8 threads).
   The waits for compressed blocks happen while `dualMtx` is held. Reopening a `.gz`
   output for each input file re-reads the whole compressed file.
5. **Serial tail.** `writeDereplDNA` and `finishMap` run on the main thread while all
   workers are idle. This is 45% of wall time in 97% coarse mode at 8 threads.
6. **Known derep lock.** It still applies to the global coarse modes
   (`-derepGlobal 1`, reassignment, 100% with retained qualities/differences, one
   thread): 26–41% of worker time.
7. **Small global mutex.** `Filters::updateMaxSeqL` takes a global mutex twice per
   pair to protect per-thread data.

## Setup

- Host: Intel Core Ultra 7 258V under WSL2 (kernel 6.18). **8 logical CPUs = 4
  performance + 4 low-power efficiency cores, no SMT**, so ideal 8-thread scaling
  is well below 8×.
- Build: `make WITH_HTS=0 STATIC=0 CXXFLAGS="-O3 -g -fno-omit-frame-pointer"
  LDFLAGS=-rdynamic`, g++ 13.3. Dynamic linking and frame pointers are required by
  the profiler; the absolute times are therefore a little slower than the default
  static build.
- Inputs: local ext4, page cache warm. Outputs go to local ext4 and are deleted
  after each run.
- Timings are medians of 3 unprofiled runs. The host got up to about 30% slower
  late in the session, so the knob experiments interleave their variants within
  each repeat and should be compared within their own table only.

| Scenario | Input | Work |
| --- | --- | --- |
| `pe_rewrite` | 1.2M synthetic 2×231 bp pairs, 1 file pair | no map/filter; parse and write paired FASTQ |
| `pe_rewrite_gz` | same | same, BGZF `.fq.gz` output |
| `pe_lotus` | 8 samples × 150k pairs | map/primers, filter, merge, 100% prefix derep, paired FASTQ out |
| `pe_lotus_gzin` / `_gzout` | same | plain-gzip input / BGZF main output |
| `se_demux` | 1M single-end reads | inline barcodes, filter, derep, FASTQ out |
| `apong_fa` | real MiSeq 2×300, 336k pairs, 4 samples | LotuS3-style: FASTA derep with merged pairs, `-suppressOutput 1` |
| `apong_fq` | same | FASTQ derep per SR, demultiplex into per-sample files |
| `apong_coarse` | same | `-derepIdentity 97 -derepStoreQuals 1` (64 partitions) |

## Scaling

| Scenario | 1 thr (s) | 2 | 4 | 8 | speedup @8 | CPU @8 / CPU @1 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| pe_rewrite | 4.20 | 4.81 | 4.34 | 3.68 | 1.14× | 3.05 |
| pe_rewrite_gz | 24.36 | 11.58 | 12.11 | 12.53 | 1.94× | 1.68 |
| pe_lotus | 18.76 | 14.13 | 8.31 | 7.89 | 2.38× | 1.94 |
| pe_lotus_gzin | 23.62 | 11.83 | 8.95 | 8.26 | 2.86× | 1.64 |
| pe_lotus_gzout (2 runs) | 30.27 | – | – | 13.96 | 2.17× | 1.86 |
| se_demux | 4.17 | 3.55 | 2.78 | 2.19 | 1.90× | 2.24 |
| apong_fa | 12.42 | 5.36 | 3.10 | 2.29 | 5.42× | 0.98 |
| apong_fq | 13.09 | 4.93 | 3.87 | 2.85 | 4.59× | 1.02 |
| apong_coarse | 13.31 | 4.31 | 3.03 | 2.52 | 5.28× | 0.65 |
| apong, `-derepGlobal 1` (1 run) | 14.21 | – | – | 13.58* | 1.05× | 1.49 |
| apong, identity 100 + quals (1 run) | 9.45 | – | – | 6.39* | 1.48× | 1.34 |

\* under the light profiler (5 ms sampling, no hold tracking). Its wall times fell
within the range of unprofiled runs made at the same time.

Notes on the table:

- The single-thread apong runs use more CPU than two threads (12.4 s vs 10.1 s) and
  vary widely (10.3–14.1 s), so the 2→8 thread ratio (2.3× for 4× the threads) is the
  fairer view of apong scaling. At one thread coarse mode also uses the global index.
- In every workload that writes FASTQ, 8 threads burn 2–3× the single-thread CPU for
  little or no speedup. `pe_rewrite` sys time rises from 1.2 s to 5.4 s, mostly futex
  and allocator traffic.

## Method

Linux `perf`, valgrind and gdb are not available in this WSL image. I wrote a small
`LD_PRELOAD` profiler instead (`.claude/archive/mt-audit/tools/mtprof.c`). It interposes
`pthread_mutex_*`, `pthread_rwlock_*`, `pthread_cond_*wait`, `syscall(SYS_futex)`
(which `std::future` uses), `read`/`write`/`readv`/`writev` and `pthread_create`. It
records:

- per call site: acquisitions, contended acquisitions, total and maximum wait,
  hold time;
- per thread: the exact time blocked in each wait class and on-CPU/run-queue time
  from `schedstat`;
- a 1 ms wall-clock stack sample of every thread, tagged with the thread's state
  (running, or which wait it is in).

`analyze.py` symbolizes the stacks with `addr2line`. libc is stripped, so its
unexported regions are labelled from disassembly: "malloc internals", "vector
string ops", and `__lll_lock_wait/wake_private` for glibc internal locks, which are
malloc arena locks here.

The heavy profile adds roughly 10–70% wall time on lock-heavy runs. The light
profile (5 ms samples, no hold timing) stays within the noise of unprofiled runs
made at the same time and shows the same state split. Lock numbers below are from the heavy 8-thread profiles unless stated.
"Summed wait" is the total over all workers; divide it by the worker thread-time
to get the share lost.

## Where worker time goes at 8 threads

| Scenario | running | mutex wait | idle (no task) | other |
| --- | ---: | ---: | ---: | --- |
| pe_rewrite (light) | 33% | **44%** | 21% | write 3% |
| pe_rewrite_gz | 16% | **61%** | 17% | 5% waiting on BGZF futures |
| pe_lotus (light) | 39% | **30%** | 29% | write 2% |
| pe_lotus_gzin | 38% | 24% | 36% | reader-bound |
| pe_lotus_gzout | 27% | **46%** | 24% | 3% BGZF futures |
| se_demux | 40% | 24% | 36% | |
| apong_fq | 46% | 25% | 27% | |
| apong_fa (light) | 64% | 2% | 34% | serial tail + reader |
| apong_coarse (light) | 33% | 1% | **66%** | serial tail |
| apong, `-derepGlobal 1` | 15% | **41%** | 44% | |

## Findings

### 1. Output is serialized per read pair

Every accepted pair goes through `OutputStreamer::writeForWriteImpl`
(`OutputStreamer.cpp:923`):

- a shared lock on `outputEpochMTX` (`:944`);
- the global `sqfqostrMTX`, taken to copy up to six `shared_ptr` stream handles
  and bump the counters (`:980`);
- `dualOfBufStream::write2`, which holds `dualMtx` (`:305`) while it appends to
  both mates' `ofbufstream`, each under its own `output_mtx`.

When a mate's buffer is full, `ofbufstream::operator<<` calls `primary->write`
on the whole buffer while holding `output_mtx`, and therefore also `dualMtx`
(`:53-69`). The primary buffers are 250/350 kB. The asynchronous writer thread
exists but is disabled: `shouldUseAsyncWriter()` returns `false` (`:247`).

| Lock site | Scenario | Acquisitions | Contended | Summed wait | Max hold |
| --- | --- | ---: | ---: | ---: | ---: |
| `dualMtx` (write2) | pe_rewrite | 1.20M | 22.6% | 18.3 s (40% of worker time) | 9.5 ms |
| | pe_lotus | 1.20M | 15.5% | 15.3 s (23%) | 8.7 ms |
| | pe_lotus_gzin | 1.20M | 14.1% | 17.7 s (18%) | – |
| | pe_rewrite_gz | 1.20M | 18.5% | **55.3 s (56%)** | **104 ms** |
| | pe_lotus_gzout | 1.20M | 13.3% | **68.0 s (41%)** | 7.4 ms |
| `sqfqostrMTX` (writeForWriteImpl) | all FASTQ runs | 0.95–1.2M | 4–11% | 0.7–1.6 s | ≤2.5 ms |
| `output_mtx` (SE `operator<<`) | se_demux | 0.95M | 7.2% | 2.75 s (11%) | – |
| `dmltMTX` (writeDemultiRecords) | apong_fq | 254k | 15.9% | 5.7 s (18%) | – |

Demultiplexing (`writeDemultiRecords`, `OutputStreamer.cpp:703`) serializes every
sample's files behind one `dmltMTX`, including the flushes of their 50 kB buffers.
`writeNonBCReads` (`:1177`) uses the same lock.

### 2. glibc malloc contention is a hidden global lock

The reader fills `multi_tmp_lines` and the workers move the header and sequence
strings into `DNA` objects (`IO.cpp:72-76`, `DNA.cpp:160-184`). On the next block
the reader therefore allocates two fresh strings per mate in the main arena. It
never frees, so its tcache stays empty and every allocation takes the arena path.
The workers free the strings later; frees that miss their tcache go back to the
main arena and take its lock. The result is visible in the stack samples as
glibc's internal lock wait/wake (not interposable, so reported as "running"):

- **Reader** (`pe_rewrite`, light profile): 23% of all its samples are in
  `__lll_lock_wake/wait_private`, plus about 10% in malloc internals. In
  `pe_lotus` (heavy profile), about half of the reader's CPU samples are
  allocator-related.
- **Workers**: 4–15% of samples are in allocator locks, plus 1–8% in malloc
  internals.

Runtime-only knobs, 8 threads, 3 interleaved repeats (medians):

| Scenario | base | `tcache_count=65535` | `arena_max=1` | 8× larger blocks |
| --- | ---: | ---: | ---: | ---: |
| pe_rewrite | 5.61 s | 5.18 s (−8%) | 5.51 s | 6.05 s |
| pe_lotus | 8.59 s | 8.26 s (−4%) | 25.61 s (3.0×) | 8.97 s |
| apong_fa | 3.08 s | 2.85 s (−8%) | 7.96 s (2.6×) | 3.09 s |
| se_demux | 2.81 s | 2.39 s (−15%) | 7.82 s (2.8×) | 3.21 s |

- A large tcache cuts CPU by 6–14%, but peak RSS grows by 160–200 MB.
- Forcing one arena is catastrophic, so allocation is on the critical path.
- Larger dispatch blocks (`-iniBlockSize 1000`) do not help. Per-task overhead is
  not the problem.

### 3. One reader thread does all parsing

`read_sequences` (`IO.cpp:387`) runs `getDNAlines` on the main thread for every
block. That covers synchronous 80 kB `read()`s (`InputStream.h:294`), gzip
inflation (`InputFileStream.cpp:110`), line splitting, and
`validateSynchronizedHeaders` (`InputStream.cpp:60`).

- **Header check:** `validateSynchronizedHeaders` copies both headers and runs
  `remove_paired_info` on each. It takes 25–37% of the reader's running time in
  `pe_rewrite`/`pe_lotus`. Workers repeat an equivalent check (`sameHeadPosFree`,
  `IO.cpp:328`).
- **Gzip input:** in `pe_lotus_gzin` the reader is running 73% of the time, and 65%
  of that is zlib inflate and crc32 for R1 and R2 in turn. It is the bottleneck
  once the output locks are out of the way.
- **Dispatch:** only `threads` blocks can be in flight (`slots(Nthreads)`,
  `IO.cpp:412`). When every slot is busy, the reader waits on the future of one
  specific slot (`IO.cpp:518-524`), even if another slot has already finished. In
  the profiles, workers are 10–19% idle while the reader is blocked this way,
  and there is no read-ahead. This is secondary for now, because most reader
  waits are really waits on workers stuck in output locks.

### 4. BGZF output

- `ThreadPool::configure` sizes the gzip pool at `threads/4`, minimum 2
  (`ThreadPool.h:154`). At 8 threads 2 compressors do all deflate work. They were
  77% busy in `pe_rewrite_gz`, while the 8 workers were mostly blocked.
- `ParallelBgzfOstream::queueCurrentBuffer` waits on the oldest compression future
  once `maxPending_` (≤8) blocks are queued (`OutputFileStream.cpp:318-366`). That
  wait runs on a worker that holds `output_mtx` and `dualMtx`: 5.1 s of futex waits,
  and the 104 ms maximum hold.
- Every input file creates a new `OutputStreamer` (`IO.cpp:1008`) and reopens the
  main outputs for append. For BGZF, `prepareAppend` (`OutputFileStream.cpp:389`)
  reads every existing compressed block to rebuild the index. `pe_lotus_gzout`
  read 1.40 GB (125,734 `read()` calls) to write 0.40 GB over 8 inputs. The total is
  about (N−1)/2 × the final size, so it grows quadratically with the number of
  samples. It cost 0.34 s from page cache here; on network storage with hundreds
  of samples it would dominate.

### 5. Serial tail: dereplication output and teardown

`writeDereplDNA` (`containers.cpp:897`), `outputGroups`, and `finishMap`
(`containers.cpp:828`) run on the main thread after all workers are idle. So does
`finishMap`'s destruction of the shard tables.

| Scenario, 8 threads | Profile | wall | serial tail |
| --- | --- | ---: | ---: |
| apong_coarse | light | 3.67 s | **1.65 s (45%)** |
| se_demux | heavy | 3.15 s | 0.45 s (14%) |
| apong_fq | heavy | 3.99 s | 0.45 s (11%) |
| apong_fa | light | 2.38 s | 0.25 s (10%) |
| pe_lotus | light | 7.08 s | 0.55 s (8%) |

In coarse mode about a third of the main thread's CPU samples are in
`diffDNA::writeFastq`/`forEachPair`. About 14% are in `CoarseDereplication::clear`,
freeing robin_hood nodes. With `-derepPerSR 1` the same pause (`LifecyclePause`)
also stops all producers mid-run, once per SequencingRun.

### 6. Dereplication locks

Status of the known issue, measured again:

- **Default exact/prefix derep**: `flush_batch` locks the 64 shards in ascending
  order (`containers.cpp:356-363`). Only 0.4–0.8% of acquisitions are contended,
  for 0.05–0.28 s of summed wait (<1% of worker time), with an average hold of
  10–18 µs. The fixed-order convoy has not shown up at 8 threads.
- **Default coarse** (97%, 64 partitions, `PartitionWorkPool`): 51 ms of mutex plus
  165 ms of condvar waits in `PartitionWorkPool::run`, 0.8% of worker time. The
  pool's `notify_all` per task (`PartitionWorkPool.h:55,95`) is fine at 8 threads.
  Coarse mode is now limited by the serial tail (finding 5).
- **Global coarse index** (`globalCoarse`, `containers.cpp:453`): chosen for
  `-derepGlobal 1`, `-derepReassign 1`, identity 100 with retained qualities or
  differences, and one thread. `flush_batch` holds one mutex for a whole
  256-entry batch.

| Global coarse run | Contended | Summed wait | 8-thread speedup |
| --- | ---: | ---: | ---: |
| `-derepGlobal 1` | 94% | 44.5 s (41% of worker time) | 1.05× |
| identity 100 + `-derepStoreQuals 1` | 90% | 13.3 s (26%) | 1.48× |

- **Latent:** all workers share `Dereplicate::merger`
  (`containers.cpp:666-673`). When `-merg_readpos`/`-qual_readpos` are set, it takes
  `statsMutex` for an O(read-length) histogram update on every merge
  (`ReadMerger.cpp:444-455`). Not exercised in these scenarios.

### 7. Smaller items

- **`Filters::updateMaxSeqL`** takes the global `updateMaxSeqMutex`
  (`Filters.cpp:1425`) twice per pair from every worker. It protects
  `maxReadLength`/`min_l`, which belong to the calling thread's own sub-filter.
  2.4M acquisitions, 2–5% contended, 0.15–0.86 s of summed wait. It is the most
  frequently contended lock outside output.
- **Uncontended per-thread mutexes:** `ReportStats::addMedian2Histo`
  (`Statistics.cpp:511`) and `PrimerVariantTracker::record` (`Statistics.cpp:10`)
  lock per-thread objects about 8 times per pair (9.6M lock operations in
  `pe_lotus`). They cost CPU but no contention.
- **Not a threading issue, but the largest worker hot spot in LotuS mode:**
  `DNA::binomialFilter`/`sum_of_binomials` take 21–28% of worker CPU in apong.
  With `n == 1`, `prob_j_errors` is zero for i > 1, yet `sum_of_binomials` loops
  i = 0..j (`DNA.cpp:424-440`). Limiting it to two terms makes the filter
  O(L·E) instead of O(L·E²).

## Proposed mitigations

Ordered by expected gain per effort. The first four need no change in output
semantics.

1. **Write output per block and outside locks.**
   - Each worker already formats records into `thread_local` strings. Have it
     format the whole block (R1, R2, singletons, FASTA/QUAL) into per-worker
     buffers instead.
   - Publish each buffer with one short critical section per block, or hand it
     to a single writer thread. That writer can also keep blocks in input order
     for deterministic output.
   - This replaces about 5 lock operations and 6 `shared_ptr` refcount updates per
     pair with one per block. Mates stay synchronized because a block's R1 and R2
     buffers are published together.
   - `maxReadsPerOFile`/`maxRdsOut` accounting moves to the writer, or to an
     atomic reservation of the block's record count.
   - Never call `primary->write` or wait on compression while holding an append
     lock; swap the full buffer out and write it after unlocking. Re-enabling the
     existing async writer path is the smallest version of this.
   - For demultiplexing, use per-sample locks or per-worker per-sample buffers, and
     reserve the `DemultiBPperSR` budget with an atomic CAS instead of `dmltMTX`.
   - Expected gain: the 23–61% of worker time now spent waiting.
2. **Compress BGZF in the workers.** A worker's block buffer can be deflated into
   64 kB BGZF blocks before hand-off, so the writer only writes bytes and index
   entries. Otherwise, size the gzip pool to the worker count when a `.gz` output
   is open. libdeflate is a 2–3× faster drop-in at level 1.
3. **Stop cross-thread string hand-off.** In `str2DNA`, copy instead of moving the
   header and sequence out of the reusable block. The reader's strings then keep
   their capacity and it stops allocating after warm-up, while each worker
   allocates from its own arena. Recycling `DNA` objects per worker, or an optional
   mimalloc/jemalloc build, are alternatives. The tcache experiment gives 4–15%
   before any other fix; the gain should grow once the output locks are gone.
4. **Make the reader's header check allocation-free.** Compare `string_view`s with
   the same `remove_paired_info` rules, or drop it in favour of the worker-side
   `sameHeadPosFree`. That frees 25–37% of the reader's CPU.
5. **Remove `updateMaxSeqMutex`.** The data is per sub-filter, so the lock can go.
   Also drop the mutexes from the per-thread `ReportStats` and
   `PrimerVariantTracker` paths, keeping them only where results are merged.
6. **Parallelize the serial tail.**
   - Partitions (coarse) and shards (exact) are independent until the final
     abundance sort. Build their output groups in parallel on the idle pool, sort
     `(count, key)` handles globally, then format records in parallel chunks
     (`diffDNA::writeFastq` is the heavy part) and write the chunks in order.
   - Skip or parallelize the teardown of the shard tables in `finishMap`, e.g.
     by clearing shards on the pool.
   - For `-derepPerSR`, swap in a fresh `Dereplicate` and write the finished one
     in the background, so producers are not paused.
7. **Feed the workers better.**
   - Inflate R1 and R2 on two reader threads that fill bounded byte queues. This
     roughly doubles gzip-input throughput; zlib-ng or libdeflate add more.
   - Longer term, let the reader hand raw byte chunks (cut at record boundaries
     counted in both mates) to workers that split the lines themselves.
   - Keep 2–3× as many block buffers as workers, and wait on any completed slot
     through a completion queue rather than a fixed future. Index per-thread
     `Filters`/`ReadMerger` by worker id rather than slot.
8. **Keep outputs open across input files.** Do not rebuild the `OutputStreamer`
   and reopen outputs for each input file. Where reopening is unavoidable, restore
   the BGZF index from the `.gzi` sidecar or scan headers with seeks instead of
   reading every block. With many small samples this also removes the per-file
   pipeline drain and re-setup (`IO.cpp:906-1112`), which was not measured here.
9. **Global coarse modes (known).**
   - Follow the plan in
     [derep-parallel-partitioning.md](derep-parallel-partitioning.md): exact-key
     sharding for fixed-length 100%, and frozen-snapshot parallel search with
     ordered commit for approximate global mode.
   - Give each worker its own merge statistics instead of the shared
     `Dereplicate::merger`.
   - If more threads are used later, replace the work pool's `notify_all` with
     per-batch completion signalling.
10. **Algorithmic, outside threading.** Restrict the `sum_of_binomials` inner loop
    to two terms. That is about a quarter of worker CPU in LotuS mode.

## Reproduction

The raw data is in [benchmarks/mt-audit](benchmarks/mt-audit/); the tools are in
[.claude/archive/mt-audit/tools](../.claude/archive/mt-audit/tools/):

- `mtprof.c`: the `LD_PRELOAD` profiler. Build it with
  `gcc -O2 -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0 -fno-omit-frame-pointer -fPIC -shared -o libmtprof.so mtprof.c -ldl -lpthread`.
- `analyze.py`, `summarize.py`, `hol.py`, `timesum.py`: analysis.
- `scenarios.sh`, `run_matrix.sh`, `run_knobs.sh`: drivers. They expect
  `~/sdm-mtprof/{src,prof,runs}` and the datasets under `~/sdm-audit`.
- `time.csv`, `knobs.csv`, `time2.csv`: raw timings.
- `profiles-summary.txt`, `profiles-light-summary.txt`, `serial-tail.txt`, and
  `reports/*.txt`: full per-scenario profile reports.

Run the profiler on `sdm` only, as
`/usr/bin/time … env LD_PRELOAD=… MTPROF_OUT=… sdm …`. If `/usr/bin/time` is
itself preloaded, it overwrites the profile when it exits.
