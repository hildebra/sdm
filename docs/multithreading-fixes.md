# Multithreading fixes and deterministic dereplication

Implemented 2026-09-24 on the working tree of `master` (c91f5f0 plus uncommitted
changes). This follows up [multithreading-audit.md](multithreading-audit.md): it
covers mitigations 1–8 from that report and makes multithreaded dereplication
deterministic again. Nothing is committed.

**Summary.**

- At 8 threads, runs that write FASTQ are 1.4–5.1× faster: `pe_rewrite` 3.59 → 0.70 s,
  `pe_lotus` 5.71 → 3.87 s, BGZF output 7.79 → 3.16 s. Derep runs are 1.1–1.5× faster,
  including 97% coarse (2.01 → 1.60 s over 5 paired repeats).
- Single-thread runs are unchanged to 9% faster.
- Dereplication output (`derep.fna`, `.rest`, `.map`, `.hq.fq`) is now byte-identical
  at every thread count, and identical to the old single-threaded output.
- Main FASTQ/FASTA output keeps the same records with mates in sync, but the record
  order still depends on thread timing, as before.
- Peak RSS at 8 threads grows by 10–45 MB.
- The release build passes 23/23 test suites. The TSan build passes 23/23 suites and
  the 9 scenario runs with no warnings.
- **Follow-up** ([below](#follow-up-global-coarse-commits-and-variant-export)):
  - Global-coarse commits no longer block the workers, and variant export (97%
    coarse with retained qualities) is formatted in parallel.
  - Against the state above, the coarse derep modes are another 1.15–1.65× faster
    at 8 threads.
  - Peak RSS of the parallel export stays within run-to-run noise (about ±12 MB) of
    the serial export at 1×, 2× and 4× the input. The naive version needed +77, +100
    and +155 MB.
  - Output is unchanged.

## What changed

| # | Audit mitigation | Implementation | Main files |
| --- | --- | --- | --- |
| 1 | Write output per block, outside locks | Per-worker `WriteBatch` (paired, single and per-sample demultiplexed buffers), published once per read block or every 4 MiB. New `AsyncOutputWriter` thread writes full buffers, so disk writes and compression no longer run under append locks. `DemultiBPperSR` budget reserved with a CAS. | `OutputStreamer.{h,cpp}`, `IO.cpp` |
| 2 | BGZF compression capacity | Gzip pool sized to the worker count (was `threads/4`). Up to `2 × pool` blocks in flight (cap 64). | `ThreadPool.h`, `OutputFileStream.cpp` |
| 3 | No cross-thread string hand-off | `str2DNA` copies from a thread-local record instead of moving the reader's strings. The reader keeps its buffers and workers allocate from their own arenas. | `IO.cpp` |
| 4 | Allocation-free header check | `validateSynchronizedHeaders` compares `string_view`s (`normalizedRecordIdView`). | `InputStream.cpp` |
| 5 | Remove `updateMaxSeqMutex` from the hot path | Relaxed atomic check first. The mutex is taken only when a new maximum length is seen. | `Filters.cpp` |
| 6 | Parallel serial tail | `writeDereplDNA` formats abundance-ordered chunks on the pool into recycled 1 MiB segments and writes them in order. Shard tables are cleared on the pool. | `containers.cpp`, `DNA.{h,cpp}` |
| 7 | Feed the workers | Read-ahead thread per plain-gzip input (4 × 256 KiB). Twice as many block slots as workers, with a completion queue instead of waiting on one fixed future. Per-thread filters and mergers indexed by worker. | `InputFileStream.cpp`, `IO.cpp`, `ThreadPool.h` |
| 8 | Keep outputs open across input files | `OutputStreamRetention` keeps main output streams open while the `OutputStreamer` is rebuilt per input file. BGZF append scans block headers with seeks instead of reading every block. | `OutputStreamer.cpp`, `IO.cpp`, `OutputFileStream.cpp` |

Scope limits:

- **Mitigation 2:** workers do not compress their own batches. The larger pool was
  enough to make BGZF runs compression-bound (see below).
- **Mitigation 5:** the per-thread `ReportStats`/`PrimerVariantTracker` mutexes were
  left alone. They are uncontended.
- **Mitigation 6:** `-derepPerSR` still pauses producers once per SequencingRun.
- **Batching exceptions:** with `-maxReadsPerOFile` or `-maxRdsOut`, output uses the
  per-record path, so file splitting and read limits stay exact.
- **Output retention:** it is off when read subselection, `-maxRdsOut` or
  optimal-cluster output is active.
- **Variant-export modes:** coarse modes that export differences or qualities first
  kept formatting the derep output serially. Chunked formatting kept 80–120 MB of
  `diffDNA` temporaries alive in worker arenas, and that amount grows with the input.
  The follow-up formats them in parallel with bounded memory.

## Deterministic dereplication

Before, a failed read joined an exact key only if a passing copy had already been
applied, and representatives depended on arrival order. The old binary's 8-thread
map totals therefore changed from run to run (282,544 vs 282,552 in one pair of runs).

Now dereplication is applied in input order per shard:

1. **Numbering.** The reader numbers each read block in input order
   (`OutputStreamer::reserveDerepBlock`).
2. **Worker side.** A worker processes the block and collects its derep entries in a
   thread-local batch. It does not flush mid-block. `endDerepBlock` seals the block,
   also when processing throws.
3. **Submission.** The thread that seals the next expected block submits every
   contiguous sealed block to the per-shard queues of `PartitionWorkPool`, in block
   order. The pool's new `afterEnqueue` hook releases the submission lock only once
   the tasks are queued.
4. **Execution.** The pool runs each partition's tasks in FIFO order, so every shard
   sees its entries in input order. Different shards still run in parallel.
5. **Global coarse index** (one partition: `-derepGlobal 1`, reassignment, 100% with
   retained qualities or differences):
   - Workers search the candidates in parallel against a read-only snapshot, under a
     shared lock.
   - The ordered commit takes the exclusive lock and searches only the anchors created
     after the snapshot.
   - It keeps the first-created (smallest-id) match rule.
   - Per-target updates are applied afterwards, in commit order per target. The
     follow-up makes both steps non-blocking.

Checks:

- **Byte-identical at 1, 4 and 8 threads** (8 threads run twice): `apong_fa`,
  `apong_fq`, `pe_lotus`, `se_demux`, `apong_global` and `apong_id100`. All also match
  the old binary at one thread.
- **97% coarse (`apong_coarse`):** identical at 2, 4 and 8 threads.
  - One thread uses the global index by design, so in general it can choose
    different representatives.
  - On apong its output still matches the 8-thread output.
- **Regression tests:** `multithreaded_derep_is_deterministic` in
  `tests/test_regressions_cli.py` compares the derep outputs at 1/4/8/8 threads, with
  failed reads ahead of their passing copies. `testSnapshotCommits` in
  `tests/TestCoarseDerep.cpp` checks that snapshot probes plus delta commits
  reproduce the serial assignments exactly.

## Output validation

`val2.sh` compared the new binary at 1 and 8 threads with the old binary at 1
thread, on all 11 scenarios:

- **New binary at 1 thread:** every output is identical. With BGZF output, only the
  block layout differs; the decompressed content is identical.
- **New binary at 8 threads:** every derep output is byte-identical. Main and
  demultiplexed FASTQ contain the same records with mates in sync.
- **One exception:** the demultiplexed SMPL92 files in `apong_fq`, which hit their
  `-DemultiBPperSR` cap (see Open items). The old binary behaves the same way.

## Benchmarks

Same host and datasets as the audit (Core Ultra 7 258V, 4 P + 4 LP cores, WSL2),
release build with frame pointers. The old and new binaries alternate within each
repeat. Values are medians of 3 wall-clock runs (s). CPU is user + sys; RSS is peak.

| Scenario | t | old | new | speedup | CPU old → new | RSS MB old → new |
| --- | ---: | ---: | ---: | ---: | --- | --- |
| pe_rewrite | 1 | 2.31 | 2.29 | 1.01× | 2.3 → 2.3 | 22 → 22 |
| | 4 | 3.03 | 0.71 | 4.27× | 7.3 → 3.6 | 24 → 33 |
| | 8 | 3.59 | 0.70 | **5.13×** | 13.5 → 4.4 | 27 → 56 |
| pe_rewrite_gz | 1 | 12.25 | 11.38 | 1.08× | 12.2 → 11.4 | 23 → 23 |
| | 4 | 7.56 | 3.69 | 2.05× | 22.4 → 19.1 | 28 → 50 |
| | 8 | 7.79 | 3.16 | **2.47×** | 25.6 → 21.3 | 33 → 72 |
| pe_lotus | 1 | 11.07 | 10.54 | 1.05× | 11.0 → 10.5 | 254 → 245 |
| | 4 | 5.72 | 4.66 | 1.23× | 18.3 → 18.3 | 259 → 269 |
| | 8 | 5.71 | 3.87 | **1.48×** | 27.4 → 25.7 | 265 → 299 |
| pe_lotus_gzin | 1 | 13.78 | 13.31 | 1.04× | 13.7 → 13.3 | 255 → 246 |
| | 8 | 6.10 | 4.31 | 1.42× | 28.3 → 28.4 | 269 → 311 |
| pe_lotus_gzout | 1 | 19.49 | 19.14 | 1.02× | 19.4 → 19.1 | 255 → 246 |
| | 8 | 9.47 | 5.76 | 1.64× | 39.8 → 39.0 | 273 → 321 |
| se_demux | 1 | 2.75 | 2.77 | 0.99× | 2.7 → 2.8 | 120 → 120 |
| | 8 | 1.91 | 1.17 | 1.63× | 8.5 → 6.9 | 126 → 157 |
| apong_fa | 1 | 5.64 | 5.26 | 1.07× | 5.6 → 5.2 | 174 → 175 |
| | 4 | 2.02 | 2.01 | 1.00× | 7.0 → 7.7 | 180 → 202 |
| | 8 | 1.97 | 1.52 | 1.30× | 10.1 → 10.4 | 187 → 221 |
| apong_fq | 1 | 5.92 | 5.70 | 1.04× | 5.9 → 5.7 | 190 → 191 |
| | 8 | 2.38 | 1.58 | 1.51× | 11.3 → 10.9 | 202 → 240 |
| apong_coarse | 1 | 6.27 | 5.76 | 1.09× | 6.3 → 5.8 | 378 → 377 |
| | 4 | 1.96 | 1.88 | 1.04× | 5.1 → 5.3 | 392 → 415 |
| | 8 | 2.02 | 2.03\* | 1.00× | 7.1 → 9.0 | 398 → 435 |
| apong_global | 1 | 6.15 | 5.80 | 1.06× | 6.1 → 5.8 | 378 → 377 |
| | 8 | 4.70 | 3.50 | 1.34× | 7.3 → 7.4 | 387 → 422 |
| apong_id100 | 1 | 4.55 | 4.28 | 1.06× | 4.6 → 4.3 | 392 → 392 |
| | 8 | 2.78 | 2.54 | 1.09× | 5.5 → 5.4 | 402 → 434 |

\* One outlier repeat. The 5-repeat rerun below gives 2.01 → 1.60 s.

All 4-thread results and per-run ranges are in
[benchmarks/mt-audit/ab.csv](benchmarks/mt-audit/ab.csv). These absolute times are
lower than those in the audit because the host was faster during this session, so
compare within this table only.

`apong_coarse` at 8 threads:

- The first new-binary repeat was a 3.62 s outlier. The `apong_global` run next to it
  was slow too (5.66 s). Both used about 1.5× their usual CPU time.
- The other two repeats took 2.03 s and 1.68 s, against 2.01–2.06 s for the old binary.

A second run with 5 alternating repeats settled this
([ab_coarse.csv](benchmarks/mt-audit/ab_coarse.csv)):

| Scenario | t | old | new | speedup | CPU old → new |
| --- | ---: | ---: | ---: | ---: | --- |
| apong_coarse | 4 | 2.01 | 1.86 | 1.08× | 5.3 → 5.2 |
| | 8 | 2.01 | 1.60 | **1.26×** | 7.1 → 6.7 |
| apong_global | 4 | 4.54 | 3.57 | 1.27× | 7.0 → 7.6 |
| | 8 | 5.24 | 3.81 | 1.38× | 8.1 → 8.0 |

- At 8 threads the new binary was faster than the old one in every repeat of
  `apong_coarse`, by 1.2–1.4×.
- The host slowed down in repeats 4–5 for both binaries: `apong_global` took 11.4 s
  (new) and 13.3 s (old) in repeat 4.

### Where worker time goes now (8 threads, light profile)

| Scenario | running (old → new) | lock wait (old → new) | Now limited by |
| --- | --- | --- | --- |
| pe_rewrite | 33% → 52% | 44% → 25% | the reader (76% busy, 28% in `read()`) |
| pe_rewrite_gz | 16% → 18% | 61% → 73% | deflate: the 8-thread gzip pool is 70% busy |
| pe_lotus | 39% → 91% | 30% → 0.8% | worker CPU |
| pe_lotus_gzin | 38% → 85% | 24% → 1.4% | worker CPU |
| apong_fa | 64% → 92% | 2% → 0.4% | worker CPU |
| apong_coarse | 33% → 55% | 1% → 0.3% | serial tail: 0.75 s of 1.81 s |
| apong_global | 15% → 23% | 41% → 6.5%, plus 45% waiting for ordered commits | the single global partition |

In `pe_rewrite_gz`, workers now wait in `flushWriteBatch` because the BGZF stream
blocks until the compressors catch up. The CPUs are saturated, so a faster deflate
(libdeflate, zlib-ng) is the next step there, not locking.

Glibc's internal lock waits still take 10–16% of worker samples in `pe_lotus`,
`apong_fa` and `apong_coarse`. The workers no longer free strings the reader
allocated, so this is ordinary arena and `mprotect` traffic. A different allocator
or a larger tcache is the remaining lever.

## Follow-up: global coarse commits and variant export

Two open items from the first round, implemented afterwards.

### Non-blocking global-coarse commits

With a single global partition (`-derepGlobal 1`, reassignment, 100% with retained
qualities or differences), every worker waited in `apply_blocks_ordered` until its
block had been committed. That was 13.1 s of 29.1 s worker time at 8 threads.

Handing the blocks off alone did not help (0.89–0.99×). A profile of that version
showed why:

- The exclusive commits took only 0.37 s in total.
- Applying their updates (representative, counts and quality evidence per cluster)
  took 2.5 s, all on the one committing thread.
- The committer waited another 1.1 s to get the write lock past the probes, because
  the lock prefers readers.

The final design in `containers.cpp`:

1. **Committer.** A worker seals its block and puts it on a FIFO commit queue in
   block order. The worker that finds no commit running becomes the committer and
   drains the queue; the others return to their next block.
2. **Exclusive section.** Per block, the committer takes the write lock once, runs
   `CoarseDereplication::add` with the probe (exact-key recheck, delta search), and
   collects the per-target updates.
3. **Update lanes.** Updates go to 64 lanes by target. Any thread may apply a ready
   lane, but only one at a time per lane, so every cluster's updates stay in commit
   order. Probes run meanwhile, since updates do not touch the index.
4. **Backlog.** Queued blocks and updates hold reads in memory. When more than
   2 × threads blocks, or 256 × threads updates, are pending, producers apply lanes
   instead of starting a new block.
5. **Probe chunks.** Probes take the read lock in chunks of 8 and step aside while a
   commit waits, so the committer waits at most one chunk.
6. **Freeing.** Committed entries and applied updates go back to the worker that
   allocated their reads, which frees them before its next block. Frees in a foreign
   malloc arena had cost 11% of worker samples in lock waits.
7. **Lifecycle.** The committer returns only when the queue and all lanes are
   empty, so a lifecycle pause (output, `reset`) waits for it as before.

### Parallel variant export with bounded memory

The naive parallel export uses more memory as the input grows:

- **Chunk buffering.** Chunks were a fixed fraction of the total output (1/(32 × pool)),
  so the in-flight chunks buffer about 4% of all output.
- **Decoding.** `diffDNA::writeFastq` sorts a cluster's variants, so it decodes every
  exact variant of the cluster (about 700 B each for paired reads). Formatting
  several large clusters at once multiplies that peak. Each worker's malloc arena
  also keeps its own peak.

The new plan in `writeDereplDNA`:

- **Small clusters:** chunks of at most 1024 reads. A cluster has at most as many
  variants as reads, so chunk output and per-worker decoding stay small.
- **Large clusters** (more than 1024 reads):
  - They are written by the output thread one at a time, in order.
  - Their variants are decoded once, into one sorted snapshot, as the serial path
    does.
  - Decoding and sorting run as pool tasks (sorted runs merged pairwise). Variants
    are unique, so the order is the serial one.
  - The snapshot is formatted in 1024-record ranges on the pool.
- **Window:** at most pool size + 2 chunks or ranges are buffered.

Peak RSS at 8 threads on apong and on copies scaled 2× and 4× with extra distinct
variants (`scale_fastq.pl`). Values are MiB from three runs per cell
([memscale.csv](benchmarks/mt-audit/memscale.csv),
[memscale-final.csv](benchmarks/mt-audit/memscale-final.csv)); differences are to the
serial export in the same run. Output was byte-identical in every run.

| Input | Serial export | Naive parallel | Bounded parallel |
| --- | ---: | ---: | ---: |
| 97% coarse ×1 | 438–441 | +74 to +78 | −1 to +5 |
| ×2 | 945–955 | +100 to +113 | −6 to +12 |
| ×4 | 1,835–1,839 | +153 to +156 | −0 to +7 |
| global coarse ×1 | 420–423 | +79 to +82 | +13 to +18 |
| ×2 | 904–909 | +105 to +111 | +4 to +12 |
| ×4 | 1,763–1,767 | +143 to +148 | +9 to +14 |

The bounded rows vary by about ±12 MB between runs and do not grow with the input.

- The 4× bounded runs are as fast as naive; both are limited by the single output
  thread writing about 1.3 GB of variant FASTQ.
- The small extra in the global modes is the commit backlog, not the export.

The record writer was reworked too:

- `writeQualityRecord` builds each FASTQ record directly into a pre-sized buffer,
  instead of a temporary `DNA` per variant.
- A snapshot decodes each variant once.
- Clusters with a single variant skip the snapshot.

A microbenchmark of the serial `writeFastq` (`bench_export.cpp`: one 20,000-variant
paired cluster, and 4,000 clusters of 1–3 variants) went from 28 to 22–25 ms and from
5.6 to 3.2 ms.

### Results

Previous state (the table above) against the final build. Hosts were noisy, so the
ranges span three A/B runs of 3–5 alternating repeats
([ab_final.csv](benchmarks/mt-audit/ab_final.csv),
[ab_final3.csv](benchmarks/mt-audit/ab_final3.csv)):

| Scenario | 1 thread | 4 threads | 8 threads |
| --- | --- | --- | --- |
| apong_coarse (97%, 64 partitions) | 1.05–1.07× | 1.09–1.19× | 1.15–1.37× |
| apong_global (`-derepGlobal 1`) | 1.02–1.09× | 1.44–1.61× | 1.43–1.65× |
| apong_id100 (100% + qualities) | 0.97–1.24× | 1.35–1.43× | 1.36–1.64× |

- The other eight scenarios do not use the changed code. Their medians differ by
  0.75–1.19×, within the run-to-run spread; `pe_rewrite`'s old runs alone ranged
  0.85–2.77 s.
- At 8 threads, unpinned, the global modes report 20–45% more CPU-seconds. This is
  mostly how this laptop accounts CPU time, not extra work; see
  [CPU time on this host](#cpu-time-on-this-host).
- Derep outputs are byte-identical to the previous build at 1 and 8 threads, and to
  the serial export.
- All suites pass in release and TSan builds (23/23, no TSan warnings), as do the
  9 TSan scenario runs.
- New tests in `tests/TestCoarseDerep.cpp`:
  - `testParallelVariantSnapshot`: parallel snapshots written in ranges equal the
    serial export, single-end and paired.
  - `testParallelGlobalDerep`: 1, 2, 4 and 8 workers sealing blocks in shuffled order
    write byte-identical output, including a cluster with more than 4096 variants.

### CPU time on this host

CPU-seconds are not a measure of work on this machine. A fixed CPU-bound loop costs
0.65 CPU-s alone, 0.80–0.89 CPU-s with 4 copies running at once, and 1.44–1.46 CPU-s
with 8 copies. Under full load, Hyper-V places vCPUs on the low-power cores and turbo
drops.

The old commit path left workers blocked for about 45% of the time, so fewer threads
ran at once and each ran faster. The new path keeps all 8 busy, so its CPU-seconds
inflate even for identical instructions. WSL2 exposes no hardware counters
(`perf_event_open` fails; `icount.c`), so the comparison below uses matched
concurrency: both builds at `-threads 8`, pinned to the same vCPUs, which they both
saturate.

| Pinned to | Scenario | Wall old → new (s) | User CPU old → new (s) |
| --- | --- | --- | --- |
| 2 vCPUs | apong_global | 4.68 → 3.92 | 6.65 → 6.72 (+1%) |
| | apong_id100 | 3.84 → 2.99 | 5.35 → 4.85 (−9%) |
| 4 vCPUs | apong_global | 4.22 → 2.81 | 7.60 → 7.93 (+4%) |
| | apong_id100 | 2.84 → 2.10 | 5.14 → 5.40 (+5%) |

Sys time was equal. At matched concurrency the new path does about the same work
(−9% to +5%) and finishes 1.2–1.5× sooner.

Three candidate overheads were tested and rejected (5 alternating repeats each,
`-threads 8`, unpinned and pinned to 4 vCPUs):

- **Probe gate spinning.**
  - `std::atomic::wait` spins and calls `sched_yield` before it sleeps. That showed
    up as 11% of worker samples when pinned.
  - Pinned CPU was equal with a blocking turnstile or a condition-variable broadcast
    (7.20 vs 7.24 s; 5.10 vs 5.11 s). The yields ran when no other thread needed the
    core.
  - Unpinned, both sleeping gates were about 10% slower (`apong_global` 1.79 → 2.04 s,
    with more CPU-seconds). A sleeping probe pays the futex, vCPU wake-up and
    clock ramp-up latency after each of the ~700 commits.
- **Read-lock traffic.** Probe chunks of 32 instead of 8 (4× fewer shared-lock
  round trips) changed nothing (7.22 vs 7.27 s; 5.22 vs 5.07 s).
- **Kept as is:** the spinning gate and chunks of 8.
- **Residual:** about 0.2 s of malloc lock waits per run. Update vectors are
  allocated by the committer and freed by the reads' owners.

## Open items

- **Probes wait for commits (global modes).** Probes pause during each exclusive
  commit, about 0.33 s in total on apong at 8 threads.
  - Deferring new-target initialization to the lanes would shorten the exclusive
    section.
  - An append-only index readable without the lock would remove the wait.
- **Output-thread bound tail.**
  - One thread writes every derep output file.
  - The clusters are sorted by abundance serially before output.
  - The single global index is freed on one thread (`CoarseDereplication::clear`,
    about 0.3 s on apong).
- **Main output order.** Still not input-ordered at more than one thread. The block
  numbers used for dereplication could also order the output batches, at the cost of
  holding back finished blocks.
- **Output caps are timing-dependent.** When a sample reaches its `-DemultiBPperSR`
  cap, which reads make it in depends on thread timing. In `apong_fq`, SMPL92 hits
  the cap: every run writes 134,284 records, but 236–320 of them differ from the
  single-threaded run. The old binary at 8 threads differs in 709. Input-ordered
  output batches would make this deterministic too.
- **Extra threads and buffers:**
  - Threads: one async writer, one read-ahead thread per plain-gzip input stream, and
    a gzip pool as large as the worker pool (idle unless BGZF is written).
  - Buffers: up to 16 MiB queued for the writer, and up to 4 MiB of batch per worker.
- **BGZF output container:** with retention, the main BGZF output of a multi-file run
  is one continuous stream rather than one segment per input file. The decompressed
  content is identical.

## Reproduction

Tools are in [.claude/archive/mt-audit/tools](../.claude/archive/mt-audit/tools/):

- `run_ab.sh "SCENARIOS" "THREADS" REPEATS`: alternating old/new timing into
  `runs/$TAG.csv`; `OLD`/`NEW` select the binaries. `PIN=0-3` pins both builds to those
  vCPUs. `absum.py` or `cpusplit.pl` (median wall/user/sys) summarizes it.
- `selfsec.pl ROLE A B`: per-function on-CPU seconds of one thread role in two
  profile reports, largest difference first.
- `icount.c`: counts user-space instructions of a command through `perf_event_open`.
  It needs hardware counters, which WSL2 does not provide.
- `cmpbin.sh "SCENARIOS" THREADS A B`: byte-compares every derep output of two
  binaries, with wall time and peak RSS.
- `scale_fastq.pl K IN OUT`: K copies of a FASTQ file, each copy with one
  substitution per read, so the clusters gain new distinct variants.
- `memscale.sh`: peak RSS of the serial, naive-parallel and bounded builds on the
  scaled inputs.
- `bench_export.cpp`: `diffDNA::writeFastq` microbenchmark, linked against a
  build's object files.
- `det.sh`: runs the derep scenarios at several thread counts and compares their
  outputs byte for byte.
- `val2.sh`: compares new outputs at 1 and 8 threads with the old binary at 1 thread,
  using `validate.py`. Derep outputs must be identical; other outputs must contain the
  same records with mates in sync.
- `tsan_runs.sh`: the scenario runs under ThreadSanitizer (run with
  `setarch x86_64 -R`).
- `restore_eol.pl REF WORK FILE...`: restores each file's original line endings and
  trailing whitespace after editing the mixed CRLF/LF sources.

Raw data in [benchmarks/mt-audit](benchmarks/mt-audit/):

- A/B timings: `ab.csv`, `ab_coarse.csv`, `ab_final.csv` and `ab_final3.csv`.
- Profiles of the first-round binary: `profiles-light-new.txt`.
- Peak-RSS scaling: `memscale-first.csv`, `memscale.csv` and `memscale-final.csv`.
- CPU at matched concurrency: `cpu-pinned.csv`.
- Gate and chunk-size experiments: `cpu-gate-experiments.csv`, whose first field
  is the experiment name.
