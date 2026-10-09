# Speed audit 2: filtering, dereplication and SIMD

Audited 2026-10-09 on the working tree after [speed-fixes.md](speed-fixes.md) and the 3.61 audits.
No production source was changed. Every proposed change exists as a prototype in
[.claude/archive/speed-audit-2/prototype.patch](../.claude/archive/speed-audit-2/prototype.patch),
each behind its own `-DXn` switch, and was measured in WSL builds of that patch.

**Summary.**

- **Biggest remaining costs** are per-read passes over qualities, the merge seed search, a header
  helper that scans the whole header four times, the binomial filter, coarse-search pointer chasing
  and malloc. Locks no longer matter at these thread counts.
- **14 exact changes** were prototyped (a 15th, X10, gave nothing and was dropped). Together, at 1 thread:
  - LotuS-style derep (`apong_fa`, `apong_fq`): **1.5–1.7×**
  - filtering only: **1.55–1.6×**
  - 97% coarse: 1.3–1.4×
  - 100% (global index): 1.25×
  - BGZF output: 1.3–1.4×
  - plain rewrite: 1.25–1.4×
- **At 6 threads** wall time improves 1.1–1.4× (user CPU 1.2–1.5×). The 100% mode is bound by its
  serial commit.
- **Outputs** are byte-identical in all 8 scenarios at 1 thread, and all 23 `make test` suites pass.
  The new kernels were fuzzed against the original scalar code (0 differences in 550,000 cases).
- **SIMD:** worthwhile places need only SSE2 (x86-64 baseline, so no dispatch). AVX2 would add
  little on 200–600 bp reads, and AVX-512 nothing measurable. POPCNT is worth enabling. The largest
  remaining hot spots are serial floating-point sums, cache misses and allocation, which SIMD cannot
  fix without changing results.

## Method

- **Data:** real apong MiSeq 2×300 (336,573 pairs, 4 samples) with its `in.map` and `sdm_miSeq.txt`
  (binomial filter on, `maxAccumulatedError 0.95`, primers).
- **Scenarios** (`tools/scen.sh`):
  - `apong_fa`, `apong_fq`, `apong_coarse` and `apong_id100`, as in
    [multithreading-audit.md](multithreading-audit.md);
  - `apong_filt`: the same filtering and FASTQ output without derep;
  - `apong_filtgz`: the same with BGZF output;
  - `apong_rewrite`: one sample, no map;
  - `apong_gzin`: the same from plain gzip.
- **Host:** WSL2 on a Core Ultra 7 258V (6 vCPUs, AVX2, no AVX-512), g++ 13 `-O3`, `WITH_HTS=0`,
  libdeflate.
- **Profiles:** the `LD_PRELOAD` stack sampler from the multithreading audit, 1 ms period. Self time
  is attributed to the innermost inlined function; `tools/callers.pl` attributes libc time to sdm
  callers.
- **Timing:**
  - Interleaved A/B, 5–7 repeats, medians (`tools/ab.sh`).
  - The host throttled late in the session. The final 1-thread comparison therefore ran both builds
    at the same time, pinned to separate vCPUs and swapped every pair (`tools/conc.sh`, 4 pairs).
- **Equivalence:**
  - Every output byte-compared, `.gz` by content (`tools/cmp.sh`).
  - `make test` on the prototype build.
  - `tools/fuzz.cpp` compares each new kernel with the original code on random inputs.

## Where the time goes (before)

Share of single-thread CPU (profile self time, %):

| Area | filt | fa | coarse 97% | id100 |
| --- | ---: | ---: | ---: | ---: |
| Accumulated error: two serial `double` sums (`qualAccumulate`, `getAccumError`) | 15.9 | 9.4 | 4.7 | 6.1 |
| Binomial filter (`errorCountDistribution`, `binomialFilter`) | 11.3 | 10.6 | 2.8 | 3.3 |
| `std::string::rfind` in `remove_paired_info` | 7.2 | 2.0 | 1.7 | 2.2 |
| Average quality, quality window, FASTQ quality characters | 12.8 | 5.5 | 2.8 | 3.9 |
| Ambiguous-base and homopolymer scans | 6.1 | 6.4 | 3.3 | 4.8 |
| Merge seed search (inclusive, incl. its hash map) | – | 17.1 | – | – |
| Consensus construction (`merge`, `getMergedQual`, complement) | – | 11.1 | – | – |
| Coarse candidate verification (`std::min` + `mateKey` = cache misses; compare) | – | – | 23.0 | 6.8 |
| `encodePatch` | – | – | 3.6 | 3.5 |
| malloc/free | 2 | 12 | 19 | 25 |
| `shared_ptr` reference counts | 3 | 4 | 4 | 4 |

At 6 threads the workers of `apong_fa` are 92% busy, so per-read savings translate into wall time.
In `apong_id100` they are 46% idle, waiting for the global index's commit.

## Findings and prototypes

Gains are 1-thread medians, interleaved, each item alone against the unmodified build.

| # | Change | Files | Gain |
| --- | --- | --- | --- |
| X1 | `remove_paired_info`: test the last two characters instead of four `rfind` over the whole header | `Common.cpp` | filt 1.10× |
| X2 | Binomial filter: skip the recurrence when a tail bound already proves the read passes | `DNA.cpp` | filt 1.11×, fa 1.10× |
| X3 | Ambiguity count and homopolymer runs from SSE2 masks | `DNA.cpp` | filt 1.05× |
| X4, X7 | `getAvgQual` with `psadbw`; FASTQ quality characters with a saturating add and clamp | `DNA.cpp` | filt 1.14× |
| X5, X6, X14 | Statistics take `const shared_ptr&`; no mutex per histogram value; primer-variant read id as a view | `Statistics.*`, `Filters.cpp` | filt 1.03×, fa 1.02× |
| X8 | Merge seed search: direct 1,024-entry table of 5-mers instead of a `string_view` hash map | `ReadMerger.*` | fa 1.13× |
| X9, X15 | Consensus built by segment (copies outside the overlap), quality table hoisted, complement table | `ReadMerger.*` | fa 1.06× |
| X10 | One `try_emplace` instead of `find` + `try_emplace` for derep keys | `containers.cpp` | none, dropped |
| X12 | Coarse search reads anchors from a flat id-indexed array of (pointer, length) | `CoarseDereplication.*` | coarse 1.11× |
| X13 | `encodePatch` finds mismatch positions 16 bytes at a time | `diffDNA.cpp` | coarse 1.03× |
| X16, X17 | SSE2 lowercase check in `seal` and N count; R1 and R2 pre-filter error sums in one interleaved loop | `DNA.*`, `IO.cpp` | filt 1.10×, fa 1.06× on top of the rest |

### X1: header pair suffix

- **Cost:** `remove_paired_info` calls `rfind("/1")`, `rfind("/2")`, `rfind(".1")` and `rfind(".2")`.
  Each scans the whole header backwards, although only a match at the last two characters counts.
- **Callers:** it runs for every read through `getPositionFreeId` (primer-variant tracking, a copy
  per read) and `normalizedRecordId`.
- **Exact:** `rfind` returns the end position exactly when the pattern sits there, so checking the
  last two characters is equivalent. X14 additionally passes `getPositionFreeIdView()`, which applies
  the same rules to a view, instead of a copied string.

### X2: binomial filter shortcut

- **Bound:** for independent per-base errors with mean μ (the sum of error probabilities), P(X ≥ m+1)
  ≤ μ^(m+1)/(m+1)!.
- **Rule:** when this bound is below α by more than the float recurrence's rounding error, the
  recurrence would stop at m = `maxErr` errors or fewer. Its result then cannot exceed the caller's
  threshold (≥ `maxErr`), so the read passes.
- **Computing μ:** an integer sum of the per-quality error probabilities, rounded up in 32-bit fixed
  point. It has no serial floating-point dependency.
- **Slack:** 8·(m+2)·(L+1)·`FLT_EPSILON`, several times the recurrence's worst-case error.
- **Effect:** on apong, 77–78% of calls skip the recurrence (487,718 of 636,889 in `apong_fa`,
  263,653 of 336,573 in `apong_filt`).
- **Checks:**
  - A cross-check build (`-DX2CHECK`) recomputed the full value for every skipped call on all
    scenarios; none exceeded `maxErr`.
  - The fuzzer did the same for 33,000 skips, with α down to 1e-4 and reads up to 600 bp.
- **Caveat:** the value returned for skipped reads is 0. Its only caller compares it with the
  threshold. A production version should become a predicate, e.g. `exceedsBinomialErrors(threshold)`,
  so that nobody later reads the value.

### X3, X4, X7, X16: per-read passes with SSE2

| Pass | Kernel | Exact because |
| --- | --- | --- |
| `scanSequenceChecks` (ambiguous bases, homopolymer runs) | Equal-neighbour and non-ACGT bitmasks per 64 bases; a run of k equal pairs ends where `eq & eq<<1 & … & eq<<(k-1)` is set, with the previous word shifted in. | Same first failing position, and the same partial ambiguity count up to it. Runs longer than 65 bases fall back to the scalar loop. |
| `getAvgQual` | `max(0,q)` by mask, then `psadbw` | Integer sum, so the order does not matter. The old loop also reloaded the member accumulator every iteration (aliasing). |
| `buildFastqQualityString` | Saturating `q + offset`, clamped to [max(33, offset), 126] in the biased unsigned domain | Equal to `sanitizeFastqQualityChar` for every input. Fuzzed with offsets −5, 0, 33, 64, 100 and 127 and random negative qualities. |
| `seal` lowercase test, `numNonCanonicalDNA` | `(c - 'a') ≤ 25` unsigned; `(c | 0x20)` in {a, c, g, t} | Same predicates as the loops and the `DNA_amb` table. |

### X8, X9: read merging

- **Seed search before:** `ReadMerger::findSeed` filled a `robin_hood::unordered_map<string_view, int>`
  with up to 92 five-base seeds of read 2 for every pair. Many are duplicates, because 13 of the 23
  seed offsets are 1.
- **Lookups before:** it then hashed one 5-mer per position of read 1, about 300 hash lookups per pair.
  That map, not the derep hash table, was most of the `robin_hood` time in `apong_fa`. This explains
  why X10 gained nothing.
- **X8:** ACGT-only 5-mers index a 1,024-entry table directly. Entries carry an epoch stamp, so the
  table is never cleared. The code rolls along read 1, so each position costs one shift and one load.
  Seeds containing other bytes (N) are kept in a short list and matched by their exact bytes.
- **X8 is exact:** same seeds, first inserted value wins as with `insert`, same scan order, same
  verification. 100,000 fuzzed pairs gave identical results, as did all merges in the apong runs.
- **X9:** the consensus loop tested coverage for every position and called `getMergedQual` out of
  line. That function guards a static local on each call.
  - The prototype copies the parts outside the overlap with `memcpy` and uses a table pointer
    inside it.
  - The table is filled from `getMergedQual` itself, so its values are identical.
  - X15 replaces the `complement` switch with a 256-byte table.

### X12, X13: coarse search

- **Before:** each candidate check followed group pointer → `shared_ptr` → string → bases, twice (R1 and
  mate). The profile shows that as `std::min` (9.3%) and `mateKey` (8.7%), where the first loads miss.
- **X12:** a flat array of anchor views (data pointer, length, mate pointer, mate length), indexed by
  id, is filled when an anchor is created. Anchor strings never change afterwards, and the array grows
  under the same exclusive lock as `groups_`.
- **After:** verification touches only the bases. The remaining misses on the bases themselves now
  show as `_mm256_movemask_epi8` (8.5%).
- **X13:** `encodePatch` compared byte by byte twice (count, then emit). It now collects mismatch
  positions with compare and movemask, then emits the same bytes.

## SIMD, AVX2, AVX-512, POPCNT

- **Already in place:** AVX2 bounded mismatch counting with runtime dispatch (speed-fixes item 3),
  used by coarse verification and merge-seed verification. libdeflate brings its own SIMD and
  carry-less-multiply CRC.
- **SSE2 is enough for the new kernels** (X3, X4, X7, X13, X16).
  - SSE2 is part of x86-64, so the static binary needs no dispatch.
  - At 200–600 bp a kernel runs 13–40 iterations, and each is now under 2% of CPU. AVX2 would halve
    the iterations for well under 1% overall.
  - AVX2 versions with the existing `__builtin_cpu_supports` dispatch are worth it only if
    kilobase reads (ONT) dominate.
  - Production versions need an `#if defined(__SSE2__)` (or x86-64) guard with the scalar loop
    as fallback, as the existing AVX2 code has.
- **POPCNT:** baseline x86-64 has no POPCNT, so `__builtin_popcount*` becomes a libgcc call
  (`__popcountdi2`, 1.3% of filter-mode CPU in the prototype). Two options:
  - Mark the few kernels `__attribute__((target("popcnt")))` behind the existing dispatch.
  - Build with `-march=x86-64-v2` (SSE4.2 and POPCNT, every x86-64 CPU since about 2009).
    [speed-audit.md](speed-audit.md#9-build-options) found `-march=x86-64-v3` alone gives nothing,
    which still holds: the gain comes from kernels written for it.
- **AVX-512:** not available on this host, so not measured. Expected to add nothing measurable:
  - per-read loops are short;
  - the remaining hot spots are serial floating-point sums, cache misses and allocation;
  - AVX-512BW mask registers would only simplify loop tails.
- **What SIMD cannot do exactly:**
  - The accumulated-error sum (sequential `double` adds) and the binomial recurrence (float) are
    latency chains with a fixed rounding order. Vectorizing along the read changes the last bits,
    which can flip a threshold decision or a representative choice.
  - Exact options:
    - skip the work (X2);
    - interleave independent reads (X17 does R1 and R2; a 4-way block pass is the next step);
    - batch the binomial recurrence across reads, one read per lane (8 AVX2 / 16 AVX-512 float
      lanes), each lane in the original operation order without FMA. After X2 the recurrence is
      only 3–4% of CPU, so this is no longer worth it.
- **Hashing:** hardware CRC32C (SSE4.2) or wyhash could replace `robin_hood::hash_bytes` for long
  derep keys. That is 1.4–2% after the prototypes, so small.

## Results

Unmodified build against the prototype build (all items except X10):

| Scenario | 1 thread, interleaved | 1 thread, concurrent pinned | 6 threads, wall | 6 threads, user CPU |
| --- | ---: | ---: | ---: | ---: |
| apong_fa | 1.57× (1.66× with X16, X17) | 1.53× | 1.38× | 1.50× |
| apong_fq | 1.56× | 1.49× | 1.28× | 1.40× |
| apong_coarse | 1.31× | 1.39× | 1.14× | 1.19× |
| apong_id100 | 1.24× | 1.26× | 1.12× | 1.22× |
| apong_filt | 1.45× (1.60× with X16, X17) | 1.55× | 1.28× | 1.43× |
| apong_filtgz | 1.28× | 1.38× | 1.24× | 1.23× |
| apong_rewrite | 1.26× | 1.41× | 1.33× | 1.50× |
| apong_gzin | – | 1.17× | – | – |

- **Columns:**
  - The interleaved and 6-thread columns come from one clean session with the prototype
    without X16/X17.
  - The concurrent column uses the full prototype while the host was throttling: 4 pairs, median
    per-pair ratio.
  - A later sequential 6-thread rerun on the throttled host varied 1.3–2.3 s for identical runs,
    and is kept only as raw data.
- **Equivalence:**
  - At 1 thread, every output of all 8 scenarios is byte-identical (`.gz` by content).
  - At 6 threads the derep outputs of `apong_fa` and `apong_fq` are identical; FASTQ files match
    as sorted sets (record order varies across threads).
  - `SMPL92` under `-DemultiBPperSR` is the known timing-dependent subset.
  - `apong_coarse` and `apong_id100` differ only as in the side finding below.
- **Memory:** peak RSS is unchanged.

## Not prototyped: next steps

1. **Allocator.** malloc and free are still 15% of single-thread CPU in `apong_fa` and 22% in
   `apong_coarse`. jemalloc is the better candidate; see [Allocators](#allocators-jemalloc-and-mimalloc).
2. **Accumulated error.** `qualAccumulate` (10–12%) and the pre-filter `getAccumError` remain serial
   `double` chains. Exact options:
   - compute the pre-filter sums for 4 reads of a block at a time;
   - fold `qualWinfloat` (4–6%) into the `qualAccumulate` pass, where its integer window sum costs
     nothing next to the floating-point chain.
3. **Record construction.** After the prototypes, about 19% of filter-mode CPU is the `DNA`
   constructor (8.7%) and libc copies (10.7%): the per-record copy in `process_read_block`, `erase` of
   the header's `@`, and a zero-filled quality vector. Converting qualities straight from the line and
   skipping the `@` by offset would cut that. Workers parsing raw chunks (speed-audit item 8) is the
   larger version.
4. **`shared_ptr` by value** remains on the per-read path, for example:
   - `Dereplicate::addDNA`;
   - `ReadMerger::merge` and `findSeedForMerge`;
   - `Filters::cutPrimer`, `findPrimer` and `checkIfPrimerHits`;
   - `DNAunique::matchedDNA`.

   Reference-count atomics are still about 5% of `apong_fa`. Passing `const shared_ptr&` where
   ownership is not taken is mechanical.
5. **Coarse candidates.**
   - The seed-count filter from speed-audit item 3 is still open: 8.8 million bounded comparisons
     per apong run at 97%, 2.4 million at 100%.
   - Prefetching the next candidate's bases while merging the posting lists would hide misses.
   - The posting-list heap is 6% after X12.
6. **Writer limits at 6 threads.**
   - In a 0.5 s filter-only run, the workers wait 106 ms in total on the output batch mutex
     (`flushWriteBatch`).
   - The writer thread of a plain rewrite is 95% in `write()`.
   - Both are I/O-bound; nothing cheap here.
7. **Not worth doing:** X10, and the derep hash table generally (1–2% after X8).

## Allocators: jemalloc and mimalloc

Measured 2026-10-09 with Ubuntu 24.04's jemalloc 5.3.0 and mimalloc 2.1.2, preloaded (`LD_PRELOAD`) into
the dynamic builds of the unmodified tree and the prototype. Other sessions kept WSL at load 2–10
throughout, so:

- 1-thread speed comes from concurrent pinned pairs (`tools/conc.sh`, 4 pairs each).
- Multithreaded effects come from profile shares and long runs on scaled inputs.

**Outputs.** With either allocator, every output of all 8 scenarios is byte-identical at 1 thread,
for both builds. At 6 threads the derep outputs of `apong_fa` and `apong_fq` are identical too.

**Speed at 1 thread** (glibc time / allocator time, median of 4 pairs):

| Scenario | Unmodified + jemalloc | Unmodified + mimalloc | Prototype + jemalloc | Prototype + mimalloc | Prototype + jemalloc, decay 0 |
| --- | ---: | ---: | ---: | ---: | ---: |
| apong_fa | 1.06× | 1.07× | 1.11× | 1.13× | 1.09× |
| apong_fq | 1.08× | 1.09× | 1.13× | 1.12× | 1.11× |
| 97% coarse | 1.09× | 1.12× | 1.09× | 1.08× | 1.07× |
| 100% | 1.09× | 1.12× | 1.12× | 1.16× | 1.06× |
| filter, `.gz`, rewrite, gzip input | 1.00–1.02× | 0.98–1.01× | 1.01–1.04× | 0.98–1.08× | – |

"Decay 0" is `MALLOC_CONF=dirty_decay_ms:0,muzzy_decay_ms:0` (see Memory below).

**At 6 threads,** allocator functions take the following share of the worker threads' running samples
(unmodified build; `tools/profalloc.sh`, `tools/allocshare.pl`):

| Scenario | glibc | jemalloc | mimalloc |
| --- | ---: | ---: | ---: |
| apong_fa | 10.5% | 2.6–3.6% | 5.9–7.2% |
| apong_fq | 10.3–11.2% | 3.0–4.0% | 5.7–8.2% |
| 97% coarse | 14.4–14.8% | 4.7–4.8% | 13.0–14.1% |

- **Coarse mode:** mimalloc does not help at 6 threads. It frees many reads on threads other than the
  one that allocated them.
- **Wall time:** the 6-thread runs on 2× and 4× inputs (`tools/memalloc.sh`) were 1.1–1.25× faster
  with either allocator, by the minimum of 2 runs. These numbers are noisy.

**Memory** (peak RSS at 6 threads, maximum of 2 runs):

| Input | Coarse: glibc | +jemalloc | +mimalloc | `apong_fa`: glibc | +jemalloc | +mimalloc |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 1× | 490 MB | +10 | +9 | 185 MB | +9 | +0 |
| 2× | 1,010 MB | +21 | +8 | 506 MB | +22 | +28 |
| 4× | 1,915 MB | +38 | +22 | 975 MB | +35 | +35 |

- **Default settings:** the overhead grows with the input, from about 1% to 4%. Both allocators keep
  freed pages for a while.
- **jemalloc with decay 0:** returning pages at once cut the 4× overhead to +10 to +13 MB, about the
  same as at 1×, so it is roughly constant. It still kept most of the 1-thread gain (1.06–1.11×).
- **mimalloc:** `MIMALLOC_PURGE_DELAY=0` did not help; one coarse run peaked at +195 MB.

**Recommendation:**

- **Choice:** jemalloc with decay 0 gives about 1.05–1.12× in the derep modes with constant memory;
  filtering-only runs barely change.
- **Settings:** jemalloc reads its options from a `malloc_conf` string compiled into the binary, so the
  settings need no environment variable.
- **Linking:** `libjemalloc.a` (Ubuntu's `libjemalloc-dev`) can link into the static build. A
  `WITH_JEMALLOC=auto` switch like the one for libdeflate fits the Makefile.
- **Before making it the default,** confirm the memory overhead on a multi-gigabyte production run.

## Side finding: 97% coarse derep is not deterministic at 6 threads

- **Observed:** with the unmodified build, three pairs of identical 6-thread `apong_coarse` runs
  differed in "non-passing reads recovered into passing groups" (1,343–1,345). With it,
  `derep.map`, `derep.fas.rest` and both `.hq.fq` files differed.
- **Conflicts with:** [multithreading-fixes.md](multithreading-fixes.md), which states that 97%
  coarse partitions are identical for any thread count ≥ 2.
- **Not caused by** the prototype; the unmodified build shows it.
- **Cause and fix:** a separate session found the cause the same day: one deferred-read byte budget
  shared by all shards. It moved the budget to each shard. That fix is described in
  `multithreading-fixes.md` ("Fix 2026-10-09") in that session's worktree, and is not yet in this tree.
- **`apong_id100`:** only the "bounded comparisons" counter varies, which is expected for snapshot
  probes.

## Reproduction

Tools are in [.claude/archive/speed-audit-2/tools](../.claude/archive/speed-audit-2/tools/), raw data
in [.claude/archive/speed-audit-2/results](../.claude/archive/speed-audit-2/results/). Paths refer to
the WSL work area `~/sdm-perf` and this session's scratch directory.

- **Setup and builds:**
  - `setup.sh`: rsync of the tree into WSL, apong inputs, the profiler, release and profiling builds.
  - `build.sh NAME "FLAGS"`: builds the patched copy with any `-DXn` set, e.g. `-DX2 -DX2CHECK` for
    the binomial cross-check.
  - `prototype.patch`: made with `diff --strip-trailing-cr` against the 2026-10-09 tree. Strip CRs
    before applying, and restore line endings afterwards with
    `.claude/archive/mt-audit/tools/restore_eol.pl`.
- **Profiles:**
  - `prof.sh` writes them; `hot.pl` lists hot functions from `analyze.py`'s JSON.
  - `callers.pl` attributes libc samples to sdm call sites.
- **Comparisons:**
  - `cmp.sh A B "SCEN" T OUT`: all outputs, paths normalized, `.gz` by content.
  - `ab.sh`: interleaved timing. `conc.sh`: concurrent pinned 1-thread pairs.
- **`fuzz.cpp`:** the prototype kernels against the original code; build with `g++ -O2 -std=c++20`.
- **Results:**
  - `profiles/p1_*`, `p6_*`: unmodified build at 1 and 6 threads. `pall_*`: prototype build.
  - `timing/*.tsv`: every run.
  - `make_test_prototype.txt`: the 23 suites on the prototype build.
