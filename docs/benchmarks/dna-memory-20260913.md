**DNA / DNAunique memory implementation and validation — 13 September 2026**

All five recommendations were implemented against a captured copy of the current working tree, including its existing quality-difference and seed-statistics work. The baseline is not Git HEAD. Raw quality values and their production representation were preserved.

**Implementation**

1. `DNAQualityEvidence` reuses the existing 64-bit sum arithmetic and overflow checks. Equal-length observations use one 64-bit coverage count. A differing observation length or nonuniform evidence source materializes the existing per-position vector. Tail evidence survives replacement and pseudo-trimming. Sample abundance is never substituted for coverage. The three averaged-quality writers share an output helper that selects coverage once per read; singleton coverage skips division by one and rounding of an already integral value.
2. `DNA` keeps its output header equal to the original implicitly until edited. An explicitly empty edited header remains distinct. All active writers and header transformations use the same lazy accessors. `saveMem()` reuses its existing header allocation and shrinks buffers only when capacity exceeds twice the used length with at least 64 bytes of slack.
3. Detached merge history is an 80-byte quality-evidence structure, instead of a complete 416-byte baseline `DNAunique`. It retains byte-sized raw evidence until promotion. Exclusive old attachments donate their buffers; externally held old attachments keep their data. Output snapshots copy history independently.
4. Aligned member grouping shrinks `DNA` from 240 to 224 bytes. `DNAunique` shrinks from 416 to 408 bytes: the base saves 16 bytes, while scalar coverage adds 8. `ElementsDetection` shrinks from 16 to 12 bytes. No packed structs, narrowed counters, or unaligned fields were introduced.
5. Parsing constructors skip the speculative 151-byte sequence reserve, string constructors move owned parameters, subsequence extraction moves its quality vector, physical trimming erases in place, and the move-aware FASTQ parser validates before consuming the header to avoid an error-only copy. Persistent-storage compaction retains the full backing sequence and raw qualities needed for later HQ output.

The sum vector remains available under its existing public name and its `uint64_t` element type. Quality history arithmetic is shared instead of duplicated. Full rebuilds are required because object layout changed.

**Before / after measurements**

Five independent runs of each version per workload, alternating version order and shuffling workload order, pinned to Linux CPU 0. These are medians, not best runs. RSS is `/proc/self/status` `VmHWM`; MiB is 1,048,576 bytes. Time-change percentages below zero mean faster. All 120 runs completed and output checksums matched across versions and repetitions.

Host: Intel Core Ultra 7 258V, x86-64 WSL2 Linux 6.18.33.2, approximately 16 GiB RAM. Both builds used GCC 13.3.0, `-O3 -std=c++20 -pthread`, `WITH_HTS=0 STATIC=0`. No build, sanitizer, or competing benchmark job ran during these measurements.

| Workload | Before MiB | After MiB | RSS change | Before s | After s | Time change |
|---|---:|---:|---:|---:|---:|---:|
| 300,000 retained DNA reads, 150 bp | 236.1 | 204.1 | -13.6% | 0.422 | 0.354 | -16.1% |
| 1,000,000 reads: parse, trim, subsequence, write | 4.2 | 4.3 | +0.4% | 1.133 | 1.014 | -10.5% |
| 200,000 unique reads, 150 bp, equal coverage | 642.1 | 401.1 | -37.5% | 0.787 | 0.548 | -30.4% |
| 200,000 unique reads, 150 bp, variable coverage | 642.2 | 633.0 | -1.4% | 0.695 | 0.666 | -4.1% |
| 10,000 unique reads, 10,000 bp, equal coverage | 1727.3 | 963.8 | -44.2% | 2.765 | 1.605 | -42.0% |
| 100,000 detached promoted merge histories | 665.0 | 327.8 | -50.7% | 0.789 | 0.579 | -26.6% |
| 100,000 detached unpromoted merge histories | 204.2 | 127.9 | -37.4% | 0.499 | 0.443 | -11.3% |
| Exact FASTQ dereplication: 200,000 × 200 bp | 11.8 | 11.3 | -4.3% | 0.305 | 0.290 | -4.9% |
| Exact FASTQ dereplication: 30,000 × ~1,500 bp | 185.8 | 168.8 | -9.1% | 0.508 | 0.490 | -3.6% |
| Exact FASTQ dereplication: 100,000 × ~1,500 bp | 598.1 | 522.0 | -12.7% | 1.708 | 1.625 | -4.8% |
| 97% FASTA dereplication: 200,000 × 200 bp | 5.7 | 5.7 | +1.0% | 0.337 | 0.309 | -8.3% |
| 97% dereplication + variants/qualities/reassignment: 20,000 × ~1,500 bp | 61.4 | 59.6 | -2.8% | 3.546 | 3.486 | -1.7% |

Object workloads retain all objects until measurement, then exercise production output snapshots/writers. Each accumulated object has the initial observation plus four updates. Detached histories use 300-base merges for 150-base representatives. The streaming case retains no read collection. Its peak RSS is primarily process overhead; small differences there are noise. Object timings include deterministic input construction and output hashing, but exclude destruction of the retained collection.

Production workloads reuse `TestCoarseBenchmark` and its existing conserved/variable sequence simulation. They include dereplication and output writing, with FASTQ enabled for the exact cases. The short exact case has low error and high duplication. The longer exact case uses ±100 bp length variation. The final coarse case retains variant qualities and performs reassignment. Output hashing occurs after the production timer and covers generated sequence, quality, map, and variant files. These production cases do not retain hundreds of thousands of promoted objects, so their RSS savings are smaller than the deliberately memory-heavy object cases.

The raw CSV includes every measurement; the summary CSV also records minimum and maximum runtime. This is a synthetic workload evaluation on one host, not a claim about every dataset or machine. Small changes near zero should be treated as run-to-run variation.

**Validation**

- Before and after: the portable `make test` suite (input helpers, regression tests, search/edit-distance oracle tests, coarse dereplication, ONT pipeline, and retained-quality pipeline).
- Final optimized code: the same suite under AddressSanitizer and UndefinedBehaviorSanitizer (`-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -fno-pie`, linked with `-no-pie`), with no reported sanitizer errors.
- Before and after, plus a final after-only run: `test_audit_pipeline.py`, `test_thorough_audit.py`, `test_merged_pipeline.py`, `test_seed_extension.py`, and `test_prefix_derep.py`.
- New regression checks: lazy and explicitly empty headers; independent header copies; physical and logical trimming; repeated storage compaction; an independent per-base sum/coverage oracle over 80 observations; uniform and nonuniform 64-bit overflow; all quality values from -5 through 127 at FASTQ offsets 33 and 64 in the singleton fast path; old merge attachments held by callers; raw-history promotion and snapshot independence.
- All benchmark outputs matched. Production checksum coverage includes HQ output and sample maps, not just cluster counts.

HTS-specific tests were not rebuilt: neither HTSlib headers/libraries nor the earlier `/tmp/sdm-hts-deps` installation are available in this environment. Portable parsing, signed quality boundaries, paired reads and merging were tested. No new quality-domain restrictions were introduced.

**Should raw `qual_` use a different element type?**

Keep `signed char` for now. It already occupies one byte. `uint8_t` also occupies one byte and would require preserving the existing handling of negative values during offset conversion. The existing code supports more than Q0–Q60: HTS decoding accepts values through Q93; FASTA/QUAL parsing accepts through `maxSAqualP` (109); regression tests exercise negative intermediate scores and values through 127. Ignoring those exceptions would change accepted inputs, coverage or emitted qualities. Merged overlap quality itself is capped at Q40 in the current merger; that does not constrain non-overlap input qualities.

Packing Q0–Q60 into six bits can theoretically reduce the raw quality payload by 25%. For a 150-base read, that is 150 → 113 bytes before representation tags/padding. It requires a packed container and changes to consumers of `getQual()`, not simply another vector element type. Four-bit storage cannot represent all 61 scores losslessly.

The isolated prototype compares signed-byte, unsigned-byte and six-bit storage on 60 million Q0–Q60 values, eight scan/reversal passes, five independent processes per representation, randomized order, pinned to CPU 0. Values above Q60 are absent, so this is the favorable case with no exception-handling cost. Final checksums match for all representations. Results are medians:

| Quality representation | Payload MiB | Peak RSS MiB | Initialization s | Scan/reverse s |
|---|---:|---:|---:|---:|
| signed_byte | 57.22 | 60.90 | 0.118 | 0.219 |
| unsigned_byte | 57.22 | 60.90 | 0.116 | 0.213 |
| packed_six | 42.92 | 46.56 | 0.190 | 1.267 |

The six-bit prototype's scan/reversal workload was **5.79×** as slow as signed-byte storage. This is a straightforward packed implementation, not a lower bound on a specialized SIMD implementation. It does not justify changing production storage under the requirement to avoid substantial runtime increases. A future lossless packed design would need a byte fallback or explicit exception encoding and benchmarks of the actual filter/merge/output paths. Skipping exceptions would be a separate behavioural change.

The 64-bit quality *sums* also cannot be narrowed merely because each individual input quality is small: accumulated depth can exceed 32-bit totals. Uniform coverage removes the larger redundant allocation while preserving those limits.

**Artifacts and reproduction**

- `docs/benchmarks/dna-memory-before-after.csv`: all before/after measurements.
- `docs/benchmarks/dna-memory-summary.csv`: medians and runtime ranges.
- `docs/benchmarks/dna-quality-storage.csv`: representation experiment.
- `docs/benchmarks/dna-memory-manifest.json`: source hashes, compiler and conditions.
- `.build/dna-memory-20260913/before-source.tar.gz`: captured baseline sources, with the same benchmark harnesses added.
- `.build/dna-memory-20260913/before` and `after`: exact benchmark executables and application binaries.
- `.build/dna-memory-20260913`: validation and benchmark logs.

From the repository root, rerun the saved binaries with:

```sh
python3 tests/benchmark_dna_memory.py --before .build/dna-memory-20260913/before --after .build/dna-memory-20260913/after --output /tmp/dna-memory-repeat.csv --repeats 5 --cpu 0
python3 tests/benchmark_quality_storage.py --binary .build/dna-memory-20260913/TestQualityStorageBenchmark --output /tmp/dna-quality-repeat.csv --cpu 0
```

To rebuild the production benchmark harnesses in either source tree:

```sh
make WITH_HTS=0 STATIC=0 -j3 .build/portable-0/TestDNAMemory .build/portable-0/TestCoarseBenchmark
```

To rebuild the separate representation prototype:

```sh
g++ -std=c++20 -O3 tests/TestQualityStorageBenchmark.cpp -o /tmp/quality-storage-benchmark
```
