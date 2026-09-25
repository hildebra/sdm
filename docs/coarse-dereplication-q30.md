# Q30, 200 bp production clustering benchmark

Measured on 2026-09-13 using **1,000,000 reads of 200 bp**, simulated Q30 substitution errors, and the production `Dereplicate` path. The observed error rate was **0.1001%** (200,200 errors in 200 million bases). No indels or length jitter were simulated.

The generator uses 50 equally sampled templates in five lineages, with conserved and variable blocks scaled from the earlier rough 16S-like model. Lineages mutate 20% of variable-block positions; ten variants per lineage mutate 2% of those positions again. This is a small synthetic community, not a biological accuracy benchmark. All reported read qualities are Q30; error positions receive no special quality penalty. The three repeats use the same templates, errors, and read order.

Times include sequence/error generation, DNA construction, production batching and clustering, the optional final pass, and writing the representative FASTA, sample map, HQ FASTQ, and enabled binary differences. They exclude template construction, FASTQ parsing, and upstream filtering. Exact-subcluster FASTA export was disabled. Every run checks that representative abundances sum to one million.

Each cell below gives the median of three separate processes. Peak memory is the executable’s Linux `/proc/self/status` `VmHWM`, in MiB; unlike `getrusage().ru_maxrss`, it avoids the launcher’s inherited pre-exec high-water floor. Memory includes the executable, retained clusters, output snapshots, and any variants/indexes used by reassignment. System load and CPU frequency were not controlled, so timing differences of a few tenths of a second should be treated cautiously.

Environment: Intel Core Ultra 7 258V, x86-64 Linux under WSL/Microsoft virtualization; GCC 13.3, `-O3 -std=c++20`, `WITH_HTS=0 STATIC=0`. One producer thread, identical small reserve settings for all modes, default seed length 20 and homopolymer cutoff 4. Private build: `.build/coarse-dev-20260912`.

## Counts only: no binary difference output

| Identity | Clusters | Streaming time (s) | Streaming peak (MiB) | With reassignment (s) | With reassignment peak (MiB) | Reads moved |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| exact | 44,728 | 1.199 | 68.61 | — | — | — |
| 99% | 4,144 | 1.956 | 11.75 | 4.488 | 42.57 | 287,873 |
| 98% | 556 | 1.546 | 6.35 | 2.683 | 22.56 | 469,235 |
| 97% | 47 | 1.363 | 5.43 | 1.770 | 17.18 | 333,957 |

## With `derepStoreDiffs=1`

| Identity | Clusters | Without reassignment (s / MiB peak) | With reassignment (s / MiB peak) | `.diff` size before / after (MiB) |
| --- | ---: | ---: | ---: | ---: |
| 99% | 4,144 | 1.993 / 16.34 | 5.243 / 47.76 | 1.210 / 1.202 |
| 98% | 556 | 1.800 / 9.89 | 2.906 / 26.25 | 0.550 / 0.526 |
| 97% | 47 | 1.703 / 9.06 | 1.994 / 20.89 | 0.530 / 0.495 |

At 200 bp, 99%, 98%, and 97% allow 2, 4, and 6 substitutions, respectively. In this dataset, 97% retained 47 representatives compared with 44,728 exact dereplicates. Higher thresholds retained more error-derived or closely related clusters and used more time and memory. The dataset is too small and homogeneous to predict memory on a highly diverse run; longer reads, error structure, and read ordering will change these results.

Reassignment preserved the number and sequences of representatives in this simulation and improved membership relative to those fixed representatives. At 97%, 333,957 reads moved as whole counted exact variants; 98% moved 469,235 and 99% moved 287,873. The move count is not a biological correctness score. The sample-count checks and exhaustive-search tests establish that counts are conserved and each move follows the requested identity rule.

Reassignment increases peak memory because it must retain exact variants and their sample counts. Optional `.diff` export adds a second set of output patches rebased against destination references. With export disabled, only output counts are rebuilt. Defaults (`derepReassign=0`, `derepStoreDiffs=0`) avoid these costs.

## Reproduce

From the repository root:

```sh
make WITH_HTS=0 STATIC=0 BUILD_DIR=.build/coarse-dev-20260912 -j4 \
    .build/coarse-dev-20260912/TestCoarseBenchmark
python3 .claude/archive/benchmarks/benchmark_coarse_q30.py \
    --binary .build/coarse-dev-20260912/TestCoarseBenchmark \
    --reads 1000000 --repeats 3
```

A single case, with identity 97%, no delta export, and final reassignment enabled:

```sh
.build/coarse-dev-20260912/TestCoarseBenchmark 97 1000000 0 0 200 30 1
```

Positional arguments are identity (`exact` or percentage), reads, delta export, length jitter, mean length, quality, and final reassignment. Default mean length and quality remain 1500 and Q20 for the earlier ONT experiment. The simulator uses 0.01% error probability increments, suitable for these Q20/Q30 runs.

Raw measurements: [coarse-production-q30-200bp.csv](benchmarks/coarse-production-q30-200bp.csv). Implementation semantics: [optional final reassignment](coarse-dereplication.md#optional-final-reassignment).

Validation included exhaustive streaming and nearest-representative oracles; variable lengths, caps, merged source selection, low-complexity fallback, sample counts, delta/FASTA exports, snapshots and resumed input, reset and per-run output; existing regression and ONT CLI checks; and address/undefined-behavior sanitizers with leak detection.
