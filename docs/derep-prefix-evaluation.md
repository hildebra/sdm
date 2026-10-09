# 100% prefix dereplication: implementation and measurements

Measured 2026-09-12. This documents the first implementation; it does not implement 99% identity matching.

## Startup policy

`derepPrefix=auto` selects the existing length-matched hash path when absolute length bounds guarantee one clustering length. In particular, an effective `minSeqLength` equal to `TruncateSequenceLength` needs no prefix index. Merged search keys share that truncation cap, so differences outside the sequence supplied to UPARSE cannot generate duplicate truncated inputs. Full merged representatives and HQ mates remain available.

With variable clustering lengths, parallel exact-key collection is followed by a sorted prefix pass when output is assembled. Identity is checked at the same start position across the shorter key. Incompatible longer sequences remain separate; an ambiguous shorter read is assigned once to the lexicographically first maximal passing key. Prefix groups are rebuilt for each snapshot without modifying the live exact-key groups.

An explicit positive `derepSrchLen` retains capped, length-matched keys in automatic mode for compatibility. `derepPrefix=1` forces prefix matching (still bounded by the clustering cap); `derepPrefix=0` preserves legacy search-length behavior. See [all option details](options.md#dereplication).

Failed and alternate-tier reads are buffered separately. A matching passing group is required to recover their counts and quality evidence. Deferred failures retain no mate representative, and retain merged quality evidence only when merged dereplication output needs it. This avoids much of the initial buffer overhead, but memory use still increases for variable-length runs.

## Measurements

The fixture contained **170,565 read pairs**. These are whole-pipeline timings, including input, filtering, merging, dereplication, and output. The same source snapshot was used for the baseline, with changes built privately using GCC, C++20, `-O3`, and the portable configuration. No shared `sdm` executable was replaced.

The fixed case used `minSeqLength=200` and `TruncateSequenceLength=200`. The variable case changed the truncation setting to `-1`; remaining filtering and abundance thresholds were preserved. Each configuration ran once for warm-up and five further times, with interleaved order. The table gives medians of those five measured runs.

| Configuration | 1 worker, seconds | 4 workers, seconds | 4 workers, peak MiB |
| --- | ---: | ---: | ---: |
| Fixed 200 bp, before | 3.73 | 1.32 | 238.8 |
| Fixed 200 bp, auto | 3.63 | 1.27 | 238.7 |
| Variable output, previous 200 bp keys | 3.93 | 1.42 | 235.1 |
| Variable output, full exact keys | 4.03 | 1.52 | 333.8 |
| Variable output, auto prefix | 4.13 | 1.57 | 401.1 |

**Timing was noisy on this host.** In the final series, individual four-worker full-exact runs ranged from 1.42 to 2.78 seconds, and prefix runs from 1.57 to 3.43 seconds. The measurements suggest modest prefix overhead, but cannot establish a precise slowdown. Fixed-length mode showed no consistent regression across the initial and final series.

The memory effect was clearer: variable prefix mode used about **20% more peak memory than full-length exact matching**, and about **71% more than the previous 200 bp key search** in this fixture. The latter also compares different matching semantics: the previous capped search ignores differences after 200 bp. Prefix mode additionally recovered 20,526 non-passing reads into passing groups. Existing abundance thresholds are applied after aggregation, so retained group counts can change.

The fixed-length path avoids sorting a prefix index and retaining deferred failures. Prefix collection remains parallel; final matching costs a sort over unique passing keys plus binary searches for unique deferred-failure keys. There is no pairwise scan over all reads and no external database to build.

## Validation

- Full C++ regression suite, including new startup selection, ambiguous-prefix, concurrent-batch, per-sample count, per-base quality, reset, and snapshot tests.
- Independent randomized pairwise-reference comparison with one and four workers, with and without truncation.
- Merged products that differ beyond the clustering cap collapse to one group while preserving a full merged representative and full HQ mates.
- A shorter failed read contributes accepted merged quality evidence only at covered positions.
- Existing abundance/quality audit (19 cases), merged-output pipeline tests, and seed-extension pipeline tests (18 cases).

The new integration test can be run with:

```sh
perl tests/test_prefix_derep.pl /path/to/sdm
```

Validation used an isolated portable build. HTS-enabled linking was not part of this test run; the concurrent HTSlib integration was left in place.
