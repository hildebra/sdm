# SDM 3.50: evaluation of 99% identity clustering

Evaluated 2026-09-12. Version **3.50 beta** marks this development stage. The candidate indexes and sensitivity probes are standalone experiments; production dereplication still uses the documented 100% modes. No 99% command-line option has been added.

**Recommendation:** retain exact-key collection, then cluster unique keys using an in-memory positional seed index, querying the rarest disjoint seeds and verifying with bounded Hamming distance. Start with packed 20-mers; the proposed sparse 15-mer grid is a useful lower-memory alternative. The current measurements favor both over the tested minimizer indexes. No external reference database is needed.

Follow-up: [1,500 bp 16S-like reads with Q20 substitution errors](99-identity-ont-simulation.md) confirm the fixed-seed advantage in a rough simulation, while showing that a 99% cutoff leaves most independently noisy reads in separate clusters.

## Identity definition and seed placement

The requested comparison starts at position zero, permits substitutions only, and covers the shorter logical clustering sequence:

```text
L = min(length(query_key), length(centroid_key))
E = floor(L / 100)
accept if mismatches(query_key[0:L], centroid_key[0:L]) <= E
```

Logical keys must obey the existing truncation/search-length policy, including merged-key caps for UPARSE. Full merged reads and HQ mates remain separate from these keys. At 199 bp the budget is one substitution; at 200 bp it is two. Floating-point threshold rounding is unnecessary. Empty keys should be rejected before clustering.

For approximately three 15-mers per 100 bp, use a shared coordinate grid at `0, 33, 66, 100, 133, 166, ...`, retaining only complete words. Do not redistribute seed positions as a fraction of each read's length: a short prefix and its longer counterpart must sample corresponding positions. Under the requested semantics, index `(position, word)` directly. Same-order or nearby-position hits are useful for a future offset/gap mode, but add candidates without improving recall here.

## The first four seeds: useful evidence, unsafe as a mandatory match

Four disjoint 15-mers cover 60 bases. Requiring all four to match rejects a legitimate 99% pair whenever a substitution falls in those bases. For a 200 bp pair with exactly two uniformly positioned substitutions, the probability of passing is:

```text
C(140, 2) / C(200, 2) = 48.89%
```

In 5,000 simulated 200 bp pairs with exactly two substitutions, only **2,488 passed** the strict first-four rule. Allowing each word to occur within five bases of its expected position gave the same result. The two substitutions usually destroy the word; relaxing its position does not repair it.

A safe early rejection instead counts broken, non-overlapping seeds. Each broken seed proves at least one mismatch, so reject only after more than `E` such seeds fail. The first four therefore require at least two matches at 200 bp and one at 300 bp; at 400 bp their failure alone cannot prove rejection. Counting actual mismatches in their covered bases gives a stronger bound than counting broken words. Overlapping seeds require care because one substitution can break several of them.

For lookup, choose **any `E+1` disjoint seeds** and union their candidate lists. At most `E` can be broken, so every acceptable centroid must occur in at least one list. Choosing the `E+1` least frequent seeds minimizes work while retaining this guarantee. A seed absent from the index has frequency zero and is useful evidence. Index all eligible seed positions for every centroid; query only the selected lists. Fall back to exhaustive verification when too few usable seeds exist, including very short or ambiguous keys.

This guarantee concerns candidate sensitivity for a fixed centroid set. It does not make the relation transitive or make greedy clustering independent of processing order.

## Measurements on captured SDM keys

Input was the passing-key capture from the earlier SDM fixture: **142,960 records**, **37,845 exact unique keys**, lengths **200–428 bp**, all A/C/G/T. Capping the same records at 200 bp produced **14,115 unique keys**. The capture SHA-256 is `5893d861a895549eca76da56ec28aaa70da292e412e5da1f44577ea92b330aa1`.

Each experiment first collapses full exact duplicates, then sorts by descending key length, descending original abundance, and lexicographic sequence. Centroids stay fixed. A query joins the first acceptable centroid in that ordering. This guarantees that an existing centroid is at least as long as the query, so the query determines both the compared length and error budget.

Timings include seed extraction, index allocation/build, lookup, candidate enumeration, and verification. They exclude input, exact duplicate collection, packing, output construction, and production count/quality aggregation. They are **single-threaded clustering-stage timings, not whole-pipeline slowdowns**. Exact collapse plus packing separately took 77–83 ms for full keys and 28 ms in the first capped-key run.

The following medians are from four interleaved runs of each method in the positional comparison. Host: Intel Core Ultra 7 258V, WSL2, GCC 13.3, C++20, `-O3`, without `-march=native`.

| Method | Full keys, seconds | Full-key candidate checks | 200 bp keys, seconds | 200 bp candidate checks |
| --- | ---: | ---: | ---: | ---: |
| Disjoint 20-mers, rarest `E+1`, packed verification | 0.253 | 2,339,265 | 0.028 | 364,640 |
| Sparse 15-mer grid, rarest `E+1`, packed verification | 0.555 | 7,944,395 | 0.075 | 1,459,241 |
| Positional 15-mer minimizers, `w=50`, byte verification | 6.191 | 34,018,658 | 0.337 | 2,568,291 |
| Positional 15-mer minimizers, `w=50`, packed verification | 4.172 | 34,018,658 | 0.202 | 2,568,291 |

All these methods produced identical assignments: **11,388 full-key clusters** and **3,717 capped-key clusters**. Positional minimizers enumerate the union of all selected-word lists; arbitrary thinning of those lists would need a separate sensitivity argument.

In the first series, byte-to-packed verification reduced the 20-mer median from 0.303 to 0.151 seconds and the sparse 15-mer median from 0.934 to 0.458 seconds. An all-seed rejection pass before byte verification was inconsistent: 0.858 seconds for full keys and 0.125 seconds for capped keys, versus 0.934 and 0.105 seconds without it. Checking extra seeds is therefore not automatically cheaper than starting the bounded verifier.

The strict first-four variant changed **9,258 full-key assignments** and produced 14,286 clusters. These are differences from the reference greedy assignments, not a direct pairwise false-negative rate; a missed assignment can also change later centroids.

The original minimizer index ignored positions and took 16.391 seconds with 81.37 million candidate checks on full keys. Including positions reduced that to 34.02 million checks. These observations apply to the tested index and enumeration strategy, not every possible minimizer algorithm.

The 20-mer index stored **228,408 centroid/seed entries**, versus **138,522** for the sparse grid and **172,474** for positional minimizers. The 65% increase over the grid is a count of postings, not a measured increase in total memory; hash-node and vector overhead also matter.

A separate synthetic fixture used 5,000 unrelated 300 bp sequences sharing their first 115 bases. Rare disjoint seeds in the variable tails rejected every candidate before full verification. Both minimizer indexes admitted all 12,497,500 centroid/query pairs through their common early words. This illustrates why early conserved matches can be weak evidence.

The host was noisy: full-key timing varied markedly across repetitions and between series. The candidate-count advantage is more robust than an exact speedup factor. These results make modest added cost look plausible on this fixture, especially with 200 bp keys, but do not establish production runtime or memory overhead on larger datasets. Repetitive data can still produce quadratic candidate counts.

Raw results: [initial comparison](benchmarks/99-identity-candidates.csv), [positional comparison](benchmarks/99-identity-position-candidates.csv), [seed sensitivity](benchmarks/99-identity-seed-sensitivity.csv), and [synthetic counterexamples](benchmarks/99-identity-seed-counterexamples.txt).

## Minimizers and an overall sketch hash

Minimizers can support lossless candidate lookup with a sufficiently small window. Let `k` be word length and `w` the number of consecutive word starts in one minimizer window. The window covers `k+w-1` bases. With at most `E` substitutions, a common full window is guaranteed when:

```text
L >= (E + 1) * (k + w - 1)
```

Derivation: the mismatches divide the matching bases into at most `E+1` exact runs. If every run were shorter than the window span, their total plus the mismatches would be smaller than this bound. Identical full windows choose identical minimizers with consistent hashing/tie-breaking. This also works when the longer read extends beyond the compared prefix.

For `k=15, w=50`, the span is 64 bases and the bound holds throughout the captured 200–428 bp range. It also holds at 99% for every length at least 128 bp. Shorter or ambiguous queries need an explicit fallback in the tested implementation. All selected minimizer lists are retained; skipping common tokens or taking only a few rare minimizers can break this proof because substitutions outside a token can change its selection.

A sparser `w=66` approaches three seeds per 100 bp on long random sequences, but its 80-base span does not guarantee a common seed at 200 bp with two errors. The probe found both random misses (2 of 5,000 such pairs) and a concrete adversarial two-substitution example. At 200 bp, measured mean token counts were 6.32 for `w=50` and 4.57 for `w=66`; finite-read density differs from the long-sequence approximation.

A whole-sketch hash is **a candidate fingerprint, not proof of 100% identity**. The probe found different 200 bp sequences with identical selected words and positions. Differences outside the selected words can be invisible before any hash collision is considered. A longer prefix extension can also add sketch entries, so whole-sketch equality is not a general prefix-equality key.

SDM already hashes exact logical sequences and checks key equality. Retain that fast path. Adding minimizer extraction and a second aggregate hash has no demonstrated benefit for exact duplicates here. Any sketch-based exact-match cache would still need full byte comparison before merging counts.

## What transfers from protal and USEARCH

Reviewed protal commit `014f4a94f7426fd09d93c040567b505a5c2e2538`.

| Technique in protal | Adaptation for SDM |
| --- | --- |
| Collect lookup ranges and process small ranges first | Use posting-list lengths to select the rarest `E+1` disjoint seeds before materializing hits. This is already represented in the prototype. [ChainAnchorFinder](https://github.com/4less/protal/blob/014f4a94f7426fd09d93c040567b505a5c2e2538/src/Core/ChainAnchorFinder.h) |
| Separate a compact exact core from packed flanking bases; count flank similarity with bit operations | Store packed flanks or a nearby packed word beside candidate IDs when it avoids fetching a whole sequence. Reject only if observed mismatches exceed the remaining budget. [Seedmap](https://github.com/4less/protal/blob/014f4a94f7426fd09d93c040567b505a5c2e2538/src/Hash/Seedmap.h) |
| Positional anchor consistency and extension | Use exact coordinate equality for this mode. Short fixed comparisons or packed extension can reject cheaply. Chaining and loose diagonals become relevant only if offsets or indels are enabled. [ChainAnchorFinder](https://github.com/4less/protal/blob/014f4a94f7426fd09d93c040567b505a5c2e2538/src/Core/ChainAnchorFinder.h) |
| Best-flank selection, ubiquity limits, early stopping | Treat these as optional heuristics, or provide a guaranteed fallback. Dropping a frequent seed or a non-best flank can discard the only valid candidate. [KmerLookup](https://github.com/4less/protal/blob/014f4a94f7426fd09d93c040567b505a5c2e2538/src/Hash/KmerLookup.h) |

The dense `SeedmapFlex` key table has `2^30` one-byte core entries plus four-byte control blocks per eight entries, about **1.5 GiB before values** in that configuration. SDM should instead index only observed centroid seeds with its existing hash containers and compact posting vectors. The lookup ideas do not require importing protal's reference database. [SeedmapFlex source](https://github.com/4less/protal/blob/014f4a94f7426fd09d93c040567b505a5c2e2538/src/Hash/SeedmapFlex.h)

This is structurally similar to USEARCH clustering: maintain an index of current centroids, prioritize plausible matches, then verify identity. USEARCH's documentation explicitly describes building that index during clustering and using word counts for prioritization rather than as an exact identity estimate. Its rejection-count limits trade sensitivity for speed; the proposed disjoint-seed rule avoids needing that tradeoff for the requested mode. [USEARCH algorithm](https://drive5.com/usearch/manual/usearch_algo.html)

## What transfers from the sylph paper

The linked paper is **Rapid species-level metagenome profiling and containment estimation with sylph**. It uses FracMinHash sampling and a statistical genome-to-metagenome containment-ANI estimator. Its duplicate handling combines alternate-base masks with a membership filter. These are useful design ideas, but the estimator does not decide exact short-read 99% membership. [Paper and Methods](https://www.nature.com/articles/s41587-024-02412-y.pdf)

My assessment: alternate-base masks could become additional candidate keys, provided enough disjoint masks survive the error budget and every candidate is verified. A cuckoo filter should be a negative lookup prefilter only; a positive result must not merge reads or drop abundance evidence. Its extra lookup is worthwhile only if measured cache savings exceed its cost.

The paper's default sampling rate, `c=200`, is too sparse for this task. For a 200 bp sequence and `k=31`, a simple independent-hash calculation gives only `170/200 = 0.85` expected retained words and approximately `(199/200)^170 = 43%` probability of no retained word. This calculation is an inference under random distinct-word assumptions, not a performance result from sylph. [Sketching Methods](https://www.nature.com/articles/s41587-024-02412-y.pdf)

The current sylph code also dispatches seed extraction to AVX2 when supported, with a scalar fallback. That is worth considering after profiling SDM's packed implementation. It is an observation from current code, not a claim about a technique described in the paper. [sketch.rs, commit cf6ee068](https://github.com/bluenote-1577/sylph/blob/cf6ee068c028f6589352c5e66014717ec5a0333f/src/sketch.rs)

## Exact verification and integration

For substitutions only, use **bounded Hamming distance**. The existing pointer overload of `ReadMerger::percentIdentity` already compares a specified length and stops after exceeding the mismatch budget. Reuse or extend that helper; the benchmark calls it as the scalar reference. The string-view overload currently requires equal lengths and is not the right entry point for shorter-prefix comparison.

For A/C/G/T, two-bit packing permits 32 bases per 64-bit word. Given XOR word `x`, nucleotide mismatches are:

```cpp
std::popcount((x | (x >> 1)) & 0x5555555555555555ULL)
```

Mask unused tail bits and stop once the total exceeds `E`. Plain XOR popcount counts differing bits, not differing bases. Preserve byte-comparison semantics for ambiguous characters with a fallback. Packing once per unique key avoids repeating this work for every candidate.

If indels become a requirement, **Edlib** is a reasonable first implementation to benchmark: bounded unit-cost edit distance with a bit-vector algorithm. Its end-gap mode must match the chosen semantics, and edit distance alone does not define the ANI denominator. **WFA2** is a candidate for longer sequences or affine gap costs. Neither is needed for the requested substitution-only check. [Edlib](https://github.com/Martinsos/edlib), [WFA2](https://github.com/smarco/WFA2-lib)

The production integration should:

1. Preserve parallel exact-key collection, shard/batch synchronization, and the existing fast path for 100% matching. Fixed-length settings still need approximate matching when 99% is explicitly selected.
2. Build the 99% stage from exact groups in a paused snapshot. Do not first collapse prefixes: their membership information can affect later approximate assignments. Reuse logical key/cap construction.
3. Freeze the centroid sequence and deterministic processing order. A read matching an existing member does not automatically match its centroid. Keep the emitted clustering prefix equal to the centroid; changing it through later quality-based representative replacement can invalidate earlier memberships.
4. Aggregate 64-bit abundances and sample counts only after verified assignment. Rejected-only evidence must not seed clusters. Keep snapshot/reset behavior and assign ambiguous reads only once.
5. Define quality handling explicitly. Exact-prefix quality aggregation cannot be blindly applied when bases disagree. Initially retain centroid qualities and matching full-read/HQ provenance; base-aware consensus would be a separate change requiring validation.
6. Measure the integrated pipeline, peak memory, and multicore behavior before choosing defaults. Compact posting storage and packed-flank rejection are reasonable next optimizations if profiling supports them.

## Reproducing the experiments

Sources: [candidate benchmark](../tests/Test99Candidates.cpp) and [sensitivity probe](../tests/Test99Seeds.cpp). They are not included in the default production test target. Build only these executables in a private directory from the repository root:

```sh
make -j4 WITH_HTS=0 STATIC=0 BUILD_DIR=/tmp/sdm-99-build \
  /tmp/sdm-99-build/Test99Candidates /tmp/sdm-99-build/Test99Seeds
/tmp/sdm-99-build/Test99Candidates --self-test
/tmp/sdm-99-build/Test99Seeds > sensitivity.csv 2> counterexamples.txt
/tmp/sdm-99-build/Test99Candidates /path/to/keys.tsv > candidates.csv
/tmp/sdm-99-build/Test99Candidates /path/to/keys.tsv \
  tile20_packed,grid15_packed,min15_w50_position,min15_w50_position_packed \
  > position-candidates.csv
```

The input is one LF-delimited record per occurrence: `1<TAB>ACGT...` for a passing logical key. Other records are ignored. The local captured reads are not bundled; the sensitivity and correctness probes generate their own deterministic sequences.

Validation passed for **1,298 unique boundary/ambiguity fixtures against exhaustive greedy clustering**, **5,000 differential packed-verifier cases**, and **35,000 generated true-99% pairs** across lengths 100–1,000 bp. The safe disjoint-seed rejection rules retained all generated pairs. The probes also assert finding the minimizer-miss and identical-sketch counterexamples. The saved experiments were compiled against the current portable sources in an isolated build; the shared SDM executable and concurrent IO work were not modified.
