# Pseudo-dereplication before Savont

Evaluated 2026-09-12. The objective is to lower SDM's live memory use and send fewer candidate reads to Savont. Some loss of closely related variants is acceptable. This changes the recommendation from strict 99% clustering to **streaming coarse grouping with capped retention of real reads**.

My initial choice would be to benchmark **97.5% grouping with 20–50 retained reads per group**, then tune the threshold and cap against actual Savont results. A more generous 32–64-read cap is a reasonable starting range for that downstream test. These are candidate settings, not validated Savont defaults. Keep groups and counts separate from biological ASV labels.

## The change that actually lowers SDM memory

SDM currently admits full exact keys into `Dereplicate::apply_pending_entry_locked`; `outputGroups` performs the prefix pass later. The earlier 99% experiments also operated on an already collected set of exact unique reads. Reducing output there does not avoid the memory peak from collecting those reads.

Compression needs to happen **before admitting a fresh full read to the permanent exact-key map**, or within bounded batches that are compacted and released. Store a compact group index, a fixed search representative, a small retained set of sequence/quality records, and sparse 64-bit sample counts. Release other read objects after recording their counts and any requested disk-backed mapping.

A cap per group alone does not impose a global memory bound: novel reads and errors can keep creating groups. Add an explicit overall memory limit with partition/batch spilling. Grouping each bounded chunk separately may produce extra representatives across chunks, but still avoids retaining the complete unique-read population in RAM.

For scale, 1,500 bp A/C/G/T requires 375 bytes when packed. A full sequence plus one-byte qualities costs approximately 3,000 payload bytes per retained read. A 50-read reservoir is therefore about 150 kB per full group, before IDs, counters, indexes, allocation overhead, and search anchors. A thousand such groups would hold about 150 MB of raw payload. These are storage calculations, not measured SDM RSS.

## Small follow-up experiment

The same 5,000-read Q20 fixture was reused: 50 artificial 1,500 bp templates, 100 independent reads per template, 1% substitutions. Templates differ by 1.35% on average within each of five lineages and 13.34% between lineages. Reads were shuffled in three deterministic arrival orders, with the first matching fixed centroid receiving the read. The positional 20-mer lookup and packed verifier were reused with a configurable mismatch budget.

The table gives medians across the three orders. Retained-read columns are `sum(min(group_size, cap))`; they count how many real reads a capped reservoir would retain. Reservoir selection and Savont execution were not benchmarked. A reservoir does not provide extra searchable centroids in this experiment.

| Raw-read identity threshold | Groups / one representative | Retained with cap 20 | Retained with cap 50 | Reads outside their group's dominant template | Grouping time |
| --- | ---: | ---: | ---: | ---: | ---: |
| 99% | 4,844 | 5,000 | 5,000 | 0.00% | 0.084 s |
| 98.5% | 3,805 | 4,777 | 4,985 | 0.00% | 0.159 s |
| 98% | 1,539 | 2,932 | 4,263 | 0.40% | 0.074 s |
| 97.5% | 338 | 1,413 | 2,934 | 3.78% | 0.026 s |
| 97% | 74 | 1,012 | 2,225 | 37.98% | 0.028 s |
| 96% | 27 | 242 | 438 | 86.42% | 0.019 s |

No groups mixed the five distant lineages in these runs. The dominant-template metric measures mixing of known simulation templates inside a group; it is not an ASV loss rate. A reservoir can retain members from several templates, and Savont may subsequently split them if enough evidence remains.

At 97.5%, one representative per group would give about **14.8-fold compression**, but retaining up to 20 real reads gives about **3.5-fold compression**. At 97%, the corresponding factors are about 67.6 and 4.9, with substantially more mixing of nearby templates. Thus simply counting group representatives overstates the reduction available when preserving independent reads for downstream inference.

At greater depth, an established group with 1,000 reads and a cap of 50 would pass only 50 onward, a 20-fold reduction for that group. Total compression will also depend on new-group creation, outliers, sample stratification, and diversity. The small fixture does not establish scaling or peak memory.

The 99% count differs slightly from the previous 4,855 because this experiment uses randomized arrival order rather than lexicographic order. The 99% default was regression-checked against the existing exhaustive fixtures; every assignment in the new sweep was separately checked to satisfy its chosen mismatch bound.

## Strategies worth considering

### 1. Loose streaming groups with capped read retention — best first implementation

Use the existing positional seed and mismatch machinery at a looser threshold, initially around 97.5–98%. Admit a read to an existing group after a full-length comparison; otherwise create a new group. Keep its search anchor fixed to limit drift and avoid union-find chaining through many individually similar reads.

Keep all reads from small groups. Once a group reaches its cap, retain a mixture of high-quality reads and an unbiased reservoir of the remaining reads. A few slots for supported sequence diversity are an optional extension. Selecting only the first reads risks poor representatives; selecting only the most unusual reads favors sequencing errors. Keep sample provenance and account for sample imbalance when pooling.

Advantages: directly reuses SDM's experimental index/verifier, predictable work, compact counts, and explicit control of output depth. Limitations: order dependence, noisy anchors can cause unnecessary groups, and different true variants can be combined. The current benchmark covers substitutions and corresponding positions; real ONT indels need separate evaluation with an appropriate verifier or positional tolerance.

A bounded quality-sorted input batch can establish cleaner anchors before processing the remaining reads. For example, a representative with 0.2% expected errors and an incoming read with 1% expected errors differ by about 1.2% under the same independent substitution model, rather than the roughly 2% expected between two Q20 reads. This should improve compression at a given threshold when quality estimates are calibrated; the current equal-error simulation did not test that selection policy.

At fixed Q20, both a query and its representative contribute errors. A quality-aware version could adjust its mismatch allowance using the sum of per-base error probabilities, with a hard divergence ceiling. This may reduce unnecessary splitting caused by poor reads, but calibration errors can also make it over-merge. It is a later optimization, not required for the first implementation.

### 2. Multiple small sequence fingerprints or MinHash LSH buckets — strongest alternative index

Use several small signatures that summarize different parts of a read. Look up exact matching signatures, rank the resulting group candidates, then check overall similarity before choosing a group. This avoids enumerating every group sharing any common minimizer, which was the expensive behavior in the previous benchmark.

For example, a signature over 32 selected positions is error-free in approximately `0.99^32 = 72.5%` of Q20 reads. Two independent copies agree at those positions with probability approximately 52.6% under the same random substitution model. Several signatures provide multiple chances to retrieve a group. Spread positions across the amplicon; a short conserved prefix alone is a poor discriminator. Missing sampled differences are biological information loss, not hash collisions.

An alternative is MinHash banding: hash small groups of sketch entries into several tables, rather than hashing the entire sketch or querying the union of all individual minimizers. A stable shared subset can retrieve a candidate even when other sketch entries differ. A limited top-candidate policy is acceptable for this compression objective: missed lookups primarily create extra groups, while a final similarity check controls overly broad assignment.

Savont itself provides a useful implementation reference. Its inspected primary clustering path uses multiple LSH signature tables, ranks candidates by bucket hits, and verifies a sketch-based similarity estimate. Its parameters include 20 tables and three selected values per signature. This is closer to the desired pseudo-dereplication lookup than the earlier exhaustive union of minimizer postings. [Savont clustering](https://github.com/bluenote-1577/savont/blob/d0f38da70961a6058bc5563326eae77aeb3c392f/src/asv_cluster.rs), [signatures](https://github.com/bluenote-1577/savont/blob/d0f38da70961a6058bc5563326eae77aeb3c392f/src/types.rs), [constants](https://github.com/bluenote-1577/savont/blob/d0f38da70961a6058bc5563326eae77aeb3c392f/src/constants.rs)

I would treat LSH as a candidate for reducing index and lookup cost. It does not by itself determine how many independent reads should be kept for Savont.

### 3. Coverage normalization using k-mer counts — attractive for a strict memory budget

Track approximate k-mer abundance, retain reads with insufficiently represented sequence, and discard further high-coverage copies. This is the digital-normalization approach implemented by khmer; its median-based command skips reads whose median k-mer abundance exceeds a chosen cutoff. Bignorm adds quality-aware filtering. [Digital normalization paper](https://arxiv.org/abs/1203.4802), [khmer command](https://khmer.readthedocs.io/en/latest/user/scripts.html#normalize-by-median-py), [Bignorm](https://arxiv.org/abs/1610.03443)

For SDM, compact approximate counters could bound the index memory and avoid full read-to-centroid comparisons. This is especially attractive when the desired output is a coverage-capped candidate set rather than an exact partition of input reads.

My adaptation would consider coverage in several separate read windows, including variable regions, rather than only a whole-read median. A dominant template can make conserved 16S words abundant even when a new template is rare. Conversely, random sequencing errors create novel words; accepting a read whenever any word is new would defeat compression. Quality filtering and requiring repeated support for novel words can help, at the cost of a warm-up or second pass.

This approach has no natural exact read-to-representative map. Abundance estimation should use preserved full-read input or explicitly accept an approximation based on group/count summaries. It also needs collision/occupancy monitoring so a saturated approximate counter does not silently suppress novel sequence.

### 4. Signatures built from repeatedly supported k-mers — promising, more state

Build sketches from words seen repeatedly and supported by adequate base quality. Independent error words tend to disappear, while shared template words remain. Index those sketches using small fingerprint bands. This addresses the reason an exact full-read hash fails at Q20: almost every read has its own error pattern.

A stable, small set of supported markers per group is preferable to accumulating the union of every observed word, which can make memory grow with error diversity. A second pass can establish stable counts before assigning signatures; an online version must manage changing markers and reindexing. isONclust3's use of high-confidence minimizers to update cluster representations is a relevant design reference, although its target is transcript-family clustering rather than ASV compression. [isONclust3 paper](https://academic.oup.com/bioinformatics/article/41/5/btaf207/8118158)

Given the accepted loss of close variants, intentionally excluding weak minority markers could improve compression. It would also make abundance and sequence diversity explicitly approximate.

### 5. Small-group consensus before collapsing — strong compression, greater implementation cost

Accumulate a handful of reads, form a consensus, and use that improved sequence to absorb later noisy copies. A cleaner representative reduces the two-noisy-read problem and can produce far fewer residual groups than a frozen raw anchor.

However, this adds alignment, consensus state, error calibration, and mixed-template handling. Per-base counts do not preserve long-range linkage when several variants share a group. Emitting one consensus to Savont also requires a deliberate downstream interface; repeating it according to abundance creates copies of one error pattern, not independent evidence. This is worth exploring if approximate template consensuses are acceptable as the final compressed object, but is a larger step than capped raw-read reservoirs.

### 6. Supporting measures

Homopolymer-compressed keys can reduce distinctions caused by run-length differences, but do little to eliminate random substitution diversity by themselves. Likewise, quality masking should affect lookup/selection first; destructively masking every low-confidence base in exported reads can remove evidence Savont needs.

Random downsampling is a useful speed baseline, but ignores local abundance and can remove rare groups. Per-group reservoirs or sequence-aware normalization target redundant coverage more selectively.

Disk-backed exact groups or representative-plus-difference encoding can lower RAM without losing reads, but do not themselves reduce Savont's input record count. They are useful complements if abundance remapping or later recovery of discarded variants matters.

## Savont interface and evidence

The inspected Savont revision is `d0f38da70961a6058bc5563326eae77aeb3c392f`. It performs read-based polymorphic-marker detection, consensus generation, and abundance refinement. Its default minimum cluster size is **12 reads**, and the inspected CLI/parser do not expose input representative multiplicities. A FASTQ record with a count in its name therefore should not be assumed to replace a set of independent reads. [CLI](https://github.com/bluenote-1577/savont/blob/d0f38da70961a6058bc5563326eae77aeb3c392f/src/cli.rs), [parser](https://github.com/bluenote-1577/savont/blob/d0f38da70961a6058bc5563326eae77aeb3c392f/src/seq_parse.rs), [pipeline](https://github.com/bluenote-1577/savont/blob/d0f38da70961a6058bc5563326eae77aeb3c392f/src/main.rs)

My proposed interface is consequently **a reduced FASTQ of retained real reads, plus SDM group/sample totals in a sidecar**. Preserve all members of small groups. The cap should leave enough reads for downstream consensus and any intended separation of mixed groups; exceeding 12 at the coarse-group level alone does not guarantee this.

SDM group totals can remain exact even while sequence representation is lossy. They cannot give an exact final ASV abundance split when a discarded set contained multiple variants. Remapping original reads from disk to the recovered ASVs can estimate abundances for recovered sequences without retaining all reads in RAM; it cannot restore an ASV that candidate reduction prevented Savont from discovering. If only the compressed set is retained, abundance allocation across later splits is necessarily approximate.

The existing ONT/LotuS documentation describes the full-read Savont route. A reduced-read route would be a separate explicit option with wrapper support and downstream testing, rather than silently changing the meaning of the existing dereplication output.

## Suggested next experiment

Compare streaming 97.5% and 97% groups with 20, 50, and optionally 100 retained reads, against LSH-based coarse groups with the same retention caps. Include a small quality-preferred subset plus reservoir sampling. Measure peak RSS and actual retained FASTQ counts, then run Savont on the original and reduced inputs. Evaluate recovered consensuses, errors in those consensuses, loss of close/rare variants, and abundance distortion. Include real ONT indels and uneven abundance; the current uniform, substitution-only fixture is insufficient to select production defaults.

## Reproduction

The existing [ONT simulation](../.claude/archive/benchmarks/Test99OntSimulation.cpp) now has a `--compression` experiment mode. It reuses [the candidate benchmark](../.claude/archive/benchmarks/Test99Candidates.cpp), whose default remains 99%:

```sh
make -j4 WITH_HTS=0 STATIC=0 BUILD_DIR=/tmp/sdm-pseudo-build \
  /tmp/sdm-pseudo-build/Test99OntSimulation /tmp/sdm-pseudo-build/Test99Candidates
/tmp/sdm-pseudo-build/Test99Candidates --self-test
/tmp/sdm-pseudo-build/Test99OntSimulation 100 --compression \
  > compression.csv 2> simulation.txt
```

Saved outputs: [threshold sweep](benchmarks/pseudo-derep-compression.csv) and [fixture diagnostics](benchmarks/pseudo-derep-simulation.txt). Timing excludes the stored synthetic fixture and diagnostic bookkeeping. The source changes are experimental tests and documentation only; production dereplication and the concurrent ONT/IO work were not modified.
