# Dereplication and representative-selection audit (sdm 3.60 beta, 2026-09-25)

Second part of the wider audit (first part: [I/O audit](io-audit.md)). Read-only review of the
dereplication path (`Dereplicate` in `containers.cpp`, `DNAunique`/`DNAuniqSet` in `DNA.cpp`,
`CoarseDereplication.cpp`, the merge step in `ReadMerger.cpp`) and of the "better DNA" selector
(`seed_selection::better` in `SeedSelection.h`, `whoIsBetter`, `betterSeedStruct::whoIsBetter2`,
`DNAuniqSet::setBest`). Nothing was changed. Items marked *reproduced* were run against the
`WITH_HTS=0 STATIC=0` build in WSL with `.claude/archive/derep-audit/`.

## Status

Fixed on 2026-09-25 (verified with `TestAuditFixes` cases
`exact_derep_recovery_order_independent` and `rejected_merge_leaves_no_marks`, the updated
`TestCoarseDerep` partition checks, and the pipeline suites):

- Finding 1: exact-key mode defers a non-passing read whose key has no passing dereplicate yet
  (`Dereplicate::apply_pending_entry`). The deferred table holds only sample counts under a
  64-bit key hash (`DeferredCounts`, about 21 bytes per key, capped at 4 million keys); the
  variant-retaining modes keep the deferred reads themselves within a 64 MiB budget, since each
  counted read must become an exact variant. Beyond either cap reads count as unmatched (as in
  3.60) and the report says so. A first version kept whole reads in every mode and added
  74 MiB (`apong_fq`) to 216 MiB (97% coarse with retained qualities) peak memory on Apong; the
  bounded version adds 1 MiB and 77 MiB. At output, deferred keys are attached to the passing
  group of the same key (`outputGroups`), counted like prefix recovery and reported as
  `Exact-key matching: N non-passing reads recovered`; keys without a passing group are
  dropped and never written to `.rest` or the map. A record set bound for the mid-quality
  output is either counted (key already has a passing dereplicate; the write is withheld) or
  written there and not counted (`Dereplicate::publish_uncounted` through the
  `OutputStreamer` callback), so LotuS3's back-mapping of that tier counts each read once.
- Finding 2: `ReadMerger::merge` no longer marks its inputs; `DNA::adoptMergeStats` is called
  only where a consensus is accepted (`Dereplicate::addDNA`, the output re-merge,
  `OutputStreamer::write2Demulti`, `DNAunique::attachMerge`).
- Finding 3: `-derepCoarseClusters 1` selects the global coarse index like `-derepReassign 1`;
  the 64-partition mode remains the memory-saving arrangement for exact-parent output.

Also on 2026-09-25, for findings 4 to 8:

- Finding 4: `seed_selection::better` now uses the overlap-error term (`overlapRatio`, weight
  0.12 in dereplication, 0.10 in seed selection, taken off `wAccErr`/`wMergeLenDev`) and
  compares expected errors per base; the dead `thScore`/`rScore` code is gone.
- Finding 5: `DNAuniqSet::setBest` chooses merge-offset subgroups by abundance; consensus length
  only breaks ties.
- Finding 7: with `-derepCoarseClusters 1` the streaming index assigns a read to the nearest
  anchor within budget (ties to the earliest anchor) instead of the first one found; exact-parent
  output keeps first fit, where the anchor is only a compression reference.
- Finding 8: `ReadMerger` requires 16 bases of overlap again (`MergeMinOverlap` in the options
  file) and an exact match below 20 bases.

Finding 6 (merged and R1 keys in one dereplication) and the items under 9 remain open.

## What holds up

- Sample counting is consistent: `diffDNA::addSample` ignores negative sample ids, so deferred
  prefix failures created with sample `-1` add no phantom counts; `.map` totals, `size=`
  headers and per-sample columns come from the same store.
- `-min_derep_copies` rules are evaluated correctly, including mixed rules such as `20,3:3`:
  `getDerepMapSort` returns every sample count (only the top ones sorted), so the bare total
  uses the full sum and `X:Y` uses the Y-th largest sample (`containers.cpp:2108-2135`).
- Multi-threaded admission is deterministic (ordered blocks), and the prefix pass assigns an
  ambiguous shorter key to the lexicographically first maximal key, independent of scheduling
  (`containers.cpp:1325-1364`).
- Coarse cluster membership is always judged against the immutable first anchor sequence, so
  replacing a cluster's representative cannot chain clusters (`CoarseDereplication.cpp:307`).

## Findings

### 1. Recovery of non-passing reads depends on input order in exact-key mode (reproduced)

In exact-key mode (`derepPrefix=auto` with a fixed clustering length, i.e. the usual LotuS
setting `minSeqLength == TruncateSequenceLength`) a read that fails quality filtering is
counted only if its key already exists: `apply_pending_entry_locked` returns without buffering
(`containers.cpp:383`). Prefix mode buffers such reads in `shard.rejected` and attaches them
at output; exact mode has no equivalent.

Reproduced with one high-quality and one identical low-quality read: order `good,bad` gives
`size=2`, order `bad,good` gives `size=1`. Shuffling an input file therefore changes
dereplicated abundances, and the effect scales with the fraction of failed reads (the Apong
run recovered thousands of low-quality reads). The log's `Low-quality reads recovered` counter
reads 0 in both runs; it is known to be incomplete (see [derep-behavior-parity.md](derep-behavior-parity.md)).

Suggested fix: reuse the existing deferred table for exact mode. Store non-passing reads whose
key is unknown in `shard.rejected` (R1 evidence only, as prefix mode does) and, in
`outputGroups`, attach them by exact key lookup instead of the prefix search. Memory is the
same class as prefix mode already pays.

### 2. A rejected merged consensus leaves the mates marked as merged (reproduced)

`ReadMerger::merge` writes the consensus length and overlap errors into both input reads
(`ReadMerger.cpp:482-486`) before the caller has decided whether the consensus is acceptable.
When `Dereplicate::addDNA` then rejects it through `checkMergedDNA` (`containers.cpp:864-868`)
only the `dna_merged` pointer is dropped; `dna->getMergeLength()` stays positive. The same
happens in `OutputStreamer::write2Demulti`/`saveForWrite_merge` and in
`UClinks::prepareSeedCandidate`. Only `takeOverDNA` (`DNA.cpp:1935`) and `attachMerge`
(`DNA.cpp:1735`) ever correct the value, and a representative created from such a read via
`addNewDNAuniq` keeps the stale length.

Consequences:

- The dereplication report counts the read as merged (`containers.cpp:1708`). Reproduced:
  two pairs whose R2 carries an `N` outside the overlap give `key sources: 0 merged; 2 R1`
  but `Dereplication: ... 1 merged)` and an empty `.merg.fa`.
- `seed_selection::better` treats the read as merged (`dMerge`/`rMerg` from
  `getMergeLength()`, `SeedSelection.h:42-51`) although `dM`/`rM` are null. A representative
  whose merge was rejected then blocks every unmerged candidate through the length gate
  `curL/refL < 0.8` (`:67`) when the rejected product was longer than 1.25 x R1, and a
  candidate whose merge was rejected receives the merge preference bonus (`:105-116`).
- `SeedReadStats` copies the stale value into seed extension (`SeedStats.cpp:9`).

Suggested fix: do not write merge metadata into the inputs inside `ReadMerger::merge`; set it
on the mates only where a consensus is accepted, or clear it (`setMergeLength(-1)`,
`setMergeErrors(-1, 0)`) at each rejection site.

### 3. Multi-threaded coarse-cluster output is effectively exact dereplication (reproduced)

Below 100% identity with more than one worker, reads are routed to 64 independent
partitions by a hash of their full search key (`containers.cpp:892-895`), and the partitions
are never consolidated (`:643-653`). Two reads one substitution apart hash to the same
partition with probability 1/64, so they cannot be grouped. Storage-only mode (the default
exact-parent output) uses a min-hash over 20-mers instead, which keeps most neighbours
together, but explicit coarse output (`-derepCoarseClusters 1`) uses whole-key hashing.

Reproduced with 13 variants (0-3 substitutions from one 200 bp sequence, 3 copies each) at
`-derepIdentity 97 -derepCoarseClusters 1`:

| Threads | Output groups | Sizes |
| --- | ---: | --- |
| 1 | 1 | 39 |
| 4 | 13 | 3 each |
| 4 with `-derepGlobal 1` | 1 | 39 |

`docs/coarse-dereplication.md` only says that counts "can differ". Since `-derepReassign 1`
already forces the global index, `-derepCoarseClusters 1` should do the same (or refuse to run
partitioned), and the report line `64 independent partitions; no final consolidation` should
say that coarse clusters are not formed across partitions.

### 4. The selector ignores merge overlap errors and un-normalized accumulated error

`seed_selection::better` computes normalized overlap-error scores from `getMergeErrors()` and
`getMergeErrorsQual()` (`SeedSelection.h:73-87`) and never uses them; `thScore`/`rScore`
(`:92-97`) are dead as well. The legacy `whoIsBetter_old` had the same dead computation. A
consensus assembled over many overlap mismatches (mis-seeded merge, chimera) therefore competes
on equal terms with a clean one.

`AccErrRatio` compares the *sums* of expected errors (`:98-99`); the legacy scorer divided by
length. For candidates of equal per-base accuracy the longer one is penalized. The effect is
small at similar lengths (<1% of the score for a 10% length difference) but systematic in
coarse mode and in merged-vs-unmerged comparisons.

Suggested: fold the overlap-error term into the product (it already has the normalized form)
or delete the computation; compare `getAccumError()/length()`.

### 5. Merge-offset subgroup choice is dominated by consensus length

`DNAuniqSet::setBest` picks the representative among merge-offset subgroups of one key by
`count x lengthRatio`, where `lengthRatio` is the ratio of merged lengths (`DNA.cpp:2262-2266`).
Different merge offsets for the same R1 mean at least one offset is wrong, and abundance is the
reliable signal; with the linear length factor, 10 reads merged at 450 bp outrank 14 reads
merged at 300 bp. Suggested: choose by count, use length only as a tie-break.

### 6. Merged and R1 keys mix in one dereplication (design note)

With `searchWithMerg` (default) a pair that merges is keyed by the consensus prefix, a pair of
the same fragment that does not merge by its R1 prefix (`containers.cpp:872-887`). Where the
search cap reaches into the overlap and the consensus took R2's base, the two keys differ and
the fragment's counts split across two dereplicates; prefix mode only rejoins them when one key
is a prefix of the other. With a cap inside the non-overlap region there is no effect. Worth a
sentence in `docs/workflows.md` next to `noSearchWithMerge`.

### 7. Coarse anchors are first-fit by arrival (design note)

`CoarseDereplication::search` with `nearest=false` accepts the first anchor within the
mismatch budget in creation order (`CoarseDereplication.cpp:173-189`, `CoarseDereplication.h:43-45`).
A read within budget of two anchors joins the older one; anchors are whichever reads arrived
first, error-containing or not. Documented in `coarse-dereplication.md` ("greedy grouping
depends on input/batch order"); `-derepReassign 1` corrects it at output. Worth recommending
reassignment whenever coarse clusters are the requested output.

### 8. Merge seeding accepts very short overlaps

`findSeed` uses 5-mer seeds and accepts the first hit with overlap >= 12 bp and at most
`overlap/10` mismatches (`ReadMerger.h:77-83`, `ReadMerger.cpp:335-341`; the minimum was 16).
A 12-19 bp overlap with one mismatch is a weak signal: for two random non-overlapping reads
the chance of a spurious 12-15 bp seed-plus-identity match is of the order of 1e-4 per pair,
which over 20 M pairs is thousands of chimeric consensus reads unless the merged-length checks
catch them. Consider restoring 16 bp, or requiring zero mismatches below 16 bp.

### 9. Smaller items

- `betterPreSeed` is unreachable (`exit(1243)`, `DNA.cpp:1428-1473`); `whoIsBetter_old`
  remains selectable via `-DNAseedSelectionClassic 1`.
- `mergeLenDevRatio` carries the largest weight in dereplication (0.38, `SeedSelection.h:163`);
  its running average is per cluster over both compared lengths, so it is symmetric for the
  first comparison. Fine, but the weight deserves a comment since it exceeds quality and error.
- The `Low-quality reads recovered` filtering counter cannot reflect deferred or exact-key
  recovery; either finalize it from the dereplicator's admission statistics or drop it.
- `-min_derep_copies` parses with `atoll`/`atoi` (`containers.cpp:700-712`): `abc` silently
  becomes `0` and accepts everything.
