# Apong dereplication count reconciliation

Investigated on 2026-09-15 using the newer supplied filtering and dereplication logs. These supersede the earlier 199,185-count comparison in the [behavioral audit](derep-behavior-parity.md). The attachment contains filtering summaries and 166 sample assignment rows per run, but no input reads or parent/sample maps. This report distinguishes arithmetic established by those logs from mechanisms reproduced on local controlled inputs.

## What changed between the runs

| Quantity | Ordinary | Retained coarse | Coarse minus ordinary |
| --- | ---: | ---: | ---: |
| Input pairs | 28,083,099 | 28,083,099 | 0 |
| High-quality R1 | 22,142,083 | 22,142,083 | 0 |
| Counts in the complete map | 23,649,182 | 23,607,760 | -41,422 |
| Counts in passing main + merged output | 22,440,309 | 22,343,347 | -96,962 |
| Counts in `.rest` | 1,208,873 | 1,264,413 | +55,540 |
| All exact output parents | 1,151,296 | 1,208,518 | +57,222 |
| Passing output parents | 186,002 | 192,047 | +6,045 |
| Parents in `.rest` | 965,294 | 1,016,471 | +51,177 |
| Search observations | 28,083,099 | 28,083,094 | -5 |
| Sum of filtering-log sample assignments | 22,148,052 | 22,142,083 | -5,969 |

Each map total equals passing plus rest. The passing deficit decomposes exactly as:

```text
96,962 = 41,422 fewer map counts + 55,540 more counts in rest
```

The second term is a **net difference in output partitions**, not proof that exactly those same fragments moved between corresponding parents. Parent identities differ between these workflows.

The coarse log also reconciles independently:

```text
22,399,242 exact-key hits + 1,208,518 newly admitted exact keys
    = 23,607,760 map counts

28,083,094 queries - 23,607,760 admitted observations
    = 4,475,334 unmatched, non-passing observations

5,683,852 candidate searches - 4,475,334 rejected new keys
    = 1,208,518 newly admitted exact keys
```

Thus the reported exact-key hits and new parents account for every coarse map count. The 770,104 internal approximate groups do not replace the 1,208,518 exact output parents or own another set of abundance counts. These numbers provide no evidence of counts lost while combining internal coarse groups. The source code likewise increments abundance once in the exact parent for each admitted observation.

## Which admission differs

Default exact output admits any eligible high-quality R1, independently of R2 quality. In these logs R1 passes a 200 bp minimum, barcodes are identified, and physically empty R2 exclusions are zero. Subtracting high-quality R1 from map abundance therefore gives admitted **mid-quality plus failed R1**:

```text
Ordinary: 23,649,182 - 22,142,083 = 1,507,099
Coarse:   23,607,760 - 22,142,083 = 1,465,677
Difference:                              41,422
```

The remaining map deficit is in non-high-quality admission/recovery, rather than the earlier R2 eligibility bug. Even if every one of the 941,846 mid-quality R1 observations were admitted, coarse mode still counted at least 523,831 fully rejected R1 observations. Its displayed zero low-quality recoveries cannot describe actual recovered abundance.

The filtering counter checks `isDereplicated()` before queued updates necessarily finish and excludes mid-quality reads. Paired retained failed repeats additionally update their variant counts without running the representative-update path that sets this flag. Consequently the ordinary filtering log's 5,969 and coarse log's zero are not comparable totals of recovered evidence. Summing the supplied sample rows confirms that they are high-quality R1 plus that partial counter; they are not `.map` sample totals. Actual per-sample map differences cannot be reconstructed from these rows.

Only five fewer observations reach coarse search. The code skips empty coarse search keys, while ordinary fixed-key submission has no equivalent guard and may also have a usable merged sequence. This is a plausible explanation for those five, but the old logs do not identify them. Five direct observations cannot explain a 41,422-count gap. New logs explicitly count empty search keys.

## Why merging matters at a 200 bp search cap

The ordinary log reports 24,829,474 preprocessing merges and 170,035 passing merged representatives. The retained run disables those merges and produces R1 representatives. Ordinary search uses a filtered merged consensus when available by default; retained-quality search currently forces R1. The old logs do not explicitly record `noSearchWithMerge`, so the effective ordinary search policy cannot be independently verified from these summaries alone.

A cap of 200 bp does **not** make merged and R1 keys equivalent:

1. Within the overlapping part of the first 200 bases, a higher-quality R2 base can correct an R1 sequencing error. Merged search then combines observations which R1-only search keeps as separate exact keys. Small R1 groups may fail `8:1,4:2,3:3` although their combined consensus group passes. This can increase both exact parent count and rest abundance after disabling merges.
2. Quality trimming can leave a failed R1 shorter than the search cap. Merging uses full available mate buffers and can produce a valid filtered search sequence that matches an existing parent. An R1-only short key cannot equal that capped parent in fixed-key mode. Identical per-mate filtering summaries therefore do not establish equal dereplication admission.
3. Fixed-key failed-read recovery is arrival-order dependent in both ordinary and retained storage: an unmatched failed key is discarded before its first passing anchor exists, but is counted after that anchor. Worker scheduling can change this ordering. Prefix mode separately buffers non-passing evidence for final recovery; this is not selected automatically when a positive search cap selects fixed keys.

The supplied aggregate logs cannot assign each of the 41,422 missing map observations, or each differing rest parent, to consensus correction, trimming, or scheduling. The mechanisms above explain how the difference can arise and have controlled reproductions. It would be too strong to attribute the complete Apong difference to merging alone or to claim that these logs establish byte-identical inputs, options, and executables.

## Controlled reproduction

The existing `tests/test_derep_behavior_parity.pl` now has `--consensus-errors`. It generates **1,000 paired fragments** in two samples. Two 30-copy families per sample contain different low-confidence R1 errors within the overlap and the first 200 bases; R2 contains the correct high-confidence bases. The R1 reads still pass filtering. With a minimum-copy cutoff of 40, merged consensus combines their capped keys into a passing 60-copy group, whereas the two R1-only groups remain in rest. The fixture also includes one-passing-mate pairs, failed and mid-quality recovery, failures preceding anchors, and varying lengths.

With `--complete-mates`, no physical-empty-mate exclusion is present. Results at **1, 4, and 12 threads**, using one binary and the same input/filter/search/cutoff settings for each comparison:

| Profile | Ordinary merged: map / passing / rest | Ordinary R1-only | Retained 100% | Retained 97% |
| --- | ---: | ---: | ---: | ---: |
| Search cap 200 | 900 / 780 / 120 | 900 / 660 / 240 | 900 / 660 / 240 | 900 / 660 / 240 |
| Search and truncation cap 200, quality trimming | 840 / 720 / 120 | 760 / 520 / 240 | 760 / 520 / 240 | 760 / 520 / 240 |
| Uncapped, automatic prefix matching | 960 / 780 / 180 | 960 / 780 / 180 | 960 / 780 / 180 | 960 / 780 / 180 |
| Search cap 200, alternate quality tier enabled | 900 / 780 / 120 | 900 / 660 / 240 | 900 / 660 / 240 | 900 / 660 / 240 |

The capped case isolates a 120-count passing/rest difference with **no map loss**. The trimmed case additionally isolates **80 fewer map counts** when merging is disabled. Prefix-mode aggregate agreement in this fixture does not imply identical merged/R1 sequences. The existing failure-before-anchor family contributes 60 counts in fixed mode versus 120 with deferred prefix recovery.

These tests verify per-sample maps, exact representative sequences and minimum-copy decisions between ordinary R1 and retained output. Retained 100% and 97% also agree on the complete HQ pair sequences, selected observed qualities, and variant counts. Every HQ sequence/quality pair is checked against the original full available mate after configured physical cuts. Optional main FASTQ retains its existing processed-output contract. The full portable suite runs the consensus-error fixture with main FASTQ and the original mixed-quality fixture with main FASTA, including seed reconstruction tests.

## Changes made during this investigation

**Finalized admission diagnostics.** Both serial/sharded and work-pool processing now use the same pending-entry helper to report:

- Search policy: merged when available or R1 only, exact-prefix/fixed/coarse output policy, and search cap.
- Search sources: merged-key observations, R1-key observations, and empty keys skipped.
- Admission before prefix recovery: high-quality, mid-quality, failed, unmatched, and deferred observations.

For default exact output, `high + mid + failed + recovered-prefix = total map abundance`. Also, `high + mid + failed + unmatched + deferred = submitted search observations`. Prefix recovery remains a separate finalized line. These counters are updated under the existing shard ownership, avoiding new per-read atomic contention. The legacy filtering recovery counter remains for compatibility and should not be used for abundance reconciliation.

**Skip candidate searches that cannot admit a read.** In default exact-output mode, a non-passing observation with no existing exact key cannot establish a parent, even if it matches an approximate anchor. Previously SDM still ran the approximate search before rejecting it. It now checks exact reuse and rejects that miss immediately. Existing exact failed repeats still count, and a later passing observation of a previously rejected key still searches and creates its parent normally. Prefix recovery and explicit coarse-parent output retain their existing policies.

In the supplied coarse run, **4,475,334 of 5,683,852 candidate-search calls (78.7%)** had this unnecessary admission outcome. The shortcut removes that class of work; this is not a prediction of 78.7% lower runtime or comparisons. Candidate-search cost varies, and filtering, exact hits, retained variants and output still require work. Logs report the shortcut separately so it is not misreported as a candidate search.

Before/after comparisons of all 48 consensus-fixture runs preserved map totals, each sample and family count, passing/rest partitions, filtering results, HQ lengths, and retained HQ sequence/quality/count hashes. At one thread, capped coarse candidate searches fell from 136 to 36 (105 to 65 comparisons); the trimmed case fell from 276 to 36 searches (785 to 65 comparisons). These small controlled measurements demonstrate the skipped work, not an HPC runtime forecast.

The portable `make WITH_HTS=0 test` suite passed. The additional native test checks failed misses, later passing anchors, preserved counts and statistics reset. Integration tests require the admission ledger to reconcile with the map and deferred recovery and check the shortcut count against unmatched observations.

## Remaining contract difference

These changes diagnose and avoid wasted work; they do not enable native merged search with retained unmerged variants. That requires separating the merged search/output reference from the full R1/R2 variant references and their observed qualities, as described in the [required implementation changes](derep-behavior-parity.md#what-remains-unsupported-and-the-required-change). Running R1-only retention is still a change in upstream dereplication behavior from ordinary merged search. Seed merging later cannot reverse that change in admission or minimum-copy filtering.
