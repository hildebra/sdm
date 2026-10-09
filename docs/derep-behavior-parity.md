# Ordinary and retained-quality dereplication: behavioral audit

Investigated on 2026-09-15. Default coarse output preserves exact dereplicates, but quality retention is not yet a behavior-preserving switch for an ordinary workflow using merged search. Two retention-specific admission bugs were corrected in this update. HQ seed exports now preserve full available mates and their observed qualities, including hidden tails. Native merged-search support and abundance-only handling of physically empty mates still require changes.

The newer supplied Apong logs, with a **96,962-count** passing difference, are reconciled in the [follow-up investigation](apong-derep-count-debug.md). That update adds finalized admission diagnostics and a 200 bp merged-consensus regression. The observations below describe the earlier comparison.

## What the supplied Apong logs establish

Both logs report 28,083,099 input pairs, identical per-mate filtering totals, and accepted R1 length quantiles of 200/200/200. Ordinary processing reports 24,829,474 merges (88.4%) and 5,855 low-quality recoveries; the coarse-retention log reports zero recoveries and no merge total. Both use the same version label, `3.53 beta`, which does not establish that their executables were byte-identical.

The reported dereplication difference is 22,440,274 minus 22,241,089 = **199,185 counts**. These filter logs alone cannot attribute that difference to particular read pairs or output partitions. The HPC inputs, effective options, dereplication reports, and map/main/merged/rest files were not accessible in this workspace. The local comparisons below use one binary for all four configurations and controlled, reproducible inputs; they are not measurements of the Apong dataset.

Three distinctions matter when interpreting the logs:

- The older `Dereplication: ... counts` total includes only parents passing the minimum-copy rules. The `.map` also includes below-cutoff parents in `.rest`. A group splitting below a cutoff can change the displayed count without losing map abundance. When merged output is enabled, passing counts include both the main file and its `.merg` sibling.
- The filtering log's low-quality-recovery counter is incomplete for batched/deferred dereplication. Filtering statistics inspect the read before its queued dereplication update may have finished; prefix recovery also happens later, at output. This counter excludes alternate-tier reads. Zero does not establish zero recovered abundance, and 5,855 is not a complete reconciliation of the 199,185 difference. This update does not replace that filtering counter with a finalized recovery ledger.
- `min(accepted R1, accepted R2)` is an **upper bound** on accepted pairs, not a measured joint count or lower bound. The old `Minimum pairs with both reads accepted` label was wrong and is now `Maximum pairs with both reads accepted`. The supplied value, 21,208,730, cannot establish how many pairs have only one accepted mate.

## Corrections made

1. **Default exact-parent admission uses R1 eligibility.** Quality retention previously required R2 to pass too, both for creating a parent and for considering a better representative. It now reuses the ordinary R1 rule and paired representative selector. Passing R2 alone cannot establish a new R1 parent. Once admitted, complete repeated pairs contribute once even if a mate fails. The selected quality witness for each exact full-pair variant remains a separate decision and still prefers an observation with both mates passing.
2. **Default exact output uses ordinary prefix consolidation and deferred recovery.** `derepPrefix=auto/1` now also applies when coarse or retained storage supplies exact output parents. Variable-length exact keys are collected during processing; the existing deterministic prefix pass attaches shorter compatible keys and failed/alternate-tier evidence at output, before minimum-copy filtering. Retained variants and counts follow that evidence through the existing `diffDNA` merge. Failed reads cannot supply the parent representative or bridge conflicting longer anchors. An explicitly positive search length in automatic mode still selects capped, equal-length keys; forced prefix mode overrides it as in ordinary dereplication.
3. **HQ seed FASTQ preserves full mates and observed qualities.** Retained variants now use the full backing sequences after physical preprocessing cuts, including tails hidden by truncation or quality trimming. Their logical lengths remain separate in the scoring metadata and are respected by final reassignment. Tail differences remain distinct variants. The two quality vectors still belong to one observed fragment. Ordinary representative HQ files already use full backing mates. This contract covers `.1.hq.fq` and `.2.hq.fq`; the optional main `-derep_format fq` retains its existing processed/averaged-quality behavior.
4. **Finalized abundance is reported explicitly:** `Dereplication abundance: X total in map; Y passing; Z in rest`. Retention additionally reports pairs whose R2 became empty and were skipped before admission. That exclusion number counts attempted fragments, not the counterfactual number that ordinary mode would admit.

Explicit `derepCoarseClusters=1` or `derepReassign=1` still requests coarse-parent output and its existing admission policy. Those modes intentionally have different parent counts and cutoff decisions. The new prefix pass consolidates exact compatible prefixes; it does not reconcile approximate search groups across partitions.

## Controlled comparison

`tests/test_derep_behavior_parity.pl` generates **1,000 paired fragments**, split equally between two samples. Most mates are 230 bp from 300 bp inserts. The fixture includes R1-only and R2-only passing observations, both-failed repeats, failures before passing anchors, different merged tails sharing R1, long/short prefixes, singletons, and 40 R2 mates removed completely by a physical cut. Nearby R1 variants exercise actual 97% grouping.

Every profile uses the same binary, inputs, filters, search settings, and minimum-copy cutoff (40) across these modes:

| Mode | Preprocessing merge flags: derep/filter/demulti | Search | Retained qualities |
| --- | --- | --- | --- |
| Ordinary merged | 1 / 1 / 1 | Merged where available | 0 |
| Ordinary R1 | 0 / 0 / 0 | Explicit `noSearchWithMerge=1` | 0 |
| Retained 100% | 0 / 0 / 0 | R1, currently forced by retention | 1 |
| Retained 97% | 0 / 0 / 0 | R1, currently forced by retention | 1 |

All three original merge-flag values were not supplied, so `1/1/1` is an explicit controlled configuration rather than a claim about the original HPC command. The runner stores every command and the binary SHA-256. The tested portable executable had SHA-256 `ab3d1b02072a950de148754e536cc31c11cc0af2361d6b49039272287fea7ced`.

Results after the corrections were identical at **1, 4, and 12 processing threads**. Each cell is **map total / passing / rest**, in fragments:

| Profile | Ordinary merged | Ordinary R1 | Retained 100% | Retained 97% |
| --- | ---: | ---: | ---: | ---: |
| Search cap 200; no trimming | 900 / 780 / 120 | 900 / 780 / 120 | 860 / 780 / 80 | 860 / 780 / 80 |
| Uncapped, automatic prefix matching | 960 / 780 / 180 | 960 / 900 / 60 | 920 / 900 / 20 | 920 / 900 / 20 |
| Search/truncation cap 200; quality trimming | 840 / 720 / 120 | 760 / 640 / 120 | 720 / 640 / 80 | 720 / 640 / 80 |

The two samples have equal totals in this fixture. The audit additionally checks each sample independently, each parent header against its map abundance, and that main + merged + rest partitions the map exactly. Retained 100% and 97% agree on parent counts, cutoff decisions, representative sequences, and linked variant sequences/counts including selected qualities. Equal-scoring duplicate observations can supply different original parent names; comparisons normalize those names after validating each output set's parent linkage. Thread count changes internal search partitioning but did not change these results.

Before the corrections, retained 100% and 97% both produced:

| Profile | Map / passing / rest before | Map / passing / rest after |
| --- | ---: | ---: |
| Capped | 680 / 540 / 140 | 860 / 780 / 80 |
| Prefix | 680 / 420 / 260 | 920 / 900 / 20 |
| Trimmed | 540 / 400 / 140 | 720 / 640 / 80 |

The causes can be isolated in the fixture:

- **180 counts** were lost through the extra R2 eligibility gate: families with only passing R1, and R1-passing observations preceding the first both-passing pair. Restoring these also moves some parents over the cutoff.
- In prefix mode, another **60 counts** from failed observations preceding their passing anchor were lost because coarse storage bypassed deferred recovery. Separately, consolidating the passing long/short prefixes restores **120 counts** to passing output from `.rest` without changing their total abundance.
- The remaining **40-count map difference** against ordinary R1 is exactly the physically empty R2 family. These groups happen to be below the cutoff; this is why passing totals agree in this fixture. They could affect passing abundance on another dataset.
- Uncapped ordinary merged search distinguishes two merged sequences that share R1. Their 30-copy parents fall below the cutoff; R1-only search combines them into a 60-copy parent. That accounts for the **120 passing/rest difference** between the ordinary modes without any map loss.
- With quality trimming, merged processing admits **80 more counts** than ordinary R1 processing in this fixture. The merger can use full backing sequences to reconstruct a valid search sequence from pairs whose failed R1 logical sequence is too short to match the capped R1 key. Disabling merges therefore changes admission even though the per-mate filtering summaries agree.

An additional alternate-quality profile exercises mid-tier recovery. `--complete-mates` replaces the 40 empty mates with complete mates while retaining 1,000 fragments; with that fixture, `--require-parity` requires ordinary R1 and retained-mode map/passing/rest equality. Neither assertion requires ordinary merged output to equal R1-only output.

## Full HQ seed contract

“Full” means the available biological read after configured **physical** barcode, primer, adapter and fixed-base cuts. HQ export restores bases hidden by search limits, `TruncateSequenceLength`, quality-window trimming and logical homopolymer trimming. It does not restore physically excised technical sequences, input read names, or discarded duplicate quality observations. Keeping those cuts is intentional: both ordinary and reconstructed seed readers use `prepareSeedCandidate` on already-processed mates and do not repeat fixed preprocessing cuts. With `seedSubclusters=1`, candidates also bypass barcode/primer filtering. Set `keepBarcodeSeq=0` during preprocessing when barcode bases should be absent from seeds; do not rely on seed extension to remove retained barcodes.

Before this fix, the trimmed comparison exported retained R1 at 200 bp and some retained R2 at 41 bp, while ordinary HQ mates were 230 bp. Both workflows now preserve 230 bp source mates at 230 bp with their observed quality strings; naturally shorter 200 bp variants remain 200 bp. Search and main representative output still use the configured logical lengths. A mate whose logical length reaches zero can still be retained if its backing sequence exists; physically empty R2 remains an explicit exclusion.

The seed-extension regression additionally removes real barcodes/primers plus fixed 5′ cuts, caps preprocessing reads at 100 bp and trims low-quality tails. It checks that ordinary and retained HQ workflows supply full mates to seed extension, including with `merge_pairs_seed=1`, without applying technical cuts a second time. Regenerate older retained HQ exports to obtain bases lost by the previous implementation; changing only the seed command cannot recover them.

## What remains unsupported, and the required change

Current `derepStoreQuals=1` forces R1 search and rejects native `merge_pairs_derep`, `merge_pairs_filter`, and `merge_pairs_demulti`. `merge_pairs_seed=1` remains available later, but it cannot undo different upstream grouping/admission or recover bases missing from older exports.

The merger already builds a separate consensus object while the original mate buffers remain available. The obstruction is in the dereplication representation and IO contracts:

1. **Separate the search anchor from the variant references.** `CoarseDereplication::Group::reference` currently serves both candidate search and R1 delta reconstruction, and admission assumes the search key is a prefix of the retained R1 sequence. Preserve the ordinary merged-or-R1 search sequence independently of full R1/R2 references. Keep ordinary merge-offset subgroup selection and the merged representative/output attached to the selected pair. Apply `noSearchWithMerge`, search caps, and ordinary merged filtering to that search/output view.
2. **Separate merge participation metadata from merged-sequence identity.** `ReadMerger` sets merge-length metadata on the original mates as well as creating the consensus. Current quality validation treats that metadata as a reason to reject the observation. Merely removing the command-line check would still fail and would not fix the shared-reference assumption.
3. **Define non-seed abundance explicitly.** Ordinary dereplication can count a passing R1 whose R2 was physically removed. Paired variant storage/seed input requires two nonempty mates, and `seedSubclusters` requires variant abundances to sum to the parent map total. Preserve such fragments as count-only evidence, with an explicit IO representation and a completeness check of `seed-variant counts + non-seed counts = parent count`; alternatively retain an appropriate pre-trim fragment with explicit eligibility. Do not assign their counts to a different full-pair variant or remove the seed reader's count checks.

These are required changes for full ordinary merged-search/output parity with retained full pairs; they are not implemented by the admission fixes above. LotuS should record the preprocessing merge policy and retention policy separately in its cache metadata and should not present the current retention switch as preserving the original merged behavior.

## Validation

The portable `make WITH_HTS=0 test` suite passed with the full-HQ changes, including the four-mode mixed-quality comparison, paired/single-end retained variants and seed extension. The complete-mate comparison also passed with both main FASTA and main FASTQ, requiring ordinary R1/retained abundance parity while checking full observed HQ bases and qualities. The coarse/diffDNA C++ suite passed AddressSanitizer, UndefinedBehaviorSanitizer and leak detection, including a regression ensuring hidden HQ tails do not change logical reassignment ties.

## Reproduce and audit real outputs

```sh
make WITH_HTS=0 test
perl tests/test_derep_behavior_parity.pl ./sdm --output /tmp/sdm-parity-new
perl tests/test_derep_behavior_parity.pl ./sdm --complete-mates --require-parity
perl tests/test_derep_behavior_parity.pl ./sdm --main-format fq --complete-mates --require-parity
# Inspect an older binary without requiring the corrected admission behavior:
perl tests/test_derep_behavior_parity.pl /path/to/old/sdm --observe
```

For the Apong reconciliation, run the four configurations above with the same executable and inputs, preserving the original merge values for the first run. Use fresh output directories, the same options/search limits, `derepPerSR=0`, and the same minimum-copy rule (`8:1,4:2,3:3` in the supplied command). Record the actual commands and executable/configuration hashes. Do not reuse an output stem that might have stale `.merg` files.

For each completed output set:

```sh
perl tests/audit_derep_counts.pl /path/to/run/derep.fas > /path/to/run/count-audit.json
```

The audit reads `derep.map`, `derep.fas`, optional `derep.merg.fas`, and `derep.fas.rest`; it reports per-sample totals in each partition and rejects missing/duplicate parents or abundance mismatches. It accepts ordinary SDM FASTA or four-line FASTQ output. The standalone audit covers parent abundance, not the completeness of retained HQ variants; the regression runner separately verifies paired variant sums and qualities. Apply it to each complete per-run output set rather than a cumulative `derepPerSR=1` map with only one run's main file.
