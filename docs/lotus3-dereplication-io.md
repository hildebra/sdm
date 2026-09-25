# LotuS3 integration: standard dereplication IO

This contract supersedes the earlier instruction to enable full-variant retention automatically for coarse LotuS runs. Coarse storage preserves ordinary external dereplicates and seed candidates; full-variant retention is an explicit optional policy. Apply these instructions in the LotuS3 repository; this SDM change does not edit LotuS3 itself. The SDM version label alone may not distinguish earlier development builds. Its `-help_flags` output describes retained qualities as `Standard HQ output: .1.hq.fq; paired additionally .2.hq.fq`.

## Files to use

Given `-o_dereplicate /scratch/run/derep.fas`, use the same paths for ordinary and coarse dereplication:

| Path | Pipeline role |
| --- | --- |
| `/scratch/run/derep.fas` | Passing representatives for the existing clustering/mapping stage. Default coarse mode preserves exact R1 dereplicates. |
| `/scratch/run/derep.merg.fas` | Additional passing merged representatives when configured dereplication merging is enabled; include these in passing-abundance audits. |
| `/scratch/run/derep.fas.rest` | Below-cutoff parents for the existing additional-count/recovery stage. |
| `/scratch/run/derep.map` | Parent sample abundances, including below-cutoff parents. |
| `/scratch/run/derep.1.hq.fq` | R1 seed candidates, also the filename for single-end input. |
| `/scratch/run/derep.2.hq.fq` | Synchronized R2 seed candidates for paired input. |

SDM derives the stem by removing the last filename extension. Keep the main FASTA output for clustering; if a consumer explicitly needs main FASTQ, supply `-derep_format fq` and an appropriate `-o_dereplicate` filename. A `.fq` extension alone does not select FASTQ. The HQ files are always FASTQ.

By default, coarse groups are internal candidate indexes and optional compression references. The main FASTA/FASTQ, merged/rest files, `.map` and representative HQ pairs preserve ordinary exact search-key dereplicates, search/truncation caps, prefix recovery and copy cutoffs. Storage-only mode preserves the configured merged/R1 search policy and ordinary merge-offset subgroup selection. With explicit full-variant quality retention, the search source is currently R1 and R2 differences do not split that parent; the existing paired better-seed selection can replace both mates. Select `-derepCoarseClusters 1` only when the output should describe coarse clusters instead. `-derepReassign 1` also explicitly selects coarse-cluster output and performs final reassignment. `-derepGlobal 1` controls search scope, not output granularity.

With `-derepStoreQuals 0`, HQ files contain one selected representative per parent, as before. With `-derepStoreQuals 1`, they contain all retained exact variants and their selected quality vectors. The paired records represent full available `(R1,R2)` variants after physical preprocessing cuts: identical R1 with different R2 can share a parent while remaining separate seed candidates. No additional parent-only HQ copy or `.subclusters.fq`/`.subclusters.1.fq`/`.subclusters.2.fq` is written. There is no file-level split between cluster reads and subreads: the representative's own 100% variant, other retained variants, and singleton clusters are all written together. Each exact variant has its own count and appears once; never append a second cluster-count record for its representative.

## Preprocessing command changes

1. Use the updated SDM executable for both preprocessing and seed extension.
2. Keep `-o_dereplicate`, `-derepPerSR`, copy thresholds, sample delimiter, output quality offset, and normal filtered-output paths. Existing `-suppressOutput` settings do not suppress dereplication files.
3. For storage-only coarse processing, change `derepIdentity` to the intended identity and preserve the configured quality-retention option (default `derepStoreQuals=0`) and preprocessing merge settings. Preserve `noSearchWithMerge` as configured; its presence selects R1 search even if given value `0`, so omit the flag for default merged search. Main FASTQ retains ordinary accumulated qualities, and HQ pairs retain the selected observed qualities without requiring full-variant retention. Only an explicit `derepStoreQuals=1` request changes HQ granularity to full variants and currently requires disabling native preprocessing merges.
4. Keep `-derepStoreDiffs 0 -derepSubclusterFasta 0 -derepReassign 0 -derepCoarseClusters 0` unless those optional outputs or reassignment are required. Binary differences are not required for retained-quality output or seed extension. Default partitioned coarse processing adds no final cluster-merging pass.
5. Replace any coarse-specific FASTQ path selection with the normal `.1.hq.fq` and `.2.hq.fq` paths. Update expected-output checks, file-size checks, cleanup/archive lists, and resume metadata together. Do not wait for `.subclusters` FASTQ files.

A storage-only paired example with preprocessing merging enabled is:

```sh
sdm -i_path INPUT_DIRECTORY -map primary/in.map -options sdm_miSeq.txt \
    -o_dereplicate /scratch/run/derep.fas -derep_format fa \
    -o_fna /scratch/run/demulti.1.fna,/scratch/run/demulti.2.fna \
    -paired 2 -derepIdentity 97 -derepStoreQuals 0 \
    -derepStoreDiffs 0 -derepSubclusterFasta 0 -derepReassign 0 -derepCoarseClusters 0 \
    -merge_pairs_derep 1 -merge_pairs_filter 1 -merge_pairs_demulti 1 \
    -min_derep_copies 8:1,4:2,3:3 -dere_size_fmt 0 \
    -sample_sep ___ -derepPerSR 0 -o_qual_offset 33 \
    -suppressOutput 1 -threads 12
```

The three merge values above are an example: preserve the actual configured values rather than forcing them to 1. Retain LotuS3's existing demultiplexed/filtered output and log arguments as needed; the example only shows the dereplication integration. Search/truncation settings continue to come from the same options and command line.

## Behavioral limits and count audit

Storage-only mode preserves ordinary merged/unmerged search and output, including merge-offset subgroup selection, main FASTQ accumulation and prefix recovery. The following limits apply specifically to explicit full-variant retention. Default retained-variant output uses ordinary R1 admission and ordinary prefix consolidation/recovery. It still cannot combine native merged search/output with paired quality retention. Both ordinary and retained HQ files now preserve full available mates and observed qualities independently of `TruncateSequenceLength` and logical quality trimming. Barcode/primer/adapter and fixed-base cuts remain applied; use `keepBarcodeSeq=0` during preprocessing when seeds must not contain barcode sequence. Reconstructed seed input does not repeat these cuts. Older retained HQ exports may lack tails; regenerate them before seed extension. Pairs with a physically empty R2 remain excluded and are reported. See the [four-mode comparison and required storage/IO changes](derep-behavior-parity.md) before treating retention as equivalent to the original pipeline.

Audit parent abundance as `map = main + merged + rest`, counting each paired fragment once. The new finalized dereplication abundance line reports this partition; the filtering log's low-quality-recovery counter is not a complete count ledger. Run `python3 tests/audit_derep_counts.py /scratch/run/derep.fas` on each fresh complete output set (`derepPerSR=0`).

Record merge flags, `noSearchWithMerge`, search/truncation limits, retention/output policy, and the executable hash in preprocessing cache metadata. Changing these requires regenerating preprocessing and downstream assignments together.

## Seed-extension command changes

Select the reader from the **effective quality-retention setting used to create the files**, not from `derepIdentity` or the filenames:

```text
if preprocessing used derepStoreQuals=1:
    seedSubclusters = 1
else:
    seedSubclusters = 0

R1 = stem + ".1.hq.fq"
R2 = stem + ".2.hq.fq"  # paired input only
```

This rule also applies at `derepIdentity=100` when quality retention is explicitly enabled. Pass the generated parent `.map`, the existing parent-to-OTU assignments, and the original fallback OTU FASTA. For paired retained variants:

```sh
sdm -i_fastq /scratch/run/derep.1.hq.fq,/scratch/run/derep.2.hq.fq \
    -paired 2 -seedSubclusters 1 \
    -derep_map /scratch/run/derep.map \
    -optimalRead2Cluster parent_assignments.uc \
    -OTU_fallback clustered_otus.fna -options sdm_miSeq.txt \
    -sample_sep ___ -i_qual_offset 33 -o_qual_offset 33 \
    -merge_pairs_seed 1 -o_fna OTU.fna -otu_matrix OTU.txt
```

Keep the pipeline's existing assignment format, reference/fallback paths, identity/coverage settings, additional-count arguments, logs, and final-output choices. `-derep_map` is the **input** abundance map for seed extension; preprocessing derives its output map from `-o_dereplicate`.

For single-end retained variants, supply only `derep.1.hq.fq` with `-paired 1 -seedSubclusters 1`, and omit pair merging. For ordinary representative HQ input, use `-seedSubclusters 0` or omit that flag.

Variant headers have the form `parent.subN;size=M;`. The native variant reader resolves the final `.subN` suffix to the parent, evaluates every variant as a seed, and applies the parent's sample vector once. Keep paired files synchronized and each parent's `.sub1`, `.sub2`, ... order intact. Do not strip suffixes, assign every variant the full parent count, independently sort mates, or add R1 and R2 abundances together. Do not filter variant records by the parent minimum-copy threshold: the map and HQ files include below-cutoff parents needed by recovery. Existing seed extension avoids repeating preprocessing cuts, truncation, or quality trimming.

## Mid-quality reads: counted by SDM or published, never both

Since 2026-09-25 a record set bound for the mid-quality tier (`-o_fna2`/`-o_fastq2`) is
handled in one of two ways, decided when the dereplicator admits it:

- Its key already has a passing dereplicate: the set is counted into that dereplicate's
  abundance (`derep.map`, `;size=N;`) and **withheld** from the mid-quality output.
- Otherwise (no passing copy of that key seen so far, or a coarse cluster rejects it): the set
  is **written** to the mid-quality output and not counted. The dereplication report shows
  these as `N published to the mid-quality output instead of counted`.

Failed reads (no output tier) are deferred instead: attached at output to the passing group of
their key when one appears (`Exact-key matching` / `Prefix matching: N non-passing reads
recovered`), otherwise dropped. The deferred tables have fixed limits (4 million keys, or
64 MiB of retained reads with `-derepStoreQuals`/`-derepStoreDiffs`); reads beyond them are
not counted, as in 3.60, and the report line `Deferred table full` gives their number. Nothing that SDM counted appears in the mid-quality files or in
`derep.fas.rest`.

**LotuS3 needs no change.** Its back-mapping of the mid-quality files (`demulti.add.fna` and
siblings, `-ucAdditionalCounts`) and of `derep.fas.rest` (`-ucAdditionalCounts1`) now counts
each read once, because those files contain only reads SDM did not count. The mid-quality
files are smaller than before, so the mapping step is faster. Reads that would only match an
OTU seed by alignment are still caught by that back-mapping. The SDM `-help_flags` output
describes this contract in the sentence starting `Non-passing reads matching a passing
dereplicate are counted in its abundance`; older builds published every mid-quality read and
counted a copy-order-dependent subset of them as well.

## Sequencing runs and resuming older jobs

With `-derepPerSR 1`, the main output receives a sanitized run label, such as `derep.runA.fas`. The normal `derep.map`, `derep.fas.rest`, `derep.1.hq.fq`, and paired `derep.2.hq.fq` accumulate across completed runs. Do not look for run-suffixed HQ files or concatenate the cumulative HQ files once per run. Optional binary differences and exact FASTA exports retain their run suffixes.

This output correction also changes coarse-parent FASTA/map rows into exact-R1 rows by default. Regenerate preprocessing and downstream clustering/assignments together when resuming outputs from earlier coarse builds.

Record the effective output granularity (`exact R1` or `coarse clusters`), quality-retention setting and the HQ record layout (`representatives` or `exact variants`) in the preprocessing completion/cache metadata. Earlier development builds used the same HQ filenames for parent-only records and separate `.subclusters` FASTQ files for variants. Regenerate preprocessing and dependent outputs when migrating those cached jobs, or explicitly keep the old matching reader/path contract for the entire old job. Never infer that an old HQ file contains variants merely because its preprocessing command included `derepStoreQuals=1`. Generated `.diff` files cannot restore retained qualities.

## Integration checks

- For ordinary exact mode without retention, standard filenames and representative seed behavior remain unchanged.
- Compare ordinary 100% (`derepStoreQuals=0`) directly with storage-only 97% (`derepStoreQuals=0`), preserving all merge/search options: main FASTA, `.rest`, `.map`, and HQ records must agree, ignoring record order when seed winners are unambiguous. The all-passing 1,000-pair regression covers this at 1, 4, and 12 workers. The mixed-quality behavioral comparison additionally reconciles failed/alternate-tier recovery and passing/rest totals; it records the remaining empty-mate difference from ordinary R1 processing.
- At 97% with quality retention, every exact full pair appears once in the two HQ files with the same ID/count and the selected observed quality vectors; parent FASTA/map outputs remain available.
- For each parent, HQ variant abundances sum to the parent map abundance. One paired fragment contributes one count, not two. Below-cutoff parents remain available.
- A parent with multiple variants can select a better later variant during seed extension without adding its sample counts again. R2-only variants remain linked to their corresponding R1 qualities.
- Validate single-end input, paired input, `derepStoreQuals=1` at 100%, reassignment if enabled, and cumulative sequencing-run output.
- A fresh run produces no `.subclusters` FASTQ files. Optional `.diff` and `.subclusters*.fna` files only appear when requested.

The storage-only acceptance test `tests/test_storage_only_parity.py` compares IDs, sequences, main/HQ quality strings, map vectors, mate relationships, all output partitions and normal seed commands/statistics/final seeds/matrices at 1, 4 and 12 workers. It includes 1,000 mixed-quality pairs, R2-only variants, search caps, prefix recovery and sequencing runs. See the [tie/order assessment and measured memory](storage-only-dereplication.md); identical totals alone are insufficient.

The additional SDM regressions for explicit variant policies are `tests/test_diff_quality_pipeline.py`, `tests/test_paired_diff_quality_pipeline.py`, `tests/test_derep_behavior_parity.py`, `tests/test_seed_subclusters.py`, and the output/run cases in `tests/TestCoarseDerep.cpp`.
