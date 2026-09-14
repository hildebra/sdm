# Filtering-option file reference

[Back to the README](../README.md) · [Sample mapping](mapping.md)

Use `-options sdm_options.txt` to select the file that controls filtering, trimming, and primer handling. Sample names, barcodes, and primer sequences belong in the separate [mapping file](mapping.md).

## File format and defaults

Write one setting and value per line, separated by a **literal tab**. Setting names are case-sensitive and have no leading dash. Lines beginning with `#` are comments. Preserve the spellings used below, including `Threshhold` and `Homonucleotide`.

```text
# Example settings: adjust for your data
minSeqLength	250
maxSeqLength	1000
minAvgQuality	25
maxAmbiguousNT	0
```

Run this example with:

```sh
sdm -i_fastq reads.fastq -o_fastq filtered.fastq -options sdm_options.txt
```

If `-options` is omitted, sdm looks for `sdm_options.txt` in the current directory. **If it cannot open the file, it continues with filtering disabled.** Check the startup messages before relying on the output as filtered data.

The defaults below apply when an options file is successfully opened and the setting is omitted. `sdm -help_options` prints these defaults as a starting configuration, with optional settings commented out and a separate example of alternate thresholds.

## Length and quality

| Setting | Default | What it controls |
| --- | --- | --- |
| `minSeqLength` | `250` | Minimum accepted length after removal and trimming. Use a value greater than `1` for an absolute number of bases. |
| `maxSeqLength` | `1000` | Maximum accepted read length. |
| `minAvgQuality` | `25` | Minimum whole-read average quality. |
| `maxAmbiguousNT` | `0` | Maximum number of ambiguous bases; a negative value disables the check. |
| `maxHomonucleotide` | `12` | Reject reads with longer homopolymer runs; `0` disables the check. |
| `trimHomonucleotide` | `12` | Trim trailing homopolymer runs at or above this length; `0` disables it. |
| `QualWindowWidth` | `50` | Window size for checking average quality along the read. |
| `QualWindowThreshhold` | `0` | Minimum window-average quality; `0` disables this check. |
| `TrimWindowWidth` | `15` | Window size for trimming low-quality 3′ ends. |
| `TrimWindowThreshhold` | `20` | Quality threshold for 3′ window trimming. |
| `maxAccumulatedError` | `-1` | Trim based on accumulated expected error; `-1` disables it. |
| `BinErrorModelMaxExpError` | `2.5` | Expected-error limit for the binomial error model. |
| `BinErrorModelAlpha` | `-1` | Binomial-model alpha, from `0` to `1`; `-1` disables the model. |

`minSeqLength` also supports a fraction from `0` to `1`, interpreted relative to the maximum read length observed. Use an absolute length when you need a fixed cutoff.

## Barcodes and primers

| Setting | Default | What it controls |
| --- | --- | --- |
| `maxBarcodeErrs` | `0` | Allowed barcode substitutions by default; edit distance (0–3, including indels) with `-ontMode 1`. |
| `maxPrimerErrs` | `1` | Allowed primer mismatches. |
| `keepBarcodeSeq` | `0` | `0` removes matched barcodes; `1` retains them. |
| `keepPrimerSeq` | `0` | `0` removes matched primers; `1` retains them. |
| `RejectSeqWithoutFwdPrim` | `F` | Set to `T` to require a forward primer. |
| `RejectSeqWithoutRevPrim` | `F` | Set to `T` to require a reverse primer. |
| `ExtensivePrimerChecks` | `F` | Set to `T` for additional primer matching checks. |

### Primer-variant reporting

These settings report substitution variants of primers. Reporting is disabled until `PrimerVariantMinCount` is positive; they do not replace `maxPrimerErrs`.

| Setting | Default | What it controls |
| --- | --- | --- |
| `PrimerVariantMaxMismatches` | `3` | Maximum substitutions considered; accepts `1` to `3`. |
| `PrimerVariantMinCount` | `0` | Minimum variant count; `0` disables reporting. |
| `PrimerVariantMinPercent` | `0` | Minimum percentage for reporting a variant; accepts `0` to `100`. The help example also uses `0`. |
| `PrimerVariantWarningPercent` | `50` | Percentage at which a variant triggers a main-log warning; accepts `0` to `100`. |

## Trimming and read formatting

| Setting | Default | What it controls |
| --- | --- | --- |
| `TrimStartNTs` | `0` | Number of leading bases to remove. |
| `TruncateSequenceLength` | `-1` | Length cap for ordinary filtered reads and unmerged dereplication output; `-1` disables it. |
| `TechnicalAdapter` | Unset | Adapter sequence to remove. Omit this setting if unused. |
| `fastqVersion` | `auto` | Input quality encoding: `auto`, `1` (offset 33), `2` (offset 64), or `3` (legacy offset 59). |
| `PEheaderPairFmt` | `0` | `0` preserves paired headers; `1` converts Illumina headers to `/1` or `/2`. Conversion to mode `2` is not implemented. |

`fastqVersion` is a version selector, while `-i_qual_offset` takes a raw offset. If both are supplied, the options-file setting is applied after the command-line offset. HTS input bypasses both and uses numeric Phred qualities.

A positive `TruncateSequenceLength` below an absolute `minSeqLength` lowers the effective minimum for ordinary reads to that cap. This does not depend on the order of the lines in the file.

Merged reads retain the configured minimum and are not capped by `TruncateSequenceLength`. High-quality mate output retains full-length reads for later merging. Seed extension ignores this cap for all outputs; merged seeds still undergo minimum- and maximum-length checks.

## Paired reads and adapters

| Setting | Default | What it controls |
| --- | --- | --- |
| `AmpliconShortPE` | `F` | Set to `T` to check for reverse primers on read 1 in short-amplicon data. |
| `CheckForMixedPairs` | `F` | Set to `T` to check for switched mates. |
| `CheckForReversedSeqs` | `F` | Set to `T` to check for reversed sequences. |
| `illuminaFwd`, `illuminaRev` | Unset | Forward and reverse paired-end adapter sequences. |
| `illuminaSngUni`, `illuminaSngIdx` | Unset | Universal and indexed single-end adapter sequences. |

Use `T` or `F` for the paired-read checks. Input identifiers are always checked for synchronization, independently of these settings.

Illumina adapter clipping also requires `-illuminaClip 1` on the command line. Supply the adapter sequences relevant to your library and omit unused adapter lines.

## Dereplication

`derepPrefix` selects 100% identity matching at the same start position. Set it in the options file, or use `-derepPrefix auto|0|1` on the command line. The command-line setting takes precedence.

| Value | Behavior |
| --- | --- |
| `auto` (default) | Use prefix matching for variable clustering lengths. Use the faster exact-key path when the length settings guarantee one clustering length, or when a positive `derepSrchLen` was explicitly supplied. |
| `1` | Force prefix matching across the shorter clustering sequence. This overrides `derepSrchLen`; `TruncateSequenceLength` still caps both keys. |
| `0` | Retain legacy length-matched keys and search-length behavior. |

For example, `minSeqLength=200` and `TruncateSequenceLength=200` select the exact-key path. A minimum above the truncation cap also qualifies because the ordinary-read minimum is reduced to the cap. An absolute equal minimum and maximum can also guarantee fixed lengths. A relative minimum, an uncapped variable length, or a minimum below the cap requires prefix matching unless an explicit positive search length requests legacy key comparison.

In automatic and forced-prefix modes, `TruncateSequenceLength` caps **both R1 and merged search keys**. Differences beyond this cap cannot create duplicate short sequences when the outputs are truncated for UPARSE. Full merged output and full HQ mates remain available for later use; their storage is independent of the clustering key.

`derepSrchLen` accepts a positive integer or `-1`. Its command-line value takes precedence over the options-file value. An explicitly positive value retains length-matched comparison of capped keys in automatic mode, bounded further by `TruncateSequenceLength`. With `derepPrefix=0`, the legacy default is the effective absolute `minSeqLength`, or `150` when unavailable; `-1` then compares complete logical sequences with equal lengths required.

Prefix matching allows no substitutions, gaps, offsets, or reverse-complement search. It operates on filtered logical sequences, not hidden backing-sequence tails. Longer sequences that disagree stay in separate groups. A shorter sequence compatible with several longer sequences is assigned once, to the lexicographically first maximal passing sequence. That maximal sequence supplies the representative; existing quality and merge-position selection still applies within its exact-key group.

The implementation retains parallel hash collection, then sorts unique passing keys at output time. Failed and alternate-tier reads are buffered separately and contribute only if they match a passing prefix group. They cannot seed a group or bridge conflicting longer sequences. Recovery is reported separately in the dereplication report, without changing filtering results. Retaining full variable-length keys and deferred failures can increase memory use; the fixed-length path avoids the prefix index and deferred-failure buffer.

See [implementation measurements and validation](derep-prefix-evaluation.md) for runtime and memory tradeoffs.

### Streaming coarse dereplication

SDM **3.52 beta** adds selected-quality retention for exact single-read and paired variants. Use the same `derepStoreQuals` flag for either input layout; paired input does not require an additional retention flag. Boolean flags below accept `0` or `1`. In commands, prepend `-`; in the options file, use the names shown.

| Option (CLI or options file) | Default | Behavior |
| --- | --- | --- |
| `derepIdentity` | `100` | Percentage in `(0,100]`, with up to six decimal places. Values below 100 select streaming substitution clustering. |
| `derepSeedHomopolymer` | `4` | Exclude a search 20-mer if it contains this many consecutive identical bases; integer `2..20`. |
| `derepReassign` | `0` | `1` enables a final pass moving exact variants and their sample counts to strictly closer qualifying high-quality representatives. Automatically retains variants/sample provenance in memory; does not enable sidecar output. |
| `derepStoreDiffs` | `0` | `1` retains counted nucleotide variants and writes `.diff`: `SDMDIFF2` for single sequences, or linked `SDMDIFF3` for pairs when `derepStoreQuals=1`. Also selects streaming mode at identity 100. Neither binary format stores qualities. |
| `derepStoreQuals` | `0` | `1` retains selected qualities per exact sequence or `(R1,R2)` variant. Writes `.subclusters.fq` for single reads or synchronized `.subclusters.1.fq`/`.2.fq` for pairs. Requires complete unmerged reads and qualities; binary export remains optional. |
| `derepSubclusterFasta` | `0` | `1` writes exact sequences with their own `;size=N;` counts. Output is `.subclusters.fna`, or synchronized `.subclusters.1.fna`/`.subclusters.2.fna` for paired quality retention. Requires `derepStoreDiffs=1`. |

For 97% clustering with selected-quality FASTQ output, an options-file configuration is:

```text
derepIdentity	97
derepStoreQuals	1
derepStoreDiffs	0
derepSubclusterFasta	0
derepReassign	0
```

Set `derepStoreDiffs=1` to additionally save binary nucleotide deltas, `derepSubclusterFasta=1` with those deltas for FASTA export, or `derepReassign=1` for the final reassignment pass. These outputs and reassignment are independently optional; quality retention itself writes subcluster FASTQ.

Command-line settings take precedence over file entries. Streaming mode takes precedence over `derepPrefix`. It uses full logical sequences unless an explicit positive `derepSrchLen` or `TruncateSequenceLength` caps the clustering keys. A minimum length alone does not cap the comparison. Identity is checked at the same start across the shorter key, allowing at most `floor(length * (100 - derepIdentity) / 100)` substitutions. No gaps, orientation search, or ignored homopolymer positions are involved.

Clusters retain their initial matching anchor and one representative selected by the existing quality rules. All assigned reads contribute to the usual sample counts and minimum-copy filtering. Per-base quality sums are omitted: coarse FASTQ output uses the selected representative's qualities. Failed reads may join an existing passing cluster, but cannot seed one; unlike deferred 100% prefix mode, unmatched failures are not retained for later recovery.

This is greedy streaming clustering: results depend on read order and worker scheduling. Use `-threads 1` for repeatability. A short initial anchor can admit longer reads whose tails disagree; membership is defined against the initial anchor, not by all-pairs identity. The output representative can differ from the anchor.

All dereplicates, including default 100% exact/prefix groups, use `diffDNA` for their sparse sample occurrences and total abundance. Without differences, quality retention, or reassignment, no nucleotide reference or variant table is allocated in this object. Enabling any of these options also retains exact variants. For reassignment, each variant's sparse sample entries provide both its multiplicity and the aggregate cluster sample map; a separate variant-count table is unnecessary.

When nucleotide storage is enabled, each distinct full processed sequence and length is counted separately inside `diffDNA`, including repeated non-reference variants. Coarse-cluster identity does not replace these exact subcluster counts. Without paired quality retention, the nucleotide-only `SDMDIFF2` sidecar covers the processed search sequence (merged when available and enabled, otherwise R1), including bases beyond a search-only cap. Paired quality retention instead writes both processed mates together in `SDMDIFF3`. Neither format retains raw pre-filter tails, per-read names, qualities, or per-read sample labels. Aggregate sample counts remain in `.map`. Per-run output receives the same run suffix as the corresponding main output. Groups below the minimum-copy cutoff are also included, as they are in `.map` and `.rest`.

For single-end input with `derepStoreQuals=1`, `.subclusters.fq` reconstructs the full processed sequence and the selected observation's corresponding quality vector for each exact variant. The shared dereplication selector compares stored scoring summaries, so it need not expand the previous sequence for each comparison. Passing observations take priority over failed observations. A failed read admitted to an existing passing cluster can supply the first quality witness for its exact variant; later passing observations replace it. Subsequent failed observations only add counts. Quality vectors are never averaged. Each record has its own `;size=N;` count; original per-read names are discarded. Output uses the configured FASTQ quality offset and includes groups below the minimum-copy cutoff, with the same per-run suffix as other dereplication outputs. At identity 100 this flag selects streaming exact-prefix mode. Quality vectors are retained in memory and exported as FASTQ; they are not serialized in the nucleotide-only binary formats.

With paired input and `derepStoreQuals=1`, exact identity includes both full processed mates and both lengths. Both selected quality vectors belong to the same observation; each pair contributes one fragment count. The initial references stay fixed for both mates. The candidate index uses R1, then verifies the identity threshold independently on R1 and R2; each search cap applies independently. A new cluster or replacement witness must pass filtering on both mates. Complete failed pairs can contribute to a qualifying existing cluster, but cannot seed one or replace a passing representative. Pairs with a mate completely removed by filtering are skipped. Missing mates and native merging (`merge_pairs_derep`, `merge_pairs_filter`, `merge_pairs_demulti`) are rejected. Paired output uses matching names, abundances, and record order. Optional binary output uses `SDMDIFF3`, with two references and linked pair patches; optional FASTA uses `.subclusters.1.fna` and `.subclusters.2.fna`. These paired semantics apply to quality retention; existing dereplication with that flag disabled keeps its previous behavior.

For paired final reassignment, both mates must qualify against the target. Closeness is the total substitutions divided by the sum of the two shorter comparison lengths; ties retain the original group. Both nucleotide patches, quality payloads, and shared sample counts move together. See [paired implementation and validation](paired-diff-quality.md).

With `derepReassign=1`, each exact variant is compared against the fixed final representatives using the same identity threshold and search caps. Its copies move together only if another representative is strictly closer and meets the threshold. A tie keeps the original group; equal best alternatives use the earliest-created target. Members with no better qualifying target stay, even when representative replacement has put them outside the configured radius. The pass does not merge entire clusters, create new clusters, iterate, or choose new representatives.

Reassignment also applies to variants from failed reads that were already admitted to a passing cluster. Per-variant sample abundances are retained internally only for this option, so `.map` remains correct after moving variants. These sample details are not written to the nucleotide-only `.diff` format. The main output, `.map`, minimum-copy decisions, and optional `.diff`/`.subclusters.fna`/`.subclusters.fq` exports agree on final membership. Output snapshots leave streaming anchors and counts unchanged.

See [implementation, delta format, and limitations](coarse-dereplication.md) and the [Q30 200 bp benchmark](coarse-dereplication-q30.md).

## A second filtering tier

A leading `*` supplies an alternate setting for reads that do not pass the primary tier. Configure the associated output with `-o_fastq2`, or `-o_fna2` and optionally `-o_qual2`.

For example, this configuration distinguishes reads passing average quality 25 from those passing a lower threshold of 20:

```text
minAvgQuality	25
*minAvgQuality	20
```

```sh
sdm -i_fastq reads.fastq -o_fastq primary.fastq -o_fastq2 additional.fastq -options sdm_options.txt
```

The `*` prefix is supported for `minSeqLength`, `minAvgQuality`, `maxAmbiguousNT`, `QualWindowThreshhold`, `TrimWindowThreshhold`, `maxPrimerErrs`, `keepPrimerSeq`, `maxAccumulatedError`, `RejectSeqWithoutRevPrim`, and `RejectSeqWithoutFwdPrim`. Put the primary setting before its alternate and specify the alternate thresholds you need explicitly.

For paired output, mates with different passing tiers are written together to the less stringent tier when that output is available.

## More help

Run `sdm -help_options` for example settings supplied with your installed version, or `sdm -help_flags` for command-line controls such as input paths, output files, merging, and worker counts.

## ONT end matching

Command-line options: `-ontMode 0|1` (default `0`), `-barcodeSearchWindow 1..10000` (ONT default `200`), and `-ontBarcodeEnds either|both` (default `either`). The latter two require ONT mode. See [supported layouts, examples and integration requirements](ont.md). Existing primer, quality, length and output options still apply.
