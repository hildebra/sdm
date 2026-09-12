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

The defaults below apply when an options file is successfully opened and the setting is omitted. The values printed by `sdm -help_options` are an example configuration, not always the defaults: notably, that example requires the forward primer, while the default does not.

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
| `maxBarcodeErrs` | `0` | Allowed barcode mismatches. |
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
| `PrimerVariantMinPercent` | `0` | Minimum percentage for reporting a variant; accepts `0` to `100`. The help example uses `1`. |
| `PrimerVariantWarningPercent` | `50` | Percentage at which a variant triggers a main-log warning; accepts `0` to `100`. |

## Trimming and read formatting

| Setting | Default | What it controls |
| --- | --- | --- |
| `TrimStartNTs` | `0` | Number of leading bases to remove. |
| `TruncateSequenceLength` | `-1` | Length cap for ordinary filtered reads and unmerged dereplication output; `-1` disables it. |
| `TechnicalAdapter` | Unset | Adapter sequence to remove. Omit this setting if unused. |
| `fastqVersion` | `auto` | Input quality encoding: automatic detection or a numeric FASTQ version recognized by sdm. |
| `PEheaderPairFmt` | `0` | Paired-end output header formatting mode. |

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

`derepSrchLen` sets the sequence length used for dereplication searches. Use a positive integer, or `-1` for the full sequence.

The command-line `-derepSrchLen` takes precedence over the options-file setting. When omitted and filtering is enabled, the value is derived from the effective positive `minSeqLength`, including any reduction by `TruncateSequenceLength`. Otherwise, the dereplicator uses its fallback of `150`.

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
