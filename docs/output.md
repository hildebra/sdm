# Output files and reports

[Back to the README](../README.md) · [Filtering options](options.md) · [Advanced workflows](workflows.md)

Output filenames and report settings are command-line options. Filtering thresholds belong in the [options file](options.md), and sample names come from the [mapping file](mapping.md).

## Choose sequence output

| Option | Output |
| --- | --- |
| `-o_fastq <file>` | Combined FASTQ file. Takes precedence over `-o_fna` and `-o_qual`. |
| `-o_fna <file>` | Combined FASTA file. |
| `-o_qual <file>` | Separate quality file accompanying FASTA output. |
| `-o_fastq2 <file>` | Additional FASTQ file for the less stringent filtering tier. |
| `-o_fna2 <file>`, `-o_qual2 <file>` | Additional FASTA and quality files for the less stringent tier. |
| `-o_demultiplex <directory>` | One set of FASTQ files per sample; the directory must already exist. |
| `-o_demultiplex_gz <0 or 1>` | Set to `1` to compress per-sample output with gzip. |
| `-o_dereplicate <file>` | Dereplicated sequences; see [dereplication](workflows.md#dereplication). |
| `-o_fastq_noBC <file>` | Reads for which no barcode was detected. |

The `2` suffix on output options means a **second filtering tier**, not read 2 of a pair. See [alternate thresholds](options.md#a-second-filtering-tier).

### Paired output

Give two comma-separated filenames to keep mates separate:

```sh
sdm -i_fastq R1.fastq,R2.fastq -paired 2 -o_fastq filtered_R1.fastq,filtered_R2.fastq -options sdm_options.txt
```

The same convention applies to `-o_fna` and `-o_qual`. Paired records remain in matching order, including when output files are split. Mates passing different quality tiers are written together to the less stringent tier when that output is available.

For paired no-barcode output, `-o_fastq_noBC` accepts either one interleaved filename or two comma-separated filenames. Single-end input uses one filename.

### Per-sample output

Create the destination directory before running sdm:

```sh
mkdir -p demultiplexed
sdm -i_fastq reads.fastq -map mapping.txt -o_demultiplex demultiplexed -o_demultiplex_gz 1 -options sdm_options.txt
```

Files are named from `SampleID`:

| Layout | Filenames without gzip |
| --- | --- |
| Single-end | `sample_1.fq` |
| Paired-end | `sample_1.1.fq`, `sample_1.2.fq` |
| Merged paired-end | `sample_1.merg.fq` |

Gzip output adds `.gz`. A merged file may be created even when no merged reads are written. Set `-pairedDemulti 1` to write only complete pairs to ordinary per-sample output.

## Output formatting

| Option | Purpose |
| --- | --- |
| `-o_qual_offset <n>` | Output FASTQ quality offset; default `33`. |
| `-oneLineFastaFormat <0 or 1>` | Set to `1` to write each FASTA sequence or quality record on one line. |
| `-sample_sep <text>` | Delimiter between sample name and read identifier. |
| `-pairedRD_HD_out <0 or 1>` | Control inclusion of paired-read identifiers in headers. |

Two specialized cleanup flags affect dereplicated output: `-onlyCanonicalNTs 1` replaces noncanonical bases with `A` and sets their qualities to zero; `-strictPositiveQualsOut 1` changes negative quality scores to zero. These transform output values rather than rejecting reads. Use filtering settings such as `maxAmbiguousNT` when you want to reject reads instead.

## Output limits and splitting

| Option | Purpose |
| --- | --- |
| `-maxReadsPerOutput <n>` | Split ordinary sequence output after the requested number of record groups. |
| `-XfirstReadsWritten <n>` | Limit the total number of output record groups. |
| `-XfirstReadsRead <n>` | Limit how much input is processed. |
| `-DemultiBPperSR <n>` | Cap bases emitted to demultiplexed output per sequencing run; `0` means unlimited. |

For `-maxReadsPerOutput` and `-XfirstReadsWritten`, one single read, one read pair, or one merged read counts as one record group. Output-file numbering and the current file's count continue across successive input files.

The demultiplexed base budget counts both mates. Merged output has its own budget of the same size. Reads and pairs are admitted whole, so output can stop below the limit. A merged read exceeding its budget is skipped without falling back to paired output. `-DemultiBPperSR` accepts nonnegative whole numbers, including scientific notation, up to `4294967295`.

These limits do not randomly subsample a dataset.

## Logs and reports

| Option | Report |
| --- | --- |
| `-log <file>` | Main filtering and demultiplexing log. Use `-log nolog` to disable it. |
| `-length_hist <file>` | Read-length histogram. |
| `-qual_hist <file>` | Read-quality histogram. |
| `-merg_readpos <file>` | Merge-position report. |
| `-qual_readpos <file>` | Quality-by-position report associated with merging. |
| `-logLvsQ <file or 1>` | Length-versus-quality data. |

Check the main log for retained reads and rejection reasons. Set report filenames explicitly when collecting results from multiple runs.

Detailed merge-position collection is enabled only when `-merg_readpos` or `-qual_readpos` is requested. Positions are zero-based in the original mate orientation; unequal mate lengths are supported.

## Dereplication output behavior

Dereplicated FASTQ qualities, including merged consensus qualities, average the observations covering each base. Replacing a representative preserves the accumulated quality evidence, and combining subgroups combines their evidence and occurrence counts.

Intermediate output snapshots include all merge-position subgroups. Each snapshot refreshes ordinary dereplication files; per-sequencing-run output remains cumulative. Final dereplication files are opened when written, so a missing input does not truncate an existing result during option setup.

With multiple workers, representative choice, some abundance counts, and record order can vary between runs. See [performance and repeatability](workflows.md#performance-and-repeatability).
