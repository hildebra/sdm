# sdm — simple demultiplexer
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sdm/badges/downloads.svg)](https://anaconda.org/bioconda/sdm)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sdm/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/sdm)


**sdm** sorts DNA sequencing reads into samples (demultiplexing) and removes reads that do not meet your quality criteria. It can also remove primers and barcodes, merge paired reads, and combine duplicate sequences (dereplication).

It accepts FASTA and FASTQ files, including gzip-compressed files, and supports single-end and paired-end data.

## Build from source

The supplied Makefile uses a C++20 compiler, GNU Make, POSIX threads, and zlib. On Linux, it links statically, so static versions of the required libraries must also be installed.

From the repository directory, run:

```sh
make -j2
./sdm -v
```

The executable is written to the repository root. The examples below assume `sdm` is on your `PATH`; otherwise, use `./sdm` from that directory. For Windows, use a Linux environment such as WSL with the same build dependencies; this checkout does not include a Visual Studio project.

## Set up filtering

Create `sdm_options.txt` with tab-separated settings, for example:

```text
minSeqLength	250
maxSeqLength	1000
minAvgQuality	25
maxAmbiguousNT	0
```

Adjust these values for your reads. An options file is not included in this checkout. **If the options file cannot be opened, sdm reports that filtering is disabled and continues without filtering.** It looks for `sdm_options.txt` in the current directory unless you supply `-options <file>`.

## Your first run

Quality-filter FASTQ reads:

```sh
sdm -i_fastq reads.fastq -o_fastq filtered.fastq -options sdm_options.txt
```

Demultiplex and filter FASTA reads using a mapping file:

```sh
sdm -i test.fna -map mapping.txt -o_fna filtered.fna -options sdm_options.txt
```

For FASTA input, sdm looks for `test.qual_` when `-i_qual` is not specified. Use `-i_qual <file>` to supply a different quality-file path. See the sample mapping example below for the structure of `mapping.txt`.

## Assign reads to samples

Use `-map mapping.txt` to tell sdm how to identify each sample. The mapping file is a text file with one sample per row and **tabs between columns**. Sample names must be unique.

For reads with a barcode at the start, a simple mapping file looks like this:

```text
#SampleID	BarcodeSequence
sample_1	ACGTACGT
sample_2	TGCATGCA
```

Replace the names and barcodes with your own. Each barcode must identify a single sample. If you use two barcodes per sample, add a `Barcode2ndPair` column; each barcode combination must be unique.

You can also add `ForwardPrimer` and `ReversePrimer` columns to describe your primers. Use `RejectSeqWithoutFwdPrim` or `RejectSeqWithoutRevPrim` with a value of `T` in your options file if reads must contain the corresponding primer.

Other layouts are supported: `SampleIDinHead` identifies samples from text in read headers, and `fastqFile` or `fnaFile` identifies samples by input filename. See the [mapping-file reference](docs/mapping.md) for all columns and examples, or run `sdm -help_map`.

Omit `-map` when you only want to filter reads without assigning them to samples.

## Common workflows

These examples assume you have created `sdm_options.txt` as described above. Replace filenames and filtering settings with those appropriate for your data.

### Write a separate file for each sample

```sh
mkdir -p demultiplexed
sdm -i_fastq reads.fastq -map mapping.txt -o_demultiplex demultiplexed -options sdm_options.txt
```

`-o_demultiplex` selects an existing directory for per-sample files. Use `-o_fastq filtered.fastq` instead when you want a combined FASTQ file, or `-o_fna filtered.fna` for FASTA output.

### Filter paired-end reads

Supply the two input filenames separated by a comma, and use `-paired 2`. Supply two output filenames to keep the mates in separate files.

```sh
sdm -i_fastq reads_R1.fastq,reads_R2.fastq -paired 2 -o_fastq filtered_R1.fastq,filtered_R2.fastq -options sdm_options.txt
```

Add `-map mapping.txt` to assign reads to samples. Paired files must have matching read identifiers in the same order.

### Use a separate barcode read file

For data with a separate MID (barcode) read file, provide it through `-i_MID_fastq`:

```sh
sdm -i_fastq reads_R1.fastq,reads_R3.fastq -i_MID_fastq reads_R2.fastq -paired 2 -map mapping.txt -o_fastq filtered_R1.fastq,filtered_R3.fastq -options sdm_options.txt
```

Here, R1 and R3 contain the paired sequences, and R2 contains the barcodes.

### Process one input file per sample

List the filenames in your mapping file:

```text
#SampleID	fastqFile
sample_1	sample_1.fastq.gz
sample_2	sample_2.fastq.gz
```

Then give the directory containing those files with `-i_path`:

```sh
sdm -i_path reads -map mapping.txt -o_fastq filtered.fastq -options sdm_options.txt
```

In this mode, include either `-o_fastq` or `-o_fna` to name the combined output.

### Combine duplicate reads

Dereplication groups repeated sequences and records their abundance. For example, to retain sequences with at least two copies:

```sh
sdm -i_fastq reads.fastq -o_dereplicate dereplicated.fna -min_derep_copies 2 -options sdm_options.txt
```

For more control over copy counts across samples, see `sdm -help_flags`.

### Try a small run first

```sh
sdm -i_fastq reads.fastq -o_fastq subset.fastq -XfirstReadsRead 10000 -XfirstReadsWritten 10000 -options sdm_options.txt
```

These limits process the first reads in the file; they do not select a random sample.

## Adjust filtering and check results

Filtering settings belong in the tab-delimited options file. Common settings are:

| Setting | What it controls |
| --- | --- |
| `minSeqLength`, `maxSeqLength` | Accepted read length after trimming. |
| `minAvgQuality` | Minimum average read quality. |
| `maxAmbiguousNT` | Number of ambiguous bases allowed. |
| `maxBarcodeErrs`, `maxPrimerErrs` | Allowed mismatches when matching barcodes and primers. |
| `keepBarcodeSeq`, `keepPrimerSeq` | Set to `1` to keep matched sequences, or `0` to remove them. |
| `RejectSeqWithoutFwdPrim`, `RejectSeqWithoutRevPrim` | Set to `T` to require a matching primer. |

See the [filtering-option reference](docs/options.md) for defaults and examples, or run `sdm -help_options`. Choose settings for your sequencing data; the examples are not suitable for every experiment.

After a run, review the log and the length and quality reports to see how many reads were retained and why others were rejected. Use `-log run.log` to choose the main log filename.

A few behaviors matter when interpreting results:

- **Filtering requires a readable options file.** Check the startup messages to confirm that filtering is enabled.
- **Length caps and merging behave differently.** `TruncateSequenceLength` caps ordinary reads but does not cap merged reads or seed-extension output. A shorter cap can also lower the effective minimum length for ordinary reads.
- **More threads can affect dereplication results.** `-threads 4` enables four workers. For repeatable dereplication runs, use `-threads 1`; with multiple workers, representative sequences, some abundance counts, and output order can vary.

## Documentation

- [Filtering-option file](docs/options.md): format, defaults, trimming, primer handling, and a second filtering tier.
- [Sample mapping file](docs/mapping.md): columns, barcode layouts, primers, and input files.
- [Output files and reports](docs/output.md): output formats, paired files, limits, logs, and quality reports.
- [Advanced workflows](docs/workflows.md): merging, dereplication, seed extension, and performance settings.

## More help

The built-in help describes the options supported by your installed version:

```sh
sdm -help_flags
sdm -help_options
sdm -help_map
sdm -v
```

These references also cover read merging, additional output formats, read selection, seed extension, and GoldenAxe processing for PacBio concatenated reads.

Except for help and version commands, command-line options take a value: `-option value`. Quote paths that contain spaces.

## Citation

If you use sdm, please cite:

Hildebrand, Falk, Raul Tadeo, Anita Voigt, Peer Bork, and Jeroen Raes. “LotuS: An Efficient and User-Friendly OTU Processing Pipeline.” Microbiome 2, no. 1 (2014): 30. https://doi.org/10.1186/2049-2618-2-30.

## License and support

sdm is free software licensed under the [GNU General Public License, version 3 or later](https://www.gnu.org/licenses/gpl-3.0.html).

Report bugs and contribute through the [sdm GitHub repository](https://github.com/hildebra/sdm). For questions, contact Falk Hildebrand at `falk.hildebrand@gmail.com`.
