# sdm — simple demultiplexer
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sdm/badges/downloads.svg)](https://anaconda.org/bioconda/sdm)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sdm/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/sdm)


**sdm** is a fast, memory-efficient command-line program for demultiplexing and quality-filtering DNA sequencing reads. It supports FASTA, FASTQ, and FASTA-plus-quality input; single-end and paired-end runs; barcode/MID-aware workflows; read merging; and dereplication.

It is intended for reproducible sequencing-data preprocessing pipelines where a compact command-line interface and detailed output control are needed.

Originally developed as a simple demultiplexer, sdm has grown into the initial filtering step of many sequencing workflows. It is written in C++11 and is designed to efficiently handle common sequence inputs, including FASTA, FASTQ, and—when gzip support is enabled—compressed input. It can detect input formats and, when configured, continue processing beyond corrupt FASTQ records.

## Highlights

- Demultiplex FASTA and FASTQ reads using a sample mapping (`-map`) file.
- Demultiplex multi-sample sequence files using DNA sequences or header keywords.
- Process single-end, paired-end, and three-file MID read layouts.
- Automatically detect FASTA or FASTQ for `-i` input.
- Convert between FASTA and FASTQ output formats.
- Read FASTA quality data from a matching `.qual_` file when available.
- Write FASTA, FASTQ, quality, per-sample demultiplexed, and dereplicated output.
- Remove sequencing and PCR primers and filter reads through an sdm options file.
- Synchronize read pairs and optionally recover from corrupt FASTQ records.
- Extract named read subsets, including from paired-end input.
- Combine multiple inputs into one filtered sequence output.
- Support externally merged paired reads and built-in pair-merging workflows.
- Produce length, quality, merge-position, and length-versus-quality reports.
- Limit reads for quick test runs or downsampling.
- Support multithreaded processing and threaded I/O.
- Include specialized GoldenAxe support for PacBio concatenated-read libraries.

## Quick start

### Prerequisites

- A C++ toolchain supported by the project configuration.
- Visual Studio is supported through the included `sdm3.vcxproj` project.
- Optional gzip support depends on the features enabled when the executable is built.

### Build on Windows with Visual Studio

1. Open `sdm3.vcxproj` in Visual Studio.
2. Select a configuration and platform, such as **Release** and **x64**.
3. Build the project.
4. Run the generated `sdm` executable from a terminal or a pipeline.

At startup, sdm reports its version and compiled capabilities, including gzip and multithreading support.

### Minimal example

```text
sdm -i test.fna -map mapping.txt -o_fna filtered.fna
```

For FASTA input, sdm looks for a matching quality file named `test.qual_` when `-i_qual` is not specified.

For FASTQ input:

```text
sdm -i_fastq reads.fastq -map mapping.txt -o_fastq filtered.fastq
```

## Input modes

Choose one primary input mode per run:

| Option | Description |
| --- | --- |
| `-i_path <path>` | Directory or path containing input files. |
| `-i <file>` | Sequence file; sdm detects FASTA or FASTQ. |
| `-i_fastq <file>` | FASTQ input. |
| `-i_fna <file>` | FASTA input. |
| `-i_qual <file>` | Quality file corresponding to FASTA input. |
| `-i_MID_fastq <file>` | FASTQ file containing MID sequences. |

Use `-paired <mode>` to describe the input layout:

- `-paired 1` — singleton input (default).
- `-paired 2` — paired input.
- `-paired 3` — R1/R3 paired input plus a separate R2 MID file.

Paired filenames supplied to a single argument are comma-delimited. Use `-onlyPair 1` or `-onlyPair 2` to process only one member of a pair.

## Mapping and filtering

Use `-map <file>` to provide the tab-delimited sample mapping file used for demultiplexing. Without `-map`, sdm performs quality filtering without assigning reads to samples.

Use `-options <file>` to select a tab-delimited filtering option file. The default is `sdm_options.txt`.

Run the built-in reference help for the current executable:

```text
sdm -help_options
sdm -help_map
```

## Mapping-file reference

### File format

The mapping file is tab-delimited. Its first row is the header and may begin with `#`; subsequent lines beginning with `#` are comments. All data rows must have the same number of columns. Header names are case-sensitive, may appear in any order, and unrecognized columns are ignored.

`SampleID` values must be unique. Barcode values must be unique for single-index data; the combination of `BarcodeSequence` and `Barcode2ndPair` must be unique for dual-index data.

### How sdm identifies samples

Use one identification scheme consistently within a mapping file:

- **Single barcode:** `BarcodeSequence`
- **Dual index:** `BarcodeSequence` and `Barcode2ndPair`
- **Header identifier:** `SampleIDinHead`, which matches a string in the FASTA or FASTQ header
- **One sample per file:** unique `fastqFile` or `fnaFile` values with `-i_path`

### Mapping columns

| Group | Column | Description |
| --- | --- | --- |
| Identity | `SampleID` | Unique output sample identifier. |
| Identity | `CombineSamples` | Combine the row's sample into the specified sample identifier. |
| Identity | `SampleIDinHead` | Header string used instead of barcode matching. |
| Identity | `SequencingRun` | Sequencing-run identifier. |
| Barcode | `BarcodeSequence` | Primary barcode/MID; IUPAC bases are supported. |
| Barcode | `Barcode2ndPair` | Second-read barcode for dual-index input. |
| Primers | `ForwardPrimer` | Forward amplification primer. |
| Primers | `LinkerPrimerSequence` | Backward-compatible alias for `ForwardPrimer`. |
| Primers | `ReversePrimer` | Reverse amplification primer. |
| Primers | `HetSpacerFwd`, `HetSpacerRev` | Forward and reverse heterogeneity spacers. Both columns are required for explicit spacer handling. |
| Input files | `fastqFile` | Relative FASTQ path, resolved from `-i_path`. |
| Input files | `fnaFile` | Relative FASTA path, resolved from `-i_path`. |
| Input files | `qualFile` | Relative FASTA-quality path. If omitted, sdm derives `<fasta-base>.qual_`. |
| Input files | `MIDfqFile` | Relative FASTQ path containing MID sequences. |
| Dereplication | `derepMin` | Empty or a non-negative, sample-specific minimum copy count. |

When `-i_path` is used, the map must include either `fastqFile` or `fnaFile`.

### Mapping examples

Single-barcode demultiplexing:

```text
#SampleID	BarcodeSequence	ForwardPrimer
sample_1	ACGTACGT	CCTACGGG
sample_2	TGCATGCA	CCTACGGG
```

Dual-indexed reads:

```text
#SampleID	BarcodeSequence	Barcode2ndPair
sample_1	ACGTACGT	GATCAAGT
```

One FASTQ file per sample with `-i_path`:

```text
#SampleID	fastqFile
sample_1	sample_1.fastq.gz
```

## Option-file reference

The `-options` file contains tab-separated `setting<TAB>value` lines. Lines beginning with `#` are comments. A leading `*` creates an alternate setting for the modified output stream; it is supported for `minSeqLength`, `minAvgQuality`, `maxAmbiguousNT`, `QualWindowThreshhold`, `TrimWindowThreshhold`, `maxPrimerErrs`, `keepPrimerSeq`, `maxAccumulatedError`, `RejectSeqWithoutRevPrim`, and `RejectSeqWithoutFwdPrim`.

### Quality and length filtering

| Setting | Default | Purpose |
| --- | --- | --- |
| `minSeqLength` / `maxSeqLength` | `250` / `1000` | Accepted length after barcode, primer, adapter removal, and trimming. |
| `minAvgQuality` | `25` | Minimum whole-read average quality. |
| `maxAmbiguousNT` | `0` | Maximum ambiguous bases. |
| `maxHomonucleotide` / `trimHomonucleotide` | `12` / `12` | Reject or trim homonucleotide runs. |
| `QualWindowWidth` / `QualWindowThreshhold` | `50` / `0` | Whole-read sliding-window check; a threshold of `0` disables it. |
| `TrimWindowWidth` / `TrimWindowThreshhold` | `15` / `20` | 3′ low-quality trimming window. |
| `maxAccumulatedError` | `-1` | Maximum accumulated expected error; `-1` disables it. |
| `BinErrorModelMaxExpError` | `2.5` | Maximum expected error for the binomial model. |
| `BinErrorModelAlpha` | `-1` | Binomial model alpha; use `-1` to disable or a value from `0` to `1`. |

### Barcode, primer, and read settings

- `maxBarcodeErrs` (default `0`) and `maxPrimerErrs` (default `1`) set accepted mismatches.
- `keepBarcodeSeq` and `keepPrimerSeq` use `0` to remove matched sequences or `1` to retain them.
- `RejectSeqWithoutFwdPrim` defaults to `T`; `RejectSeqWithoutRevPrim` defaults to `F`.
- `ExtensivePrimerChecks` enables additional primer matching checks when set to `T`.
- `TechnicalAdapter` supplies an adapter sequence to remove.
- `TrimStartNTs` removes leading nucleotides; `TruncateSequenceLength` limits ordinary filtered reads and unmerged dereplication output (`-1` disables it). As in the historical workflow, a smaller positive cap lowers the effective absolute `minSeqLength` for those ordinary reads; this adjustment is now independent of option-file line order. Merged sequences retain the configured minimum and are not capped. HQ mate files retain full length for subsequent merging, and seed extension ignores the truncation setting for all its outputs.
- `fastqVersion` accepts `auto` or a numeric version recognized by sdm.
- `PEheaderPairFmt` controls paired-end output header formatting.

### Paired reads, dereplication, and adapters

- `AmpliconShortPE`, `CheckForMixedPairs`, and `CheckForReversedSeqs` accept `T` or `F`.
- `-derepSrchLen` on the command line takes precedence over `derepSrchLen` in the option file. If neither is set, the search length defaults to a positive `minSeqLength`, otherwise 150. A value of `-1` uses the full sequence.
- With `-illuminaClip 1`, provide adapter sequences through `illuminaFwd`, `illuminaRev`, `illuminaSngUni`, and `illuminaSngIdx` as needed.

## Output modes

| Option | Description |
| --- | --- |
| `-o_fna <file>` | FASTA output. |
| `-o_qual <file>` | Quality output paired with FASTA output. |
| `-o_fastq <file>` | FASTQ output; takes precedence over `-o_fna` and `-o_qual`. |
| `-o_fna2`, `-o_qual2`, `-o_fastq2` | Additional output for reads passing the less stringent filter tier. |
| `-o_demultiplex <path>` | Write individual demultiplexed files. |
| `-o_demultiplex_gz <0|1>` | Gzip demultiplexed output. |
| `-o_dereplicate <file>` | Write dereplicated DNA reads. |
| `-o_fastq_noBC <file>` | Write reads for which no barcode was detected. |
| `-log <file|nolog>` | Set the demultiplexing log destination or disable logging. |

Additional output controls include `-o_qual_offset`, `-oneLineFastaFormat`, `-maxReadsPerOutput`, `-pairedDemulti`, and `-pairedRD_HD_out`.

`-maxReadsPerOutput` counts a single read, a read pair, or a merged read as one output record group. The file generation and its count persist across successive input files. `-XfirstReadsWritten` uses the same counting rule for the total output limit.

`-DemultiBPperSR` caps emitted bases for demultiplexed output; `0` means unlimited. Both mates count toward the ordinary-output budget, and merged output has an independent budget of the same size. Reads and pairs are admitted whole, so a budget can retain fewer bases than its limit. A merged read exceeding its budget is skipped without falling back to paired output. Counts must be nonnegative whole numbers; this option accepts scientific notation within the unsigned 32-bit range.

Dereplicated FASTQ qualities, including merged consensus qualities, average the observations covering each base. Replacing a representative preserves this evidence, and absorbing subgroups combines their quality evidence with their occurrence counts. Final dereplication output is opened only when it is written, so a missing input cannot truncate an existing result during option setup.

Dereplication batches updates into 64 independently locked hash tables. Updates for one sequence key remain serialized, including sample counts, quality accumulation, and representative selection; unrelated keys can proceed concurrently. This particularly helps workloads with many actively repeated sequence groups. One dominant sequence still shares one lock. Worker batches retain pending reads until the dereplicator flushes them, including when a worker exits. Output snapshots include all merge-position subgroups, preserve the live groups, and refresh ordinary dereplication files on each snapshot; per-sequencing-run output remains cumulative. Buffered output failures are checked through close, and a replacement abundance map is promoted only after successful writing.

Multithreaded dereplication does not guarantee insertion order or byte-identical results between runs. As before, scheduling can change the chosen representative and whether a failed read arrives before the first passing read for its group; the latter can affect abundance counts. Equal-abundance output records may also change order. Use one worker when repeatability is required. `-derep_reserve` remains the total requested reservation, divided across the tables.


For paired output, give two comma-separated filenames to `-o_fna` or `-o_fastq`, and to `-o_qual` when requesting separate qualities. Concurrent writers keep R1/R2 and their QUAL records in matching order, including across output-file rotation. Mates with different passing quality tiers are written together to the less stringent tier when that output is available. No-barcode output accepts one filename for single reads, or one interleaved filename or two separate filenames for paired reads.

If input parsing or a worker fails, all submitted workers finish before the error is returned and their buffers and output objects are released.

## Common workflows

### Quality-filter FASTQ reads

```text
sdm -i_fastq reads.fastq -o_fastq filtered.fastq -options sdm_options.txt
```

### Demultiplex paired reads

```text
sdm -i_fastq reads_R1.fastq,reads_R2.fastq -paired 2 -map mapping.txt -o_demultiplex demultiplexed
```

### Use a separate MID read file

```text
sdm -i_fastq reads_R1.fastq,reads_R3.fastq -i_MID_fastq reads_R2.fastq -paired 2 -map mapping.txt -o_fastq filtered.fastq
```

### Dereplicate reads

```text
sdm -i_fastq reads.fastq -o_dereplicate dereplicated.fna -min_derep_copies 2
```

`-min_derep_copies` accepts rules such as `10:1,3:3`, meaning at least 10 copies in one sample or three copies in three samples. Select `-dere_size_fmt 0` for `size=X;` headers or `-dere_size_fmt 1` for `_X` headers.

### Test a subset of a dataset

```text
sdm -i_fastq reads.fastq -o_fastq subset.fastq -XfirstReadsRead 10000 -XfirstReadsWritten 10000
```

## Useful command-line options

### Read handling

- `-i_qual_offset <0|auto|n>` — input FASTQ quality offset.
- `-o_qual_offset <n>` — output FASTQ quality offset; default: `33`.
- `-ignore_IO_errors <0|1>` — continue after corrupted FASTQ records; default: `0`.
- Paired-read and MID identifiers are always checked for synchronized input order.
- `-mergedPairs <0|1>` — indicate reads were merged externally.
- `-5PR1cut <n>` and `-5PR2cut <n>` — remove leading bases from read 1 or 2.

### Filtering and reporting

- `-onlyCanonicalNTs <0|1>` — retain only canonical nucleotides.
- `-strictPositiveQualsOut <0|1>` — require positive output quality values.
- `-specificReads <file>` and `-excludeFile <file>` — select named reads.
- `-length_hist <file>` and `-qual_hist <file>` — write length and quality histograms.
- `-merg_readpos <file>` and `-qual_readpos <file>` — request optional merge and quality position reports. Collection is disabled unless at least one report is requested. Positions are zero-based in the original mate orientation; reports support unequal mate lengths and concurrent merging.
- `-logLvsQ <file|1>` — write length-versus-quality data.

### Performance

- `-threads <n>` — worker thread count; default: `1`.
- `-threadIO <0|1>` — enable threaded I/O; default: `1`.
- `-iniBlockSize <n>` — initial read block size; default: `120`.

### Merging, dereplication, and OTU support

- `-merge_pairs_filter`, `-merge_pairs_seed`, and `-merge_pairs_demulti` — enable pair merging in the corresponding workflow.
- `-merge_pairs_derep <0|1>` — merge pairs while dereplicating.

- `-derep_reserve <n>`, `-derepPerSR <0|1>`, `-derep_format <text>`, and `-derepSrchLen <n>` — tune dereplication.
- `-OTU_fallback <file>` and `-OTU_fallback_refclust <file>` — fallback OTU FASTA sequences.
- `-otu_matrix <file>` — OTU matrix output.

Merging uses the retained overlap evidence from both mates, after physical primer/barcode removal, without reversing or shortening the input reads in place. The resulting consensus is filtered with the configured quality and sequence checks. Primer/barcode removal is preserved across read-through overlaps; keeping primers or barcodes retains them in the merged amplicon. `TruncateSequenceLength` does not cap merged output, including merged dereplication and demultiplexed output. The configured `minSeqLength` is checked again after quality trimming before merged output is written. A rejected merge falls back to the existing eligible mate-output path. Seed extension consumes already processed HQ mates and preserves the full consensus of the selected pair without repeating quality trimming or 5′/barcode/primer cuts. `TruncateSequenceLength` is ignored for all seed-extension outputs, including FASTQ, FASTA, unmerged mates and fallback seeds; minimum/maximum length checks still apply to merged seeds. Dereplication HQ mates retain their full length so that a clustering length cap cannot remove their overlap. If HQ files were generated with the truncating behavior in 3.44–3.45, regenerate them from the original reads before rerunning seed extension. With `-merge_pairs_seed 1 -o_fna OTU.fna`, paired input produces merged seeds directly in the main `OTU.fna` (and `-o_qual`, if requested), using the OTU identifier without a mate suffix. A seed without an eligible merge uses its selected R1 or fallback sequence in that file. With two FASTA output paths, merged seeds go to the first path; FASTQ seed output retains its separate merged file. Representative scoring is unchanged.


### Specialized processing

- `-GoldenAxe <0|1>` — enable GoldenAxe processing for PacBio concatenated reads.
- `-GoldenAxeMinAmpli <n>` and `-GoldenAxeMaxAmpli <n>` — constrain accepted GoldenAxe amplicon counts; `-1` disables the corresponding bound. GoldenAxe requires forward and reverse primers for each sample; their index tables may differ.
- `-illuminaClip <0|1>` — enable Illumina adapter clipping.

## Complete command reference

The executable is the authoritative reference for the version you are running:

```text
sdm -help_flags
sdm -help_options
sdm -help_map
sdm -v
```

Options are supplied as pairs: `-option value`. Quote paths or values that contain whitespace.

## License

sdm is free software licensed under the [GNU General Public License, version 3 or later](https://www.gnu.org/licenses/gpl-3.0.html).

## Support and contributions

Report issues and contribute through the [sdm GitHub repository](https://github.com/hildebra/sdm). For questions, contact Falk Hildebrand at `falk.hildebrand@gmail.com`.
