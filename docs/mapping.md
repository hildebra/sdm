# Mapping-file reference

[Back to the README](../README.md) · [Filtering options](options.md)

Use `-map mapping.txt` to assign reads to samples. This file describes sample names, barcodes, primers, and optionally input filenames. Filtering settings belong in the separate [options file](options.md).

## File format

The mapping file is tab-delimited. The header may begin with `#`. Blank lines and comments are allowed before it; subsequent lines beginning with `#` are comments. All data rows must have the same number of columns. Header names are case-sensitive, may appear in any order, and unrecognized columns are ignored.

`SampleID` values must be unique. Barcode values must be unique for single-index data; the combination of `BarcodeSequence` and `Barcode2ndPair` must be unique for dual-index data.

## How sdm identifies samples

Use one identification scheme consistently within a mapping file:

- **Single barcode:** `BarcodeSequence`
- **Dual index:** `BarcodeSequence` and `Barcode2ndPair`
- **Header identifier:** `SampleIDinHead`, which matches a string in the FASTA or FASTQ header
- **One sample per file:** unique `fastqFile` or `fnaFile` values with `-i_path`

## Mapping columns

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

When `-i_path` is used, the map must include either `fastqFile` or `fnaFile`, and the command must supply `-o_fna` or `-o_fastq`. Paths in the file columns are resolved relative to `-i_path`. For example, `-i_path reads` and `fastqFile` value `sample_1.fastq.gz` refer to `reads/sample_1.fastq.gz`.

## Mapping examples

These examples contain literal tabs. Preserve the tabs when saving the file, and replace the example sequences with those used in your experiment.

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

Run the single-barcode example with:

```sh
sdm -i_fastq reads.fastq -map mapping.txt -o_fastq filtered.fastq -options sdm_options.txt
```

For the one-file-per-sample example, use:

```sh
sdm -i_path reads -map mapping.txt -o_fastq filtered.fastq -options sdm_options.txt
```

## Paired reads and separate MID files

For paired FASTQ input, use `-paired 2` and comma-separated R1,R2 filenames, either with `-i_fastq` or in the `fastqFile` column. The mates must have matching read identifiers in the same order.

For a separate MID (barcode) file, use `-i_MID_fastq` on the command line or `MIDfqFile` in the map. Alternatively, `-paired 3` accepts three comma-separated FASTQ filenames in R1,R2,R3 order, with R2 containing the MID reads.

## Primer handling

The map supplies primer sequences; the options file controls how they are handled. `keepPrimerSeq` controls whether matched primers remain in output, and `maxPrimerErrs` sets the allowed mismatches. Set `RejectSeqWithoutFwdPrim` or `RejectSeqWithoutRevPrim` to `T` to require the corresponding primer. Both requirements default to `F`.

GoldenAxe processing requires both forward and reverse primers for each sample.

## Troubleshooting

- If sdm reports a column-count mismatch, check that every row has the same number of tabs as the header, including tabs for empty cells.
- If reads are not assigned as expected, check column capitalization, barcode sequences, and the chosen identification scheme. For header matching, use strings that distinguish samples unambiguously.
- If input files cannot be found, check the filenames relative to `-i_path`.

Run `sdm -help_map` for the reference supplied with your installed version.
