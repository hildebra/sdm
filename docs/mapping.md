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
- **One sample per file:** unique `fastqFile`, `alignmentFile` or `fnaFile` values with `-i_path`
- **SAM/BAM/CRAM sample barcode:** `BarcodeSequence` (and optionally `Barcode2ndPair`) with `-htsBarcode BC`; see [alignment input](hts.md)

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
| Input files | `fastqFile` | Relative FASTQ path, resolved from `-i_path`; use comma-separated R1,R2 paths for paired FASTQ. Leave empty for a sample using `alignmentFile`. |
| Input files | `alignmentFile` | One relative SAM, BAM or CRAM path, resolved from `-i_path`. This is a separate column, mutually exclusive with a nonempty `fastqFile` entry in the same row. |
| Input files | `fnaFile` | Relative FASTA path, resolved from `-i_path`. |
| Input files | `qualFile` | Relative FASTA-quality path. If omitted, sdm derives `<fasta-base>.qual_`. |
| Input files | `MIDfqFile` | Relative FASTQ path containing MID sequences. |
| Dereplication | `derepMin` | Empty or a non-negative, sample-specific minimum copy count. |

When `-i_path` is used, the map must include `fastqFile`, `alignmentFile` or `fnaFile`, and the command must supply `-o_fna` or `-o_fastq`. Paths in the file columns are resolved relative to `-i_path`. For example, `-i_path reads` and `fastqFile` value `sample_1.fastq.gz` refer to `reads/sample_1.fastq.gz`.

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

## SAM, BAM and CRAM file columns

Use the fixed, case-sensitive header `alignmentFile` for all three alignment formats. An alignment-only map can omit `fastqFile` entirely:

```text
#SampleID	alignmentFile
sample_sam	sample_sam.sam
sample_bam	sample_bam.bam
sample_cram	sample_cram.cram
```

```sh
sdm -i_path reads -map alignments.tsv -o_fastq filtered.fastq -options sdm_options.txt
```

For a mixed FASTQ/alignment run, include both columns and fill exactly one per sample:

```text
#SampleID	alignmentFile	fastqFile
sample_fq		sample_fq.fastq.gz
sample_bam	sample_bam.bam	
sample_cram	sample_cram.cram	
```

The last two rows end with a tab for the empty `fastqFile` cell. Keep that tab when saving the map; do not write a dash or `NA` for an empty cell. Either column order is supported. The FASTQ header is `fastqFile`, not `fqFile`.

With `alignmentFile` present, every row must fill exactly one of `alignmentFile` or `fastqFile`. SDM rejects rows with both or neither, non-alignment input in `alignmentFile`, and comma-separated alignment paths. Each alignment file contains either singleton reads or both mates of every pair. Mixed FASTQ and alignment inputs must have the same singleton/paired layout. Use a separate mapping run for FASTA inputs. `qualFile` is not supported in these maps, and `MIDfqFile` must be empty for alignment rows.

SDM detects the alignment format as it does for command-line input, including extensionless BAM/CRAM files. Existing maps that put alignment paths in `fastqFile` alone remain readable for compatibility; use `alignmentFile` in new maps.

For externally referenced CRAM, add `-cramRef reference.fa`; this option applies to the whole run and defaults to no explicitly supplied reference. Supply the matching reference with its indexes as described in [CRAM reference handling](hts.md#cram-references). The default HTS-enabled build is required for alignment input.

To demultiplex a shared alignment file using `BC` tags, repeat its path for each sample, supply barcode columns, and enable `-htsBarcode BC`:

```text
#SampleID	alignmentFile	BarcodeSequence
sample_A	shared.bam	ACGT
sample_B	shared.bam	TTTT
```

```sh
sdm -i_path reads -map barcodes.tsv -htsBarcode BC \
  -o_fastq filtered.fastq -options sdm_options.txt
```

Repeated paths are grouped and read once. Without barcode or header assignment, paths must be unique for one-sample-per-file assignment. Inline barcode matching remains the default; see [alignment barcode modes](hts.md#standard-sample-barcode-tags).

## Paired reads and separate MID files

For paired FASTQ input, use `-paired 2` and comma-separated R1,R2 filenames, either with `-i_fastq` or in the `fastqFile` column. The mates must have matching read identifiers in the same order.

For paired SAM/BAM/CRAM input, `alignmentFile` contains one query-name-sorted file holding both mates; use `-paired 2` and separate R1,R2 output paths. Singleton alignment files may use any order. See [pairing requirements](hts.md#input-and-pairing).

For a separate MID (barcode) file with FASTQ input, use `-i_MID_fastq` on the command line or `MIDfqFile` in the map. Alternatively, `-paired 3` accepts three comma-separated FASTQ filenames in R1,R2,R3 order, with R2 containing the MID reads.

## Primer handling

The map supplies primer sequences; the options file controls how they are handled. `keepPrimerSeq` controls whether matched primers remain in output, and `maxPrimerErrs` sets the allowed mismatches. Set `RejectSeqWithoutFwdPrim` or `RejectSeqWithoutRevPrim` to `T` to require the corresponding primer. Both requirements default to `F`.

GoldenAxe processing requires both forward and reverse primers for each sample.

## Troubleshooting

- If sdm reports a column-count mismatch, check that every row has the same number of tabs as the header, including tabs for empty cells.
- If reads are not assigned as expected, check column capitalization, barcode sequences, and the chosen identification scheme. For header matching, use strings that distinguish samples unambiguously.
- If input files cannot be found, check the filenames relative to `-i_path`.

Run `sdm -help_map` for the reference supplied with your installed version.
