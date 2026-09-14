# SAM, BAM and CRAM input

[README](../README.md) · [Mapping files](mapping.md)

The default `sdm` build includes HTSlib and reads alignment files. It uses HTSlib to recover sequencing reads and sends them through the same filtering, demultiplexing, merging and dereplication pipeline as FASTQ input. Output remains in SDM's existing formats; alignments, CIGAR strings and auxiliary tags are not preserved as SAM/BAM/CRAM output.

## Input and pairing

`-i` detects SAM/BAM/CRAM using extensions and file content, including extensionless BAM/CRAM. Explicit `-i_sam`, `-i_bam` and `-i_cram` aliases are also accepted. Each input set is one alignment file containing either singleton reads or complete paired-end reads. SDM infers the layout from the first primary record; an explicit `-paired 1` or `-paired 2` must agree. Do not pass a comma-separated pair of BAMs.

```sh
./sdm -i reads.bam -o_fastq filtered.fq -options sdm_options.txt

samtools sort -n -o paired.by-name.bam paired.bam
./sdm -i paired.by-name.bam -paired 2 \
  -o_fastq filtered.1.fq,filtered.2.fq -options sdm_options.txt
```

Paired files must be query-name sorted, with primary R1 and R2 records adjacent after secondary and supplementary records are excluded. Either mate may appear first. Both natural (`samtools sort -n`) and lexicographic (`samtools sort -N`) query-name order are accepted. SDM checks order while streaming, using bounded memory; it does not sort or retain a cache of distant mates. A coordinate-sort header is rejected for paired input. Name collation alone does not guarantee the required sorted order. If distinct names differ only in leading zeroes in numeric fields, use lexicographic sorting (`-N`) to keep the order unambiguous.

Singleton files may have any order. Mixing singleton and paired records in one file, orphan pairs, duplicate primary records within a query-name group, and ambiguous READ1/READ2 flags are rejected. `-onlyPair 1` or `-onlyPair 2` selects one mate from an otherwise valid, complete, name-sorted paired file.

Both mapped and unmapped primary records are retained, including reads flagged as duplicates or QC failures. Secondary and supplementary alignments are excluded. Reverse-strand sequence and quality are restored to sequencing orientation using SDM's existing reverse-complement helper. Soft-clipped bases are retained; hard-clipped bases were not stored and cannot be recovered.

Only regular local input files are supported in this first implementation. Input probing and restart reopen the file. Stdin, FIFOs, remote URLs and a separate MID FASTQ combined with alignment input are unsupported.

## Quality scores and errors

HTS qualities are decoded as numeric Phred scores and passed to SDM as Phred+33. FASTQ quality autodetection and a forced input offset do not apply to alignment records. Output quality-offset settings still apply.

Missing sequence, missing qualities and qualities outside the FASTQ-representable range (0–93) are errors. SDM does not invent quality values. The unusual SAM `=` sequence symbol (reference-equality encoding) is rejected; provide explicit sequence bases. Corrupt records, reference failures and invalid pairs are fatal even with `-ignore_IO_errors 1`, because skipping individual records could change read pairing.

Order and corruption checks are incremental. An error late in a file can leave partial output, which must not be treated as a completed run.

## CRAM references

`-cramRef` defaults to no explicitly supplied reference. Reference-free and embedded-reference CRAMs can be read without this option. Externally referenced CRAMs need the matching reference, either supplied explicitly or available through HTSlib's normal lookup mechanisms.

```sh
./sdm -i reads.cram -cramRef reference.fa \
  -o_fastq filtered.fq -options sdm_options.txt
```

Use a reference FASTA with its `.fai` index (`samtools faidx reference.fa`); BGZF-compressed reference FASTA also needs its `.gzi`. HTSlib's `REF_PATH` and `REF_CACHE` conventions remain available. Missing or mismatched references cause an error. No BAM/CRAM alignment index is required for sequential input. See [HTSlib reference handling](https://www.htslib.org/doc/reference_seqs.html).

## Inline barcodes and one file per sample

The default `-htsBarcode inline` uses SDM's normal inline barcode logic. An existing `BC` tag does not override that default.

For already separated samples, use the dedicated `alignmentFile` column for SAM, BAM and CRAM paths. Omit `fastqFile` when all inputs are alignment files:

```text
#SampleID	alignmentFile
sample_A	sample_A.bam
sample_B	sample_B.bam
```

```sh
./sdm -i_path /data/reads -map samples.tsv \
  -o_fastq filtered.fq -options sdm_options.txt
```

For mixed FASTQ/alignment maps, include both `alignmentFile` and `fastqFile`, filling exactly one per row. Leave the FASTQ cell empty for alignment samples, including its trailing tab when it is the final column. SDM rejects conflicting or missing entries. `MIDfqFile` and `qualFile` cannot supply data for an alignment row; FASTA inputs require a separate mapping run. See [mapping-file examples and validation](mapping.md#sam-bam-and-cram-file-columns).

All files in a run must use a compatible singleton/paired layout. Paired sample alignment files must also be query-name sorted. `-cramRef` works identically for command-line and mapping-file inputs and applies to all CRAM files in the run.

## Standard sample-barcode tags

Use `-htsBarcode BC` to assign samples from the standard `BC:Z` sample/library barcode tag. This convention is defined by the [SAM tag specification](https://samtools.github.io/hts-specs/SAMtags.pdf) and used by [Picard's IlluminaBasecallsToSam](https://gatk.broadinstitute.org/hc/en-us/articles/360056968872-IlluminaBasecallsToSam-Picard).

```sh
./sdm -i reads.by-name.bam -htsBarcode BC -map barcodes.tsv \
  -o_fastq filtered.1.fq,filtered.2.fq -options sdm_options.txt
```

For single indexes, `BC:Z:ACGT` matches `BarcodeSequence=ACGT`. For dual indexes, `BC:Z:ACGT-TGCA` matches a mapping row with `BarcodeSequence=ACGT` and `Barcode2ndPair=TGCA`. Match the recorded index orientation in the mapping file; SDM does not reverse-complement barcode tags based on alignment strand.

```text
#SampleID	BarcodeSequence	Barcode2ndPair
sample_A	ACGT	TGCA
sample_B	TTTT	CCCC
```

Tag matching is exact against the complete barcode or index pair, including for singleton records. It does not use the inline mismatch allowance or IUPAC wildcard matching. Unknown barcodes fail sample assignment through the normal filtering path. Barcode bases are not removed from the biological sequence, regardless of `keepBarcodeSeq`.

For paired reads, a tag present on either primary mate applies to both; tags present on both must agree. A missing tag on both mates, a non-string tag, or conflicting tags is an error. `QT` barcode qualities are not used for matching. `CB`/`CR` cell barcodes, `RX` molecular identifiers, and `RG` read-group assignments are separate conventions and are not used by this mode. Use a barcode mapping rather than `SampleIDinHead` or per-file assignment when enabling `-htsBarcode BC`.
