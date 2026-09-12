# Advanced workflows and command-line options

[Back to the README](../README.md) · [Filtering options](options.md) · [Output files](output.md)

The README covers common starting points. This reference describes additional controls for input, merging, dereplication, and specialized workflows. All command-line options below take a value; run `sdm -help_flags` for the reference supplied with your installed version.

## Input and read handling

| Option | Purpose |
| --- | --- |
| `-i <file>` | Detect FASTA or FASTQ automatically. |
| `-i_fastq <file>` | Supply FASTQ input explicitly. |
| `-i_fna <file>` | Supply FASTA input explicitly. |
| `-i_qual <file>` | Supply separate FASTA quality data; for `test.fna`, the inferred filename is `test.qual_`. |
| `-i_path <directory>` | Resolve input filenames from the [mapping file](mapping.md). |
| `-i_MID_fastq <file>` | Supply a separate barcode/MID read file. |
| `-paired <1, 2, or 3>` | Single-end, paired-end, or R1/R3 with MID reads in R2. Default: `1`. |
| `-onlyPair <1 or 2>` | Process only the selected mate. |
| `-i_qual_offset <0, auto, or n>` | Input FASTQ quality offset; `0` or `auto` requests detection. |
| `-mergedPairs <0 or 1>` | Indicate that input reads were merged externally. |
| `-5PR1cut <n>`, `-5PR2cut <n>` | Remove leading bases from read 1 or read 2. |
| `-ignore_IO_errors <0 or 1>` | Attempt to continue after corrupt FASTQ records; default `0`. |

Choose one primary input mode: `-i`, `-i_fastq`, `-i_fna`, or `-i_path`. Quality and MID files are companion inputs. Paired filenames are comma-separated; with `-paired 3`, list R1,R2,R3 in that order.

Paired and MID read identifiers are always checked for synchronized order. Error recovery does not disable these checks.

For named-read extraction, `-specificReads <file>` selects reads and can be combined with `-excludeFile <file>` in the extraction workflow. `-number T` with an optional `-prefix <text>` renumbers FASTA headers and exits. Use these as separate processing workflows rather than ordinary filtering thresholds.

## Pair merging

Enable merging for the output you need:

| Option | Workflow |
| --- | --- |
| `-merge_pairs_filter 1` | Ordinary filtered paired reads. |
| `-merge_pairs_demulti 1` | Per-sample demultiplexed output. |
| `-merge_pairs_derep 1` | Dereplication. |
| `-merge_pairs_seed 1` | Seed extension. |

For example, merge paired reads while demultiplexing:

```sh
mkdir -p demultiplexed
sdm -i_fastq R1.fastq,R2.fastq -paired 2 -map mapping.txt -options sdm_options.txt -o_demultiplex demultiplexed -merge_pairs_demulti 1
```

Merging uses overlap evidence after primer and barcode removal. If you choose to retain those sequences, they remain in the merged amplicon. The consensus is checked against the configured quality and sequence filters, including minimum length after quality trimming. Rejected merges can fall back to eligible mate output; exhausted merged-output base budgets do not trigger that fallback.

`TruncateSequenceLength` does not cap merged output. High-quality (HQ) mate files retain their full length so that later merging still has the required overlap. See [length settings](options.md#trimming-and-read-formatting).

## Dereplication

```sh
sdm -i_fastq reads.fastq -options sdm_options.txt -o_dereplicate dereplicated.fna -min_derep_copies 2
```

| Option | Purpose |
| --- | --- |
| `-o_dereplicate <file>` | Dereplicated sequence output. |
| `-derep_format <fa or fq>` | FASTA (default) or FASTQ output. |
| `-derep_map <file>` | Dereplication abundance-map output. |
| `-min_derep_copies <rule>` | Minimum abundance required to retain a sequence. |
| `-dere_size_fmt <0 or 1>` | `0` uses `size=X;` abundance headers; `1` uses `_X`. |
| `-derepSrchLen <n>` | Search length; `-1` uses the full sequence. Overrides the options-file value. |
| `-derepPerSR <0 or 1>` | Enable dereplication by sequencing run, using the map's `SequencingRun` column. |
| `-derep_reserve <n>` | Reserve capacity for sequence groups. |

A copy rule such as `10:1,3:3` retains sequences with at least 10 copies in one sample **or** at least three copies in each of three samples. The mapping file also accepts a sample-specific `derepMin` value.

The default search length depends on filtering settings; see [dereplication options](options.md#dereplication). For quality averaging and snapshots, see [output behavior](output.md#dereplication-output-behavior).

## Seed extension and OTU output

Seed extension works with selected reads and read-to-cluster assignments. Relevant command-line inputs include:

| Option | Purpose |
| --- | --- |
| `-optimalRead2Cluster <file>`, `-optimalRead2Cluster_ref <file>` | Read-to-cluster assignments. |
| `-OTU_fallback <file>`, `-OTU_fallback_refclust <file>` | Fallback OTU FASTA sequences. |
| `-otu_matrix <file>` | OTU matrix output. |
| `-merge_pairs_seed 1` | Enable pair merging for seed extension. |

Seed extension uses already processed HQ mates without repeating quality trimming or 5′, barcode, or primer cuts. It ignores `TruncateSequenceLength` for all outputs, including unmerged mates and fallback seeds. Merged seeds still undergo minimum- and maximum-length checks.

With `-merge_pairs_seed 1 -o_fna OTU.fna`, merged seeds are written directly to `OTU.fna`, using the OTU identifier without a mate suffix. Corresponding qualities go to `-o_qual` when requested. A seed without an eligible merge uses the selected R1 or fallback sequence. With two FASTA output paths, merged seeds go to the first path; FASTQ seed output uses a separate merged file.

If you have HQ mate files generated with the truncating behavior in versions 3.44–3.45, regenerate them from the original reads before rerunning seed extension.

## Performance and repeatability

| Option | Default | Purpose |
| --- | --- | --- |
| `-threads <n>` | `1` | Number of processing workers. |
| `-threadIO <0 or 1>` | `1` | Enable threaded input/output. |
| `-iniBlockSize <n>` | `120` | Initial number of reads per processing block. |

Try increasing `-threads` for large datasets and compare elapsed time on your machine. `-derep_reserve` can reserve capacity for dereplication; it specifies the total reservation, not an amount per worker.

Multithreaded dereplication does not guarantee byte-identical results. Scheduling can affect representative selection, counts when a failed read arrives before the first passing read in a group, and the order of equally abundant records. Use `-threads 1` when repeatability is required.

## Specialized library processing

- `-GoldenAxe 1` enables processing of PacBio concatenated-read libraries. Supply forward and reverse primers for every sample in the mapping file.
- `-GoldenAxeMinAmpli <n>` and `-GoldenAxeMaxAmpli <n>` set accepted amplicon-count bounds; `-1` disables the corresponding bound.
- `-illuminaClip 1` enables Illumina adapter clipping. Set the relevant adapter sequences in the [options file](options.md#paired-reads-and-adapters).
