# Advanced workflows and command-line options

[Back to the README](../README.md) · [Filtering options](options.md) · [Output files](output.md)

The README covers common starting points. This reference describes additional controls for input, merging, dereplication, and specialized workflows. All command-line options below take a value; run `sdm -help_flags` for the reference supplied with your installed version.

## Input and read handling

| Option | Purpose |
| --- | --- |
| `-i <file>` | Detect FASTA, FASTQ, SAM, BAM or CRAM automatically. |
| `-i_fastq <file>` | Supply FASTQ input explicitly. |
| `-i_sam <file>`, `-i_bam <file>`, `-i_cram <file>` | Supply alignment input; requires an HTS-enabled build. |
| `-cramRef <reference.fa>` | Matching indexed reference for CRAM; default: none explicitly supplied. |
| `-htsBarcode <inline or BC>` | Inline barcodes (default), or exact sample assignment from `BC:Z` tags. |
| `-i_fna <file>` | Supply FASTA input explicitly. |
| `-i_qual <file>` | Supply separate FASTA quality data; for `test.fna`, the inferred filename is `test.qual_`. |
| `-i_path <directory>` | Resolve input filenames from the [mapping file](mapping.md). |
| `-i_MID_fastq <file>` | Supply a separate barcode/MID read file. |
| `-paired <1, 2, or 3>` | Single-end, paired-end, or R1/R3 with MID reads in R2. Default: infer the input layout. |
| `-onlyPair <1 or 2>` | Process only the selected mate. |
| `-i_qual_offset <0, auto, or n>` | Input FASTQ quality offset; `0` or `auto` requests detection. |
| `-mergedPairs <0 or 1>` | Indicate that input reads were merged externally. |
| `-5PR1cut <n>`, `-5PR2cut <n>` | Remove leading bases from read 1 or read 2. |
| `-ignore_IO_errors <0 or 1>` | Attempt to continue after corrupt FASTQ records; default `0`. |

Choose one primary input mode: `-i`, `-i_fastq`, `-i_fna`, an explicit alignment input, or `-i_path`. Quality and MID files are companion inputs. Paired filenames are comma-separated; with `-paired 3`, list R1,R2,R3 in that order.

SAM/BAM/CRAM input uses one file per input set, containing singletons or complete query-name-sorted pairs. Separate MID FASTQ files cannot accompany alignment input. See [alignment input](hts.md) for orientation, qualities, references and error handling.

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

Merging uses overlap evidence after primer and barcode removal; the overlap must span at least `MergeMinOverlap` bases (default 16, exact match below 20 bases). If you choose to retain those sequences, they remain in the merged amplicon. The consensus is checked against the configured quality and sequence filters, including minimum length after quality trimming. Rejected merges can fall back to eligible mate output; exhausted merged-output base budgets do not trigger that fallback.

`TruncateSequenceLength` does not cap merged output. High-quality (HQ) mate files retain their full length so that later merging still has the required overlap. See [length settings](options.md#trimming-and-read-formatting).

## Dereplication

```sh
sdm -i_fastq reads.fastq -options sdm_options.txt -o_dereplicate dereplicated.fna -min_derep_copies 2
```

| Option | Purpose |
| --- | --- |
| `-o_dereplicate <file>` | Dereplicated sequence output. |
| `-derep_format <fa or fq>` | FASTA (default) or FASTQ output. |
| `-derep_map <file>` | Read the dereplication abundance map during seed extension; preprocessing writes it beside `-o_dereplicate`. |
| `-min_derep_copies <rule>` | Minimum abundance required to retain a sequence. |
| `-dere_size_fmt <0 or 1>` | `0` uses `size=X;` abundance headers; `1` uses `_X`. |
| `-derepSrchLen <n>` | Explicit search length; see the interaction with prefix mode below. Overrides the options-file value. |
| `-derepPrefix <auto or 0 or 1>` | Automatic (default), legacy exact keys, or forced 100% prefix matching. |
| `-derepPerSR <0 or 1>` | Enable dereplication by sequencing run, using the map's `SequencingRun` column. |
| `-derep_reserve <n>` | Reserve capacity for sequence groups. |

A bare number such as `8` retains sequences with at least 8 copies in total across all samples. A per-sample rule `X:Y` requires at least X copies in each of at least Y samples: `10:1,3:3` retains sequences with at least 10 copies in one sample **or** at least three copies in each of three samples. Separate rules with `,` (or `;`); a sequence passes if any rule holds. The mapping file also accepts a sample-specific `derepMin` value.

The default search length depends on filtering settings; see [dereplication options](options.md#dereplication). For quality averaging and snapshots, see [output behavior](output.md#dereplication-output-behavior).

## Seed extension and OTU output

Seed extension works with selected reads and read-to-cluster assignments. Relevant command-line inputs include:

| Option | Purpose |
| --- | --- |
| `-optimalRead2Cluster <file>`, `-optimalRead2Cluster_ref <file>` | Read-to-cluster assignments. |
| `-OTU_fallback <file>`, `-OTU_fallback_refclust <file>` | Fallback OTU FASTA sequences. |
| `-otu_matrix <file>` | OTU matrix output. |
| `-minQueryCov <fraction>` | Minimum seed/OTU query coverage; default `0.8`. |
| `-id <percent>` | Minimum seed/OTU identity; default `97`. |
| `-merge_pairs_seed 1` | Enable pair merging for seed extension. |
| `-seedSubclusters <0 or 1>` | Use reconstructed exact-variant FASTQ from `-derepStoreQuals 1`; default `0` reads one ordinary representative per parent. |

Seed extension uses already processed HQ mates without repeating quality trimming or 5′, barcode, or primer cuts. It ignores `TruncateSequenceLength` for all outputs, including unmerged mates and fallback seeds. Merged seeds still undergo minimum- and maximum-length checks.

With `-merge_pairs_seed 1 -o_fna OTU.fna`, merged seeds are written directly to `OTU.fna`, using the OTU identifier without a mate suffix. Corresponding qualities go to `-o_qual` when requested. A seed without an eligible merge uses the selected R1 or fallback sequence. With two FASTA output paths, merged seeds go to the first path; FASTQ seed output uses a separate merged file.

If you have HQ mate files generated with the truncating behavior in versions 3.44–3.45, regenerate them from the original reads before rerunning seed extension.

### Reconstructed exact variants

Use `-seedSubclusters 1` with the original parent assignments, parent `derep.map`, and fallback OTU FASTA. Retained variants now use the standard HQ filenames, replacing parent-only HQ records when `-derepStoreQuals 1` is enabled. Single-read input is `clusters.1.hq.fq`. There are no separate subcluster FASTQ files. See the [LotuS3 integration checklist](lotus3-dereplication-io.md). Paired input uses both synchronized exports:

```sh
sdm -i_fastq clusters.1.hq.fq,clusters.2.hq.fq -paired 2 \
    -seedSubclusters 1 -optimalRead2Cluster parent_assignments.uc \
    -derep_map clusters.map -OTU_fallback clustered_otus.fna \
    -options sdm_options.txt -merge_pairs_seed 1 \
    -o_fastq seeds.1.fq,seeds.2.fq -otu_matrix OTU.txt
```

The final generated `.subN` suffix resolves to the parent query; the original UC, PAF, UPARSE, or CD-HIT parser supplies its OTU assignment and identity. Every exact variant of a mapped parent competes through the existing seed selector in export order. Its `;size=N;` is a candidate multiplicity, not a sample assignment. A parent with variants of sizes two and one contributes its three parent counts once. Both mates of a fragment count together; identical R1 with different R2 remains separate candidates, with the corresponding quality vectors kept together.

Parent map rows and mapper assignments may be in a different order from the FASTQ exports. Complete groups are read even when below the preprocessing copy cutoff. Unassigned parents remain available to the existing additional-count and `.rest`/reference recovery paths; these paths retain their count-only behavior and fallback rules. Multi-target chimera assignments retain the existing split-count rules and do not supply seed candidates.

The reader validates synchronized mate IDs, abundances, record completeness, generated subcluster order, and each included parent's summed variant abundance. Keep each parent's generated `.sub1`, `.sub2`, … order intact and supply both mate files for paired exports. Quality trimming, preprocessing cuts, and truncation are not repeated, including with `-merge_pairs_seed 1`. Binary `.diff` output is optional and is not read by seed extension.

Run the preprocessing-to-seed regressions with:

```sh
python3 tests/test_seed_subclusters.py ./sdm
```


## Performance and repeatability

| Option | Default | Purpose |
| --- | --- | --- |
| `-threads <n>` | `1` | Number of processing workers. |
| `-threadIO <0 or 1>` | `1` | Enable threaded input/output. |
| `-iniBlockSize <n>` | `120` | Initial number of reads per processing block. |

Multithreading is included in current SDM builds. Default help and `-v` report threading, gzip and HTSlib availability. `-threads` controls processing workers; HTS decoding is currently sequential, with decoded reads processed by those workers. `-threadIO` controls the FASTA/FASTQ I/O backend independently.

Try increasing `-threads` for large datasets and compare elapsed time on your machine. `-derep_reserve` can reserve capacity for dereplication; it specifies the total reservation, not an amount per worker.

Multithreaded dereplication does not guarantee byte-identical results. Scheduling can affect representative selection, counts when a failed read arrives before the first passing read in a group, and the order of equally abundant records. Use `-threads 1` when repeatability is required.

## Specialized library processing

- `-GoldenAxe 1` enables processing of PacBio concatenated-read libraries. Supply forward and reverse primers for every sample in the mapping file.
- `-GoldenAxeMinAmpli <n>` and `-GoldenAxeMaxAmpli <n>` set accepted amplicon-count bounds; `-1` disables the corresponding bound.
- `-illuminaClip 1` enables Illumina adapter clipping. Set the relevant adapter sequences in the [options file](options.md#paired-reads-and-adapters).
