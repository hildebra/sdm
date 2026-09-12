# sdm — simple demultiplexer
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sdm/badges/downloads.svg)](https://anaconda.org/bioconda/sdm)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sdm/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/sdm)


**sdm** sorts DNA sequencing reads into samples (demultiplexing) and removes reads that do not meet your quality criteria. It can also remove primers and barcodes, merge paired reads, and combine duplicate sequences (dereplication).

It accepts FASTA and FASTQ files, including gzip-compressed files, and supports single-end and paired-end data. The default build also reads SAM, BAM and CRAM; see [alignment input](docs/hts.md).

## Build from source

The supplied Makefile uses a C++20 compiler, GNU Make, POSIX threads, zlib and HTSlib. By default, `make` produces **one fully static executable, `sdm`, with HTSlib enabled**. It reads FASTA/FASTQ (including gzip), SAM, BAM and CRAM.

| Command | Linking | Input support |
| --- | --- | --- |
| `make` | Static | FASTA/FASTQ plus SAM/BAM/CRAM |
| `make STATIC=0` | Shared libraries | FASTA/FASTQ plus SAM/BAM/CRAM |
| `make WITH_HTS=0` | Static | FASTA/FASTQ only |
| `make WITH_HTS=0 STATIC=0` | Shared libraries | FASTA/FASTQ only |

Every configuration writes `sdm`; switching configurations replaces that executable. `make hts` and `make portable` are aliases for selecting `WITH_HTS=1` and `WITH_HTS=0`. The former `sdm-hts` executable and `HTS_STATIC` setting are no longer used.

The default build needs HTSlib headers, `libhts.a`, and static archives for its enabled dependencies as well as the compiler runtime and system libraries. `pkg-config --static --libs htslib` supplies the dependency flags. On Debian/Ubuntu, start with `build-essential zlib1g-dev libhts-dev pkg-config`; additional development packages depend on how HTSlib was built (for example libbz2, liblzma, libdeflate or libcurl). Some distributions do not supply every required static archive; build those dependencies from source or use `STATIC=0` in that case. Missing dependencies produce a build error with setup instructions.

```sh
make -j2
./sdm -v
```

For a custom HTSlib installation, set `PKG_CONFIG_PATH` to its `lib/pkgconfig` directory, or pass `HTS_CFLAGS` and `HTS_LIBS`. When specifying `HTS_LIBS` for a static build, include HTSlib's transitive dependencies too. Use `LDFLAGS` for additional library search/runtime paths. The existing `isa_FLAGS` and `isa_FLAGS2` overrides remain available. Each configuration has separate objects and automatic header dependencies under `.build/`. Run `make clean` when changing compiler flags or the HTSlib installation within a configuration.

```sh
make test                  # helper, regression and HTS input/pipeline tests
make test WITH_HTS=0        # helper and FASTA/FASTQ regression tests
```

### Make options

| Variable | Default | Purpose |
| --- | --- | --- |
| `WITH_HTS` | `1` | Enable HTSlib (`1`); use `0` to build without SAM/BAM/CRAM input. |
| `STATIC` | `1` | Fully static linking (`1`); use `0` for shared linking. Applies to `sdm` and test executables. |
| `CXX`, `CC` | GNU Make/compiler defaults | C++ and C compilers; C++20 support is required. |
| `CXXFLAGS` | `-O3` | C++ optimization/compiler flags; C++20 and thread flags are added separately. |
| `CPPFLAGS` | `-I. -D__STDC_CONSTANT_MACROS` | Preprocessor/include flags. |
| `CFLAGS` | Empty | C compiler flags. |
| `PKG_CONFIG` | `pkg-config` | Program used to query HTSlib compiler and link flags. |
| `PKG_CONFIG_PATH` | Environment | Directories containing custom `.pc` files, including `htslib.pc`. |
| `HTS_CFLAGS` | From pkg-config | HTSlib include/compiler flags. |
| `HTS_LIBS` | From pkg-config; fallback `-lhts` | HTSlib link flags; static linking needs all transitive dependencies. |
| `LDFLAGS`, `LDLIBS` | Empty | Additional linker flags and libraries. |
| `isa_FLAGS`, `isa_FLAGS2` | Empty | Existing CPU/compiler and linker overrides. |

For a custom installation with complete pkg-config metadata:

```sh
PKG_CONFIG_PATH=/opt/htslib/lib/pkgconfig make -j2
```

Alternatively, pass explicit flags. This example assumes HTSlib was built with zlib, bzip2 and liblzma, with no additional dependencies:

```sh
make -j2 HTS_CFLAGS='-I/opt/htslib/include' \
  HTS_LIBS='-L/opt/htslib/lib -lhts -llzma -lbz2 -lz -lm'
```

Use the libraries required by your own HTSlib configuration. Shared builds with a custom library directory may also need a runtime search path in `LDFLAGS`. No HTSlib source is downloaded or copied by `make`.

`make all` and `make sdm` select the default executable target. `make hts` / `make portable` select the same target with HTS enabled/disabled. `make test` tests the selected configuration; `make test-hts` / `make test-portable` select the corresponding test configuration. `make clean` removes executables and `.build/`; `make distclean` is an alias. Reapply custom build variables when invoking test targets. Do not build different configurations concurrently in the same checkout, since they write the same `sdm`.

A fully static executable does **not** need `libhts.a`, `libhts.so`, or an HTSlib installation on the recipient's machine. Build it for the recipient's operating system and CPU. Externally referenced CRAM files can still need their matching reference FASTA and indexes; see [CRAM references](docs/hts.md#cram-references).

See [license and redistribution information](#license-and-support) for SDM and the linked HTS dependencies. The Makefile builds executables and tests; it does not assemble a redistribution package.

Examples below assume `sdm` is on your `PATH`; otherwise, use `./sdm`. For Windows, use a Linux environment such as WSL with the same build dependencies; this checkout does not include a Visual Studio project.

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

Other layouts are supported: `SampleIDinHead` identifies samples from text in read headers, and `fastqFile`, `alignmentFile` or `fnaFile` identifies samples by input filename. See the [mapping-file reference](docs/mapping.md) for all columns and examples, or run `sdm -help_map`.

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

### Read SAM, BAM or CRAM files listed in a mapping file

Use the fixed `alignmentFile` header for all three formats, and omit the `fastqFile` column when it is not needed:

```text
#SampleID	alignmentFile
sample_A	sample_A.sam
sample_B	sample_B.bam
sample_C	sample_C.cram
```

```sh
sdm -i_path reads -map alignments.tsv -o_fastq filtered.fastq -options sdm_options.txt
```

For mixed FASTQ/alignment inputs, both columns may be present, but each row must fill exactly one: leave `fastqFile` empty for alignment samples and `alignmentFile` empty for FASTQ samples. Preserve tabs for empty cells. All inputs must use a compatible singleton/paired layout; paired alignment files must be query-name sorted. For paired output, use `-paired 2 -o_fastq filtered.1.fq,filtered.2.fq`. Add `-cramRef reference.fa` when an external CRAM reference is needed. See [mapping examples](docs/mapping.md#sam-bam-and-cram-file-columns) and [alignment input](docs/hts.md) for barcode modes and reference requirements.

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

Matching is selected at startup using `derepPrefix` (default: `auto`). Set it with `-derepPrefix auto|0|1` on the command line, or with a tab-separated `derepPrefix` entry in the options file. The command-line value takes precedence.

| Value | Matching behavior |
| --- | --- |
| `auto` | Use 100% prefix matching for variable clustering lengths, or the fast exact-key path when length settings guarantee a fixed length. An explicitly positive `derepSrchLen` also selects length-matched keys. |
| `1` | Force 100% prefix matching from the same start across the shorter clustering sequence. |
| `0` | Use legacy length-matched keys and search-length behavior. |

For fixed 200 bp clustering, these options keep the fast exact-key path:

```text
minSeqLength	200
TruncateSequenceLength	200
derepPrefix	auto
```

In automatic and forced-prefix modes, `TruncateSequenceLength` caps **both R1 and merged search keys**. Differences beyond the cap cannot create duplicate truncated inputs for UPARSE. Full merged output and full HQ mates remain available for later use.

`-derepSrchLen <n>` (or `derepSrchLen` in the options file) accepts a positive integer or `-1`; the command-line value takes precedence. An explicit positive length keeps capped, length-matched comparison in automatic mode. Forced prefix mode overrides this search length but still respects the truncation cap. Use `-derepPrefix 0 -derepSrchLen -1` for legacy full-length exact matching.

Prefix matching allows no substitutions, gaps, or offsets. Longer sequences that disagree stay separate; a short read compatible with several longer sequences is assigned once, deterministically. See [dereplication options](docs/options.md#dereplication) for matching and count-recovery details, and [runtime and memory measurements](docs/derep-prefix-evaluation.md) for the tradeoffs.

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
- **Length caps apply to clustering keys too.** `TruncateSequenceLength` caps ordinary reads and, in automatic or forced-prefix dereplication mode, merged search keys. Full merged output and seed-extension output retain their full length. A shorter cap can also lower the effective minimum length for ordinary reads.
- **More threads can affect dereplication results.** `-threads 4` enables four workers. For repeatable dereplication runs, use `-threads 1`; with multiple workers, representative sequences, some abundance counts, and output order can vary.

## Documentation

- [Filtering-option file](docs/options.md): format, defaults, trimming, primer handling, and a second filtering tier.
- [Sample mapping file](docs/mapping.md): columns, barcode layouts, primers, and input files.
- [SAM/BAM/CRAM input](docs/hts.md): pairing, CRAM references, and barcode tags.
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

sdm is free software licensed under the [GNU General Public License, version 3 or later](LICENSE).

The default build links HTSlib, copyright Genome Research Limited. Its code outside `cram/` uses the MIT/Expat license; `cram/` uses the modified 3-clause BSD license. Its included htscodecs library uses BSD terms, with public-domain/CC0 exceptions identified in its notice. Full upstream texts are included here:

- [HTSlib MIT/Expat and modified BSD licenses](licenses/HTSlib-LICENSE).
- [htscodecs license and public-domain notices](licenses/htscodecs-LICENSE.md).
- [Versions and notice provenance](licenses/README.md).

These copies come from HTSlib 1.24 and its included htscodecs 1.6.7. The Makefile uses the installed HTSlib selected by your build flags; include matching notices if you distribute a different release, plus required notices/materials for the actual linked compression, networking and runtime libraries.

For binary distribution under SDM's GPL license, include the license and provide the corresponding source and build scripts as required by [GPLv3 section 6](https://www.gnu.org/licenses/gpl-3.0.html#section6), including non-system dependencies. A matching source archive alongside the binary is a practical release layout. Other linked libraries retain their own license obligations. A fully static binary needs no accompanying `libhts.a` or `libhts.so` at runtime.

Report bugs and contribute through the [sdm GitHub repository](https://github.com/hildebra/sdm). For questions, contact Falk Hildebrand at `falk.hildebrand@gmail.com`.
