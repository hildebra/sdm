# sdm — simple demultiplexer
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sdm/badges/downloads.svg)](https://anaconda.org/bioconda/sdm)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/sdm/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/sdm)


**sdm** sorts DNA sequencing reads into samples (demultiplexing) and removes reads that do not meet your quality criteria. It can also remove primers and barcodes, merge paired reads, and combine duplicate sequences (dereplication).

It accepts FASTA and FASTQ files, including gzip-compressed files, and supports single-end and paired-end data. The default build also reads SAM, BAM and CRAM; see [alignment input](docs/hts.md).

**3.60 beta:** audit bug fixes (see the changelog in `DNAconsts.h`), lower dereplication memory, mate tags (`/1` `/2` or Illumina `1:N`/`2:N`) that match the output file after orientation swaps, and a single `make test` entry point for all tests.

**3.52 beta:** retain selected qualities for exact single-read and `(R1, R2)` variants, count fragments once, and export synchronized paired subclusters. Optional reassignment preserves mate linkage and sample counts; paired nucleotide sidecars use `SDMDIFF3`. See the [dereplication flags](docs/options.md#streaming-coarse-dereplication) and [paired-mode guide](docs/paired-diff-quality.md).

## Build from source

The supplied Makefile uses a C++20 compiler, GNU Make, POSIX threads, ISA-L, zstd and HTSlib. By default, `make` produces **one fully static executable, `sdm`, with HTSlib enabled**. It reads FASTA/FASTQ (plain, gzip/BGZF or zstd), SAM, BAM and CRAM.

| Command | Linking | Input support |
| --- | --- | --- |
| `make` | Static | FASTA/FASTQ plus SAM/BAM/CRAM |
| `make STATIC=0` | Shared libraries | FASTA/FASTQ plus SAM/BAM/CRAM |
| `make WITH_HTS=0` | Static | FASTA/FASTQ only |
| `make WITH_HTS=0 STATIC=0` | Shared libraries | FASTA/FASTQ only |

Every configuration writes `sdm`; switching configurations replaces that executable. `make hts` and `make portable` are aliases for selecting `WITH_HTS=1` and `WITH_HTS=0`. The former `sdm-hts` executable and `HTS_STATIC` setting are no longer used.

The default build needs HTSlib headers, `libhts.a`, and static archives for its enabled dependencies as well as the compiler runtime and system libraries. `pkg-config --static --libs htslib` supplies the dependency flags. On Debian/Ubuntu, start with `build-essential libisal-dev libzstd-dev libhts-dev pkg-config`; additional development packages depend on how HTSlib was built (always zlib, and for example libbz2, liblzma, libdeflate or libcurl). sdm itself uses neither zlib nor libdeflate; a portable build (`WITH_HTS=0`) and its tests need neither. Some distributions do not supply every required static archive; build those dependencies from source or use `STATIC=0` in that case. Missing dependencies produce a build error with setup instructions.

ISA-L (igzip) reads every gzip input, BGZF or not, and compresses the BGZF of `.gz` output; zstd reads and writes `.zst`. A static build needs `libisal.a` and `libzstd.a`, and Debian/Ubuntu's `libisal-dev` has no `libisal.a`: `make deps` (`./build_deps.sh`) builds static ISA-L 2.32.1 and zstd 1.5.7 from their checksummed release tarballs into `.deps/` (needs cmake, curl or wget, and nasm 2.14.01 or later for x86-64, which it builds itself if missing), and `make` uses `.deps/` from then on. `make distclean` removes it; `make clean` keeps it. See [compressed input and output](docs/compression.md) for the formats and how they are read and written.

To build the default static `sdm` with micromamba (no root access needed), run `./build_static.sh -j8`. The first run creates the micromamba env `sdm-static` from conda-forge (`gcc_linux-64 gxx_linux-64 sysroot_linux-64 make cmake pkg-config curl nasm zlib bzip2 liblzma-devel liblzma-static`). It then builds HTSlib 1.24 from source into that env, because bioconda's `htslib` needs libcurl and OpenSSL, which conda does not provide as static archives, and static ISA-L and zstd into `<env>/sdm-deps` (`build_deps.sh`). This HTSlib uses zlib, bzip2 and liblzma (no libdeflate) and has no remote-URL (libcurl, S3, GCS) or plugin support. Later runs only rebuild sdm, with objects in `.build/conda-*`. Arguments are passed to make, e.g. `./build_static.sh -j8 test`. Set `SDM_ENV`, `SDM_MICROMAMBA`, `HTSLIB_VERSION`, `ISAL_VERSION` or `ZSTD_VERSION` to change the defaults; `micromamba env remove -y -n sdm-static` starts over. A shared build (`make STATIC=0`) in a conda env needs only the compilers plus `htslib isa-l zstd pkg-config`.

```sh
make -j2
./sdm -v
```

For a custom HTSlib installation, set `PKG_CONFIG_PATH` to its `lib/pkgconfig` directory, or pass `HTS_CFLAGS` and `HTS_LIBS`. When specifying `HTS_LIBS` for a static build, include HTSlib's transitive dependencies too. Use `LDFLAGS` for additional library search/runtime paths. The existing `isa_FLAGS` and `isa_FLAGS2` overrides remain available. Each configuration has separate objects and automatic header dependencies under `.build/`. Run `make clean` when changing compiler flags or the HTSlib installation within a configuration.

```sh
make test                  # all unit, regression and pipeline suites (HTS suites included)
make test WITH_HTS=0        # the same without SAM/BAM/CRAM suites
make test-quick             # skip the slow parity/seed-extension suites
```

`make test` builds the C++ test programs and runs every suite through `tests/run_tests.pl` (Perl 5 with core modules; no Python needed). It continues after a failure, prints `PASS`/`FAIL` per suite with the last lines of any failing output, and exits non-zero if anything failed. Full logs are written to `.build/<configuration>/test-logs/<suite>.log`. Run selected suites with `make test TEST_ARGS="test_regressions_cli TestAuditFixes"`, or list them with `perl tests/run_tests.pl --list`. Individual suites can also be run directly: `.build/portable-1/TestAuditFixes [test-name]` or `perl tests/test_regressions_cli.pl ./sdm [case-name ...]`.

GitHub Actions ([.github/workflows/tests.yml](.github/workflows/tests.yml)) runs `make test` on every push and pull request in five configurations: with HTSlib, portable, static portable (ISA-L and zstd from `make deps`), ASan + UBSan, and TSan. Failing jobs upload their test logs.

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
| `DEPS_PREFIX` | `.deps` | Where `make deps` builds static ISA-L and zstd; used when it holds `lib/libisal.a` (and `lib/libzstd.a`). |
| `ISAL_CFLAGS`, `ISAL_LIBS` | `DEPS_PREFIX`, else pkg-config `libisal`; fallback `-lisal` | ISA-L include and link flags. |
| `ZSTD_CFLAGS`, `ZSTD_LIBS` | `DEPS_PREFIX`, else pkg-config `libzstd`; fallback `-lzstd` | zstd include and link flags. |

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

`make all` and `make sdm` select the default executable target. `make hts` / `make portable` select the same target with HTS enabled/disabled. `make test` tests the selected configuration; `make test-hts` / `make test-portable` select the corresponding test configuration. `make clean` removes executables and `.build/`; `make distclean` also removes `.deps/`. Reapply custom build variables when invoking test targets. Do not build different configurations concurrently in the same checkout, since they write the same `sdm`.

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

In automatic and forced-prefix modes, `TruncateSequenceLength` caps **both R1 and merged search keys**. Differences beyond the cap cannot create duplicate truncated inputs for UPARSE. In ordinary dereplication, full merged output and full HQ mates remain available for later use.

`-derepSrchLen <n>` (or `derepSrchLen` in the options file) accepts a positive integer or `-1`; the command-line value takes precedence. An explicit positive length keeps capped, length-matched comparison in automatic mode. Forced prefix mode overrides this search length but still respects the truncation cap. Use `-derepPrefix 0 -derepSrchLen -1` for legacy full-length exact matching.

Prefix matching allows no substitutions, gaps, or offsets. Longer sequences that disagree stay separate; a short read compatible with several longer sequences is assigned once, deterministically. See [dereplication options](docs/options.md#dereplication) for matching and count-recovery details, and [runtime and memory measurements](docs/derep-prefix-evaluation.md) for the tradeoffs.

For noisy reads, **streaming coarse dereplication** shares sequence references across nearby reads while preserving exact R1 dereplicates in standard outputs. Set `-derepIdentity` below `100` to enable it:

```sh
sdm -i_fastq ont.fastq -o_dereplicate clusters.fna -options sdm_options.txt \
    -derepIdentity 97.5 -threads 1
```

This example uses a 97.5% substitution identity threshold; choose the threshold for the compression and variant resolution you need. Each exact effective R1 key keeps its high-quality representative and counts, chosen by the existing dereplication rules. Coarse groups supply shared immutable search/compression references. Matching uses an immutable initial sequence, so replacing the output representative does not move the cluster boundary. Comparisons start at position zero, cover the shorter clustering key, and respect explicit search and truncation caps. There is no indel alignment or reverse-complement search.

With identity below 100 and multiple workers, coarse dereplication uses **64 independent partitions**. For default compact storage, an eligible positioned 20-mer from the effective search key determines routing (whole-key hashing is the fallback); identical keys reuse one dereplicate even when R2 differs, and the existing better-seed selection can replace both representative mates. Exact-key reuse avoids repeating coarse candidate searches. A shared work pool schedules ready partitions across the existing workers. Partitions run concurrently without a final merge. Internal coarse boundaries may differ from global grouping, but default exact R1 counts and minimum-copy decisions do not. Exact retained variants and their selected qualities remain available. Use `-derepGlobal 1` when one global clustering is required; one worker, identity 100, and `-derepReassign 1` also retain the global path. The startup log reports which mode is active.

Search uses non-overlapping 20-mers and excludes seeds containing homopolymers of four or more bases, or ambiguous bases. `-derepSeedHomopolymer <2..20>` changes the run cutoff. Homopolymer positions still count in the final identity calculation. Reads without enough usable seeds use direct comparison.

Every dereplicate stores its sparse sample counts in a base `diffDNA` object, including ordinary 100% exact/prefix dereplication. With nucleotide storage disabled, its reference and variant entries remain unallocated. The same counts supply abundance headers, sample maps, and minimum-copy filtering.

With only `derepIdentity=97` enabled, ordinary outputs are the contract: retain the configured merged/R1 search policy, exact keys, prefix recovery, admission, copy cutoffs and sample counts. Main FASTQ uses ordinary per-base quality accumulation; standard HQ files contain one full representative pair with its observed qualities per ordinary parent. Representative qualities do not require `derepStoreQuals=1`, and the normal seed command needs no `seedSubclusters`. Explicit full-variant retention is a separate policy. See the [storage-only parity and memory assessment](docs/storage-only-dereplication.md).

By default, coarse groups are internal candidate indexes and immutable compression references. The main FASTA/FASTQ, merged/rest files, `.map` and representative HQ pairs preserve ordinary exact search-key dereplicates, search/truncation caps, prefix recovery and copy cutoffs. Storage-only mode preserves the configured merged/R1 search policy and ordinary merge-offset subgroup selection. With explicit full-variant quality retention, the search source is currently R1 and R2 differences do not split that parent; the existing paired better-seed selection can replace both mates. Select `-derepCoarseClusters 1` only when the output should describe coarse clusters instead. `-derepReassign 1` also explicitly selects coarse-cluster output and performs final reassignment. `-derepGlobal 1` controls search scope, not output granularity.

Default 97% storage packs exact representatives instead of keeping a full DNA object for each key. Each coarse group shares immutable nucleotide references. Full unmerged fragment variants and their occurrence counts are delta encoded; each ordinary merge-offset subgroup keeps one selected pair's observed quality vectors and, for main FASTQ, lossless quality sums/coverage. Quality arrays choose raw, run-length or base-offset bit packing without changing scores. Better seeds replace the selected pair through the ordinary rules. Export reconstructs one exact parent at a time, keeping the same ordinary seed-candidate granularity. Memory still depends on distinct keys, variant diversity, qualities and index overhead; see the [peak-RSS benchmark](docs/storage-only-dereplication.md#peak-memory-measured-separately).

Add `-derepStoreDiffs 1` to retain nucleotide-only observations in a compact **`clusters.diff`** sidecar. This enables the nucleotide portion of `diffDNA`, which records each distinct nucleotide-change pattern and processed length once, with its exact copy count, relative to the fixed initial sequence. It discards individual read names and qualities. High-quality representative replacement leaves the reference and all stored differences intact. Memory grows with distinct exact R1 representatives, coarse references, and sample counts rather than every repeated observation. This retains more representative pairs than earlier coarse-parent-only builds. With differences enabled, memory also grows with the distinct encoded nucleotide variants. Exact repeats of any variant only increment its count.

To export every exact subcluster within the coarse groups, also add `-derepSubclusterFasta 1`. SDM then writes `clusters.subclusters.fna`, with headers such as `>read_17.sub2;size=42;`. Each record reconstructs one distinct full processed sequence and uses that sequence's own count; the normal output contains one representative per exact effective R1 key by default. The export is off by default and requires `-derepStoreDiffs 1`.

Add **`-derepStoreQuals 1`** to retain each exact variant's selected quality vector and export `clusters.1.hq.fq`. The existing representative-selection rules choose one complete observed quality vector for each exact sequence; these are not averaged qualities. Nucleotides remain encoded as differences from the immutable initial reference. Names are generated as `@read_17.sub2;size=42;`, with the exact variant's own count. The main exact R1 representative output and `.map` remain available. Retained variants use the normal HQ FASTQ filenames, with no duplicate subcluster FASTQ files; seed extension reads them with `-seedSubclusters 1`.

```sh
sdm -i_fastq ont.fastq -o_dereplicate clusters.fna -options sdm_options.txt \
    -derepIdentity 97 -derepStoreQuals 1 -threads 1
```

Quality retention defaults to `0` and requires complete qualities for unmerged single-end reads or complete read pairs. It automatically retains nucleotide variants in memory; `-derepStoreDiffs 1` separately enables binary delta output. Binary `.diff` output remains nucleotide-only, using `SDMDIFF2` for single reads and `SDMDIFF3` for linked pairs. Selected qualities and counts also follow variants through optional final reassignment. See [single-end implementation and ONT simulation benchmarks](docs/diff-quality-retention.md).

For paired input, the same flag retains **exact `(R1, R2)` variants**. Pairs sharing R1 but differing in R2 remain separate exact subclusters. SDM selects both quality vectors from one observed fragment and counts each fragment once. It writes synchronized `clusters.1.hq.fq` and `clusters.2.hq.fq` files with matching IDs and `;size=N;` values.

```sh
sdm -i_fastq reads.1.fq,reads.2.fq -paired 2 \
    -o_dereplicate clusters.fna -options sdm_options.txt \
    -derepIdentity 97 -derepStoreQuals 1 -threads 1
```

An identical effective R1 key reuses its parent even if R2 differs. For a new R1 key, R1 supplies candidate seeds; both mates must independently meet the configured same-start identity threshold against their corresponding initial references. Search caps apply to each mate independently. Final reassignment compares both mates and moves complete variants and sample counts together. Native merging options must be disabled for this mode. Default parent admission uses R1 eligibility, while variant quality-witness selection separately prefers both mates passing. HQ variants preserve full available mates and observed qualities, including logically trimmed tails, after physical technical cuts. Physically empty mates remain excluded; enabling retention does not preserve ordinary merged search/output. See the [four-mode behavioral audit and required changes](docs/derep-behavior-parity.md). The [Apong count investigation](docs/apong-derep-count-debug.md) reconciles the newer logs and reproduces merged/R1 differences with a 200 bp search cap. Optional paired `.diff` files use the linked `SDMDIFF3` format; paired exact-subcluster FASTA produces `.subclusters.1.fna` and `.subclusters.2.fna`. See [paired behavior and validation](docs/paired-diff-quality.md) and the [LotuS3 IO integration instructions](docs/lotus3-dereplication-io.md).

For explicitly requested coarse-cluster output with final reassignment, add `-derepReassign 1` (default `0`). SDM compares each exact variant against the final high-quality representatives and moves all its copies and sample counts to a strictly closer representative that meets `derepIdentity`. Equal identities keep the original cluster. This is one pass with fixed representatives; it does not select new representatives or iterate. It works at any configured identity, including 97%, 98%, and 99%.

Reassignment automatically retains nucleotide variants and per-variant sample counts in memory, even with `-derepStoreDiffs 0`, so it uses more memory than ordinary coarse clustering. Each variant stores one sparse sample-count vector; its abundance and the cluster sample map are derived from those entries. Writing `.diff` or exact-subcluster FASTA remains separately controlled by the existing options. Output counts, sample maps, minimum-copy filtering, and enabled variant exports all use the final assignments. See the [Q30 200 bp benchmark](docs/coarse-dereplication-q30.md) for speed and memory measurements.

All of these options can be set as tab-separated entries in `sdm_options.txt`; command-line values take precedence. Default `derepIdentity=100`, `derepStoreDiffs=0`, `derepStoreQuals=0`, and `derepReassign=0` retain the existing automatic exact/prefix modes. Enabling differences or quality retention at `100` selects streaming storage with exact R1 output keys; reassignment selects coarse-cluster output. See [coarse dereplication and the delta format](docs/coarse-dereplication.md) for semantics, limits, and validation.

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

For single-end ONT amplicons, see [ONT matching and LotuS integration](docs/ont.md) and the [primer/barcode search audit](docs/search-audit.md). ONT matching requires `-ontMode 1`; the default remains off.

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
- [ISA-L BSD 3-clause license](licenses/isa-l-LICENSE) and [zstd BSD license](licenses/zstd-LICENSE), linked by every build.
- [Versions and notice provenance](licenses/README.md).

These copies come from HTSlib 1.24 and its included htscodecs 1.6.7. The Makefile uses the installed HTSlib selected by your build flags; include matching notices if you distribute a different release, plus required notices/materials for the actual linked compression, networking and runtime libraries.

For binary distribution under SDM's GPL license, include the license and provide the corresponding source and build scripts as required by [GPLv3 section 6](https://www.gnu.org/licenses/gpl-3.0.html#section6), including non-system dependencies. A matching source archive alongside the binary is a practical release layout. Other linked libraries retain their own license obligations. A fully static binary needs no accompanying `libhts.a` or `libhts.so` at runtime.

Report bugs and contribute through the [sdm GitHub repository](https://github.com/hildebra/sdm). For questions, contact Falk Hildebrand at `falk.hildebrand@gmail.com`.
