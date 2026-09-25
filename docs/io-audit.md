# I/O audit (sdm 3.60 beta, 2026-09-25)

Read-only review of the input and output layer: `InputFileStream.cpp`, `OutputFileStream.cpp`,
`InputStream.{h,cpp}` (`ifbufstream`, `InputStreamer`), `OutputStreamer.{h,cpp}` (`ofbufstream`,
`AsyncOutputWriter`, retention, demultiplexed output), `HtsInput.cpp`, the reader loop in `IO.cpp`
and the helpers in `Common.cpp`. Nothing was changed. The items marked *reproduced* were run
against a `WITH_HTS=0 STATIC=0` build in WSL with the scripts in `.claude/archive/io-audit/`.

## Status

Fixed on 2026-09-25 (same day, verified with `.claude/archive/io-audit/exp_fix.sh` and the
new `TestAuditFixes` cases `fastq_offset_autodetect` and `output_finalize_failure_fatal`):

- Finding 1: `InputStreamer::auto_fq_version` now accumulates the ASCII range over records.
  Any value above 105 or below 59 settles on Phred+33 at once; the range 59..105 stays
  undecided until the probe block is exhausted (`finalizeFastqVersion`, called from
  `read_sequences` and before immediate conversions), and only then does the legacy
  "above 74 means +64" rule apply, with a message naming the decision.
- Finding 2: `closeOutStreams_locked` and `closeOutFilesDemulti_locked` finalize every
  stream explicitly and rethrow the first failure after closing all of them; the subset and
  seed-writing paths in `IO.cpp` close explicitly before deleting their `OutputStreamer`.
  `-o_demultiplex` to a full disk now exits 1 with `Failed to finalize output file: <file>`.
- Finding 3: the `.gzi` index is opt-in with `-o_gzi 1` (`setBgzfIndexOutput`); by default
  `finalize()` removes a stale index next to the file instead of writing one.

Finding 4 was fixed on 2026-09-25 as part of the general audit: `applyFileIT` keeps `.gz`
last and ignores dots in directory names (`Common.cpp`, test `apply_file_it_naming`).
Findings 5 to 9 are still open.

## What is in good shape

- BGZF input: the `.gzi` index is validated block by block against the file, every mismatch
  falls back to sequential gzip with a warning, a missing EOF marker warns and continues, and a
  truncated gzip member returns the recovered data (`InputFileStream.cpp:334-410`, `:105-133`).
- BGZF output: parallel block compression keeps block order, delayed write failures are
  re-raised from `finalize()`, appends to an existing BGZF file rebuild the index by seeking
  over payloads (`OutputFileStream.cpp:250-291`, `:417-511`).
- Record synchronization: R1/R2/MID ids are compared per record on the reader thread, a
  stream ending early is an error, and `-ignore_IO_errors 1` drops the whole record set instead
  of desynchronizing mates (`InputStream.cpp:592-657`).
- Output ordering: one writer thread receives buffers in append order per stream, appenders
  block only on a 16 MB global cap, and per-worker batches publish R1/R2 together
  (`OutputStreamer.cpp:14-75`, `:444-584`).
- `TestGzipInput`, `TestBgzfOutput`, `TestIfBufStream` and `TestHtsInput` cover the above.

## Findings

### 1. Default FASTQ offset detection turns Phred+33 data with Q > 41 into Phred+64 (reproduced)

`InputStreamer::getDNAlines` runs `auto_fq_version` on the **first record only**
(`InputStream.cpp:726-729`); the decision is `min ASCII < 59 -> 33, else max ASCII > 74 -> 64`
(`:1128-1134`). Nothing re-checks later records: `minmaxQscore`, which could switch back to
+33 on impossible values, is never called (`DNA.cpp:156` has it commented out). Illumina
headers containing ` 1:N:` are exempt (`:1147-1156`); anything else is not.

A PacBio HiFi read (or any read) whose first record has no base below Q26 and any base above
Q41 (ASCII `K`..`~`) therefore selects +64 for the whole file. Every quality is then 31 too
low, which changes filtering, trimming, merging and dereplication seed selection silently.
This is the default (`fastqVersion auto`, `-i_qual_offset` unset).

Reproduced: a two-read FASTQ whose first read has qualities `~`x30 `I`x30 and second `5`x60
is written back with `-o_qual_offset 33` as `_`x30 `*`x30 and `!`x60; `-i_qual_offset 33`
(without the option-file override) writes it back unchanged.

Suggested rule: any ASCII value above 105 (`i`, Phred+64 Q41) is impossible in +64, so it
should force 33 immediately; and the decision should use the whole probe block that
`read_sequences` already reads (`IO.cpp:482-519`, 120 records), not one record.

### 2. Write failures on demultiplexed (and other non-retained) outputs do not fail the run (reproduced)

Demultiplexed streams are closed by dropping the `shared_ptr` (`OutputStreamer.cpp:897-904`,
`:912-919`); `~ofbufstream` calls `deactivate()` inside `try/catch` and only prints
`Failed to finalize <file>` (`:92-97`). The same applies to `dualOfBufStream`'s destructor
(`:354-363`) and to the main outputs whenever retention is off (`-maxReadsPerOutput`, read
subsets, seed extension), because `closeOutStreams_locked` only calls `emptyStream()` before
`reset()` (`:1597-1647`). Only retained streams go through `closeAll()`, which rethrows.

Reproduced: with `dm/S1.fq` symlinked to `/dev/full`, `-o_demultiplex dm` prints the message
and exits 0; `-o_fastq /dev/full` exits 1 with `Error: Failed to finalize output file`.
A pipeline consuming the per-sample files would continue on truncated data.

Suggested fix: in `closeOutFilesDemulti_locked` and `closeOutStreams_locked(true)` call
`deactivate()` explicitly on every stream, remember the first exception, and rethrow it after
all streams are closed (the pattern `closeAll()` already uses).

### 3. Every gzip output gets a `.gzi` sidecar, including per-sample demultiplexed files (reproduced)

`ParallelBgzfOstream::finalize()` always writes `<file>.gzi` (`OutputFileStream.cpp:283`,
`:513-530`); `-o_demultiplex_gz 1` therefore leaves `S1.fq.gz` **and** `S1.fq.gz.gzi` (and
`.1.fq.gz.gzi`, `.2.fq.gz.gzi`, `.merg.fq.gz.gzi` for paired input) in the demultiplex
directory. Downstream code that lists that directory with an unanchored pattern such as
`fq.gz` (R `list.files`, shell globs like `*fq.gz*`) will pick the index files up; the code
itself notes that these files feed dada2 (`OutputStreamer.cpp:1065`). The index only helps
sdm re-read its own output in parallel. Consider not writing it for demultiplexed outputs, or
making it opt-in (`-o_gzi 0|1`), and documenting it in `docs/output.md`.

Related: a crash before `finalize()` leaves the previous run's `.gzi` next to a new `.gz`.
The reader detects the mismatch and falls back, but removing a stale `.gzi` when the file is
opened for truncation (`OutputFileStream.cpp:234-240`) would avoid the warning and the wasted
validation pass.

### 4. `.gz` outputs without an inner extension lose compression when rotated or suffixed (reproduced)

`applyFileIT` looks for the dot before `.gz` (`Common.cpp:263-281`); when there is none
(`-o_fastq rot.gz`) it appends the iteration after the suffix. `-maxReadsPerOutput 1` then
produces `rot.gz` (BGZF) and `rot.gz.1` (plain text), and the singleton files become
`rot.gz.singl`-style names that are also written uncompressed. `rot2.fq.gz` rotates
correctly to `rot2.1.fq.gz`. Fix: when `isGZfile(x)` and there is no inner dot, insert the
iteration/suffix before the `.gz`.

### 5. BGZF write errors lose the file name

`ParallelBgzfOstream` sets `badbit|failbit` exceptions on its `std::ofstream`
(`OutputFileStream.cpp:239`), so an `ENOSPC` surfaces as `Error: basic_ios::clear: iostream
error` (reproduced with `-o_fastq full.fq.gz -> /dev/full`). `writeBlock`/`finalize` should
wrap the failure with the path, as the plain-file path does (`:594-609`).

### 6. Per-stream memory of the asynchronous writer scales with the number of output files

Each `ofbufstream` keeps its buffer plus up to two spares once buffers have been handed off
(`OutputStreamer.cpp:155-193`), and each BGZF stream adds a 64 KB block buffer. For paired
demultiplexed output that is about 3 x (40 + 60 + 50) KB + 3 x 64 KB ~= 640 KB per sample,
roughly double the pre-async footprint; 1000 samples ~= 640 MB. Main outputs use 5 MB
buffers, so the four singleton and the merged stream, which are always opened (the
`onlyPrep` flag never prevents activation: `:78-90` vs `:1733`), cost 25 MB even when empty.
Within the accepted constant-overhead budget for typical sample counts, but worth knowing.

### 7. `ParallelBgzfOstream` passes `this` to `std::ostream` before its `std::streambuf` base exists

`class ParallelBgzfOstream : public std::ostream, private std::streambuf` with
`std::ostream(this)` (`OutputFileStream.cpp:213-222`) converts `this` to a base class whose
construction has not started, which [class.cdtor]/3 leaves undefined. GCC/MSVC accept it.
Declaring the `streambuf` base first (or holding it as a member and calling `rdbuf()` in the
body) removes the question.

### 8. HTS paired input rejects collate-grouped BAM

`loadSet` requires query names to be strictly increasing, lexically or naturally
(`HtsInput.cpp:195-198`), so `samtools collate` output (mates adjacent, names unordered) is
refused although adjacency is all the reader needs. The check exists to catch a name that
reappears later; if that guarantee is worth keeping, a documented `-htsUnsorted` style opt-out
or a bounded set of recently seen names would let grouped input through. Documented behavior,
so a design note rather than a bug.

### 9. Dead or vestigial code in the input layer

- `InputStreamer::read_fastq_entry_fast` (`InputStream.cpp:421-454`) is unreferenced and
  writes past `Iqual` when a quality line is longer than the sequence.
- `InputStreamer::read_fasta_entry` and `jumpToNextDNA` (`:353-403`, `:523-533`) are
  unreferenced; `read_fasta_entry` never fills `tseq`/`tqual`, so it would produce empty reads.
- `read_fastq_entry`, `jmp_fastq`, `_measure`, `_drawbar`, `_print` are declared or defined
  but unused (`InputStream.h:417-441`).
- `ifbufstream` takes `Instr::TimedSharedLockGuard` on a per-instance `input_mtx` in its
  constructor, destructor and `reset()` (`InputStream.h:119-133`, `:138-143`, `:150`); nothing
  else uses that mutex. The `if (!primary)` branch after `openInputFile` is unreachable, since
  `openInputFile` throws instead of returning null, so the "Could not open or empty" messages
  in `setupFastq_2`/`setupFastaQual2`/`openMIDseqs` only fire for empty files.
- `prepareAppend` accepts only BGZF headers with `XLEN == 6` (`OutputFileStream.cpp:451`),
  while the reader parses arbitrary extra subfields (`InputFileStream.cpp:291-313`); BGZF
  written by other tools with extra subfields cannot be appended to.
- `printBuildFeatures` reports `Gzip support: enabled (zlib)` in libdeflate builds
  (`IO.cpp:1450`).
- 27 `exit()` calls remain in the six I/O source files. `main` deliberately uses `_Exit` for
  exceptions so worker and writer threads are not torn down under static destructors
  (`Demultipl.cpp:160-162`); an `exit()` reached while workers run does not get that
  protection, and retained BGZF streams finalized during static destruction would find the
  gzip pool already stopped. Converting the remaining calls to exceptions (as most of
  `readCmdArgs` already does) closes that gap.

## Scope not covered

HTSlib paths were reviewed but not executed (no HTSlib in the WSL build). Dereplication output
files (`containers.cpp:1436-1500`) were only skimmed: they use plain `std::ofstream` with
`finalizeOutputFile`, so the delayed-failure handling matches the main outputs.
