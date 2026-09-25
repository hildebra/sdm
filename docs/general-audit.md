# General audit: filtering, seed extension, records, statistics (sdm 3.60 beta, 2026-09-25)

Third part of the audit, after the [I/O audit](io-audit.md) and the [dereplication audit](derep-audit.md).
Five reviewers (read-only) covered `Filters.cpp`/`Filters.h`, the seed-extension and OTU-matrix code
(`UClinks`, `ReadSubset`, `OptContainer`), `DNA.cpp`/`CompactDNA`/`diffDNA`, and `Statistics`,
`Benchmark.h`, the thread pools, `Common.cpp` and the output batching. Every finding listed here was
re-verified in the code before it was accepted; fixed items were changed the same day and are covered
by the full test suite (`make test`, 23 suites) plus the new `TestAuditFixes` cases.

## Fixed

| Area | Defect | Change |
| --- | --- | --- |
| Filtering | A second-tier options file without `*minSeqLength` left `alt_min_l` at 0, so every read below `minSeqLength` was rescued into the mid-quality tier (`Filters.h` `check_lengthXtra`). | The tier keeps the primary minimum unless `*minSeqLength` is given (`alt_min_l_set`). |
| Filtering | A read failing the primary quality window but passing the alternate one stayed green (`Filters.cpp`, quality-window block). | It is marked mid-quality like the other rescues. |
| Filtering | The 3' quality-window trim (`qualWinPos`) ran after the length check, so a read trimmed below `minSeqLength`, even to zero bases, was still written as passing. | Length is re-checked after the trim; zero-length reads fail (`MinLengthAfterTrim`). |
| Filtering | With `ExtensivePrimerChecks T`, read 2 was tested with `pair` hard-coded to 0 and reverse-complemented whenever normal read-through contained the forward primer's reverse complement (`checkIfRevPrimerHits` call). | The real mate index is passed. |
| Filtering | `Forward/Reverse primer detected` totals also counted variants beyond `PrimerErrs` recorded on the failure path. | `PrimerVariantTracker::record` takes a `detected` flag; failure-path variants are reported only. |
| Filtering | One empty `BarcodeSequence` cell set the minimum barcode length to 0 and silently failed every read of that file group as "barcode not identified". | `BarcodePreStats` rejects an empty barcode when others are present, naming the sample. |
| Filtering | `SampleIDinHead` matching was first-wins by row order, so tag `S1` claimed `S10` reads. | Longest matching tag wins (`findTag`, `detectCutBC`). |
| Filtering | Dual-barcode headers named the other mate's barcode in `orig_bc=` (`dblBCeval` passed the `pair1` flags swapped). | Flags corrected. |
| Filtering | `setBCfixed` indexed `BCdFWDREV[fwd]` while every other site uses `[!isPair1]`. | Indexes `[!fwd]`. |
| Filtering | `-5PR1cut`/`-5PR2cut` are unsigned; a negative value wrapped to a huge cut. | Negative values are rejected. |
| Seed extension | UPARSE 9+ chimera lines: `find("parents=") + 4` pointed into the keyword, and the line then fell into the `top=` parser, which overwrote the parents. With `-count_chimeras T` no chimera ever added counts. | Offset 8; chimera lines skip the `top=` branch. |
| Seed extension | The OTU matrix streamed `double` cells with six significant digits: 1234567 became `1.23457e+06`. | Cells and the total print as integers. |
| Seed extension | PAF query coverage used the alignment block length (column 11), which includes target-relative gaps and can exceed the query. | Coverage is `(qend - qstart) / qlen`. |
| Seed extension | Secondary/inversion PAF records (`tp:A:S`, `tp:A:I`) counted a read once per record in header-sample mode. | Such records are skipped. |
| Options | `parseSizeT` accepted `-5` (wrapped by `stoull`). | Negative values fall back to the default. |
| Output naming | `applyFileIT` found the extension dot anywhere in the path (`/data/run.2024/out` became `/data/run.1.2024/out`) and put rotation numbers after `.gz` (`rot.gz.1`, written uncompressed; I/O audit finding 4). | Only the basename is inspected and `.gz` stays last. |
| Statistics | A mapping file with more than 10 000 samples aborted in `addStats` after the first input file had been processed. | Cap removed (the loop is bounded anyway); `Statistics::merge` no longer returns early on an out-of-range index. |
| Helpers | `intwithcommas(int)` printed `-123` as `-,123`; `Benchmark::printResults` consumed the total it printed; the DEBUG instrumentation arrays were per-translation-unit statics. | Forwarded to the 64-bit overload; local copy; `inline` variables. |
| Records | `DNA::cutSeq(start, stop)` had no negative-start guard. | Rejected. |

## Fixed in 3.61 (second round, same day)

| Area | Defect | Change |
| --- | --- | --- |
| Barcodes | Orientation auto-detection counted hits per worker thread: with 4 threads no thread reached the 5 000-hit decision, so reverse-oriented reads were accepted run-dependently (12 000 or 11 976 of 12 000 test reads at 4 threads versus 11 788 at 1 thread); the restart path had no caller, the read that triggered the decision was dropped, and the run log never said barcodes were reversed. | `Filters::SharedRunState` holds hits, the decision and the reversal for the per-file filter and all its sub-filters; each copy swaps its own tables in `syncBarcodeOrientation`; dead restart code removed; the main filter reports the decision. Now 11 789 at 1 and 4 threads, log line present. |
| Filtering | A fractional `minSeqLength` derived its minimum from a per-thread running maximum. | The running maximum is in the shared state; each copy re-derives its minimums when it grows. |
| Barcodes | `matchSeq_tot` returned the first offset within the error budget, so a one-error match at offset 0 shadowed an exact barcode further in (wrong cut, wrong `bc_diffs`, spurious ambiguity). | The lowest-error offset wins (ties: first). |
| Filtering | `TrimStartNTs` applied only without a mapping file; `*keepPrimerSeq` was parsed into an unused value. | Trimming runs after barcode/primer removal in both modes; `*keepPrimerSeq` is ignored with a warning and removed from help and docs. |
| Seed extension | `-dere_size_fmt 1` (`_N` headers) failed the derep.map abundance check. | `;size=N;` and trailing `_N` are both parsed. |
| Seed extension | cd-hit `.clstr` parsing looped forever on CRLF files or at EOF inside a cluster, used unchecked positions on malformed lines, and dropped members buffered before a final representative line; `getMAPPERline` returned `true` on a failed stream; a missing additional `.uc` file was parsed as a dead stream. | `\r` stripped for all mapper files, EOF and malformed lines raise clear errors, buffered members are drained, a failed stream ends parsing, missing files return. |
| Seed extension | Perfect-chimera keys were not normalized like target keys; chimera counts were split over all targets in header-sample mode but over known targets in derep mode. | One `clusterKey` normalization; both paths split over known targets. |
| Subsets | `ReadSubset` validated columns for 1000 rows only (and rejected blank lines the parser skips); single-end subsets were never renamed. | All rows validated with the parser's skip rules; single-end reads renamed. |
| Statistics | `collectstats::reset()` skipped eight counters; worker statistics were not reset; mean quality was truncated in `ReportStats`; RSQS divided by a disabled `min_q`. | Full reset including sub-filters; rounded mean; RSQS uses only configured thresholds. |
| GoldenAxe | Reads below 100 bp, unassigned reads and reads outside the amplicon-count limits vanished from all statistics; with `keepBarcodeSeq 1` every sub-read failed barcode detection. | Such reads are counted as processed and rejected with a reason; sub-reads carry the resolved barcode. |
| Thread pools | A worker submitting into its own full `ThreadPool` queue could deadlock; `PartitionWorkPool::run` could unwind with published stack jobs if `afterEnqueue` threw. | Workers run the task inline when the queue is full; the batch completes before the error is rethrown. |

## Open (verified, not changed)

- **`cutPrimerRev`'s "extensive" windows are now unconditional** (`if (true)`), so the 5' reverse-primer
  search on read 2 covers up to `length - 100` positions and a chance internal match can truncate a
  read lacking the primer.
- **Records:** `matchSeq`/`matchSeqRev` carry dead non-`_NEWMATCH` branches that would not compile;
  `CompactDNA` can serialise a partial-pair block (`SDMDIFF3` with an empty mate reference) that
  `diffDNA::read` refuses; `.1`/`.2` header suffixes are stripped as mate tags, so SRA-style spot ids
  compare equal; `getDNAsubseq`/`getQual(start, end)` are unchecked (all callers verified safe).
- **Barcode orientation, residual timing:** the decision is now shared, but reads processed by other
  workers while the 5 000th hit is being counted can still be scanned in both orientations, so the
  accepted count can differ by a few reads between runs (11 789 versus 11 790 in the test above).
  A fully deterministic decision would need a pre-pass over the first reads.

## Looked solid

Mapping-file parsing (BOM, CRLF, comments, column enforcement), ONT end matching, primer cut
coordinates, merged-read checks, per-thread statistics and their merge, barcode index alignment
including `CombineSamples`, IUPAC and complement tables, quality-probability table bounds, the
`ThreadPool` wake protocol, `PartitionWorkPool` state transitions, the `OutputStreamer` lock order,
input-record validation, quality statistics and trimming bounds, the banded edit distance, FASTQ/FASTA
writers, and the `CompactDNA`/`diffDNA` codecs with their overflow checks.
