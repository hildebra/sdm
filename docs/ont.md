# Single-end ONT amplicons

SDM **3.51 beta** adds an opt-in extension for ordinary single-end FASTQ and FASTQ.gz amplicons. This is a bounded demultiplexer, not an ONT kit catalogue. The layouts below are tested with deterministic synthetic reads; real-kit sensitivity and specificity have not been established.

## Options and assignment

| Command-line option | Default | Meaning |
| --- | --- | --- |
| `-ontMode 0|1` | `0` | Enable ONT matching with `1`; requires a mapping file. |
| `-barcodeSearchWindow <bp>` | `200` in ONT mode | Integer 1–10000. Test barcode start offsets `0..bp-1`, measured inward from either end. A complete barcode alignment can extend beyond the window. Reads shorter than the window are bounded by their actual length. |
| `-ontBarcodeEnds either|both` | `either` | `either` accepts unique evidence at one or both ends; `both` needs two separate, consistent end hits. Ignored for assignment by filename. |

The last two options require `-ontMode 1`. Invalid values are errors. ONT mode rejects paired/comma-separated inputs, MID streams, `-onlyPair`, GoldenAxe, FASTA and alignment inputs. Mapping columns `Barcode2ndPair`, `MIDfqFile`, `SampleIDinHead` and `alignmentFile` are rejected in ONT mode, even when empty. `Barcode2ndPair` retains its paired-read meaning outside ONT mode.

In the existing options file, `maxBarcodeErrs` is the **edit-distance budget per end**, including substitutions, insertions and deletions. ONT accepts integers 0–3, default 0. A barcode must contain only A/C/G/T and be longer than this budget. Outside ONT mode the existing substitution-only policy remains. Choose a small budget appropriate to the separation between the supplied barcodes.

At each end, exact matches take priority; approximate matching is restricted to the window and error budget. A unique lowest-error sample is required. Tied samples at either end are rejected, even if the other end could resolve the tie. Confident end assignments to different samples are rejected in both modes. Equally good barcode alignments that imply different trim endpoints, overlapping end hits and contradictory primer orientation evidence are also rejected as ambiguous. Internal or terminal damage that exceeds the error budget gives no evidence at that end; it is not itself a conflicting assignment.

## Supported layout and trimming

For a map barcode `B`, forward primer `F` and reverse primer `R`, the tested layout is:

```text
optional adapter — B — spacer/flank — F — amplicon — revcomp(R)
                 — optional spacer — optional revcomp(B) — optional adapter
```

The rear barcode is `revcomp(B)` on the displayed strand. Both ends use the same sample barcode. Front-only, rear-only and dual-ended forms and the reverse complement of each are supported. Arbitrary asymmetric barcode pairs, concatemer splitting and ligation rescue are outside this implementation.

When `keepBarcodeSeq=0`, SDM removes the barcode and everything exterior to it at each detected end. The barcode alignment determines how many read bases are removed, including indels. This also removes surrounding adapters without needing an adapter catalogue. Existing `TechnicalAdapter` processing remains available for per-file reads. Existing primer matching then removes spacers/flanks and configured primers; the reverse primer is checked on the same read automatically in ONT mode. Primer start offsets use the same window, measured from the ends remaining after barcode removal. Missing primers do not justify guessing unknown adapter or flank lengths.

`maxPrimerErrs` retains its existing **substitution/IUPAC** semantics; primer indels are not supported. Existing `keepPrimerSeq`, primer-presence rejection, length and quality settings apply. Strand selection uses configured primers, compares both orientations without editing the read, and reverses sequence and quality together once. Conflicting strand evidence is rejected. If neither orientation has primer evidence, input orientation is retained and configured primer-presence rules decide whether the read passes.

Symmetric barcodes alone cannot determine biological strand. With primer columns absent, barcode matching/trimming works in either input orientation, but reads retain their input strand. Reverse-complement inserts therefore need not dereplicate together. Already oriented per-sample input is appropriate when strand normalization is required after primer removal.

## Example and existing outputs

The checked-in [example map](../tests/fixtures/ont/map.tsv), [options](../tests/fixtures/ont/options.txt) and [six synthetic reads](../tests/fixtures/ont/reads.fastq) are runnable from the repository root:

```sh
mkdir -p /tmp/sdm-ont-example/samples
./sdm -i_path tests/fixtures/ont -map tests/fixtures/ont/map.tsv \
  -options tests/fixtures/ont/options.txt -paired 1 -ontMode 1 \
  -barcodeSearchWindow 200 -ontBarcodeEnds either -threads 4 \
  -i_qual_offset 33 -sample_sep ___ \
  -o_fastq /tmp/sdm-ont-example/full.fq \
  -o_demultiplex /tmp/sdm-ont-example/samples/ -DemultiBPperSR 0 \
  -o_dereplicate /tmp/sdm-ont-example/derep.fas \
  -derep_format fa -min_derep_copies 1 \
  -log /tmp/sdm-ont-example/filter.log
```

The map uses literal tabs:

```text
#SampleID	BarcodeSequence	LinkerPrimerSequence	ReversePrimer	fastqFile
S0	ACGTCAGTGCTAGACG	GATCTGACCGTACGTA	CTAGGCTACGACTGAC	reads.fastq
S1	TGCATCGACAGTTCGA	GATCTGACCGTACGTA	CTAGGCTACGACTGAC	reads.fastq
```

Expected outputs: six full processed reads, three reads in each sample FASTQ, two dereplicated sequences, and two full HQ representatives in `derep.1.hq.fq`. One insert is shared across samples, with map counts `0:3` and `1:2`; the other has `1:1`. The fixture options disable quality trimming and truncation to isolate matching behavior; production runs should retain the intended existing filtering preset.

Formats are unchanged. `-o_demultiplex` emits full processed per-sample reads independently of dereplication; `derep.fas` and `derep.map` retain abundance and sample information. `CombineSamples` and `SequencingRun` keep their existing behavior. The emitted map and HQ representatives are tested through SDM's existing `-optimalRead2Cluster` / `-otu_matrix` path.

The normal summary gains four mutually exclusive assignment counters: assigned, unassigned, ambiguous and conflicting ends. These describe assignment before subsequent primer/quality rejection, so assigned is not a count of output reads. Per-file sample assignment is included. Failed primer/quality reads still follow existing output and dereplication policy.

## One FASTQ per sample and removed primers

Omit barcode columns to retain filename-based assignment:

```text
#SampleID	fastqFile	LinkerPrimerSequence	ReversePrimer
S0	S0.fastq.gz	GATCTGACCGTACGTA	CTAGGCTACGACTGAC
S1	S1.fastq.gz	GATCTGACCGTACGTA	CTAGGCTACGACTGAC
```

Each input filename must identify one sample. This workflow performs primer processing and count tracking without barcode matching; `both` does not impose a barcode requirement. For already-trimmed reads, omit both primer columns and use `RejectSeqWithoutFwdPrim=F` and `RejectSeqWithoutRevPrim=F`. No primer-state flag was added. Use a valid `-options` file when primer processing is wanted: SDM's existing missing-options behavior disables filtering/primer removal.

## LotuS integration requirements

No LotuS files, Barbell files or bundled executables are modified by this change.

- Require SDM **3.51** or a build that explicitly advertises `ONT amplicon end matching: enabled` in `-version`. Earlier SDM command parsing may silently accept unknown flags, so successfully invoking an old binary with new flags is not a capability test.
- For the SDM ONT preprocessing path, pass `-ontMode 1 -paired 1`. Pass the final `-barcodeSearchWindow` and `-ontBarcodeEnds` options if exposing overrides; otherwise their ONT defaults apply. Leave `-ontMode` off for other technologies and later processing of HQ seeds/representatives.
- For pooled barcode input, supply `BarcodeSequence` and the actual primers. For one-file-per-sample input, omit barcode columns. Remove incompatible mapping columns entirely from the ONT map.
- Keep LotuS's existing `-ontPrimerState present|removed` map/preset preparation. In `removed` state, omit primer columns and disable primer-presence rejection. No additional SDM primer-state flag is necessary.
- Keep sending **full per-sample FASTQ** to Savont and retain `derep.map` for backmapping/counts. Do not substitute dereplicated representatives for the full-read input. Preserve the intended quality policy, avoid an unintended `TruncateSequenceLength` cap, and keep `-DemultiBPperSR 0` when all passing reads are required.
- LotuS currently defaults to `-ontDemux 0` and does not pass the new SDM options. Wiring that SDM route and installing the validated executable are separate wrapper changes; the external Barbell route remains optional.

## Validation

`make test` includes `tests/TestSearch.cpp` and `perl tests/test_ont_pipeline.pl ./sdm`. The search tests compare 2,000 deterministic cases with an independent full edit-distance reference. The pipeline tests exercise exact/substitution/indel matching, reverse qualities, 17–50 bp adapter offsets, 32–59 bp barcode-to-primer spacers, front/rear/missing/conflicting ends, ties, window boundaries, short reads, primer retention/removal, absent primers, per-file assignment, combined samples, sequencing runs, full-read/derep/HQ counts, and matrix generation at one and four threads. The runtime fixture repeats 4,000 reads against a 24-barcode panel; it is a sanity check, not a performance guarantee.
