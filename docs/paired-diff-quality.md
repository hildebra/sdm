# Exact paired variants with retained qualities

Available in **SDM 3.52 beta**: `-derepStoreQuals 1` supports complete, unmerged paired reads. Input layout selects single-end or paired storage; the flag remains off by default. Both mates must have complete input qualities.

```sh
sdm -i_fastq reads.1.fq,reads.2.fq -paired 2 \
    -o_dereplicate clusters.fna -options sdm_options.txt \
    -derepIdentity 97 -derepStoreQuals 1 -threads 1
```

With identity below 100 and multiple workers, the default uses 64 independent
coarse partitions routed by the effective R1 key. Identical R1 keys stay in
one partition even when R2 differs; a shared work pool schedules ready partitions; the final output preserves exact R1 dereplicates without merging search groups. `-derepGlobal 1`
or `-derepReassign 1` retains a global index. See the
[parallel-mode semantics](options.md#streaming-coarse-dereplication).

This selects **exact `(R1, R2)` variants inside exact effective R1 dereplicates**, sharing coarse references for compression. Different full backing sequences or lengths in either mate after physical preprocessing cuts produce distinct exact variants, even when their clustering prefixes match. A pair is counted once as a fragment. In particular, identical R1 sequences with different R2 sequences do not collapse into one exact variant.

By default, coarse groups are internal search and compression references. The main FASTA and `.map` retain exact effective search-key dereplicates, including the configured search/truncation cap and ordinary prefix consolidation when selected. With quality retention, the search source is currently R1. R2 differences do not split that parent; the existing paired better-seed selection can replace both mates. Minimum-copy rules apply to these exact dereplicates. Select `-derepCoarseClusters 1` only when the output should describe coarse clusters instead. `-derepReassign 1` also explicitly selects coarse-cluster output and performs final reassignment. `-derepGlobal 1` controls search scope, not output granularity.

## Storage and representative selection

Each paired cluster has two immutable initial nucleotide references. `diffDNA` uses one length-delimited pair of delta keys in its existing variant map, one fragment count or sparse sample vector, and one selected fragment payload. It does not put two independently counting `diffDNA` objects underneath a pair.

The selected payload contains both processed quality vectors and both compact scoring summaries. Selection reuses the existing paired `whoIsBetter` formula through the shared statistical helper. **Both quality vectors are replaced together**, so exported mates belong to one observation. Independently selecting the best R1 and best R2 is not allowed. Count/sample data and comparison history survive replacement, and both nucleotide references stay fixed.

Default exact parents use **R1 passing status** for admission and representative replacement, just like ordinary dereplication. An R1-passing pair is not discarded because R2 failed. Selection of the observed quality witness for each full-pair variant is independent and prefers both mates passing. With explicit coarse-parent output, creating a coarse parent also requires both mates to pass. In fixed-key mode, complete failed repeats count only after an exact R1 parent exists; in prefix mode, failures are deferred and attached to a compatible passing parent at output, independently of arrival order. Failed reads cannot supply the parent representative.

A mate physically reduced to zero length cannot participate in paired variant storage; the fragment is skipped and reported as an empty-R2 retention exclusion. Missing mates or missing qualities are rejected. Native `merge_pairs_derep`, `merge_pairs_filter`, and `merge_pairs_demulti` must be disabled. HQ variants retain full available mates and their observed quality vectors after physical technical cuts, including tails hidden by logical truncation or quality trimming. Logical lengths remain separate in the scoring metadata. A logically empty mate is still retainable when its backing sequence exists. See the [behavioral audit and required merged/full-pair changes](derep-behavior-parity.md).

Snapshots share immutable quality payloads, and subsequent replacements leave earlier snapshots unchanged. A merge/reassignment rebases both nucleotide deltas and carries the quality payload plus the shared fragment/sample counts. Merging compressed witnesses compares the selected surviving fragments, rather than replaying discarded observations.

## Search-group admission and final reassignment

**An identical effective R1 search key reuses its established parent dereplicate, regardless of R2 identity.** The existing paired `whoIsBetter` selection evaluates passing observations and can replace both representative mates together, including their full sequences beyond search caps and logical trimming. Different full pairs remain distinct retained variants with their own counts and quality witnesses. Exact assignments persist for the lifetime of the index and bypass seed lookup and candidate verification.

For a previously unseen R1 key, the existing R1 seed index supplies candidates. Its homopolymer exclusions and fallback sensitivity remain unchanged. Each candidate is then checked against **both** corresponding initial references:

- Same-start, substitutions-only comparison, independently across the shorter R1 and shorter R2 lengths.
- Each mate independently meets the configured `derepIdentity` threshold. A very good R1 cannot compensate for an R2 below the threshold.
- `derepSrchLen` and the effective truncation cap apply to each mate independently. A naturally short R1 does not shorten the R2 comparison. Full processed sequences beyond a search-only cap remain in the exact variant key.

There is no indel alignment, mate swapping, or reverse-complement search in this comparison. Supply consistently oriented reads. At identity 100, a previously unseen R1 key uses exact-prefix comparison on both mates; an identical known R1 key takes the reuse path above. Exact variant identity uses both full backing mate sequences and lengths after physical cuts; logical trimming does not merge distinct tails.

With `-derepReassign 1`, output explicitly describes coarse clusters, and both mates must qualify against a candidate's final representative pair. Closeness is:

```text
(R1 substitutions + R2 substitutions) / (R1 compared bases + R2 compared bases)
```

Only a strictly closer qualifying pair receives the variant and all its counts. Ties retain the original group; equal alternative targets use the earliest cluster. Reassignment remains one pass with fixed representatives. The ordinary quality-disabled dereplication paths retain their previous matching behavior.

## Output

| Option | Paired output |
| --- | --- |
| `derepStoreQuals=1` | `clusters.1.hq.fq` and `clusters.2.hq.fq` |
| Additionally `derepStoreDiffs=1` | `clusters.diff`, using linked `SDMDIFF3` blocks |
| Additionally `derepSubclusterFasta=1` | `clusters.subclusters.1.fna` and `clusters.subclusters.2.fna`; requires binary delta export as before |

Mate files have identical generated IDs, abundances, and record order. Headers look like `@parent.sub7;size=42;`; the count is 42 fragments in both files, not 84 reads. Pair exports include groups below the parent minimum-copy cutoff, and the standard HQ files accumulate across sequencing runs, like normal dereplication. Optional `.diff` and exact FASTA exports keep their run suffixes. Original read names are not archived. Quality vectors are selected observations, not per-base averages or a consensus.

The standard HQ files contain these variants instead of an additional parent-only copy. This includes the representative's exact pair, other retained pairs, and singleton clusters in the same files, each once with its own fragment count. No separate `.subclusters` FASTQ files are written. The normal parent representative output and `.map` remain available. Single-end input continues to write `.1.hq.fq` and nucleotide-only `SDMDIFF2`. Both binary versions omit qualities and sample IDs; retaining qualities in memory does not change that serialization policy.

## Native seed extension

Pass these reconstructed FASTQ exports directly to seed extension with `-seedSubclusters 1` (default: `0`). Supply both mate files, the original parent mapper assignments, `-derep_map`, and `-OTU_fallback`. Every mapped exact pair competes as a seed, while the parent's sample vector is added once. `-merge_pairs_seed 1` can merge the exported mates without repeating preprocessing cuts or truncation, and receives full available mate tails from the updated HQ exports. Barcode/primer cuts must be completed during preprocessing; the reconstructed seed reader does not repeat them. Binary `.diff` export remains optional. See the [seed-extension command and count semantics](workflows.md#reconstructed-exact-variants) and the [LotuS3 integration instructions](lotus3-dereplication-io.md).

## Paired binary format and API

`SDMDIFF3` is an explicit paired nucleotide format, so older readers will reject it rather than silently treating R2 variants as R1 duplicates. The updated `diffDNA::read` accepts single-end versions 1/2 and paired version 3.

A version 3 block contains:

1. Eight bytes `SDMDIFF3`.
2. Cluster label, R1 reference, and R2 reference, each preceded by its unsigned LEB128 byte length.
3. Total fragment count, exact-reference-pair count, and payload byte length as unsigned LEB128 values.
4. Each non-reference pair: fragment count, pair-key byte length, and pair-key bytes.

A pair key contains the R1 patch byte length, R1 patch bytes, and R2 patch bytes. Each nonempty mate patch reuses the existing processed-length/change-count/position-gap/base encoding. An empty mate patch denotes that mate's reference; the empty overall key denotes the complete reference pair. Length delimiters preserve boundaries even when two different mate pairs have the same concatenated nucleotide string.

Use `diffDNA::forEachPair` to obtain both reconstructed sequences, a fragment count, optional sample provenance, and the selected quality record when retained. Single-sequence visitors/exports reject paired objects, preventing accidental loss of mate linkage. Paired `writeFastq`/`writeFasta` overloads take two streams and export each fragment synchronously. Reloaded binary data can produce paired FASTA but cannot recover the selected qualities.

## Validation

The mixed-quality 1,000-fragment comparison at 1, 4, and 12 threads is documented in the [behavioral audit](derep-behavior-parity.md). It exercises ordinary merged search, ordinary R1 search, retained 100%, and retained 97%, including deferred recovery, empty mates and passing/rest reconciliation.

The optimized portable build validates:

- 1,600 paired observations against uncompressed paired `DNAunique` representative selection, with changes in mate eligibility, quality, primer metadata, hidden tails, and mate lengths.
- 1,200 candidate queries at 97%, 99%, and 100% against an independent exhaustive two-mate oracle, including short, ambiguous and homopolymer-rich sequences.
- Exact R1 reuse despite unrelated R2, better-pair replacement with different R1 tails beyond a search cap, persistent assignments after index growth, and R2-driven final reassignment, exact pair boundaries, one shared sample total per fragment, snapshot independence, rebasing, self-merges, binary round trips and truncated-file rejection.
- Independent mate caps, full sequences beyond those caps, per-run names, below-cutoff exports, resumed input, missing mates/qualities, and incompatible native merging flags.
- The full paired FASTQ pipeline: 240 distinct pairs and 720 observed fragments, each with three quality profiles. Exact pair sequences, both selected source quality vectors, synchronized names/counts, and sample-map totals agree at 97%, 99%, and 100%, with reassignment, and with four threads plus a 100 bp search cap. Enabled FASTA output agrees record-for-record with FASTQ.
- Existing single-end retained-quality, ONT, input-helper, search, and targeted regression suites.
- The final coarse/diffDNA suite also passes AddressSanitizer, UndefinedBehaviorSanitizer, and leak detection.

Run the paired pipeline test directly, or through `make test`:

```sh
python3 tests/test_paired_diff_quality_pipeline.py /path/to/sdm
```

The paired path adds a second reference, patch, quality vector and scoring summary per selected pair, plus R2 candidate verification. Memory depends on unique pair count and both processed lengths; R2 diversity can increase variant counts even when all R1 sequences are identical. The earlier single-end benchmark is not a measured paired performance result.
