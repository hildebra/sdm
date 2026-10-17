#!/usr/bin/env python3
"""Check full merged seed output in the requested OTU FASTA and FASTQ.

Run: python3 tests/test_seed_extension.py /path/to/sdm
"""
from pathlib import Path
import gzip
import random
import shutil
import subprocess
import sys
import tempfile

from test_merged_pipeline import fasta_records, fastq_sequences, reverse_complement


def run_case(binary, root, name, output_format, cap=-1, stringent=False,
             assignments=True, merge=True):
    folder = root / name
    folder.mkdir()
    rng = random.Random(73521)
    consensus = ''.join(rng.choice('ACGT') for _ in range(180))
    for mate, sequence in [(1, consensus[:120]), (2, reverse_complement(consensus[60:]))]:
        with gzip.open(folder / f'input.{mate}.fq.gz', 'wt') as stream:
            stream.write(f'@pair/{mate}\n{sequence}\n+\n' + 'I' * len(sequence) + '\n')
    (folder / 'fallback.fna').write_text('>OTU_1\n' + consensus[:120] + '\n')
    (folder / 'assign.paf').write_text(
        'pair/1\t120\t0\t120\t+\tOTU_1\t120\t0\t120\t120\t120\t60\n' if assignments else '')
    settings = dict(minSeqLength=50, maxSeqLength=1000, TruncateSequenceLength=cap,
                    minAvgQuality=45 if stringent else 0,
                    QualWindowWidth=10, QualWindowThreshhold=45 if stringent else 0,
                    TrimWindowWidth=10, TrimWindowThreshhold=45 if stringent else 0,
                    maxAccumulatedError=0.005 if stringent else -1,
                    BinErrorModelAlpha=-1, maxAmbiguousNT=0, maxHomonucleotide=0,
                    RejectSeqWithoutFwdPrim='T', RejectSeqWithoutRevPrim='T')
    (folder / 'options.txt').write_text(''.join(f'{k}\t{v}\n' for k, v in settings.items()))
    command = [str(binary), '-i_fastq', f'{folder}/input.1.fq.gz,{folder}/input.2.fq.gz',
               '-paired', '2', '-optimalRead2Cluster', str(folder / 'assign.paf'),
               '-OTU_fallback', str(folder / 'fallback.fna'), '-options', str(folder / 'options.txt'),
               '-merge_pairs_seed', '1' if merge else '0', '-i_qual_offset', '33',
               '-o_qual_offset', '33', '-sample_sep', '___', '-log', str(folder / 'log.txt'),
               '-threads', '4', '-oneLineFastaFormat', '1', '-5PR1cut', '100', '-5PR2cut', '100']
    write_qualities = output_format == 'paired_fasta' or stringent
    if output_format == 'fastq':
        command += ['-o_fastq', f'{folder}/OTU.1.fq,{folder}/OTU.2.fq']
    else:
        paths = str(folder / 'OTU.fna')
        quals = str(folder / 'OTU.qual')
        if output_format == 'paired_fasta':
            paths += ',' + str(folder / 'OTU.2.fna')
            quals += ',' + str(folder / 'OTU.2.qual')
        command += ['-o_fna', paths]
        if write_qualities:
            command += ['-o_qual', quals]
    result = subprocess.run(command, cwd=folder, capture_output=True, text=True, timeout=60)
    (folder / 'console.log').write_text(result.stdout + result.stderr)
    assert result.returncode == 0, f'{name}: exit {result.returncode}\n{result.stderr}'
    expected = consensus if merge and assignments else consensus[:120]
    if output_format == 'fastq':
        assert fastq_sequences(folder / 'OTU.merg.fq') == [expected], name
        assert (folder / 'OTU.merg.fq').read_text().startswith('@OTU_1\n'), name
        assert all(not (folder / f'OTU.{suffix}.fq').read_text()
                   for suffix in ['1', '2', '1.singl', '2.singl']), name
    else:
        records = fasta_records(folder / 'OTU.fna')
        assert records == [['>OTU_1', expected]], f'{name}: expected {len(expected)} bases in OTU.fna, got {records}'
        if write_qualities:
            qualities = fasta_records(folder / 'OTU.qual')
            assert len(qualities) == 1 and qualities[0][0] == '>OTU_1', name
            assert len(qualities[0][1].split()) == len(expected), name
        if output_format == 'paired_fasta':
            assert not (folder / 'OTU.2.fna').read_text(), name
            assert not (folder / 'OTU.2.qual').read_text(), name
    print(f'PASS {name}', flush=True)


def run_derep_to_seed(binary, root, threads, merge_derep):
    """A 200-base clustering cap must not destroy the overlap in HQ mates."""
    folder = root / f'derep_to_seed_{threads}_merge{int(merge_derep)}'
    folder.mkdir()
    # A representative amplicon from the reported regression. Two 300-base
    # reads overlap; their 200-base prefixes cannot reconstruct this sequence.
    consensus = (
        'AGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGGCCCTGCTTTTGT'
        'GGGGTGCTCGAGTGGCGAACGGGTGAGTAACACGTGAGTAACCTGCCCTTGACTTTGGGATAACTTCAGGAAACTGGGG'
        'CTAATACCGGATAGGAGCTCCTGCTGCATGGTGGGGGTTGGAAAGTTTCGGCGGTTGGGGATGGACTCGCGGCTTATCA'
        'GCTTGTTGGTGGGGTAGTGGCTTACCAAGGCTTTGACGGGTAGCCGGCCTGAGAGGGTGACCGGCCACATTGGGACTGA'
        'GATACGGCCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGGAAGCCTGATGCAGCAACGCCGCG'
        'TGCGGGATGACGGCCTTCGGGTTGTAAACCGCTTTCGCCTGTGACGAAGCGTGAGTGACGGTAATGGGTAAAGAAGCACC'
        'GGCTAACTACGTGCCAGCAGCCGCGGTAA')
    assert 400 < len(consensus) < 600
    mates = [consensus[:300], reverse_complement(consensus[-300:])]
    for mate, sequence in enumerate(mates, 1):
        with gzip.open(folder / f'raw.{mate}.fq.gz', 'wt') as stream:
            for i in range(300):
                stream.write(f'@read{i}/{mate}\n{sequence}\n+\n' + 'G' * len(sequence) + '\n')
    (folder / 'map.tsv').write_text('#SampleID\tfastqFile\nS\traw.1.fq.gz,raw.2.fq.gz\n')
    (folder / 'options.txt').write_text(
        'minSeqLength\t200\nmaxSeqLength\t1000\nTruncateSequenceLength\t200\n'
        'minAvgQuality\t20\nQualWindowThreshhold\t0\nTrimWindowThreshhold\t25\n'
        'TrimWindowWidth\t20\nmaxAccumulatedError\t0.95\nBinErrorModelAlpha\t-1\n'
        'maxAmbiguousNT\t0\nmaxHomonucleotide\t0\nkeepPrimerSeq\t1\n')
    common = ['-paired', '2', '-sample_sep', '___', '-i_qual_offset', '33',
              '-o_qual_offset', '33', '-options', str(folder / 'options.txt')]
    def execute(command, label):
        result = subprocess.run([str(binary), *common, *command], cwd=folder,
                                capture_output=True, text=True, timeout=60)
        (folder / f'{label}.log').write_text(result.stdout + result.stderr)
        assert result.returncode == 0, result.stderr
    execute(['-i_path', str(folder), '-map', str(folder / 'map.tsv'), '-threads', str(threads),
             '-o_dereplicate', str(folder / 'derep.fna'), '-min_derep_copies', '1',
             '-merge_pairs_derep', str(int(merge_derep)), '-derep_format', 'fa',
             '-o_fastq', f'{folder}/filtered.1.fq,{folder}/filtered.2.fq',
             '-log', str(folder / 'filter.stats')], 'derep')
    clustered = fasta_records(folder / ('derep.merg.fna' if merge_derep else 'derep.fna'))
    assert len(clustered) == 1 and clustered[0][1] == (consensus if merge_derep else consensus[:200]), clustered
    for mate, expected in enumerate(mates, 1):
        assert fastq_sequences(folder / f'filtered.{mate}.fq') == [expected[:200]] * 300
        hq = fastq_sequences(folder / f'derep.{mate}.hq.fq')
        assert hq == [expected], f'HQ R{mate}: expected 300 bases, got {[len(x) for x in hq]}'
    seed_id = (folder / 'derep.1.hq.fq').read_text().splitlines()[0][1:]
    (folder / 'tmp_otu.fa').write_text(f'>{seed_id}\n{clustered[0][1]}\n')
    (folder / 'finalOTU.uc').write_text(f'{seed_id}\totu1\t*\n')
    (folder / 'finalOTU.ADD.paf').write_text('')
    (folder / 'finalOTU.REST.paf').write_text('')
    execute(['-i_fastq', f'{folder}/derep.1.hq.fq,{folder}/derep.2.hq.fq',
             '-o_fastq', f'{folder}/otu_seeds.1.fq,{folder}/otu_seeds.2.fq',
             '-uparseVer', '12', '-optimalRead2Cluster', str(folder / 'finalOTU.uc'),
             '-derep_map', str(folder / 'derep.map'), '-merge_pairs_seed', '1',
             '-mergedPairs', '0', '-OTU_fallback', str(folder / 'tmp_otu.fa'),
             '-ucAdditionalCounts', str(folder / 'finalOTU.ADD.paf'),
             '-ucAdditionalCounts1', str(folder / 'finalOTU.REST.paf'),
             '-otu_matrix', str(folder / 'OTU.txt'), '-id', '97', '-minQueryCov', '0.8',
             '-count_chimeras', 'T', '-log', str(folder / 'seed.stats')], 'seed')
    assert fastq_sequences(folder / 'otu_seeds.merg.fq') == [consensus]
    for suffix in ['1', '2', '1.singl', '2.singl']:
        assert not (folder / f'otu_seeds.{suffix}.fq').read_text(), suffix
    matrix = (folder / 'OTU.txt').read_text().splitlines()
    header = (folder / 'otu_seeds.merg.fq').read_text().splitlines()[0][1:]
    assert matrix[1] == f'{header}\t300', matrix
    print(f'PASS {folder.name}: {len(clustered[0][1])}-base clustering, 300-base HQ mates, '
          f'{len(consensus)}-base seed', flush=True)


def main():
    binary = Path(sys.argv[1] if len(sys.argv) > 1 else './sdm').resolve()
    root = Path(tempfile.mkdtemp(prefix='sdm-seed-extension-'))
    try:
        for threads in (4, 12):
            for merge_derep in (False, True):
                run_derep_to_seed(binary, root, threads, merge_derep)
        # The paired-FASTA case previously succeeded but wrote just 120-base R1
        # even though the log reported a 180-base merged representative.
        for fmt in ['paired_fasta', 'fasta', 'fastq']:
            for stringent in [False, True]:
                for cap in [-1, 150]:
                    run_case(binary, root, f'{fmt}_strict{int(stringent)}_cap{cap}',
                             fmt, cap=cap, stringent=stringent)
        run_case(binary, root, 'fallback_in_main_fasta', 'fasta', cap=60, assignments=False)
        run_case(binary, root, 'merging_disabled', 'fasta', cap=60, merge=False)
    except BaseException:
        print(f'Regression inputs and logs retained in {root}', file=sys.stderr)
        raise
    shutil.rmtree(root)
    print('Seed-extension output regressions passed (18 cases)')


if __name__ == '__main__':
    main()
