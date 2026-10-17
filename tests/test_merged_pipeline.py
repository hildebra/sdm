#!/usr/bin/env python3
"""Exercise merged demultiplexing/dereplication with paired gzip reads.

Run: python3 tests/test_merged_pipeline.py /path/to/sdm
"""
from pathlib import Path
import gzip
import shutil
import subprocess
import sys
import tempfile


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def fastq_sequences(path):
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0, path
    for offset in range(0, len(lines), 4):
        assert lines[offset].startswith('@') and lines[offset + 2] == '+', path
        assert len(lines[offset + 1]) == len(lines[offset + 3]), path
    return lines[1::4]


def fasta_records(path):
    records = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            records.append([line, ''])
        elif line:
            records[-1][1] += line
    return records


def main():
    binary = Path(sys.argv[1] if len(sys.argv) > 1 else './sdm').resolve()
    root = Path(tempfile.mkdtemp(prefix='sdm-merged-pipeline-'))
    insert = 'ACGTTGCAAGTCGATCGGATCTAGCTAACGGTTCAGTACCGATGCTAGGCTAACGTAGCTTGACCGTACAGGTTCGATACCTGA'
    forward, reverse = 'AGGTCAGTACCGTAAC', 'TCCAGATGCTACGTCA'
    barcode1, barcode2 = 'GACTGA', 'CTGTAC'
    raw1 = barcode1 + forward + insert + reverse_complement(reverse) + reverse_complement(barcode2) + 'GATCTAGTCGAT'
    raw2 = barcode2 + reverse + reverse_complement(insert) + reverse_complement(forward) + reverse_complement(barcode1) + 'CGTAGTCATGCA'
    copies = 300
    for mate, sequence in [(1, raw1), (2, raw2)]:
        with gzip.open(root / f'input.{mate}.fastq.gz', 'wt') as stream:
            for index in range(copies):
                stream.write(f'@pair_{index}/{mate}\n{sequence}\n+\n' + 'I' * len(sequence) + '\n')
    mapping = root / 'input.map'
    mapping.write_text('#SampleID\tBarcodeSequence\tBarcode2ndPair\tLinkerPrimerSequence\tReversePrimer\tfastqFile\n'
                       f'sample\t{barcode1}\t{barcode2}\t{forward}\t{reverse}\tinput.1.fastq.gz,input.2.fastq.gz\n')
    cases = [('remove', 0, 0, -1, 50), ('cap', 0, 0, 70, 50),
             ('keep_primers', 0, 1, -1, 50), ('keep_both', 1, 1, -1, 50),
             ('short_merge', 0, 0, -1, 100), ('cap_short_merge', 0, 0, 70, 100),
             ('cap_fasta', 0, 0, 70, 50)]
    try:
        for threads in (4, 12):
            for name, keep_barcode, keep_primer, truncate, minimum in cases:
                derep_format = 'fa' if name.endswith('_fasta') else 'fq'
                folder = root / f'{name}_{threads}'
                folder.mkdir()
                demulti = folder / 'demulti'
                demulti.mkdir()
                opts = folder / 'options.txt'
                settings = dict(minSeqLength=minimum, maxSeqLength=1000,
                                TruncateSequenceLength=truncate, minAvgQuality=20,
                                QualWindowWidth=10, QualWindowThreshhold=0,
                                TrimWindowWidth=10, TrimWindowThreshhold=0,
                                maxAccumulatedError=-1, BinErrorModelAlpha=-1,
                                maxAmbiguousNT=0, maxHomonucleotide=0,
                                keepBarcodeSeq=keep_barcode, keepPrimerSeq=keep_primer,
                                maxPrimerErrs=0, maxBarcodeErrs=0,
                                RejectSeqWithoutFwdPrim='T', RejectSeqWithoutRevPrim='T',
                                AmpliconShortPE='F', CheckForMixedPairs='F', CheckForReversedSeqs='F')
                opts.write_text(''.join(f'{key}\t{value}\n' for key, value in settings.items()))
                command = [str(binary), '-i_path', str(root), '-map', str(mapping),
                           '-options', str(opts), '-paired', '2', '-i_qual_offset', '33',
                           '-threads', str(threads), '-o_fastq', f'{folder}/out.1.fq,{folder}/out.2.fq',
                           '-o_demultiplex', str(demulti), '-merge_pairs_demulti', '1',
                           '-o_dereplicate', str(folder / 'derep.fna'), '-derep_format', derep_format,
                           '-merge_pairs_derep', '1', '-min_derep_copies', '1',
                           '-derepSrchLen', '-1', '-sample_sep', '___', '-log', str(folder / 'filter.log')]
                result = subprocess.run(command, cwd=root, capture_output=True, text=True, timeout=60)
                (folder / 'console.log').write_text(result.stdout + result.stderr)
                assert result.returncode == 0, f'{name}: exit {result.returncode}\n{result.stderr}'
                merged = fastq_sequences(demulti / 'sample.merg.fq')
                derep = (fastq_sequences(folder / 'derep.merg.fna') if derep_format == 'fq'
                         else [sequence for _, sequence in fasta_records(folder / 'derep.merg.fna')])
                if name in ('short_merge', 'cap_short_merge'):
                    assert not merged and not derep, 'merged product below minSeqLength was emitted'
                else:
                    expected = insert
                    if keep_primer:
                        expected = forward + expected + reverse_complement(reverse)
                    if keep_barcode:
                        expected = barcode1 + expected + reverse_complement(barcode2)
                    assert merged == [expected] * copies, f'{name}: incorrect demultiplexed merge: {merged[:1]}'
                    assert derep == [expected], f'{name}: incorrect dereplicated merge: {derep}'
                    header = (folder / 'derep.merg.fna').read_text().splitlines()[0]
                    assert f';size={copies};' in header, f'{name}: incorrect abundance: {header}'
                if truncate > 0:
                    for output in [folder / 'out.1.fq', folder / 'out.2.fq']:
                        assert all(len(seq) <= truncate for seq in fastq_sequences(output)), f'{output}: length cap bypassed'
                    for output in [folder / 'derep.1.hq.fq', folder / 'derep.2.hq.fq']:
                        assert all(len(seq) > truncate for seq in fastq_sequences(output)), f'{output}: HQ overlap was truncated'
                print(f'PASS {name}: {threads} workers', flush=True)
    except BaseException:
        print(f'Regression inputs and logs retained in {root}', file=sys.stderr)
        raise
    shutil.rmtree(root)
    print('Merged paired-gzip pipeline regressions passed')


if __name__ == '__main__':
    main()
