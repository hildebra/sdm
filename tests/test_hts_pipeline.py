#!/usr/bin/env python3
"""End-to-end HTS regressions. Run via make test; fixtures stay in /tmp."""
from pathlib import Path
import subprocess
import sys
import tempfile

SDM = str(Path(sys.argv[1] if len(sys.argv) > 1 else './sdm').resolve())
CONVERTER = str(Path(sys.argv[2] if len(sys.argv) > 2 else '.build/hts-1/TestHtsInput').resolve())
HEADER = '@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n'
SEQ = 'ACGTCAGTACGATCGATGCATGCTACGTAGCTAGTCAGTC'
QUAL = 'FGHI' * (len(SEQ) // 4)
assert len(QUAL) == len(SEQ)


def record(name, flag=0, seq=SEQ, qual=QUAL, tag=''):
    return f'{name}\t{flag}\tchr1\t1\t60\t{len(seq)}M\t*\t0\t0\t{seq}\t{qual}{tag}\n'


def rc(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def reads(path):
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0, (path, lines)
    return sorted((lines[i].split()[0], lines[i+1], lines[i+3]) for i in range(0, len(lines), 4))


def run(root, name, source=None, paired=1, extra=(), valid=True, quality='auto', threads=4, demultiplex=False):
    directory = root / name
    directory.mkdir()
    options = directory / 'options.txt'
    options.write_text('minSeqLength\t2\nmaxSeqLength\t1000\nminAvgQuality\t0\n'
                       'maxAmbiguousNT\t10\nkeepBarcodeSeq\t0\nfastqVersion\t' +
                       {'33': '1', '64': '2'}.get(quality, quality) + '\n')
    outputs = [directory / f'out.{mate}.fq' for mate in range(1, paired + 1)]
    command = [SDM, '-options', str(options), '-threads', str(threads),
               '-iniBlockSize', '1', '-o_fastq', ','.join(map(str, outputs)),
               '-log', str(directory / 'run.log')]
    if source is not None:
        command += ['-i', str(source)]
    command += list(map(str, extra))
    if demultiplex:
        (directory / 'demux').mkdir()
        command += ['-o_demultiplex', str(directory / 'demux') + '/']
    result = subprocess.run(command, cwd=directory, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True, timeout=90)
    (directory / 'console.log').write_text(result.stdout)
    assert (result.returncode == 0) == valid, (command, result.returncode, result.stdout)
    return [reads(p) for p in outputs], result.stdout


def convert(source, target, fmt, reference=None):
    subprocess.run([CONVERTER, '--convert', str(source), str(target), fmt] + ([str(reference)] if reference else []), check=True,
                   stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def test_mapping(root):
    # The dedicated alignment column can contain every HTS format in one run.
    files = [('A', 'a.sam'), ('B', 'b.bam'), ('C', 'c.cram'), ('D', 'd.fq')]
    for sample, filename in files:
        sam = root / (sample + '_mapping.sam')
        sam.write_text(HEADER + record(sample))
        target = root / filename
        if filename.endswith('.sam'):
            target.write_text(sam.read_text())
        elif filename.endswith('.fq'):
            target.write_text(f'@{sample}\n{SEQ}\n+\n{QUAL}\n')
        else:
            convert(sam, target, target.suffix[1:])
    for order in (('fastqFile', 'alignmentFile'), ('alignmentFile', 'fastqFile')):
        mapping = root / 'mixed.tsv'
        lines = ['#SampleID\t' + '\t'.join(order)]
        for sample, filename in reversed(files):
            fields = {'fastqFile': filename if sample == 'D' else '',
                      'alignmentFile': filename if sample != 'D' else ''}
            lines.append(sample + '\t' + '\t'.join(fields[col] for col in order))
        mapping.write_text('\n'.join(lines) + '\n')
        actual, _ = run(root, 'mixed_' + order[0], quality='33', demultiplex=True,
                        extra=['-i_path', root, '-map', mapping])
        assert actual[0] == [(f'@{sample}__{sample}', SEQ, QUAL) for sample, _ in files], actual
        for sample, _ in files:
            output = root / ('mixed_' + order[0]) / 'demux' / (sample + '.fq')
            assert reads(output) == [(f'@{sample}', SEQ, QUAL)], (output, reads(output))

    # Dedicated column alone, extensionless BAM, and blank optional file cells.
    (root / 'extensionless').write_bytes((root / 'b.bam').read_bytes())
    mapping = root / 'alignment_only.tsv'
    mapping.write_text('#SampleID\talignmentFile\tfnaFile\tqualFile\tMIDfqFile\n'
                       'A\ta.sam\t\t\t\nB\textensionless\t\t\t\nC\tc.cram\t\t\t\n')
    actual, _ = run(root, 'alignment_only', extra=['-i_path', root, '-map', mapping])
    assert actual[0] == [(f'@{sample}__{sample}', SEQ, QUAL) for sample in 'ABC'], actual

    mapping.write_text('#SampleID\talignmentFile\tfnaFile\nA\ta.sam\t\n')
    actual, _ = run(root, 'empty_fasta_column', extra=['-i_path', root, '-map', mapping])
    assert actual[0] == [('@A__A', SEQ, QUAL)], actual

    mapping.write_text('#SampleID\talignmentFile\tfastqFile\n'
                       'F\t\tinput.1.fq,input.2.fq\nA\tpairs.bam\t\n')
    actual, _ = run(root, 'mixed_pairs', paired=2, quality='33',
                    extra=['-i_path', root, '-map', mapping])
    for mate, output in enumerate(actual, 1):
        # Existing per-file assignment prefixes R1; R2 retains its read header.
        expected = sorted((f'@{sample}__r{i:04}/1' if mate == 1 else f'@r{i:04}/2', SEQ, QUAL)
                          for sample in ('F', 'A') for i in range(140))
        assert output == expected, (mate, output[:2], len(output))

    for fmt in ('sam', 'bam', 'cram'):
        mapping.write_text('#SampleID\talignmentFile\nA\tpairs.' + fmt + '\n')
        actual, _ = run(root, 'mapped_pairs_' + fmt, paired=2, quality='64',
                        extra=['-i_path', root, '-map', mapping])
        assert all(len(mate) == 140 for mate in actual), actual
        assert all(seq == SEQ and qual == QUAL for mate in actual for _, seq, qual in mate), actual
    mapping.write_text('#SampleID\talignmentFile\nA\texternal.cram\n')
    actual, _ = run(root, 'mapped_reference', paired=2,
                    extra=['-i_path', root, '-map', mapping, '-cramRef', root / 'reference.fa'])
    assert all(len(mate) == 140 for mate in actual), actual

    # One alignment input shared by two barcode mapping rows retains grouping.
    tagged = root / 'mapped_tags.sam'
    tagged.write_text(HEADER + record('a', tag='\tBC:Z:ACGT') + record('b', tag='\tBC:Z:TTTT'))
    mapping.write_text('#SampleID\talignmentFile\tBarcodeSequence\n'
                       'A\tmapped_tags.sam\tACGT\nB\tmapped_tags.sam\tTTTT\n')
    actual, _ = run(root, 'mapped_tags', extra=['-i_path', root, '-map', mapping, '-htsBarcode', 'BC'])
    assert len(actual[0]) == 2 and {h[:2] for h, _, _ in actual[0]} == {'@A', '@B'}, actual

    # Short relative input directories must be resolved just like absolute ones.
    mapping.write_text('#SampleID\talignmentFile\nA\ta.sam\n')
    actual, _ = run(root, 'mapped_relative', extra=['-i_path', '..', '-map', mapping])
    assert actual[0] == [('@A__A', SEQ, QUAL)], actual

    invalid = [
        ('both', '#SampleID\tfastqFile\talignmentFile\nA\td.fq\ta.sam\n', 'exactly one'),
        ('missing', '#SampleID\tfastqFile\talignmentFile\nA\t\t\n', 'exactly one'),
        ('missing_after_valid', '#SampleID\talignmentFile\nA\ta.sam\nB\t\n', 'exactly one'),
        ('wrong_format', '#SampleID\talignmentFile\nA\td.fq\n', 'requires SAM/BAM/CRAM'),
        ('two_files', '#SampleID\talignmentFile\nA\ta.sam,b.bam\n', 'accepts one'),
        ('mid', '#SampleID\talignmentFile\tMIDfqFile\nA\ta.sam\td.fq\n', 'cannot be combined'),
        ('quality', '#SampleID\talignmentFile\tqualFile\nA\ta.sam\treads.qual\n', 'qualFile'),
        ('fasta', '#SampleID\talignmentFile\tfnaFile\nA\ta.sam\treads.fa\n', 'fnaFile'),
        ('duplicate_header', '#SampleID\talignmentFile\talignmentFile\nA\ta.sam\tb.bam\n', 'ambiguous entries'),
    ]
    for name, text, error in invalid:
        mapping.write_text(text)
        actual, log = run(root, 'invalid_map_' + name, valid=False,
                          extra=['-i_path', root, '-map', mapping])
        assert error in log and actual == [[]], (name, actual, log)


def tests(root):
    sam = root / 'pairs.sam'
    sam.write_text(HEADER + ''.join(
        record(f'r{i:04}', 145, rc(SEQ), QUAL[::-1]) + record(f'r{i:04}', 65)
        for i in range(140)))
    expected = [[(f'@r{i:04}/{mate}', SEQ, QUAL) for i in range(140)] for mate in (1, 2)]
    for fmt in ('sam', 'bam', 'cram'):
        source = sam if fmt == 'sam' else root / ('pairs.' + fmt)
        if fmt != 'sam':
            convert(sam, source, fmt)
        actual, _ = run(root, 'pairs_' + fmt, source, paired=2, quality='64')
        assert actual == expected, (fmt, actual[:1])
        limited, _ = run(root, 'limit_' + fmt, source, paired=2, extra=['-XfirstReadsRead', '6'])
        assert limited == [r[:3] for r in expected], (fmt, limited)
    reference = root / 'reference.fa'
    reference.write_text('>chr1\n' + 'A' * 1000 + '\n')
    (root / 'reference.fa.fai').write_text('chr1\t1000\t6\t1000\t1001\n')
    external_cram = root / 'external.cram'
    convert(sam, external_cram, 'cram', reference)
    actual, _ = run(root, 'external_reference', external_cram, paired=2,
                    extra=['-cramRef', reference])
    assert actual == expected
    # Equivalent FASTQ runs through the existing reader to detect HTS-only drift.
    fqs = [root / 'input.1.fq', root / 'input.2.fq']
    for mate, path in enumerate(fqs, 1):
        path.write_text(''.join(f'@r{i:04}/{mate}\n{SEQ}\n+\n{QUAL}\n' for i in range(140)))
    actual, _ = run(root, 'fastq', ','.join(map(str, fqs)), paired=2, quality='33')
    assert actual == expected, [(len(a), len(e), a[:2], e[:2]) for a, e in zip(actual, expected)]
    actual, _ = run(root, 'only_r2', root / 'pairs.bam', extra=['-onlyPair', '2'])
    assert actual[0] == expected[1]

    for dual in (False, True):
        mapping = root / f'barcode_{dual}.tsv'
        mapping.write_text('#SampleID\tBarcodeSequence' + ('\tBarcode2ndPair' if dual else '') +
                           '\nA\tACGT' + ('\tTGCA' if dual else '') +
                           ('\nB\tACGT' if dual else '\nB\tTTTT') + ('\tCCCC' if dual else '') + '\n')
        tagged = root / f'tagged_{dual}.sam'
        bc_a, bc_b = ('ACGT-TGCA', 'ACGT-CCCC') if dual else ('ACGT', 'TTTT')
        tagged.write_text(HEADER + record('a', 65, tag='\tBC:Z:' + bc_a) + record('a', 129) +
                          record('b', 65) + record('b', 129, tag='\tBC:Z:' + bc_b) +
                          record('c', 65, tag='\tBC:Z:GGGG') + record('c', 129))
        for fmt in ('sam', 'bam', 'cram'):
            source = tagged if fmt == 'sam' else root / f'tagged_{dual}.{fmt}'
            if fmt != 'sam':
                convert(tagged, source, fmt)
            actual, _ = run(root, f'tags_{dual}_{fmt}', source, paired=2,
                            extra=['-map', mapping, '-htsBarcode', 'BC'], demultiplex=True)
            for mate, output in enumerate(actual, 1):
                assert len(output) == 2, (dual, fmt, actual)
                assert all(s == SEQ and q == QUAL for _, s, q in output), actual
                assert any(h.startswith('@A') and h.endswith(f'a/{mate}') for h, _, _ in output), actual
                assert any(h.startswith('@B') and h.endswith(f'b/{mate}') for h, _, _ in output), actual
                for sample in ('A', 'B'):
                    split = reads(root / f'tags_{dual}_{fmt}' / 'demux' / f'{sample}.{mate}.fq')
                    assert len(split) == 1 and split[0][0] == f'@{sample.lower()}/{mate}' and split[0][1:] == (SEQ, QUAL), split
        # Dual-index BC tags can also identify singleton records.
        tagged.write_text(HEADER + record('single', tag='\tBC:Z:' + bc_a))
        actual, _ = run(root, f'tags_single_{dual}', tagged,
                        extra=['-map', mapping, '-htsBarcode', 'BC'])
        assert len(actual[0]) == 1 and actual[0][0][1:] == (SEQ, QUAL), actual

    # Inline barcode mode remains the default and physically removes the barcode.
    mapping = root / 'inline.tsv'
    mapping.write_text('#SampleID\tBarcodeSequence\nA\tAAAA\nB\tCCCC\n')
    inline = root / 'inline.sam'
    inline.write_text(HEADER + record('inline', seq='AAAA' + SEQ, qual='IIII' + QUAL,
                                      tag='\tBC:Z:CCCC'))
    actual, _ = run(root, 'inline', inline, extra=['-map', mapping])
    assert len(actual[0]) == 1 and actual[0][0][0].startswith('@A') and actual[0][0][1:] == (SEQ, QUAL), actual

    # Keep compatibility with alignment paths in the legacy fastqFile column.
    for sample in ('A', 'B'):
        source = root / f'{sample}.sam'
        source.write_text(HEADER + record(sample))
        convert(source, root / f'{sample}.bam', 'bam')
    for column in ('fastqFile',):
        mapping = root / f'{column}.tsv'
        mapping.write_text(f'#SampleID\t{column}\nA\tA.bam\nB\tB.bam\n')
        actual, _ = run(root, 'files_' + column, extra=['-i_path', str(root) + '/', '-map', mapping])
        assert len(actual[0]) == 2 and all(s == SEQ for _, s, _ in actual[0]), actual
        assert {h[:2] for h, _, _ in actual[0]} == {'@A', '@B'}, actual

    test_mapping(root)

    bad = root / 'bad.sam'
    bad.write_text(HEADER.replace('queryname', 'coordinate') + record('a', 65) + record('a', 129))
    _, log = run(root, 'bad_sort', bad, paired=2, valid=False)
    assert 'query-name sorted' in log
    bad.write_text(HEADER + record('a', qual='*'))
    _, log = run(root, 'bad_quality', bad, valid=False, extra=['-ignore_IO_errors', '1'])
    assert 'Missing or out-of-range QUAL' in log
    # Explicit input aliases also route through the common format detector.
    actual, _ = run(root, 'explicit', extra=['-i_bam', root / 'A.bam'])
    assert actual[0] == [('@A', SEQ, QUAL)]


if __name__ == '__main__':
    with tempfile.TemporaryDirectory(prefix='sdm-hts-pipeline-') as tmp:
        tests(Path(tmp))
    print('All HTS pipeline tests passed')
