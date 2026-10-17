#!/usr/bin/env python3
"""Regressions for subgroup counts, merged qualities, options and output failures.

Usage: python3 tests/test_thorough_audit.py /path/to/sdm
All fixtures are synthetic and use gzip input. /dev/full cases run on Linux.
"""
import gzip
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

SDM = str(Path(sys.argv[1] if len(sys.argv) > 1 else './sdm').resolve())
RNG = random.Random(73521)
REFERENCE = ''.join(RNG.choice('ACGT') for _ in range(130))
RC = lambda seq: seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def fastq(path, reads, mate=None):
    with gzip.open(path, 'wt') as out:
        for i, (sequence, quality) in enumerate(reads):
            name = f'r{i}' + (f'/{mate}' if mate else '')
            out.write(f'@{name}\n{sequence}\n+\n' + chr(33 + quality) * len(sequence) + '\n')


def records(path):
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0, path
    result = []
    for i in range(0, len(lines), 4):
        assert lines[i].startswith('@') and lines[i + 2] == '+', path
        assert len(lines[i + 1]) == len(lines[i + 3]), path
        result.append((lines[i], lines[i + 1], [ord(q) - 33 for q in lines[i + 3]]))
    return result


def run(d, extra=(), paired=1, threads=4, valid=True):
    output = str(d / 'out.fq') if paired == 1 else f'{d}/out.1.fq,{d}/out.2.fq'
    args = [SDM,
            '-paired', str(paired), '-threads', str(threads), '-i_qual_offset', '33',
            '-o_fastq', output, '-log', str(d / 'filter.log'), *map(str, extra)]
    if (d / 'map.tsv').exists():
        args += ['-i_path', str(d) + '/', '-map', str(d / 'map.tsv')]
    env = dict(os.environ)
    env.setdefault('ASAN_OPTIONS', 'detect_leaks=1:halt_on_error=1')
    env.setdefault('UBSAN_OPTIONS', 'halt_on_error=1:print_stacktrace=1')
    result = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, timeout=90, env=env)
    (d / 'command.json').write_text(json.dumps(args, indent=2))
    (d / 'console.log').write_text(result.stdout)
    assert (result.returncode == 0) == valid, (args, result.returncode, result.stdout[-8000:])
    assert 'ERROR: AddressSanitizer' not in result.stdout and 'runtime error:' not in result.stdout, result.stdout
    return result.stdout


def mapping(d, paired=False):
    files = 'in.1.fq.gz,in.2.fq.gz' if paired else 'in.fq.gz'
    (d / 'map.tsv').write_text('#SampleID\tfastqFile\nS\t' + files + '\n')


def derep_args(d, merged=False, filename='derep.fq'):
    args = ['-o_dereplicate', d / filename, '-derep_format', 'fq',
            '-derepSrchLen', '20', '-min_derep_copies', '1', '-suppressOutput', '3']
    if merged:
        args += ['-merge_pairs_derep', '1', '-merge_pairs_filter', '1']
    return args


def merged_quality(d, reverse=False, threads=4, merge_key=True):
    quals = [40, 20] if reverse else [20, 40]
    fastq(d / 'in.1.fq.gz', [(REFERENCE[:40], q) for q in quals], 1)
    fastq(d / 'in.2.fq.gz', [(RC(REFERENCE[10:]), q) for q in quals], 2)
    mapping(d, True)
    args = derep_args(d, True)
    if not merge_key:
        args += ['-noSearchWithMerge', '1']
    run(d, args, 2, threads)
    output = records(d / 'derep.merg.fq')
    assert len(output) == 1 and ';size=2;' in output[0][0], output
    assert output[0][1] == REFERENCE, output
    # Outside the overlap, consensus qualities equal the input qualities.
    assert output[0][2][:10] == [30] * 10 and output[0][2][40:] == [30] * 90, output


def subgroup_counts(d, threads):
    unmergeable = 'TTCGTTGAGGTGAGGATAGTGACCTATCCATTTTTACGGGACATCGCATTCGGAAAGTGAGTGAACTACCTGCATTACATGATTTTTAGTTCGATATTAGTTTCTCCTACTAACGGTTGAC'
    fastq(d / 'in.1.fq.gz', [(REFERENCE[:40], 40)] * 11, 1)
    fastq(d / 'in.2.fq.gz', [(unmergeable, 40)] * 10 + [(RC(REFERENCE[10:]), 40)], 2)
    mapping(d, True)
    run(d, derep_args(d, True), 2, threads)
    output = records(d / 'derep.fq') + records(d / 'derep.merg.fq')
    assert len(output) == 1 and ';size=11;' in output[0][0], output
    assert '\t0:11' in (d / 'derep.map').read_text()


def minimum_and_cap(d, reversed_order):
    fastq(d / 'in.fq.gz', [('ACGT' * 40, 40)])
    mapping(d)
    settings = ['minSeqLength\t100\n', 'TruncateSequenceLength\t70\n']
    (d / 'options.txt').write_text(''.join(reversed(settings) if reversed_order else settings))
    run(d, ['-options', d / 'options.txt'])
    output = records(d / 'out.fq')
    assert len(output) == 1 and output[0][1] == ('ACGT' * 40)[:70], output
    assert 'Sequence length below 70.00 bp' in (d / 'filter.log').read_text()


def start_trim(d):
    fastq(d / 'in.fq.gz', [('ACGT' * 40, 40)])
    # Omit -map for no-map mode, where TrimStartNTs applies.
    (d / 'options.txt').write_text('minSeqLength\t100\nTruncateSequenceLength\t-1\nTrimStartNTs\t70\n')
    run(d, ['-i_fastq', d / 'in.fq.gz', '-options', d / 'options.txt'])
    assert records(d / 'out.fq') == []


def option_precedence(d, mode):
    fastq(d / 'in.fq.gz', [('ACGT' * 25, 40), ('ACGT' * 5 + 'TGCA' * 20, 40)])
    mapping(d)
    settings = 'minSeqLength\t100\nTruncateSequenceLength\t-1\n'
    if mode in ('both', 'file'):
        settings += 'derepSrchLen\t30\n'
    (d / 'options.txt').write_text(settings)
    args = derep_args(d)
    if mode in ('file', 'default'):
        i = args.index('-derepSrchLen'); del args[i:i + 2]
    if mode == 'negative':
        args[args.index('-derepSrchLen') + 1] = '-1'
    run(d, args + ['-options', d / 'options.txt'])
    output = records(d / 'derep.fq')
    assert len(output) == (1 if mode in ('cli', 'both') else 2), output


def finalization(d, target):
    if not Path('/dev/full').exists():
        return
    fastq(d / 'in.fq.gz', [('ACGT' * 20, 40)] * 2)
    mapping(d)
    path = d / ('derep.map.tmp' if target == 'map' else 'derep.fq')
    path.symlink_to('/dev/full')
    console = run(d, derep_args(d), valid=False)
    assert 'finalize output file' in console, console
    if target == 'map':
        output = d / 'derep.map'
        assert output.is_file() and not output.is_symlink()
        assert '\t0:2' in output.read_text(), output.read_text()


def extensionless(d):
    fastq(d / 'in.fq.gz', [('ACGT' * 20, 40)] * 2)
    mapping(d)
    run(d, derep_args(d, filename='derep'))
    assert len(records(d / 'derep')) == 1
    assert (d / 'derep.map').is_file()
    assert (d / 'derep.1.hq.fq').is_file()


def goldenaxe(d, mode):
    forwards = ['AGGTCAGTACCGTAAC', 'CATGTCGACATAGTCA']
    reverse = 'TCCAGATGCTACGTCA'
    barcodes = ['GACTGA', 'CTCGTG']
    rows = ['#SampleID\tfastqFile\tBarcodeSequence\tLinkerPrimerSequence\tReversePrimer']
    combined = []
    for i in range(2):
        bc = '' if mode == 'multifile' else barcodes[i]
        sequence = bc + forwards[i] + 'ACGT' * 30 + RC(reverse)
        file = f'in{i}.fq.gz' if mode == 'multifile' else 'in.fq.gz'
        rows.append(f'{chr(65+i)}\t{file}\t{bc}\t{forwards[i]}\t{reverse}')
        if mode == 'multifile':
            fastq(d / file, [(sequence, 40)])
        else:
            combined.append((sequence, 40))
    if combined:
        fastq(d / 'in.fq.gz', combined)
    (d / 'map.tsv').write_text('\n'.join(rows) + '\n')
    (d / 'options.txt').write_text('minSeqLength\t100\nTruncateSequenceLength\t-1\n')
    args = ['-options', d / 'options.txt', '-GoldenAxe', '1',
            '-GoldenAxeMinAmpli', '-1', '-GoldenAxeMaxAmpli', '-1']
    if mode == 'multifile':
        args += derep_args(d)
    run(d, args)
    if mode == 'multifile':
        text = (d / 'derep.map').read_text()
        assert '#SMPLS\t0:A\t1:B' in text and '\t0:1' in text and '\t1:1' in text, text
    else:
        assert len(records(d / 'out.fq')) == 2


def main():
    cases = [('subgroups_4', lambda d: subgroup_counts(d, 4)),
             ('subgroups_12', lambda d: subgroup_counts(d, 12)),
             ('merged_qualities_4', lambda d: merged_quality(d)),
             ('merged_qualities_12', lambda d: merged_quality(d, True, 12)),
             ('merged_r1_key', lambda d: merged_quality(d, merge_key=False)),
             ('minimum_then_cap', lambda d: minimum_and_cap(d, False)),
             ('cap_then_minimum', lambda d: minimum_and_cap(d, True)),
             ('start_trim', start_trim), ('close_failure', lambda d: finalization(d, 'reads')),
             ('map_failure', lambda d: finalization(d, 'map')),
             ('no_extension', extensionless), ('dotted.directory', extensionless),
             ('goldenaxe_multifile', lambda d: goldenaxe(d, 'multifile')),
             ('goldenaxe_primers', lambda d: goldenaxe(d, 'primers'))]
    cases += [(f'precedence_{mode}', lambda d, m=mode: option_precedence(d, m))
              for mode in ('cli', 'both', 'file', 'default', 'negative')]
    with tempfile.TemporaryDirectory(prefix='sdm-thorough-regressions-') as root:
        for name, case in cases:
            d = Path(root) / name
            d.mkdir()
            case(d)
            print('PASS', name, flush=True)
    print(f'Passed {len(cases)} thorough audit pipeline cases')


if __name__ == '__main__':
    main()
