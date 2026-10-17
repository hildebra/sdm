#!/usr/bin/env python3
"""End-to-end regressions for output accounting and optional merge statistics.

Usage: python3 tests/test_audit_pipeline.py /path/to/sdm
The same cases can run against an ASan/UBSan build; inputs and outputs are temporary.
"""
import gzip
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile

SDM = str(Path(sys.argv[1] if len(sys.argv) > 1 else './sdm').resolve())
RNG = random.Random(73521)
REFERENCE = ''.join(RNG.choice('ACGT') for _ in range(130))


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def fastq(path, sequence, count, prefix='read', mate=None):
    text = ''.join(f'@{prefix}_{i}' + (f'/{mate}' if mate else '') +
                   f'\n{sequence}\n+\n' + 'I' * len(sequence) + '\n' for i in range(count))
    with gzip.open(path, 'wt') as output:
        output.write(text)


def records(path):
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0, path
    result = []
    for i in range(0, len(lines), 4):
        assert lines[i].startswith('@') and lines[i + 2] == '+', path
        assert len(lines[i + 1]) == len(lines[i + 3]), path
        result.append((lines[i], lines[i + 1], lines[i + 3]))
    return result


def normalized_id(record):
    return record[0].split()[0].lstrip('@').split('___')[-1].removesuffix('/1').removesuffix('/2')


def run(directory, extra, paired=1, threads=4, valid=True):
    output = directory / 'out.fq' if paired == 1 else f'{directory}/out.1.fq,{directory}/out.2.fq'
    command = [SDM, '-o_fastq', str(output), '-paired', str(paired),
               '-threads', str(threads), '-i_qual_offset', '33', '-sample_sep', '___',
               '-log', str(directory / 'filter.log'), *map(str, extra)]
    environment = dict(os.environ)
    environment.setdefault('ASAN_OPTIONS', 'detect_leaks=1:halt_on_error=1')
    environment.setdefault('UBSAN_OPTIONS', 'halt_on_error=1:print_stacktrace=1')
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, env=environment, timeout=120)
    (directory / 'command.log').write_text(result.stdout)
    assert (result.returncode == 0) == valid, (command, result.returncode, result.stdout[-8000:])
    assert 'ERROR: AddressSanitizer' not in result.stdout and 'runtime error:' not in result.stdout, result.stdout
    return result


def rotation(directory, paired, threads):
    rows = ['#SampleID\tfastqFile']
    for sample in ('A', 'B'):
        paths = []
        for mate in range(1, paired + 1):
            path = directory / f'{sample}.{mate}.fq.gz'
            fastq(path, 'ACGT', 4, sample, mate if paired == 2 else None)
            paths.append(path.name)
        rows.append(sample + '\t' + ','.join(paths))
    mapping = directory / 'map.tsv'
    mapping.write_text('\n'.join(rows) + '\n')
    run(directory, ['-map', mapping, '-i_path', str(directory) + '/',
                    '-maxReadsPerOutput', '3'], paired, threads)
    all_records = []
    for generation, expected in enumerate((3, 3, 2)):
        suffix = '' if generation == 0 else f'.{generation}'
        first = records(directory / (f'out{suffix}.fq' if paired == 1 else f'out.1{suffix}.fq'))
        assert len(first) == expected, (generation, expected, first)
        if paired == 2:
            second = records(directory / f'out.2{suffix}.fq')
            assert list(map(normalized_id, first)) == list(map(normalized_id, second))
        all_records.extend(first)
    assert sorted(map(normalized_id, all_records)) == [f'{s}_{i}' for s in ('A', 'B') for i in range(4)]


def demultiplex(directory, paired, limit, expected, merged=False, threads=12):
    paths = []
    for mate in range(1, paired + 1):
        path = directory / f'input.{mate}.fq.gz'
        sequence = (REFERENCE[:40] if mate == 1 else reverse_complement(REFERENCE[10:])) if merged else 'ACGT'
        fastq(path, sequence, 20, mate=mate if paired == 2 else None)
        paths.append(path.name)
    mapping = directory / 'map.tsv'
    mapping.write_text('#SampleID\tfastqFile\nA\t' + ','.join(paths) + '\n')
    extra = ['-map', mapping, '-i_path', str(directory) + '/',
             '-o_demultiplex', str(directory) + '/', '-suppressOutput', '3',
             '-DemultiBPperSR', str(limit)]
    if merged:
        extra += ['-merge_pairs_demulti', '1']
    run(directory, extra, paired, threads)
    if merged:
        output = records(directory / 'A.merg.fq')
        assert len(output) == expected and all(r[1] == REFERENCE for r in output), output
        assert not records(directory / 'A.1.fq') and not records(directory / 'A.2.fq')
    else:
        first = records(directory / ('A.fq' if paired == 1 else 'A.1.fq'))
        assert len(first) == expected, first
        if paired == 2:
            second = records(directory / 'A.2.fq')
            assert list(map(normalized_id, first)) == list(map(normalized_id, second))


def statistics(directory, derep, threads):
    paths = []
    for mate, sequence in ((1, REFERENCE[:40]), (2, reverse_complement(REFERENCE[10:]))):
        path = directory / f'input.{mate}.fq.gz'
        fastq(path, sequence, 20, mate=mate)
        paths.append(path.name)
    mapping = directory / 'map.tsv'
    mapping.write_text('#SampleID\tfastqFile\nA\t' + ','.join(paths) + '\n')
    merge_path, qual_path = directory / 'requested_merge.tsv', directory / 'requested_quality.tsv'
    extra = ['-map', mapping, '-i_path', str(directory) + '/',
             '-merg_readpos', merge_path, '-qual_readpos', qual_path]
    if derep:
        extra += ['-o_dereplicate', directory / 'derep.fq', '-derep_format', 'fq',
                  '-merge_pairs_derep', '1', '-merge_pairs_filter', '1',
                  '-min_derep_copies', '1', '-derepSrchLen', '20']
    else:
        extra += ['-o_demultiplex', str(directory) + '/', '-merge_pairs_demulti', '1']
    run(directory, extra, 2, threads)
    assert merge_path.exists() and qual_path.exists(), list(directory.iterdir())
    merge_rows = [row.split('\t') for row in merge_path.read_text().splitlines()[1:]]
    observations = sum(float(row[i]) for row in merge_rows for i in (2, 3, 5, 6) if row[i])
    assert observations == 20 * 30, observations
    assert all(not row[1] and not row[2] and not row[3] for row in merge_rows), merge_rows
    quality_rows = [row.split('\t') for row in qual_path.read_text().splitlines()[1:]]
    assert len(quality_rows) == 120
    assert all(float(value) == 40 for row in quality_rows for value in row[1:] if value)
    if derep:
        merged = records(directory / 'derep.merg.fq')
        assert len(merged) == 1 and merged[0][1] == REFERENCE and 'size=20;' in merged[0][0], merged


def statistics_run_reset(directory):
    rows = ['#SampleID\tfastqFile\tSequencingRun']
    for run_id, count in (('SR1', 5), ('SR2', 7)):
        paths = []
        for mate, sequence in ((1, REFERENCE[:40]), (2, reverse_complement(REFERENCE[10:]))):
            path = directory / f'{run_id}.{mate}.fq.gz'
            fastq(path, sequence, count, prefix=run_id, mate=mate)
            paths.append(path.name)
        rows.append(run_id + '\t' + ','.join(paths) + '\t' + run_id)
    mapping = directory / 'map.tsv'
    mapping.write_text('\n'.join(rows) + '\n')
    run(directory, ['-map', mapping, '-i_path', str(directory) + '/',
                    '-o_dereplicate', directory / 'derep.fq', '-derep_format', 'fq',
                    '-derepPerSR', '1', '-merge_pairs_derep', '1', '-merge_pairs_filter', '1',
                    '-derepSrchLen', '20', '-min_derep_copies', '1',
                    '-merg_readpos', directory / 'merge.tsv',
                    '-qual_readpos', directory / 'quality.tsv'], 2, 12)
    for run_id, count in (('SR1', 5), ('SR2', 7)):
        rows = [row.split('\t') for row in (directory / f'merge.{run_id}.tsv').read_text().splitlines()[1:]]
        observations = sum(float(row[i]) for row in rows for i in (2, 3, 5, 6) if row[i])
        assert observations == count * 30, (run_id, observations)
        assert (directory / f'quality.{run_id}.tsv').exists()


def preflight(directory):
    previous = directory / 'derep.fq'
    sentinel = b'PREVIOUS DEREPLICATION RESULT\n'
    previous.write_bytes(sentinel)
    run(directory, ['-i_fastq', directory / 'missing.fq', '-o_dereplicate', previous,
                    '-derep_format', 'fq'], valid=False)
    assert previous.read_bytes() == sentinel


def invalid_option(directory, option, value):
    input_path = directory / 'input.fq.gz'
    fastq(input_path, 'ACGT', 1)
    result = run(directory, ['-i_fastq', input_path, option, value], valid=False)
    assert option in result.stdout and 'nonnegative whole number' in result.stdout, result.stdout


def main():
    # Invalid-input subprocesses must not leave platform-dependent core dumps.
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_CORE, (0, 0))
    except ImportError:
        pass
    cases = [('preflight', preflight, ()), ('statistics_run_reset', statistics_run_reset, ())]
    cases += [(f'rotation_{p}_{t}', rotation, (p, t)) for p in (1, 2) for t in (1, 4, 12)]
    cases += [(f'demulti_{p}_{limit}', demultiplex, (p, limit, expected)) for p, limit, expected in
              ((1, '8e0', 2), (1, 0, 20), (1, 8, 2), (1, 9, 2), (2, 0, 20), (2, 16, 2), (2, 17, 2), (2, 7, 0))]
    cases += [(f'merged_quota_{limit}', demultiplex, (2, limit, expected, True)) for limit, expected in ((0, 20), (260, 2), (259, 1))]
    cases += [(f'statistics_{d}_{t}', statistics, (d, t)) for d in (False, True) for t in (4, 12)]
    cases += [(f'invalid_{i}', invalid_option, pair) for i, pair in enumerate([
        ('-DemultiBPperSR', '-1'), ('-DemultiBPperSR', 'nan'), ('-DemultiBPperSR', 'inf'),
        ('-DemultiBPperSR', '1.5'), ('-DemultiBPperSR', '4294967296'), ('-DemultiBPperSR', '4junk'),
        ('-maxReadsPerOutput', '-1'), ('-maxReadsPerOutput', '18446744073709551616'),
        ('-maxReadsPerOutput', '3junk'), ('-XfirstReadsWritten', '-1'), ('-XfirstReadsRead', '-1')])]
    with tempfile.TemporaryDirectory(prefix='sdm-audit-pipeline-') as temp:
        for name, function, args in cases:
            directory = Path(temp) / name
            directory.mkdir()
            function(directory, *args)
            print(f'PASS {name}', flush=True)
    print(f'Passed {len(cases)} audit pipeline cases')


if __name__ == '__main__':
    main()
