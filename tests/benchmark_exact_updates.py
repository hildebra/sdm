#!/usr/bin/env python3
"""Compare real 100% CLI paths with alternating binaries and validated abundances.

Input generation and output validation are outside process timings. Default
cases omit every derepIdentity/derepPrefix/derepStore* override. Optional quality
and diff cases explicitly enable those features, which select a different backend.
"""
import argparse
from collections import Counter
import csv
import hashlib
from itertools import zip_longest
import json
import os
from pathlib import Path
import random
import re
import resource
import shutil
import subprocess
import tempfile
import time

from test_diff_quality_pipeline import fastq, retained_fastq_paths
from test_paired_diff_quality_pipeline import fasta

CASES = ('single-exact', 'paired-exact', 'paired-prefix', 'paired-quals', 'paired-diffs')


def fixture(root, case, reads, templates, length):
    root.mkdir()
    rng = random.Random(2026091408)
    prefix = case == 'paired-prefix'
    paired = case != 'single-exact'
    maximum = length + 55 if prefix else length
    sequences = [tuple(''.join(rng.choices('ACGT', k=maximum)) for _ in range(2))
                 for _ in range(templates)]
    assert len({pair[0][:length] for pair in sequences}) == templates
    paths = [root / f'input.{mate}.fq' for mate in range(1, 3 if paired else 2)]
    outputs = [path.open('w') for path in paths]
    expected = Counter()
    expected_pairs = Counter()
    try:
        for i in range(reads):
            # First include every full-length anchor; later prefix observations
            # have exactly one compatible anchor, independent of worker order.
            pair = sequences[i] if i < templates else rng.choice(sequences)
            n = rng.randint(length, maximum) if prefix and i >= templates else maximum
            expected[pair[0]] += 1
            expected_pairs[pair] += 1
            for mate, output in enumerate(outputs):
                sequence = pair[mate][:n]
                output.write(f'@read{i}/{mate + 1}\n{sequence}\n+\n' + '?' * n + '\n')
    finally:
        for output in outputs:
            output.close()
    settings = {'minSeqLength': length, 'maxSeqLength': maximum,
                'TruncateSequenceLength': -1 if prefix else length,
                'minAvgQuality': 0, 'QualWindowThreshhold': 0, 'TrimWindowThreshhold': 0,
                'maxAccumulatedError': -1, 'BinErrorModelAlpha': -1,
                'maxAmbiguousNT': 0, 'maxHomonucleotide': 0, 'trimHomonucleotide': 0,
                'fastqVersion': 1}
    config = root / 'options.txt'
    config.write_text(''.join(f'{key}\t{value}\n' for key, value in settings.items()))
    mapping = root / 'map.tsv'
    mapping.write_text('#SampleID\tfastqFile\nS0\t' + ','.join(path.name for path in paths) + '\n')
    stem = root / 'out'
    command = ['-i_path', str(root), '-map', str(mapping), '-options', str(config),
               '-paired', '2' if paired else '1', '-o_dereplicate', str(stem) + '.fna',
               '-o_fna', ','.join(str(stem) + f'.{mate}.fna' for mate in range(1, len(paths) + 1)),
               '-sample_sep', '___', '-log', str(stem) + '.log',
               '-merge_pairs_derep', '0', '-merge_pairs_filter', '0', '-merge_pairs_demulti', '0',
               '-min_derep_copies', '1', '-derepPerSR', '0', '-suppressOutput', '1',
               '-i_qual_offset', '33', '-o_qual_offset', '33', '-oneLineFastaFormat', '1']
    if case == 'paired-quals':
        command += ['-derepStoreQuals', '1']
    elif case == 'paired-diffs':
        command += ['-derepStoreDiffs', '1', '-derepSubclusterFasta', '1']
    return command, expected, expected_pairs


def validate(root, case, expected, expected_pairs):
    stem = str(root / 'out')
    records = fasta(Path(stem + '.fna'))
    observed = Counter()
    for name, sequence in records:
        observed[sequence] += int(name.split(';size=')[1].split(';')[0])
    assert len(records) == len(expected) and observed == expected, 'Representative sequence/count mismatch'
    lines = Path(stem + '.map').read_text().splitlines()
    total = sum(int(entry.split(':')[1]) for line in lines if not line.startswith('#')
                for entry in line.split('\t')[1:] if ':' in entry)
    assert total == sum(expected.values()), 'Sample map total differs'
    if case == 'paired-quals':
        actual_pairs = Counter()
        first, second = retained_fastq_paths(stem, paired=True)
        for a, b in zip_longest(fastq(first), fastq(second)):
            assert a is not None and b is not None and a[0] == b[0] and a[3] == b[3]
            assert a[2] == '?' * len(a[1]) and b[2] == '?' * len(b[1]), 'Quality witness changed'
            actual_pairs[a[1], b[1]] += a[3]
        assert actual_pairs == expected_pairs, 'Retained exact pair/count mismatch'
    elif case == 'paired-diffs':
        variants = fasta(Path(stem + '.subclusters.fna'))
        counts = Counter()
        for name, sequence in variants:
            counts[sequence] += int(name.split(';size=')[1].split(';')[0])
        assert len(variants) == len(expected) and counts == expected, 'Diff variant counts differ'
        assert Path(stem + '.diff').read_bytes().startswith(b'SDMDIFF')
    return hashlib.sha256(json.dumps(sorted(observed.items())).encode()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', action='append', required=True, help='LABEL=/path/to/sdm; supply at least twice')
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES[:3]))
    parser.add_argument('--reads', type=int, default=500000)
    parser.add_argument('--templates', type=int, default=10000)
    parser.add_argument('--length', type=int, default=200)
    parser.add_argument('--threads', nargs='+', type=int, default=[1, 12])
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not 1 <= args.templates <= args.reads or args.length < 50 or args.repeats < 1 or min(args.threads) < 1:
        parser.error('require 1 <= templates <= reads, length >= 50, and positive repeats/threads')
    binaries = [(label, str(Path(path).resolve())) for label, path in (item.split('=', 1) for item in args.binary)]
    if len(binaries) < 2 or len({label for label, _ in binaries}) != len(binaries):
        parser.error('supply at least two binaries with unique labels')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    log_dir = args.output.with_suffix('')
    log_dir.mkdir(exist_ok=True)
    cpus = sorted(os.sched_getaffinity(0))
    metadata = dict(binaries={label: dict(path=path, sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest())
                              for label, path in binaries}, arguments=vars(args), cpus=cpus)
    args.output.with_suffix('.json').write_text(json.dumps(metadata, indent=2, default=str) + '\n')
    fields = ['repeat', 'case', 'version', 'threads', 'reads', 'templates', 'length', 'mode',
              'wall_seconds', 'cpu_seconds', 'cpu_percent', 'peak_rss_bytes', 'canonical_sha256']
    with tempfile.TemporaryDirectory(prefix='sdm-exact-updates-') as temp, args.output.open('w', newline='') as output:
        root = Path(temp)
        # Keep executable loading on the same filesystem (e.g. WSL /tmp versus
        # a Windows checkout), and retain original paths/hashes in metadata.
        staged = []
        for index, (label, binary) in enumerate(binaries):
            target = root / f'sdm-{index}'
            shutil.copy2(binary, target)
            staged.append((label, str(target)))
        binaries = staged
        fixtures = {case: fixture(root / case, case, args.reads, args.templates, args.length) for case in args.cases}
        writer = csv.DictWriter(output, fieldnames=fields)
        writer.writeheader()
        for repeat in range(1, args.repeats + 1):
            workloads = [(case, threads) for case in args.cases for threads in args.threads]
            random.Random(20260914 + repeat).shuffle(workloads)
            versions = binaries if repeat % 2 else list(reversed(binaries))
            for case, threads in workloads:
                command, expected, expected_pairs = fixtures[case]
                affinity = str(cpus[0]) if threads == 1 else ','.join(map(str, cpus))
                for label, binary in versions:
                    invocation = ['taskset', '-c', affinity, binary, *command, '-threads', str(threads)]
                    before = resource.getrusage(resource.RUSAGE_CHILDREN)
                    start = time.perf_counter()
                    process = subprocess.run(invocation, text=True, capture_output=True, timeout=300)
                    wall = time.perf_counter() - start
                    after = resource.getrusage(resource.RUSAGE_CHILDREN)
                    cpu = after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
                    log = log_dir / f'{repeat}-{case}-{threads}-{label}.log'
                    log.write_text(json.dumps(invocation) + '\n' + process.stdout + process.stderr)
                    assert process.returncode == 0, f'Process failed: {log}'
                    mode = re.search(r'Dereplication matching: (.+)', process.stderr)
                    assert mode, f'Missing mode: {log}'
                    required = 'streaming substitution' if case in CASES[3:] else ('100% prefix' if case == 'paired-prefix' else 'exact length')
                    assert required in mode[1], f'Unexpected backend: {mode[1]}'
                    checksum = validate(root / case, case, expected, expected_pairs)
                    memory = re.search(r'Peak process memory: (\d+) bytes', process.stderr)
                    row = dict(repeat=repeat, case=case, version=label, threads=threads,
                               reads=args.reads, templates=args.templates, length=args.length, mode=mode[1],
                               wall_seconds=wall, cpu_seconds=cpu, cpu_percent=100 * cpu / wall,
                               peak_rss_bytes=int(memory[1]) if memory else 0, canonical_sha256=checksum)
                    writer.writerow(row)
                    output.flush()
                    print(f'{repeat} {case} {threads}t {label}: {wall:.4f}s, {100*cpu/wall:.1f}% CPU', flush=True)
                    for path in (root / case).glob('out*'):
                        path.unlink()


if __name__ == '__main__':
    main()
