#!/usr/bin/env python3
"""Prefix dereplication integration tests; run with /path/to/sdm.

Fixtures and outputs live in a temporary directory. Expected clusters are
computed by a deliberately simple pairwise reference, independent of the C++
sorted-prefix implementation.
"""
from pathlib import Path
import random
import re
import shutil
import tempfile

from test_thorough_audit import fastq, mapping, records, run, RC


def cluster_counts(path):
    return {seq: int(re.search(r';size=(\d+);', name).group(1))
            for name, seq, _ in records(path)}


def reference(reads):
    unique = set(reads)
    anchors = sorted(s for s in unique if not any(t != s and t.startswith(s) for t in unique))
    expected = dict.fromkeys(anchors, 0)
    for sequence in reads:
        expected[next(s for s in anchors if s.startswith(sequence))] += 1
    return expected


def options(d, minimum, truncate):
    path = d / 'options.txt'
    path.write_text(f'minSeqLength\t{minimum}\nTruncateSequenceLength\t{truncate}\n'
                    'maxSeqLength\t1000\nminAvgQuality\t0\nQualWindowThreshhold\t0\n'
                    'TrimWindowThreshhold\t0\nmaxAccumulatedError\t-1\nBinErrorModelAlpha\t-1\n'
                    'maxAmbiguousNT\t0\nmaxHomonucleotide\t0\ntrimHomonucleotide\t0\n')
    return ['-options', path]


def arguments(d):
    return ['-o_dereplicate', d / 'derep.fq', '-derep_format', 'fq', '-min_derep_copies', '1']


def randomized(d, threads, cap):
    rng = random.Random(461)
    common = ''.join(rng.choices('ACGT', k=50))
    templates = [common + ''.join(rng.choices('ACGT', k=70)) for _ in range(12)]
    reads = [s[:n] for s in templates for n in [40, 50, 60, 80, 110, 120]] * 9
    rng.shuffle(reads)
    fastq(d / 'in.fq.gz', [(s, 40) for s in reads]); mapping(d)
    args = arguments(d)
    if cap:
        args += options(d, 40, 80)
    console = run(d, args, threads=threads)
    expected = reference([s[:80] if cap else s for s in reads])
    observed = cluster_counts(d / 'derep.fq')
    assert observed == expected, (observed, expected)
    assert '100% prefix' in console, console
    assert sum(observed.values()) == len(reads)
    # An explicit compatibility override retains length-matched keys.
    run(d, args + ['-derepPrefix', '0', '-derepSrchLen', '-1'], threads=threads)
    assert len(cluster_counts(d / 'derep.fq')) > len(expected)


def merged_clustering_cap(d, threads):
    rng = random.Random(91)
    shared = ''.join(rng.choices('ACGT', k=120))
    inserts = [shared + ''.join(rng.choices('ACGT', k=40)) for _ in range(2)]
    # Identical R1 clustering sequences, distinct full merged products.
    fastq(d / 'in.1.fq.gz', [(s[:100], 40) for s in inserts], 1)
    fastq(d / 'in.2.fq.gz', [(RC(s[60:]), 40) for s in inserts], 2)
    mapping(d, True)
    args = arguments(d) + options(d, 80, 80) + ['-merge_pairs_derep', '1', '-merge_pairs_filter', '1', '-derepSrchLen', '-1']
    console = run(d, args, paired=2, threads=threads)
    assert 'exact length (search cap 80)' in console, console
    merged = records(d / 'derep.merg.fq')
    assert len(merged) == 1 and ';size=2;' in merged[0][0], merged
    assert merged[0][1] in inserts, 'full merged representative must remain available'
    hq = records(d / 'derep.1.hq.fq')
    assert len(hq) == 1 and len(hq[0][1]) == 100, hq
    # Legacy full-key mode demonstrates the duplicate downstream prefixes.
    run(d, args + ['-derepPrefix', '0'], paired=2, threads=threads)
    legacy = records(d / 'derep.merg.fq')
    assert len(legacy) == 2 and len({s[:80] for _, s, _ in legacy}) == 1, legacy



def merged_prefix_recovery(d):
    rng = random.Random(612)
    sequence = ''.join(rng.choices('ACGT', k=160))
    fastq(d / 'in.1.fq.gz', [(sequence[:100], 40), (sequence[:60], 20)], 1)
    fastq(d / 'in.2.fq.gz', [(RC(sequence[60:]), 40), (RC(sequence[20:120]), 40)], 2)
    mapping(d, True)
    args = arguments(d) + options(d, 50, -1) + ['-merge_pairs_derep', '1', '-merge_pairs_filter', '1']
    path = d / 'options.txt'
    path.write_text(path.read_text().replace('minAvgQuality\t0', 'minAvgQuality\t27'))
    console = run(d, args, paired=2)
    merged = records(d / 'derep.merg.fq')
    assert len(merged) == 1 and merged[0][1] == sequence and ';size=2;' in merged[0][0], merged
    assert '1 non-passing reads recovered' in console, console
    assert merged[0][2][:20] == [30] * 20, merged
    assert merged[0][2][120:] == [40] * 40, 'uncovered tails must retain their original quality'


def main():
    root = Path(tempfile.mkdtemp(prefix='sdm-prefix-pipeline-'))
    try:
        for threads in [1, 4]:
            for cap in [False, True]:
                d = root / f'random_{threads}_{cap}'; d.mkdir()
                randomized(d, threads, cap)
            d = root / f'merged_{threads}'; d.mkdir()
            merged_clustering_cap(d, threads)
        d = root / 'merged_recovery'; d.mkdir()
        merged_prefix_recovery(d)
        d = root / 'invalid'; d.mkdir()
        fastq(d / 'in.fq.gz', [('ACGT' * 20, 40)]); mapping(d)
        console = run(d, arguments(d) + ['-derepPrefix', 'invalid'], valid=False)
        assert '-derepPrefix must be auto, 0, or 1' in console, console
    except BaseException:
        print(f'Fixtures retained at {root}')
        raise
    shutil.rmtree(root)
    print('Prefix dereplication pipeline regressions passed')


if __name__ == '__main__':
    main()
