#!/usr/bin/env python3
"""Validate exact paired variants and selected fragment qualities through real FASTQ IO."""
from collections import Counter
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile
from test_diff_quality_pipeline import fastq


def fasta(path):
    lines = path.read_text().splitlines()
    assert len(lines) % 2 == 0
    return [(lines[i][1:], lines[i + 1]) for i in range(0, len(lines), 2)]



def assert_identity_output_equivalence(exact_stem, coarse_stem):
    """Only record order may differ; retain IDs, counts, sequences and qualities.

    Comparing multisets also catches duplicate records that a set or a total
    abundance check would hide. Diagnose ID-only differences separately without
    accepting them as equal files.
    """
    differences = []
    for suffix, reader in (('.fna', fasta), ('.fna.rest', fasta),
                           ('.1.hq.fq', fastq), ('.2.hq.fq', fastq)):
        kwargs = {'allow_unannotated': True} if reader is fastq else {}
        exact = list(reader(Path(str(exact_stem) + suffix), **kwargs))
        coarse = list(reader(Path(str(coarse_stem) + suffix), **kwargs))
        if Counter(exact) == Counter(coarse):
            print(f'PASS identity equivalence {suffix}: {len(exact)} records', flush=True)
            continue
        def content(record):
            if reader is fastq:
                return record[1:]  # sequence, quality, abundance
            return record[1], int(record[0].split(';size=')[1].split(';')[0])
        same_content = Counter(map(content, exact)) == Counter(map(content, coarse))
        difference = 'read content/counts match but identifiers differ' if same_content else 'read content or abundances differ'
        differences.append(f'{suffix}: 100% has {len(exact)} records; 97% has {len(coarse)} records; {difference}')
    exact_map = Counter(Path(str(exact_stem) + '.map').read_text().splitlines())
    coarse_map = Counter(Path(str(coarse_stem) + '.map').read_text().splitlines())
    if exact_map != coarse_map:
        differences.append('.map: exact dereplicate sample abundances or identifiers differ')
    assert not differences, '100% versus 97% output equivalence failed:\n' + '\n'.join(differences)


def thousand_pair_equivalence(root, inputs, run):
    """1,000 fragments: 250 nearby R1 keys, repeats, divergent R2s and tails.

    Abundances span both sides of the copy cutoff. Each seed winner has a
    uniquely better observed quality, so worker scheduling cannot choose a tied
    read ID. Compare whole records, then check independent input-derived truth.
    """
    rng = random.Random(2026091407)
    observations = []
    for family in range(50):
        length = (200, 201, 220, 254, 255)[family % 5]
        first = ''.join(rng.choices('ACGT', k=length))
        second = ''.join(rng.choices('ACGT', k=230))
        for variant, copies in enumerate((1, 2, 3, 4, 10)):
            r1 = first if variant == 0 else first[:30 + variant] + next(
                b for b in 'ACGT' if b != first[30 + variant]) + first[31 + variant:]
            for copy in range(copies):
                pair = r1, second
                if copy >= max(1, copies // 2):
                    # An identical effective R1 must reuse its dereplicate even
                    # though the second mate falls well outside 97% identity.
                    r2 = ''.join(next(b for b in 'ACGT' if b != base) if p < 20 else base
                                 for p, base in enumerate(second))
                    tail = r1 if length == 200 else r1[:-1] + next(b for b in 'ACGT' if b != r1[-1])
                    pair = tail, r2
                quality = 10 + 3 * copy
                observations.append((pair, quality))
    assert len(observations) == 1000
    rng.shuffle(observations)
    for mate, path in enumerate(inputs):
        path.write_text(''.join(f'@equivalence{i}/{mate + 1}\n{pair[mate]}\n+\n' +
                               chr(33 + quality) * len(pair[mate]) + '\n'
                               for i, (pair, quality) in enumerate(observations)))
    expected = Counter(pair for pair, _ in observations)
    best_quality = {pair: max(q for other, q in observations if pair == other) for pair in expected}
    assert len(expected) == 450
    for cap in (200, -1):
        key = lambda pair: pair[0][:cap] if cap > 0 else pair[0]
        counts = Counter(key(pair) for pair, _ in observations)
        winners = {}
        for pair, quality in observations:
            if key(pair) not in winners or quality > winners[key(pair)][1]:
                winners[key(pair)] = pair, quality
        for retain in (0, 1):
            baseline = root / f'equivalence-100-{cap}-{retain}'
            common = ('-derepSrchLen', str(cap), '-derepStoreQuals', str(retain),
                      '-min_derep_copies', '3', '-suppressOutput', '3')
            process = run(baseline, '-threads', '1', '-derepIdentity', '100', *common)
            assert process.returncode == 0, process.stderr
            for threads in (1, 4, 12):
                stem = root / f'equivalence-97-{cap}-{retain}-{threads}'
                process = run(stem, '-threads', str(threads), '-derepIdentity', '97', *common)
                assert process.returncode == 0, process.stderr
                assert_identity_output_equivalence(baseline, stem)
                first = list(fastq(Path(str(stem) + '.1.hq.fq')))
                second = list(fastq(Path(str(stem) + '.2.hq.fq'), allow_unannotated=not retain))
                assert len(first) == len(second) == (len(expected) if retain else len(counts))
                actual = Counter()
                for a, b in zip(first, second):
                    if retain:
                        assert a[0] == b[0] and a[3] == b[3], 'Unsynchronized paired FASTQ'
                    else:
                        assert b[0].removesuffix('/2') in a[0], 'Representative mates came from different fragments'
                    pair = a[1], b[1]
                    actual[pair] += a[3]
                    quality = best_quality[pair] if retain else winners[key(pair)][1]
                    assert (a[2], b[2]) == tuple(chr(33 + quality) * len(seq) for seq in pair)
                truth = expected if retain else Counter({winners[k][0]: count for k, count in counts.items()})
                assert actual == truth and sum(actual.values()) == 1000, 'Lost or double-counted fragments'
                for suffix, passing in (('.fna', True), ('.fna.rest', False)):
                    observed = Counter((seq, int(header.split(';size=')[1].split(';')[0]))
                                       for header, seq in fasta(Path(str(stem) + suffix)))
                    truth = Counter((winners[k][0][0], count) for k, count in counts.items()
                                    if (count >= 3) == passing)
                    assert observed == truth, 'Wrong exact R1 seed or abundance cutoff'
                assert not any(root.glob(stem.name + '.subclusters*.fq'))
                if threads == 1:
                    match = re.search(r'Coarse search groups: ([\d,]+); output exact (?:R1|search-key) dereplicates', process.stderr)
                    assert match and int(match[1].replace(',', '')) < len(counts), 'Fixture did not exercise coarse grouping'
                print(f'PASS 1,000 paired fragments: cap={cap}, retention={retain}, threads={threads}; '
                      f'{len(counts)} exact R1 dereplicates, {len(first)} synchronized FASTQ pairs', flush=True)


def main(binary):
    binary = str(Path(binary).resolve())
    rng = random.Random(2026091504)
    templates = [(''.join(rng.choices('ACGT', k=230)), ''.join(rng.choices('ACGT', k=260))) for _ in range(12)]
    pairs = []
    for i in range(240):
        first, second = templates[i % 12]
        # Many exact R1 repeats; variants differ solely in R2, including length.
        first = first[:200 + (i % 3) * 10]
        second = list(second[:210 + (i % 5) * 10])
        for _ in range(1 + i % 4):
            p = rng.randrange(len(second))
            second[p] = rng.choice([b for b in 'ACGT' if b != second[p]])
        pairs.append((first, ''.join(second)))
    observations = []
    for pair in pairs:
        for q1, q2 in ((10, 20), (30, 40), (20, 30)):
            qualities = tuple(''.join(chr(33 + q + p % 5 - 2) for p in range(len(s))) for s, q in zip(pair, (q1, q2)))
            observations.append((pair, qualities, q1))
    rng.shuffle(observations)
    expected = Counter(pair for pair, _, _ in observations)
    winners = {pair: quality for pair, quality, q in observations if q == 30}
    with tempfile.TemporaryDirectory(prefix='sdm-paired-diff-cli-') as temp:
        root = Path(temp)
        inputs = [root / f'reads.{mate}.fq' for mate in (1, 2)]
        for mate, path in enumerate(inputs):
            path.write_text(''.join(f'@fragment{i}/{mate + 1}\n{pair[mate]}\n+\n{quality[mate]}\n'
                                    for i, (pair, quality, _) in enumerate(observations)))
        settings = {'minSeqLength': 50, 'maxSeqLength': 500, 'minAvgQuality': 0,
                    'QualWindowThreshhold': 0, 'TrimWindowThreshhold': 0,
                    'maxAccumulatedError': -1, 'BinErrorModelAlpha': -1,
                    'maxAmbiguousNT': 0, 'maxHomonucleotide': 0, 'trimHomonucleotide': 0,
                    'TruncateSequenceLength': -1, 'fastqVersion': 1, 'derepStoreQuals': 1}
        config = root / 'options.txt'
        config.write_text(''.join(f'{key}\t{value}\n' for key, value in settings.items()))
        mapping = root / 'map.tsv'
        mapping.write_text(f'#SampleID\tfastqFile\nS0\t{inputs[0]},{inputs[1]}\n')
        def run(stem, *extra):
            args = [binary, '-i_fastq', ','.join(map(str, inputs)), '-paired', '2',
                    '-map', str(mapping), '-options', str(config), '-o_dereplicate', str(stem) + '.fna',
                    '-o_fastq', str(stem) + '.1.fq,' + str(stem) + '.2.fq',
                    '-threadIO', '0', '-derep_reserve', '1', '-log', str(stem) + '.log', *extra]
            return subprocess.run(args, cwd=root, text=True, capture_output=True)
        for identity, reassign, threads, cap, binary_output in [(97, 0, 1, -1, 0), (99, 0, 1, -1, 1),
                (100, 0, 1, -1, 0), (97, 1, 1, -1, 1), (97, 1, 4, 100, 0),
                (97, 0, 4, -1, 0), (97, 0, 12, 100, 1)]:
            stem = root / f'out-{identity}-{reassign}-{threads}'
            process = run(stem, '-threads', str(threads), '-derepIdentity', str(identity),
                          '-derepReassign', str(reassign), '-derepSrchLen', str(cap),
                          '-derepStoreDiffs', str(binary_output), '-derepSubclusterFasta', str(binary_output))
            assert process.returncode == 0, process.stderr
            r1 = list(fastq(Path(str(stem) + '.1.hq.fq')))
            r2 = list(fastq(Path(str(stem) + '.2.hq.fq')))
            assert len(r1) == len(r2) == len(expected)
            actual = Counter()
            for a, b in zip(r1, r2):
                assert a[0] == b[0] and a[3] == b[3], 'Mate names/counts are not synchronized'
                pair = a[1], b[1]
                actual[pair] += a[3]
                assert (a[2], b[2]) == winners[pair], 'Qualities do not match the selected source fragment'
            assert actual == expected, 'Pair identities or fragment counts differ'
            if not reassign:
                parents = {}
                for record in r1:
                    key = record[1][:cap] if cap > 0 else record[1]
                    parent = record[0].rsplit('.sub', 1)[0]
                    assert parents.setdefault(key, parent) == parent, 'An exact R1 key was split across parent dereplicates'

            map_total = sum(int(entry.split(':')[1]) for line in Path(str(stem) + '.map').read_text().splitlines()
                            if not line.startswith('#') for entry in line.split('\t')[1:] if ':' in entry)
            assert map_total == len(observations), 'Counted mates rather than fragments'
            assert not any(root.glob(stem.name + '.subclusters*.fq')), 'Unexpected separate subcluster FASTQ'
            assert Path(str(stem) + '.diff').exists() == bool(binary_output)
            if binary_output:
                assert Path(str(stem) + '.diff').read_bytes().startswith(b'SDMDIFF3')
                fa1 = fasta(Path(str(stem) + '.subclusters.1.fna'))
                fa2 = fasta(Path(str(stem) + '.subclusters.2.fna'))
                assert fa1 == [(a[0], a[1]) for a in r1] and fa2 == [(b[0], b[1]) for b in r2]
            print(f'PASS paired: {identity}% reassign={reassign} threads={threads} cap={cap}: '
                  f'{len(actual)} exact pairs / {sum(actual.values())} fragments, selected source qualities')
        # One shared FASTQ output includes the original anchor, its better
        # replacement, and a singleton parent below the abundance cutoff.
        anchor = tuple(sequence[:200] for sequence in templates[0])
        variant = tuple(sequence[:31] + next(base for base in 'ACGT' if base != sequence[31]) + sequence[32:]
                        for sequence in anchor)
        singleton = tuple(sequence[:200] for sequence in templates[1])
        mixed = [(anchor, 20)] * 2 + [(variant, 40)] * 3 + [(singleton, 30)]
        for mate, path in enumerate(inputs):
            path.write_text(''.join(f'@mixed{i}/{mate + 1}\n{pair[mate]}\n+\n' +
                                    chr(33 + quality) * 200 + '\n'
                                    for i, (pair, quality) in enumerate(mixed)))
        stem = root / 'all-dereplicates'
        process = run(stem, '-threads', '1', '-derepIdentity', '97', '-derepSrchLen', '200',
                      '-min_derep_copies', '2', '-suppressOutput', '3')
        assert process.returncode == 0, process.stderr
        first = list(fastq(Path(str(stem) + '.1.hq.fq')))
        second = list(fastq(Path(str(stem) + '.2.hq.fq')))
        assert len(first) == len(second) == 3, 'Anchor, representative, and singleton must each appear once'
        actual, parents, parent_totals = Counter(), {}, Counter()
        qualities = {pair: chr(33 + quality) * 200 for pair, quality in mixed}
        for a, b in zip(first, second):
            assert a[0] == b[0] and a[3] == b[3], 'Mixed output lost mate synchronization'
            pair = a[1], b[1]
            actual[pair] += a[3]
            assert a[2] == b[2] == qualities[pair], 'Mixed output changed the selected pair qualities'
            parents[pair] = a[0].rsplit('.sub', 1)[0]
            parent_totals[parents[pair]] += a[3]
        assert actual == Counter({anchor: 2, variant: 3, singleton: 1}), 'Omitted or duplicated exact dereplicates'
        assert len(set(parents.values())) == 3, 'Nearby R1 dereplicates must retain separate output parents'
        assert parent_totals == Counter({parents[anchor]: 2, parents[variant]: 3, parents[singleton]: 1})
        main_records = fasta(Path(str(stem) + '.fna'))
        rest_records = fasta(Path(str(stem) + '.fna.rest'))
        assert Counter(seq for _, seq in main_records) == Counter((anchor[0], variant[0])), 'Both exact R1 dereplicates must pass'
        assert len(rest_records) == 1 and rest_records[0][1] == singleton[0], 'Singleton must be below the parent cutoff'
        assert not any(root.glob(stem.name + '.subclusters*.fq')), 'Separate subread FASTQ was emitted'
        assert not Path(str(stem) + '.diff').exists(), 'Unified FASTQ must not require binary export'
        print('PASS unified HQ FASTQ: anchor, selected representative, and below-cutoff singleton; 3 exact pairs / 6 fragments')
        for flag in ('-merge_pairs_derep', '-merge_pairs_filter', '-merge_pairs_demulti'):
            failed = run(root / ('invalid' + flag), '-threads', '1', flag, '1')
            assert failed.returncode != 0 and 'unmerged' in failed.stderr, failed.stderr
        # The same six observations were just run at 97% above. Change only
        # identity, retaining the same caps, quality retention, and copy cutoff.
        # Both near R1 variants must pass the cutoff separately at either identity.
        exact_stem = root / 'all-dereplicates-100'
        process = run(exact_stem, '-threads', '1', '-derepIdentity', '100', '-derepSrchLen', '200',
                      '-min_derep_copies', '2', '-suppressOutput', '3')
        assert process.returncode == 0, process.stderr
        assert_identity_output_equivalence(exact_stem, stem)
        thousand_pair_equivalence(root, inputs, run)
    print('Paired diffDNA pipeline checks passed')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else './sdm')
