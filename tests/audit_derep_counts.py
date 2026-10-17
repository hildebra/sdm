#!/usr/bin/env python3
"""Reconcile finalized SDM map, main, merged and rest counts (derepPerSR=0)."""
import argparse
from collections import Counter
import json
from pathlib import Path


def records(path):
    """Stream FASTA (including wrapped sequences) or SDM four-line FASTQ."""
    with Path(path).open() as stream:
        first = stream.readline()
        if not first:
            return
        if first.startswith('>'):
            name, length = first[1:].rstrip(), 0
            for line in stream:
                if line.startswith('>'):
                    yield name, length
                    name, length = line[1:].rstrip(), 0
                else:
                    length += len(line.strip())
            yield name, length
        elif first.startswith('@'):
            name = first
            while name:
                seq, plus, qual = (stream.readline() for _ in range(3))
                assert name.startswith('@') and plus.startswith('+') and qual, f'Malformed FASTQ: {path}'
                assert len(seq.rstrip()) == len(qual.rstrip()), f'Quality length: {path}'
                yield name[1:].rstrip(), len(seq.rstrip())
                name = stream.readline()
        else:
            raise ValueError(f'Expected FASTA or FASTQ: {path}')


def audit(main):
    main = Path(main)
    stem = main.with_suffix('')
    parents, sample_names, total = {}, {}, Counter()
    with Path(str(stem) + '.map').open() as stream:
        for line in stream:
            fields = line.rstrip().split('\t')
            if fields[0] == '#SMPLS':
                sample_names.update((int(k), v) for entry in fields[1:] for k, v in [entry.split(':', 1)])
            elif not line.startswith('#') and line.strip():
                name, values = fields[0], Counter()
                assert name not in parents, f'Duplicate map ID: {name}'
                for entry in fields[1:]:
                    if ':' not in entry:
                        continue
                    key, count = entry.split(':', 1)
                    assert int(count) >= 0 and int(key) in sample_names, entry
                    values[sample_names[int(key)]] += int(count)
                parents[name] = values
                total.update(values)
    files = {'main': main, 'merged': Path(str(stem) + '.merg' + main.suffix), 'rest': Path(str(main) + '.rest')}
    partitions, seen = {}, set()
    for label, path in files.items():
        counts, lengths, number = Counter(), Counter(), 0
        if path.exists():
            for name, length in records(path):
                assert name in parents and name not in seen, f'Unknown/duplicate output parent: {name}'
                seen.add(name)
                abundance = sum(parents[name].values())
                if ';size=' in name:
                    assert int(name.split(';size=')[1].split(';')[0]) == abundance, name
                counts.update(parents[name]); lengths[length] += 1; number += 1
        partitions[label] = dict(records=number, counts=sum(counts.values()), samples=dict(counts), lengths=dict(lengths))
    assert seen == set(parents), f'{len(set(parents) - seen)} map parents absent from passing/rest outputs'
    reconciled = Counter()
    for part in partitions.values():
        reconciled.update(part['samples'])
    assert reconciled == total, 'Per-sample map != passing + rest'
    return dict(main=str(main), map_total=sum(total.values()), samples=dict(total), parents=len(parents),
                passing=partitions['main']['counts'] + partitions['merged']['counts'],
                rest=partitions['rest']['counts'], outputs=partitions)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path, help='Literal -o_dereplicate path; derepPerSR=0')
    args = parser.parse_args()
    print(json.dumps(audit(args.output), indent=2))
