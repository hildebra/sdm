#!/usr/bin/env python3
"""Compare ordinary and 97% storage outputs, including IDs, qualities and normal seeds.

Uses the shared 1,000-pair mixed-quality fixture. No variant-retention flags are
implied by identity. --observe records differences for a historical executable.
"""
import argparse
import hashlib
import json
from pathlib import Path
import random
import re
import subprocess
import tempfile

from test_derep_behavior_parity import create_fixture
from test_merged_pipeline import fasta_records


def fq(path):
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0, path
    result = []
    for i in range(0, len(lines), 4):
        name, sequence, plus, quality = lines[i:i+4]
        assert name.startswith('@') and plus.startswith('+') and len(sequence) == len(quality), path
        result.append((name[1:], sequence, quality))
    return result


def snapshot(folder, pattern, fmt):
    result = {}
    for path in sorted(folder.glob(pattern)):
        if path.suffix in ('.log', '.json'): continue
        if path.suffix == '.map':
            # IDs and every sample entry are retained; only row/column order is irrelevant.
            value = sorted(tuple([line.split('\t')[0], *sorted(line.split('\t')[1:])])
                           for line in path.read_text().splitlines())
        elif path.name.endswith('.hq.fq') or fmt == 'fq':
            value = sorted(fq(path))
        else:
            value = sorted(map(tuple, fasta_records(path)))
        result[path.name] = value
    return result


def run(binary, folder, name, args):
    command = [str(binary), *map(str, args)]
    result = subprocess.run(command, capture_output=True, text=True, timeout=120)
    (folder / (name + '.command.json')).write_text(json.dumps(command, indent=2) + '\n')
    (folder / (name + '.console.log')).write_text(result.stdout + result.stderr)
    assert result.returncode == 0, (command, result.stderr[-3000:])
    return command


def prepare(root):
    settings, _, _ = create_fixture(root, complete_mates=True, consensus_errors=True)
    rng = random.Random(2026091523)
    raw = set()
    for sample in range(2):
        first, second = (fq(root / f'S{sample}.{mate}.fq') for mate in (1, 2))
        output = [[], []]
        unmergeable = ''.join(rng.choice('ACGT') for _ in range(230))
        for a, b in zip(first, second):
            pair = []
            for mate, (name, seq, quality) in enumerate((a, b)):
                # Variable observed qualities, with a unique strong observation
                # in each family/subgroup, avoiding accidental ID ties here.
                if min(map(ord, quality)) >= 53:
                    q = 50 if name.endswith('_0/' + str(mate + 1)) else 30 + rng.randrange(10)
                    quality = ''.join(chr(33 + q - (i % 3)) for i in range(len(seq)))
                if mate == 1 and '_fr1first_' in name and '_both_' in name:
                    # Same R1, a distinct R2-only full-fragment variant.
                    pos = 23
                    seq = seq[:pos] + next(c for c in 'ACGT' if c != seq[pos]) + seq[pos+1:]
                if mate == 1 and '_fsplit_b_' in name:
                    # Different overlap offset inside a shared capped search key.
                    seq, quality = seq[:20] + seq[40:], quality[:20] + quality[40:]
                if '_fzero_' in name:
                    copy = int(name.rsplit('_', 1)[1].split('/')[0])
                    if mate == 1 and copy < 19: seq = seq[:20] + unmergeable
                    quality = chr(33 + (60 if copy == 19 else 40 if copy == 0 else 30 + copy % 5)) * len(seq)
                pair.append((name, seq, quality))
                output[mate].append(f'@{name}\n{seq}\n+\n{quality}\n')
            raw.add((pair[0][1], pair[0][2], pair[1][1][20:], pair[1][2][20:]))
        for mate in (1, 2):
            (root / f'S{sample}.{mate}.fq').write_text(''.join(output[mate-1]))
    rows = (root / 'input.map').read_text().splitlines()
    (root / 'runs.map').write_text(rows[0] + '\tSequencingRun\n' +
        '\n'.join(row + '\trun' + str(i) for i, row in enumerate(rows[1:])) + '\n')
    return settings, raw


def validate_parents(folder, raw, fmt):
    first, second = (fq(folder / f'derep.{mate}.hq.fq') for mate in (1, 2))
    assert len(first) == len(second)
    headers = set()
    for a, b in zip(first, second):
        assert (a[1], a[2], b[1], b[2]) in raw, 'HQ bases/qualities no longer belong to one full observed pair'
        norm = lambda name: re.sub(r'/[12]$', '', name.split(';size=')[0])
        assert norm(a[0]) == norm(b[0]), 'HQ mate IDs are unlinked'
        assert '.sub' not in a[0], 'Storage-only mode exported variant candidates'
        assert a[0] not in headers, 'Duplicate representative HQ ID'
        headers.add(a[0])
    rows = [line.split('\t') for line in (folder / 'derep.map').read_text().splitlines() if not line.startswith('#')]
    assert headers == {row[0] for row in rows}, 'HQ parents and map IDs differ'
    counts = {row[0]: sum(int(entry.split(':')[1]) for entry in row[1:] if ':' in entry) for row in rows}
    emitted = {}
    for path in folder.glob('derep*.fna*'):
        records = [(name[1:], seq) for name, seq in fasta_records(path)] if fmt == 'fa' else [(name, seq) for name, seq, _ in fq(path)]
        for name, _ in records:
            assert name not in emitted, 'Parent emitted in more than one output partition'
            emitted[name] = int(re.search(r';size=(\d+);', name)[1])
    assert emitted == counts, 'Main + merged + rest records do not partition map abundance'
    return sum(counts.values())


def seed(binary, folder, config, threads, merge):
    out = folder / 'seed'; out.mkdir()
    parents = fq(folder / 'derep.1.hq.fq')
    targets = {}
    assignments = []
    for name, seq, _ in parents:
        family = re.search(r'_f([^_]+)_', name)[1]
        target = targets.setdefault(family, (f'OTU{len(targets)}', seq))
        # PAF query IDs are ordinary parent IDs, without any coarse membership.
        assignments.append(f'{name}\t{len(seq)}\t0\t{len(seq)}\t+\t{target[0]}\t{len(seq)}\t0\t{len(seq)}\t{len(seq)}\t{len(seq)}\t60\n')
    (out / 'assign.paf').write_text(''.join(assignments))
    (out / 'fallback.fna').write_text(''.join(f'>{name}\n{seq}\n' for name, seq in targets.values()))
    command = run(binary, out, 'seed', ['-i_fastq', f'{folder}/derep.1.hq.fq,{folder}/derep.2.hq.fq',
        '-paired', 2, '-options', config, '-sample_sep', '___', '-threads', threads,
        '-i_qual_offset', 33, '-o_qual_offset', 33, '-merge_pairs_seed', merge,
        '-optimalRead2Cluster', out/'assign.paf', '-OTU_fallback', out/'fallback.fna',
        '-derep_map', folder/'derep.map', '-o_fastq', f'{out}/OTU.1.fq,{out}/OTU.2.fq',
        '-otu_matrix', out/'matrix.tsv', '-log', out/'stats.log'])
    assert '-seedSubclusters' not in command
    outputs = {p.name: p.read_text().replace(str(folder), '<run>') for p in out.iterdir()
               if p.name.startswith('OTU.') or p.name in ('matrix.tsv', 'stats.log', 'assign.paf', 'fallback.fna')}
    return [arg.replace(str(folder), '<run>') for arg in command], outputs


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary', nargs='?', type=Path, default=Path('./sdm'))
    p.add_argument('--output', type=Path)
    p.add_argument('--observe', action='store_true')
    p.add_argument('--threads', type=int, nargs='+', default=[1, 4, 12])
    p.add_argument('--case', action='append')
    args = p.parse_args(); binary = args.binary.resolve()
    root = args.output or Path(tempfile.mkdtemp(prefix='sdm-storage-parity-'))
    if args.output: root.mkdir(parents=True, exist_ok=False)
    settings, raw = prepare(root)
    # name, cap, truncate, quality trim, preprocessing merge, R1-only search, per-run
    cases = [('capped-r1', 200, -1, 0, 0, 1, 0), ('capped-merged', 200, -1, 0, 1, 0, 0),
             ('merged-output-r1-search', 200, -1, 0, 1, 1, 0),
             ('prefix-r1', -1, -1, 0, 0, 1, 0), ('prefix-merged', -1, -1, 0, 1, 0, 0),
             ('trimmed-r1', 200, 200, 25, 0, 1, 0), ('trimmed-merged', 200, 200, 25, 1, 0, 0),
             ('runs-merged', 200, -1, 0, 1, 0, 1)]
    results = []; across_threads = {}
    for name, cap, truncate, trim, merge, r1, per_run in cases:
      if args.case and name not in args.case: continue
      cfg = dict(settings, TruncateSequenceLength=truncate, TrimWindowThreshhold=trim, TrimWindowWidth=20)
      if name == 'capped-r1': cfg['*minAvgQuality'] = 8
      config = root / (name + '.options');config.write_text(''.join(f'{k}\t{v}\n' for k, v in cfg.items()))
      for fmt in ('fa', 'fq'):
       for threads in args.threads:
        pair = []; seeds = []; commands = []
        for identity in ('100', '97'):
            folder = root / f'{name}-{fmt}-{threads}t-{identity}'; folder.mkdir()
            command = run(binary, folder, 'primary', ['-i_path', root, '-map', root/('runs.map' if per_run else 'input.map'),
                '-options', config, '-paired', 2, '-i_qual_offset', 33, '-o_qual_offset', 33,
                '-sample_sep', '___', '-threads', threads, '-threadIO', 0, '-derepPerSR', per_run,
                '-o_dereplicate', folder/'derep.fna', '-derep_format', fmt, '-derepIdentity', identity,
                '-derepStoreQuals', 0, '-derepStoreDiffs', 0, '-derepCoarseClusters', 0, '-derepReassign', 0,
                '-min_derep_copies', '40', '-derep_reserve', 1, '-derepSrchLen', cap,
                '-merge_pairs_derep', merge, '-merge_pairs_filter', merge, '-merge_pairs_demulti', merge,
                *(['-noSearchWithMerge', 1] if r1 else []), '-5PR2cut', 20,
                '-o_fna', f'{folder}/filtered.1.fna,{folder}/filtered.2.fna',
                '-o_fna2', f'{folder}/filtered.1.add.fna,{folder}/filtered.2.add.fna', '-log', folder/'filter.log'])
            console = (folder/'primary.console.log').read_text()
            if not args.observe:
                policy = 'R1 only' if r1 else 'merged when available'
                assert 'Dereplication search policy: ' + policy in console
                sources = re.search(r'Dereplication key sources: ([\d,]+) merged', console)
                assert (int(sources[1].replace(',','')) > 0) == bool(merge and not r1), 'Fixture used the wrong search source'
            total = validate_parents(folder, raw, fmt)
            state = snapshot(folder, 'derep*', fmt)
            assert any(state[k] for k in state if k.endswith('.rest')), 'Fixture lacks below-cutoff parents'
            assert any(state[k] for k in state if k.endswith('.map'))
            pair.append(state)
            commands.append([arg.replace(str(folder), '<run>') for arg in command])
            seeds.append(seed(binary, folder, config, threads, 0 if name == 'capped-r1' else 1))
        differing = [k for k in pair[0].keys() | pair[1].keys() if pair[0].get(k) != pair[1].get(k)]
        seed_diff = [k for k in seeds[0][1].keys() | seeds[1][1].keys() if seeds[0][1].get(k) != seeds[1][1].get(k)]
        assert seeds[0][0] == seeds[1][0], 'Coarse storage changed the normal seed command'
        for command in commands: command[command.index('-derepIdentity')+1] = '<identity>'
        assert commands[0] == commands[1], 'Preprocessing comparison changed another option'
        row = dict(case=name, format=fmt, threads=threads, map_total=total, output_differences=differing, seed_differences=seed_diff)
        results.append(row)
        for identity, state in zip(('100', '97'), pair):
            key = name, fmt, identity
            baseline = across_threads.setdefault(key, state)
            row.setdefault('thread_differences', {})[identity] = [k for k in state if state[k] != baseline[k]]
        print(json.dumps(row), flush=True)
        if not args.observe:
            assert not differing and not seed_diff, row
            assert not any(row['thread_differences'].values()), row
    (root/'results.json').write_text(json.dumps(dict(binary=str(binary), sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        fragments=1000, runs=results), indent=2)+'\n')
    print(f'Completed {len(results)} ordinary/97% comparisons; artifacts: {root}', flush=True)


if __name__ == '__main__':
    main()
