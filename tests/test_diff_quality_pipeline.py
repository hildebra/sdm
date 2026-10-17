#!/usr/bin/env python3
"""Full FASTQ pipeline checks for selected single-end diffDNA qualities.

Rough ONT-like fixture: ~1,500 bp, 1% substitutions, variable terminal length,
homopolymer insertion/deletion stress, and repeated observations with Q10/20/30
profiles. The error and quality generators are intentionally independent.
Usage: python3 tests/test_diff_quality_pipeline.py /path/to/sdm
"""
from collections import Counter
from pathlib import Path
import random
import subprocess
import sys
import tempfile


def fastq(path, allow_unannotated=False):
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0, path
    for i in range(0, len(lines), 4):
        name, sequence, plus, quality = lines[i:i + 4]
        assert name.startswith('@') and plus == '+' and len(sequence) == len(quality), name
        if allow_unannotated and ';size=' not in name:
            count = None  # Ordinary representative R2 preserves the input ID.
        else:
            assert name.count(';size=') == 1, name
            count = int(name.split(';size=')[1].split(';')[0])
        yield name[1:], sequence, quality, count



def retained_fastq_paths(stem, paired=False):
    """Benchmark compatibility with saved pre-update binaries in fresh fixtures.

    Production contract tests use the HQ paths explicitly; only comparisons
    against older binaries may read the former separate variant export.
    """
    stem = str(stem)
    old = [Path(stem + suffix) for suffix in
           (('.subclusters.1.fq', '.subclusters.2.fq') if paired else ('.subclusters.fq',))]
    if all(path.exists() for path in old):
        return old
    return [Path(stem + suffix) for suffix in
            (('.1.hq.fq', '.2.hq.fq') if paired else ('.1.hq.fq',))]


def main(binary):
    binary = str(Path(binary).resolve())
    rng = random.Random(2026091307)
    templates = []
    for _ in range(12):
        sequence = ''.join(rng.choices('ACGT', k=1600))
        templates.append(sequence[:450] + 'AAAAAA' + sequence[456:])
    observations = []
    for i in range(600):
        sequence = list(templates[i % len(templates)][:rng.randint(1400, 1600)])
        for pos, base in enumerate(sequence):
            if rng.random() < .01:
                sequence[pos] = rng.choice([b for b in 'ACGT' if b != base])
        # Preserve internal indels exactly even though grouping is substitutions-only.
        if i % 10 == 0:
            sequence.insert(453, 'A')
        elif i % 10 == 1:
            del sequence[453]
        sequence = ''.join(sequence)
        for q in (10, 30, 20):
            quality = ''.join(chr(33 + q + (p % 7) - 3) for p in range(len(sequence)))
            observations.append((sequence, quality, q))
    rng.shuffle(observations)
    counts = Counter(seq for seq, _, _ in observations)
    winners = {seq: qual for seq, qual, q in observations if q == 30}
    with tempfile.TemporaryDirectory(prefix='sdm-diff-quality-cli-') as temp:
        root = Path(temp)
        reads = root / 'reads.fq'
        reads.write_text(''.join(f'@r{i}\n{seq}\n+\n{qual}\n' for i, (seq, qual, _) in enumerate(observations)))
        settings = {'minSeqLength': 50, 'maxSeqLength': 2000, 'minAvgQuality': 0,
                    'QualWindowThreshhold': 0, 'TrimWindowThreshhold': 0,
                    'maxAccumulatedError': -1, 'BinErrorModelAlpha': -1,
                    'maxAmbiguousNT': 0, 'maxHomonucleotide': 0, 'trimHomonucleotide': 0,
                    'TruncateSequenceLength': -1, 'fastqVersion': 1, 'derepStoreQuals': 1}
        config = root / 'options.txt'
        config.write_text(''.join(f'{key}\t{value}\n' for key, value in settings.items()))
        mapping = root / 'map.tsv'
        mapping.write_text(f'#SampleID\tfastqFile\nS0\t{reads}\n')
        cases = [(97, 0, 1, -1), (98, 0, 1, -1), (99, 0, 1, -1),
                 (97, 1, 1, -1), (100, 0, 1, -1), (97, 0, 4, 700)]
        for identity, reassign, threads, cap in cases:
            stem = root / f'out-{identity}-{reassign}-{threads}-{cap}'
            args = [binary, '-i_fastq', str(reads), '-map', str(mapping), '-options', str(config),
                    '-o_fastq', str(stem) + '.fq', '-o_dereplicate', str(stem) + '.fna',
                    '-threads', str(threads), '-threadIO', '0', '-derepIdentity', str(identity),
                    '-derepReassign', str(reassign), '-derepSrchLen', str(cap), '-derep_reserve', '1',
                    '-log', str(stem) + '.log']
            process = subprocess.run(args, cwd=root, text=True, capture_output=True)
            assert process.returncode == 0, process.stderr
            observed = Counter()
            records = list(fastq(Path(str(stem) + '.1.hq.fq')))
            for name, seq, qual, count in records:
                observed[seq] += count
                assert qual == winners[seq], (name, 'Selected source quality differs')
            assert observed == counts, (identity, 'Full processed sequence/count mismatch')
            assert len(records) == len(counts), 'Exact repeats must collapse into one subcluster'
            assert not Path(str(stem) + '.diff').exists(), 'Quality flag must not enable binary output'
            assert not any(root.glob(stem.name + '.subclusters*.fq')), 'Duplicate FASTQ export'
            assert not Path(str(stem) + '.2.hq.fq').exists(), 'Single-end output has no mate file'
            map_total = sum(int(entry.split(':')[1]) for line in Path(str(stem) + '.map').read_text().splitlines()
                            if not line.startswith('#') for entry in line.split('\t')[1:] if ':' in entry)
            assert map_total == len(observations), 'Sample map counts differ'
            print(f'PASS: {identity}% reassign={reassign} threads={threads} cap={cap}: '
                  f'{len(records)} exact variants, {sum(observed.values())} reads, exact selected qualities')
        # CLI precedence: disabling retained qualities keeps normal representative HQ output.
        stem = root / 'disabled'
        args = [binary, '-i_fastq', str(reads), '-map', str(mapping), '-options', str(config),
                '-o_dereplicate', str(stem) + '.fna', '-derepIdentity', '97', '-derepStoreQuals', '0',
                '-threads', '1', '-threadIO', '0', '-log', str(stem) + '.log']
        process = subprocess.run(args, cwd=root, text=True, capture_output=True)
        assert process.returncode == 0, process.stderr
        assert Path(str(stem) + '.1.hq.fq').exists(), 'Normal HQ output must still exist without quality retention'
        assert not any(root.glob(stem.name + '.subclusters*.fq'))
        assert all('.sub' not in name for name, *_ in fastq(Path(str(stem) + '.1.hq.fq'))), 'Disabled retention must emit parent representatives'
        help_text = subprocess.run([binary, '-help_flags'], capture_output=True, text=True, check=True).stdout
        assert '-derepStoreQuals' in help_text
    print('Single-end retained-quality FASTQ pipeline checks passed')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else './sdm')
