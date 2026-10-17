#!/usr/bin/env python3
"""Measure the paired coarse CLI at different worker counts (Linux/Unix).

Reuses the retained-quality pipeline's FASTQ reader and checks exact fragment
counts and paired quality witnesses after each timed process. Synthetic input
generation and output validation are excluded from the reported timings.
"""
import argparse
from collections import Counter
import csv
import hashlib
from itertools import zip_longest
from pathlib import Path
import random
import re
import resource
import shutil
import os
import threading
import subprocess
import tempfile
import time

from test_diff_quality_pipeline import fastq, retained_fastq_paths


def timed_run(command, cwd, trace_path, sample_cpu):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    start = time.perf_counter()
    samples = []
    stop = threading.Event()
    with subprocess.Popen(command, cwd=cwd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) as process:
        def monitor():
            ticks = os.sysconf('SC_CLK_TCK')
            previous_time, previous_cpu = start, 0.0
            while not stop.wait(0.1):
                try:
                    fields = Path(f'/proc/{process.pid}/stat').read_text().rsplit(')', 1)[1].split()
                except FileNotFoundError:
                    return
                now = time.perf_counter()
                cpu = (int(fields[11]) + int(fields[12])) / ticks
                samples.append((now - start, 100 * (cpu - previous_cpu) / (now - previous_time)))
                previous_time, previous_cpu = now, cpu
        observer = threading.Thread(target=monitor) if sample_cpu else None
        if observer:
            observer.start()
        try:
            stdout, stderr = process.communicate()
        except BaseException:
            process.kill(); process.communicate()
            raise
        finally:
            stop.set()
            if observer:
                observer.join()
        result = subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
    wall = time.perf_counter() - start
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu = after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime
    if sample_cpu:
        with trace_path.open('w', newline='') as output:
            writer = csv.writer(output); writer.writerow(['seconds', 'cpu_percent']); writer.writerows(samples)
    return result, wall, cpu


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--compare-binary', type=Path, help='Alternate this baseline with --binary on identical inputs')
    parser.add_argument('--sample-cpu', action='store_true', help='Record Linux CPU use at 100 ms intervals')
    parser.add_argument('--reads', type=int, default=200000)
    parser.add_argument('--templates', type=int, default=50)
    parser.add_argument('--variants', type=int, default=10000, help='Number of generated variant pairs before resampling observations')
    parser.add_argument('--length', type=int, default=250)
    parser.add_argument('--shared-r1', action='store_true', help='Stress the R1-only candidate index')
    parser.add_argument('--global-mode', type=int, choices=[0, 1], default=0)
    parser.add_argument('--threads', type=int, nargs='+', default=[1, 4, 12])
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.reads < 1 or args.variants < 1 or args.repeats < 1 or any(t < 1 for t in args.threads) or not 1 <= args.templates <= 10000:
        parser.error('reads, variants, repeats, and threads must be positive; templates must be in 1..10000')
    if not 50 <= args.length <= 500:
        parser.error('length must be in 50..500')
    binary = str(args.binary.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rng = random.Random(2026091401)
    length = args.length
    templates = [tuple(''.join(rng.choices('ACGT', k=length)) for _ in range(2))
                 for _ in range(args.templates)]
    if args.shared_r1:
        templates = [(templates[0][0], mate) for _, mate in templates]
    variants = []
    for i in range(args.variants):
        pair = []
        for template in templates[i % len(templates)]:
            sequence = list(template)
            for position in rng.sample(range(length), i % 4):
                sequence[position] = rng.choice([b for b in 'ACGT' if b != sequence[position]])
            pair.append(''.join(sequence))
        variants.append(tuple(pair))
    expected = Counter()
    with tempfile.TemporaryDirectory(prefix='sdm-paired-thread-bench-') as temp:
        root = Path(temp)
        # Equalize executable loading when the checkout is on a Windows mount.
        versions = [('current', Path(binary))]
        if args.compare_binary:
            versions.insert(0, ('before', args.compare_binary.resolve()))
        binaries = []
        for version, source in versions:
            target = root / ('sdm-' + version)
            shutil.copy2(source, target)
            binaries.append((version, str(target)))
        inputs = [root / f'reads.{mate}.fq' for mate in (1, 2)]
        with inputs[0].open('w') as first, inputs[1].open('w') as second:
            for i in range(args.reads):
                pair = variants[rng.randrange(len(variants))]
                expected[pair] += 1
                for mate, output in enumerate((first, second)):
                    output.write(f'@fragment{i}/{mate + 1}\n{pair[mate]}\n+\n' + '?' * length + '\n')
        settings = {'minSeqLength': 50, 'maxSeqLength': 500, 'minAvgQuality': 0,
                    'QualWindowThreshhold': 0, 'TrimWindowThreshhold': 0,
                    'maxAccumulatedError': -1, 'BinErrorModelAlpha': -1,
                    'maxAmbiguousNT': 0, 'maxHomonucleotide': 0, 'trimHomonucleotide': 0,
                    'TruncateSequenceLength': -1, 'fastqVersion': 1}
        config = root / 'options.txt'
        config.write_text(''.join(f'{key}\t{value}\n' for key, value in settings.items()))
        mapping = root / 'map.tsv'
        mapping.write_text(f'#SampleID\tfastqFile\nS0\t{inputs[0].name},{inputs[1].name}\n')
        with args.output.open('w', newline='') as output:
            writer = csv.DictWriter(output, fieldnames=['repeat', 'version', 'threads', 'fragments', 'templates', 'input_variants', 'shared_r1', 'global_mode', 'partitions', 'largest_partition_queries',
                'mate_length', 'wall_seconds', 'cpu_seconds', 'cpu_percent', 'exact_pairs',
                'coarse_clusters', 'subcluster_fastq_bytes', 'peak_rss_bytes', 'comparisons', 'fallback_queries', 'exact_reuse', 'candidate_searches', 'pool_tasks', 'pool_helped', 'pool_peak_active', 'output_sha256'])
            writer.writeheader()
            for repeat in range(1, args.repeats + 1):
                for threads in args.threads:
                    for version, binary in (binaries if repeat % 2 else list(reversed(binaries))):
                        stem = root / f'out-{repeat}-{threads}-{version}'
                        command = [binary, '-i_path', str(root), '-map', str(mapping),
                            '-paired', '2', '-options', str(config), '-o_dereplicate', str(stem) + '.fna',
                            '-o_fna', str(stem) + '.1.fna,' + str(stem) + '.2.fna',
                            '-o_fna2', str(stem) + '.1.add.fna,' + str(stem) + '.2.add.fna',
                            '-sample_sep', '___', '-log', str(stem) + '.log',
                            '-merge_pairs_derep', '0', '-merge_pairs_filter', '0', '-merge_pairs_demulti', '0',
                            '-dere_size_fmt', '0', '-min_derep_copies', '8:1,4:2,3:3',
                            '-derepIdentity', '97', '-derepStoreQuals', '1', '-derepPerSR', '0',
                            '-suppressOutput', '1', '-o_qual_offset', '33', '-oneLineFastaFormat', '1',
                            '-threads', str(threads), '-derepGlobal', str(args.global_mode)]
                        log = args.output.with_name(args.output.stem + f'-{repeat}-{threads}-{version}.log')
                        result, wall, cpu = timed_run(command, root, log.with_suffix('.cpu.csv'), args.sample_cpu)
                        log.write_text(result.stdout + result.stderr)
                        if result.returncode:
                            raise RuntimeError(f'Benchmark failed: see {log}')
                        actual = Counter()
                        emitted_pairs = 0
                        paths = retained_fastq_paths(stem, paired=True)
                        for first, second in zip_longest(fastq(paths[0]), fastq(paths[1])):
                            assert first is not None and second is not None, 'Missing mate'
                            assert first[0] == second[0] and first[3] == second[3], 'Unsynchronized mates'
                            assert first[2] == second[2] == '?' * length, 'Changed quality witness'
                            actual[first[1], second[1]] += first[3]
                            emitted_pairs += 1
                        assert actual == expected, 'Exact paired fragment counts changed'
                        assert emitted_pairs == len(expected), 'Exact pairs were duplicated across partitions'
                        map_lines = Path(str(stem) + '.map').read_text().splitlines()
                        clusters = [line for line in map_lines if line and not line.startswith('#')]
                        map_total = sum(int(entry.split(':')[1]) for line in clusters
                                        for entry in line.split('\t')[1:] if ':' in entry)
                        assert map_total == args.reads, 'Sample map lost fragment counts'
                        matching = re.search(r'Coarse matching: ([0-9,]+) bounded comparisons; ([0-9,]+) queries requiring seed fallback', result.stderr)
                        assert matching, 'Missing coarse search diagnostics'
                        partitions = re.search(r'Coarse partitions: ([0-9]+); active [0-9]+; largest ([0-9,]+) of', result.stderr)
                        memory = re.search(r'Peak process memory: ([0-9]+) bytes', result.stderr)
                        checksum = hashlib.sha256()
                        for path in [Path(str(stem) + suffix) for suffix in ('.fna', '.fna.rest', '.map')] + paths:
                            checksum.update(path.read_bytes())
                        reuse = re.search(r'Coarse exact reuse: ([0-9,]+) of [0-9,]+ queries; ([0-9,]+) candidate searches', result.stderr)
                        pool = re.search(r'Coarse work pool: [0-9,]+ batches; ([0-9,]+) partition tasks; ([0-9,]+) helped by another worker; peak active ([0-9]+)', result.stderr)
                        number = lambda match, group: int(match[group].replace(',', '')) if match else 0
                        row = dict(repeat=repeat, version=version, threads=threads, fragments=args.reads,
                                   templates=args.templates, input_variants=args.variants, shared_r1=int(args.shared_r1), global_mode=args.global_mode,
                                   partitions=int(partitions[1]) if partitions else 1,
                                   largest_partition_queries=int(partitions[2].replace(',', '')) if partitions else args.reads,
                                   mate_length=length,
                                   wall_seconds=wall, cpu_seconds=cpu, cpu_percent=100 * cpu / wall,
                                   exact_pairs=len(actual), coarse_clusters=len(clusters),
                                   subcluster_fastq_bytes=sum(path.stat().st_size for path in paths),
                                   peak_rss_bytes=int(memory[1]) if memory else 0,
                                   comparisons=int(matching[1].replace(',', '')),
                                   fallback_queries=int(matching[2].replace(',', '')),
                                   exact_reuse=number(reuse, 1), candidate_searches=number(reuse, 2) if reuse else args.reads,
                                   pool_tasks=number(pool, 1), pool_helped=number(pool, 2), pool_peak_active=number(pool, 3),
                                   output_sha256=checksum.hexdigest())
                        writer.writerow(row)
                        output.flush()
                        print(row, flush=True)
                        for path in root.glob(stem.name + '*'):
                            path.unlink()


if __name__ == '__main__':
    main()
