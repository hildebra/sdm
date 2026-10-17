#!/usr/bin/env python3
"""Alternate unchanged/optimized binaries; check output identity and record RSS/time.

Build TestDNAMemory and TestCoarseBenchmark with the same release flags in both
source trees before running. Memory is each process's Linux VmHWM, in KiB.
"""
import argparse
import csv
import random
import os
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--before', type=Path, required=True, help='baseline build directory')
parser.add_argument('--after', type=Path, required=True, help='optimized build directory')
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--repeats', type=int, default=5)
parser.add_argument('--cpu', type=int, help='pin the runner and child processes to one Linux CPU')
args = parser.parse_args()
if args.cpu is not None:
    os.sched_setaffinity(0, {args.cpu})
cases = [('objects', ['dna', '300000', '150']),
         ('objects', ['stream', '1000000', '150']),
         ('objects', ['equal', '200000', '150']),
         ('objects', ['variable', '200000', '150']),
         ('objects', ['equal', '10000', '10000']),
         ('objects', ['history', '100000', '150']),
         ('objects', ['history_raw', '100000', '150']),
         ('pipeline', ['exact', '200000', '0', '0', '200', '60', '0', '0', 'fq']),
         ('pipeline', ['exact', '30000', '0', '100', '1500', '30', '0', '0', 'fq']),
         ('pipeline', ['exact', '100000', '0', '100', '1500', '30', '0', '0', 'fq']),
         ('pipeline', ['97', '200000', '0', '0', '200', '30', '0', '0', 'fa']),
         ('pipeline', ['97', '20000', '1', '100', '1500', '20', '1', '1', 'fa'])]
args.output.parent.mkdir(parents=True, exist_ok=True)
fields = ['case', 'kind', 'version', 'repeat', 'seconds', 'peak_rss_kib', 'checksum', 'dna_bytes', 'unique_bytes']
expected = {}
with args.output.open('w', newline='') as output, args.output.with_suffix('.log').open('w') as log:
    writer = csv.DictWriter(output, fieldnames=fields); writer.writeheader()
    for repeat in range(1, args.repeats + 1):
        order = cases.copy(); random.Random(1907 + repeat).shuffle(order)
        for kind, parameters in order:
            case = ':'.join(parameters)
            versions = [('before', args.before), ('after', args.after)]
            if repeat % 2 == 0: versions.reverse()
            for version, directory in versions:
                binary = directory / ('TestDNAMemory' if kind == 'objects' else 'TestCoarseBenchmark')
                command = [str(binary.resolve()), *parameters]
                process = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=log, check=True)
                row = list(csv.DictReader(process.stdout.splitlines()))[0]
                checksum = row['checksum' if kind == 'objects' else 'output_checksum']
                if case in expected and expected[case] != checksum:
                    raise AssertionError(f'Output differs for {case}: {expected[case]} vs {checksum}')
                expected[case] = checksum
                result = dict(case=case, kind=kind, version=version, repeat=repeat,
                              seconds=row['seconds'], peak_rss_kib=row['peak_rss_kib'], checksum=checksum,
                              dna_bytes=row.get('dna_bytes', ''), unique_bytes=row.get('unique_bytes', ''))
                writer.writerow(result); output.flush()
                print(f'{repeat} {version} {case}: {row["seconds"]} s, {row["peak_rss_kib"]} KiB, checksum {checksum}', flush=True)
