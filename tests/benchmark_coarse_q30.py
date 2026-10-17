#!/usr/bin/env python3
"""Repeat the production 200 bp Q30 benchmark in independent processes."""
import argparse
import csv
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--reads', type=int, default=1_000_000)
parser.add_argument('--repeats', type=int, default=3)
parser.add_argument('--output', type=Path, default=Path('docs/benchmarks/coarse-production-q30-200bp.csv'))
args = parser.parse_args()
if args.repeats < 1 or not 0 < args.reads <= 1_000_000:
    parser.error('reads must be 1..1000000; repeats must be positive')
rows = []
cases = [('exact', 0, 0)] + [(identity, diffs, reassign)
    for reassign in (0, 1) for diffs in (0, 1) for identity in ('99', '98', '97')]
args.output.parent.mkdir(parents=True, exist_ok=True)
with args.output.with_suffix('.log').open('w') as log:
    for repeat in range(1, args.repeats + 1):
        for identity, diffs, reassign in cases:
            command = [str(args.binary.resolve()), identity, str(args.reads), str(diffs), '0', '200', '30', str(reassign)]
            result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=log, check=True)
            row = list(csv.DictReader(result.stdout.splitlines()))[0]
            row['repeat'] = repeat
            rows.append(row)
            print(f"repeat={repeat} identity={identity} diffs={diffs} reassign={reassign} "
                  f"clusters={row['clusters']} moved={row['moved_reads']} "
                  f"seconds={row['seconds']} peak_KiB={row['peak_rss_kib']}", flush=True)
with args.output.open('w', newline='') as output:
    writer = csv.DictWriter(output, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
