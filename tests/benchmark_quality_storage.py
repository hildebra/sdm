#!/usr/bin/env python3
"""Compare byte qualities with a six-bit prototype on Q0..Q60 data only."""
import argparse
import csv
import os
import random
import subprocess
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--output', type=Path, required=True)
parser.add_argument('--repeats', type=int, default=5)
parser.add_argument('--bases', type=int, default=60_000_000)
parser.add_argument('--passes', type=int, default=8)
parser.add_argument('--cpu', type=int)
args = parser.parse_args()
if args.cpu is not None:
    os.sched_setaffinity(0, {args.cpu})
args.output.parent.mkdir(parents=True, exist_ok=True)
expected = None
with args.output.open('w', newline='') as output:
    writer = None
    for repeat in range(1, args.repeats + 1):
        modes = ['signed_byte', 'unsigned_byte', 'packed_six']
        random.Random(1841 + repeat).shuffle(modes)
        for mode in modes:
            result = subprocess.run([str(args.binary.resolve()), mode, str(args.bases), str(args.passes)],
                                    text=True, capture_output=True, check=True)
            row = list(csv.DictReader(result.stdout.splitlines()))[0]
            if expected is not None and row['checksum'] != expected:
                raise AssertionError('Quality contents changed across representations')
            expected = row['checksum']
            row['repeat'] = repeat
            if writer is None:
                writer = csv.DictWriter(output, fieldnames=list(row)); writer.writeheader()
            writer.writerow(row); output.flush()
            print(repeat, mode, row['scan_reverse_seconds'], row['peak_rss_kib'], flush=True)
