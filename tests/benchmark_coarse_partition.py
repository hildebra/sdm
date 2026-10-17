#!/usr/bin/env python3
"""Experimental routing, local clustering, and centroid-merge measurements."""
import argparse
import csv
import subprocess
from pathlib import Path

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--binary', type=Path, required=True)
parser.add_argument('--repeats', type=int, default=3)
parser.add_argument('--output', type=Path, default=Path('docs/benchmarks/coarse-partition.csv'))
args=parser.parse_args()
if args.repeats < 1: parser.error('repeats must be positive')
cases=[]
for mode in ('exact','exact_diff'):
    for threads in (1,4,8): cases.append((mode,threads,64,200000,200,30,100))
for reads,length,quality,identity in [(100000,200,30,97),(20000,1500,20,97)]:
    for threads in (1,4,8): cases.append(('global',threads,1,reads,length,quality,identity))
    for mode in ('roundrobin','fullhash','min20','min8','rare20','supported20'):
        for threads in (4,8): cases.append((mode,threads,threads,reads,length,quality,identity))
    for mode in ('min20','supported20'): cases.append((mode,8,64,reads,length,quality,identity))
for mode,threads in [('global',1),('global',8),('fullhash',8),('min20',8),('rare20',8)]:
    cases.append((mode,threads,1 if mode=='global' else threads,100000,200,30,99))
rows=[];args.output.parent.mkdir(parents=True,exist_ok=True)
with args.output.with_suffix('.log').open('w') as log:
    for repeat in range(1,args.repeats+1):
        for case in cases:
            result=subprocess.run([str(args.binary.resolve()),*map(str,case)],text=True,stdout=subprocess.PIPE,stderr=log,check=True)
            row=list(csv.DictReader(result.stdout.splitlines()))[0];row['repeat']=repeat;rows.append(row)
            print('repeat',repeat,case,'seconds',row['total_seconds'],'groups',row['local_clusters'],row['final_clusters'],'outside',row['outside_final_anchor'],flush=True)
with args.output.open('w',newline='') as output:
    writer=csv.DictWriter(output,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
