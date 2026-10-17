#!/usr/bin/env python3
"""Audit four dereplication configurations on 1,000 mixed-quality read pairs.

Conservation and retained-100/coarse agreement are required. Only the documented empty-mate exclusion may differ from ordinary R1 counts.
Use --observe for historical binaries, or --complete-mates --require-parity to
require full abundance parity on a fixture without that exclusion.
"""
import argparse
import hashlib
import tempfile
import json, random, re, subprocess
from collections import Counter
from pathlib import Path
from test_merged_pipeline import reverse_complement, fasta_records
from test_diff_quality_pipeline import fastq
from audit_derep_counts import audit, records as fastx_records


def create_fixture(root, complete_mates=False, consensus_errors=False):
    rng=random.Random(2026091501)
    # Each family has an independent 300 bp insert, sampled by overlapping 230 bp mates.
    settings=dict(minSeqLength=200,maxSeqLength=1000,TruncateSequenceLength=-1,minAvgQuality=27,
        QualWindowThreshhold=0,TrimWindowThreshhold=0,maxAccumulatedError=-1,BinErrorModelAlpha=-1,
        maxAmbiguousNT=0,maxHomonucleotide=0,trimHomonucleotide=0,fastqVersion=1)
    maptext='#SampleID\tfastqFile\n'
    raw_mates=[set(),set()]
    raw_pairs=set()
    for sample in range(2):
        templates={f:''.join(rng.choices('ACGT',k=300)) for f in ('pass','r1only','r2onlynew','r1first','failfirst','prefix','zero','split')}
        specs=[('pass','both',40,40,40,230),('r1only','r1only',60,40,10,230),
            ('pass','r2only',40,10,40,230),('pass','r1only',40,40,10,230),('pass','neither',30,10,10,230),
            ('r2onlynew','r2only',20,10,40,230),('r1first','r1only',30,40,10,230),('r1first','both',30,45,45,230),
            ('failfirst','neither',30,10,10,230),('failfirst','both',30,45,45,230),
            ('split','a',30,40,40,230),('split','b',30,40,40,230),
            ('prefix','long',30,40,40,230),('prefix','short',30,40,40,200),('zero','empty_r2',20,40,40,230)]
        records=[[],[]]
        n=0
        for family,kind,copies,q1,q2,length in specs:
            sequence=templates[family]
            if family=='split' and kind=='b':
                sequence=sequence[:270]+next(base for base in 'ACGT' if base!=sequence[270])+sequence[271:]
            first=sequence[:length];second=reverse_complement(sequence[-230:]) if family!='zero' or complete_mates else ''
            error_position=None
            if consensus_errors and family=='split':
                error_position=120 if kind=='a' else 140
                first=first[:error_position]+next(base for base in 'ACGT' if base!=first[error_position])+first[error_position+1:]
            for copy in range(copies):
                for mate,(seq,q) in enumerate(((first,q1),('ACGT'*5+second,q2))):
                    quality=chr(33+q)*len(seq)
                    if mate==0 and error_position is not None:
                        quality=quality[:error_position]+'5'+quality[error_position+1:]  # Q20 mismatch, Q40 correct mate.
                    records[mate].append(f'@S{sample}_f{family}_{kind}_{copy}/{mate+1}\n{seq}\n+\n'+quality+'\n')
                n+=1
        for i in range(10):
            sequence=''.join(rng.choices('ACGT',k=300))
            if i<2:
                sequence=templates['pass'];p=31+i
                sequence=sequence[:p]+next(base for base in 'ACGT' if base!=sequence[p])+sequence[p+1:]
            for mate,seq in enumerate((sequence[:230],'ACGT'*5+reverse_complement(sequence[-230:]))):
                records[mate].append(f'@S{sample}_fsingleton{i}_both_0/{mate+1}\n{seq}\n+\n'+'I'*len(seq)+'\n')
            n+=1
        assert n==500,n
        for first_record,second_record in zip(*records):
            first=first_record.splitlines();second=second_record.splitlines()
            a=(first[1],first[3]);b=(second[1][20:],second[3][20:])
            raw_mates[0].add(a);raw_mates[1].add(b);raw_pairs.add((*a,*b))
        for mate in (1,2):(root/f'S{sample}.{mate}.fq').write_text(''.join(records[mate-1]))
        maptext+=f'S{sample}\tS{sample}.1.fq,S{sample}.2.fq\n'
    (root/'input.map').write_text(maptext)
    return settings, raw_mates, raw_pairs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('binary',nargs='?',type=Path,default=Path('./sdm'))
    parser.add_argument('--output',type=Path,help='Fresh directory for commands, fixtures and results')
    parser.add_argument('--identity',default='97')
    parser.add_argument('--consensus-errors',action='store_true',help='Introduce R1-only overlap errors that merging can correct')
    parser.add_argument('--main-format',choices=('fa','fq'),default='fa')
    parser.add_argument('--require-parity',action='store_true')
    parser.add_argument('--observe',action='store_true',help='Report historical admission gaps without enforcing the fixes')
    parser.add_argument('--complete-mates',action='store_true',help='Keep the 40 otherwise empty R2 mates')
    args=parser.parse_args()
    binary=args.binary.resolve()
    root=args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix='sdm-derep-parity-'))
    if args.output: root.mkdir(parents=True,exist_ok=False)
    settings,raw_mates,raw_pairs=create_fixture(root,args.complete_mates,args.consensus_errors)
    mode_options={
     'ordinary_merged':['-merge_pairs_derep','1','-merge_pairs_filter','1','-merge_pairs_demulti','1'],
     'ordinary_r1':['-merge_pairs_derep','0','-merge_pairs_filter','0','-merge_pairs_demulti','0','-noSearchWithMerge','1'],
     'retained_100':['-merge_pairs_derep','0','-merge_pairs_filter','0','-merge_pairs_demulti','0','-derepStoreQuals','1','-derepIdentity','100'],
     'retained_coarse':['-merge_pairs_derep','0','-merge_pairs_filter','0','-merge_pairs_demulti','0','-derepStoreQuals','1','-derepIdentity',args.identity]}
    results=[]
    parent_sequences={}
    for profile,cap,truncate,trim in (('capped',200,-1,0),('prefix',-1,-1,0),('trimmed',200,200,25),('alternate',200,-1,0)):
     if profile=='alternate': settings['*minAvgQuality']=8
     settings.update(TruncateSequenceLength=truncate,TrimWindowThreshhold=trim,TrimWindowWidth=20)
     config=root/(profile+'.options')
     config.write_text(''.join(f'{k}\t{v}\n' for k,v in settings.items()))
     for threads in (1,4,12):
      for mode,options in mode_options.items():
        folder=root/f'{profile}_{threads}_{mode}';folder.mkdir(exist_ok=True)
        command=[str(binary),'-i_path',str(root),'-map',str(root/'input.map'),'-options',str(config),
         '-paired','2','-o_dereplicate',str(folder/'derep.fna'),'-o_fna',f'{folder}/filtered.1.fna,{folder}/filtered.2.fna',
         '-o_fna2',f'{folder}/filtered.1.add.fna,{folder}/filtered.2.add.fna','-threads',str(threads),'-threadIO','0',
         '-derep_format',args.main_format,'-derepSrchLen',str(cap),'-min_derep_copies','40','-derep_reserve','1','-sample_sep','___',
         '-5PR2cut','20','-log',str(folder/'filter.log'),*options]
        result=subprocess.run(command,text=True,capture_output=True,timeout=120)
        (folder/'command.json').write_text(json.dumps(command,indent=2)+'\n');(folder/'console.log').write_text(result.stdout+result.stderr)
        assert result.returncode==0,(mode,result.stderr[-2000:])
        counts=Counter();categories=Counter();mapids={}
        for line in (folder/'derep.map').read_text().splitlines():
            if line.startswith('#'): continue
            header,*entries=line.split('\t');pairs=[entry.split(':') for entry in entries if ':' in entry]
            per={int(k):int(v) for k,v in pairs};total=sum(per.values());mapids[header]=total
            family=re.search(r'_f([^_]+)_',header)[1]
            categories[family]+=total;counts.update(per)
        outputs={};seen={}
        sequences={}
        for suffix in ('.fna','.merg.fna','.fna.rest'):
            p=folder/('derep'+suffix)
            records=[]
            if p.exists():
                records=fasta_records(p) if args.main_format=='fa' else [('>'+name,seq) for name,seq,_,_ in fastq(p)]
            outputs[suffix]=sum(int(h.split(';size=')[1].split(';')[0]) for h,_ in records)
            sequences[suffix]=sorted((seq,int(h.split(';size=')[1].split(';')[0]))
                for h,seq in records if args.complete_mates or '_fzero_' not in h)
            for h,seq in records:
                assert h[1:] in mapids and h not in seen
                seen[h]=True
        parent_sequences[profile,threads,mode]=sequences
        assert len(seen)==len(mapids),(mode,len(seen),len(mapids))
        assert sum(outputs.values())==sum(counts.values())
        filt=(folder/'filter.log').read_text();recovery=re.search(r'Low-quality reads recovered \(dereplication\):\s*([\d,]+)',filt)
        row=dict(profile=profile,cap=cap,threads=threads,mode=mode,map_total=sum(counts.values()),samples=dict(counts),categories=dict(categories),outputs=outputs,
           recovered=int(recovery[1].replace(',','')) if recovery else None)
        check=audit(folder/'derep.fna')
        assert check['map_total']==row['map_total'] and check['passing']==sum(v for k,v in outputs.items() if k!='.fna.rest')
        row['passing']=check['passing'];row['rest']=check['rest']
        admission=re.search(r'Dereplication admission before prefix recovery: ([\d,]+) high-quality R1; ([\d,]+) mid-quality R1; ([\d,]+) failed R1; ([\d,]+) unmatched; ([\d,]+) deferred; ([\d,]+) published',result.stderr)
        if admission:
            admitted=list(map(lambda value:int(value.replace(',','')),admission.groups()))
            recovered_prefix=re.search(r'(?:Prefix|Exact-key) matching: ([\d,]+) non-passing',result.stderr)
            recovered_prefix=int(recovered_prefix[1].replace(',','')) if recovered_prefix else 0
            assert sum(admitted[:3])+recovered_prefix==row['map_total'], 'Final admission ledger disagrees with map'
            queries=re.search(r'Search sequence length:.*total over ([\d,]+) reads',result.stderr)
            assert sum(admitted)==int(queries[1].replace(',','')), 'Admission outcomes do not partition submitted reads'
            row['admission']=dict(zip(('high','mid','failed','unmatched','deferred','published'),admitted))
        elif not args.observe:
            raise AssertionError('Missing finalized admission diagnostics')
        if not args.observe:
            for mate in (1,2):
                for _,sequence,quality,_ in fastq(folder/f'derep.{mate}.hq.fq',allow_unannotated=True):
                    assert (sequence,quality) in raw_mates[mate-1], 'HQ output lost full observed bases or qualities'
        row['hq_lengths']={mate:dict(Counter(length for _,length in fastx_records(folder/f'derep.{mate}.hq.fq'))) for mate in (1,2)}
        if not args.observe:
            assert f"{row['map_total']:,} total in map; {row['passing']:,} passing; {row['rest']:,} in rest" in result.stderr
        row['filtering']={name:re.search(re.escape(name)+r':\s*(.*)',filt)[1] for name in
           ('Reads processed','Rejected','Accepted (high quality)','Accepted (mid quality)','Reads trimmed at the end')}
        found=re.search(r'Coarse search groups: ([\d,]+);',result.stderr)
        row['search_groups']=int(found[1].replace(',','')) if found else None
        if mode.startswith('retained'):
            reuse=re.search(r'Coarse exact reuse: ([\d,]+) of ([\d,]+) queries; ([\d,]+) candidate searches',result.stderr)
            exact_hits,queries,candidates=(int(value.replace(',','')) for value in reuse.groups())
            shortcut=re.search(r'Coarse admission shortcut: ([\d,]+)',result.stderr)
            rejected=int(shortcut[1].replace(',','')) if shortcut else 0
            assert queries==exact_hits+candidates+rejected
            row['coarse_work']=dict(exact_hits=exact_hits,queries=queries,candidates=candidates,rejected=rejected)
            if not args.observe:
                assert rejected==row['admission']['unmatched'], 'Shortcut changed exact-parent admission'
                assert queries==sum(admitted[:4]), 'Deferred reads unexpectedly entered coarse search'
            r1=list(fastq(folder/'derep.1.hq.fq'));r2=list(fastq(folder/'derep.2.hq.fq'))
            assert len(r1)==len(r2)
            observed=Counter()
            for a,b in zip(r1,r2):
                assert a[0]==b[0] and a[3]==b[3] and a[1] and b[1]
                if not args.observe: assert (a[1],a[2],b[1],b[2]) in raw_pairs, 'HQ mates no longer belong to one full observation'
                observed[a[0].split(';size=')[0].rsplit('.sub',1)[0]]+=a[3]
            assert observed==Counter({name.split(';size=')[0]:n for name,n in mapids.items()})
            row['retained_pairs']=len(r1)
            # Equal-scoring copies may supply different source IDs under scheduling.
            # Compare linked sequences, observed qualities and counts; validate IDs above.
            row['hq_sha256']=hashlib.sha256(json.dumps(sorted((a[1:],b[1:]) for a,b in zip(r1,r2))).encode()).hexdigest()
            if not args.observe:
                assert f'Retention exclusions: {0 if args.complete_mates else 40} pairs with empty R2' in result.stderr
        results.append(row)
        print(f'{profile:7} {threads:2}t {mode:17} map={row["map_total"]:4} passing={row["passing"]:4} rest={row["rest"]:3}',flush=True)
    # Counts must be independent of scheduling as well as internal identity.
    by_profile_mode={}
    for row in results:
        key=(row['profile'],row['mode'])
        baseline=by_profile_mode.setdefault(key,row)
        for field in ('map_total','passing','rest','samples','categories'):
            assert row[field]==baseline[field], f'Thread count changed {key}: {field}'
    gaps=[]
    for i in range(0,len(results),4):
        merged,ordinary,exact,coarse=results[i:i+4]
        assert all(row['filtering']==ordinary['filtering'] for row in (merged,exact,coarse)), 'Filtering outcomes changed'
        for field in ('map_total','samples','categories','outputs','retained_pairs','hq_sha256'):
            assert exact[field]==coarse[field], f'Internal coarse grouping changed {field}'
        if args.consensus_errors and ordinary['profile'] in ('capped','trimmed','alternate'):
            assert merged['passing'] > ordinary['passing'], 'Merged consensus did not repair the overlap error families'
        if not args.observe:
            expected=dict(ordinary['categories'])
            if not args.complete_mates: expected.pop('zero',None)
            assert exact['categories']==expected, (exact['profile'],'Unexpected retention admission difference',expected,exact['categories'])
            assert ordinary['passing']==exact['passing'], 'Retention changed minimum-copy decisions'
            assert parent_sequences[ordinary['profile'],ordinary['threads'],'ordinary_r1']==parent_sequences[exact['profile'],exact['threads'],'retained_100'], 'Retention changed representative sequences or counts'
        if exact['threads']==1 and float(args.identity)<=97:
            assert coarse['search_groups']<exact['search_groups'], 'Fixture did not exercise approximate grouping'
        if any(ordinary[k]!=exact[k] for k in ('map_total','passing','rest','categories')):
            gaps.append(dict(profile=exact['profile'],threads=exact['threads'],
                map_difference=ordinary['map_total']-exact['map_total'],
                passing_difference=ordinary['passing']-exact['passing'],rest_difference=ordinary['rest']-exact['rest']))
    report=dict(binary=str(binary),sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),input_fragments=1000,
        samples=2,consensus_errors=args.consensus_errors,complete_mates=args.complete_mates,main_format=args.main_format,runs=results,ordinary_retention_gaps=gaps)
    (root/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS count reconciliation and retained identity comparison; {len(gaps)} ordinary/retention gaps')
    print(f'Inputs, commands, logs and results: {root}')
    if args.require_parity: assert not gaps, 'Ordinary/retained behavior differs; see results.json' 


if __name__ == '__main__':
    main()
