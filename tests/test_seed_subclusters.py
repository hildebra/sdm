#!/usr/bin/env python3
"""Preprocess real multi-sample reads, then extend seeds from exact FASTQ exports.

Run: python3 tests/test_seed_subclusters.py /path/to/sdm [legacy-binary]
"""
from collections import Counter, defaultdict
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile

from test_diff_quality_pipeline import fastq
from test_merged_pipeline import fasta_records, reverse_complement


def run(binary, folder, label, args, error=None):
    result = subprocess.run([str(binary), *map(str, args)], cwd=folder, text=True,
                            capture_output=True, timeout=60)
    (folder / f'{label}.console').write_text(result.stdout + result.stderr)
    if error is None:
        assert result.returncode == 0, f'{label}: {result.stderr}'
    else:
        assert result.returncode != 0 and error in result.stderr, f'{label}: {result.stderr}'
    return result


def fq_records(path):
    lines = path.read_text().splitlines()
    assert len(lines) % 4 == 0, path
    return [(lines[i][1:], lines[i+1], lines[i+3]) for i in range(0, len(lines), 4)]


def write_fq(path, records):
    path.write_text(''.join(f'@{name}\n{seq}\n+\n{quality}\n' for name, seq, quality in records))


def parent_id(name):
    return name.split(';size=')[0].rsplit('.sub', 1)[0]


def paf(query, target, identity=100, coverage=1):
    length = round(200 * coverage)
    return f'{query}\t200\t0\t{length}\t+\t{target}\t200\t0\t{length}\t{round(length*identity/100)}\t{length}\t60\n'


def uc(query, target, identity=100, kind='H'):
    return f'{kind}\t0\t200\t{identity}\t+\t0\t0\t*\t{query}\t{target}\n'


def fixture(binary, root, paired, identity, cap=100, physical_tags=False):
    folder = root / (f'{"paired" if paired else "single"}-{identity}-cap{cap}' + ('-tags-trimmed' if physical_tags else ''))
    folder.mkdir()
    rng = random.Random(9814)
    templates = {label: ''.join(rng.choices('ACGT', k=350)) for label in 'ABCDEF'}
    seqs = {}
    for label, sequence in templates.items():
        seqs[label] = (sequence[:200], reverse_complement(sequence[-200:])) if paired else (sequence[:200],)
    variant = list(seqs['A'])
    changing_mate = int(paired)
    mutated = list(variant[changing_mate])
    # Past the search cap, and outside the paired overlap; 100% prefix parents
    # still contain two exact full-length variants.
    # Text export orders variants by sequence. Keep the better seed later so
    # this still tests replacement after the first candidate was evaluated.
    position = next(p for p in range(120, 150) if mutated[p] != 'T')
    mutated[position] = 'T'
    variant[changing_mate] = ''.join(mutated)
    seqs['a'] = tuple(variant)
    # Two samples, a 3-fragment parent with 2 exact variants, another parent,
    # a below-cutoff parent, chimera, unassigned fallback, and count-only recovery.
    observations = [
        [('A', 35), ('a', 45), ('B', 39), ('B', 39), ('D', 40), ('D', 40), ('E', 38), ('F', 40)],
        [('A', 35), ('B', 39), ('C', 42), ('D', 40), ('D', 40), ('E', 38), ('F', 40)],
    ]
    raw_expected = Counter()
    primers=('AGGTCAGTACCGTAAC','TCCAGATGCTACGTCA')
    barcodes=(('GACTGA','CTGTAC'),('TACCTG','AGTCGA'))
    mapping = '#SampleID\tfastqFile' + ('\tBarcodeSequence\tBarcode2ndPair\tLinkerPrimerSequence\tReversePrimer' if physical_tags else '') + '\n'
    for sample, reads in enumerate(observations):
        paths = []
        for mate in range(2 if paired else 1):
            path = folder / f'S{sample}.{mate+1}.fq'
            records = []
            for i, (label, q) in enumerate(reads):
                # Include names ending in .subN to check only the generated suffix is removed.
                name = f'S{sample}r{i}.sub9' + (f'/{mate+1}' if paired else '')
                quality = ''.join(chr(33 + q + j % 3) for j in range(200))
                prefix='ACGTA' + (barcodes[sample][mate]+primers[mate] if physical_tags else '')
                if physical_tags: quality=quality[:-20]+'-'*20  # Q12 tail must survive logical quality trimming.
                records.append((name, prefix + seqs[label][mate], 'I'*len(prefix) + quality))
                if mate == 0:
                    raw_expected[seqs[label]] += 1
            write_fq(path, records)
            paths.append(path)
        mapping += f'S{sample}\t' + ','.join(p.name for p in paths)
        if physical_tags: mapping += '\t'+'\t'.join((*barcodes[sample],*primers))
        mapping += '\n'
    (folder / 'samples.tsv').write_text(mapping)
    settings = dict(minSeqLength=50, maxSeqLength=1000, TruncateSequenceLength=-1,
                    minAvgQuality=0, QualWindowThreshhold=0, TrimWindowThreshhold=0,
                    maxAccumulatedError=-1, BinErrorModelAlpha=-1, maxAmbiguousNT=0,
                    maxHomonucleotide=0, trimHomonucleotide=0, fastqVersion=1,
                    keepPrimerSeq=1)
    if physical_tags:
        settings.update(TruncateSequenceLength=100,keepPrimerSeq=0,keepBarcodeSeq=0,
                        TrimWindowThreshhold=25,TrimWindowWidth=10)
    (folder / 'pre.options').write_text(''.join(f'{k}\t{v}\n' for k, v in settings.items()))
    common = ['-paired', 2 if paired else 1, '-i_qual_offset', 33, '-o_qual_offset', 33,
              '-sample_sep', '___', '-threads', 1, '-threadIO', 0]
    run(binary, folder, 'preprocess', [*common, '-i_path', folder, '-map', folder/'samples.tsv',
        '-options', folder/'pre.options', '-o_dereplicate', folder/'derep.fna',
        '-o_fastq', ','.join(str(folder/f'filtered.{m+1}.fq') for m in range(2 if paired else 1)),
        '-derepStoreQuals', 1, '-derepStoreDiffs', 0, '-derepIdentity', identity,
        '-derepSrchLen', cap, '-min_derep_copies', 2, '-5PR1cut', 5, '-5PR2cut', 5,
        '-derep_reserve', 1, '-log', folder/'pre.log'])
    assert not (folder/'derep.diff').exists()
    paths = [folder/f'derep.{m+1}.hq.fq' for m in range(2)] if paired else [folder/'derep.1.hq.fq']
    records = [list(fastq(path)) for path in paths]
    actual = Counter()
    variants = defaultdict(list)
    for row in zip(*records):
        assert len({r[0] for r in row}) == 1
        pair = tuple(r[1] for r in row)
        actual[pair] += row[0][3]
        variants[parent_id(row[0][0])].append(row)
    assert actual == raw_expected, (actual, raw_expected)
    map_lines = (folder/'derep.map').read_text().splitlines()
    sample_names = {int(x.split(':')[0]): x.split(':')[1] for x in map_lines[0].split('\t')[1:]}
    parents = {}
    counts = {}
    for line in map_lines[1:]:
        header, *entries = line.split('\t')
        key = header.split(';size=')[0]
        counts[key] = Counter({sample_names[int(s)]: int(n) for x in entries if ':' in x for s, n in [x.split(':')]})
        for label, sequences in seqs.items():
            if any(tuple(r[1] for r in row) == sequences for row in variants[key]):
                parents[label] = header
    assert sum(sum(c.values()) for c in counts.values()) == 15
    # R2-only variants share the ordinary R1 dereplicate at every identity.
    if paired or cap == 100:
        assert parents['A'] == parents['a']
        a = variants[parents['A'].split(';size=')[0]]
        assert len(a) == 2 and a[0][0][3] == 2 and a[1][0][3] == 1
        assert tuple(r[1] for r in a[1]) == seqs['a'], 'Better variant must be later'
    else:
        assert parents['A'] != parents['a']
    assert parents['C'] in (folder/'derep.fna.rest').read_text(), 'Below-cutoff parent missing'
    settings['TruncateSequenceLength'] = 60
    settings['TrimWindowThreshhold'] = 50
    settings['maxAccumulatedError'] = 0.001
    (folder/'seed.options').write_text(''.join(f'{k}\t{v}\n' for k, v in settings.items()))
    return dict(folder=folder, paired=paired, identity=identity, cap=cap, common=common, paths=paths,
                variants=variants, parents=parents, counts=counts, sequences=seqs)


def seed(binary, f, fmt='paf', merge=False, reordered=False, chimeras=True, extra=(), error=None, tag=None, reference=False):
    folder, parents, seqs = f['folder'], f['parents'], f['sequences']
    name = tag or f'{fmt}-merge{int(merge)}-reorder{int(reordered)}-chim{int(chimeras)}'
    out = folder/name
    out.mkdir(exist_ok=True)
    # Existing parent representatives are mapper targets, exactly as in LotuS.
    target_a, target_b, target_e = [parents[x].split(';size=')[0] for x in 'ABE']
    (out/'fallback.fna').write_text(''.join(f'>{target}\n{seqs[label][0][:100]}\n'
                                           for target, label in zip((target_a,target_b,target_e), 'ABE')))
    assignments = [(parents['A'], target_a), (parents['B'], target_b)]
    if parents['a'] != parents['A']:
        assignments.append((parents['a'], target_a))
    if reordered:
        assignments.reverse()
    flags = []
    suffix = '.paf' if fmt == 'paf' else '.uc'
    if fmt == 'paf':
        text = ''.join(paf(q, t) for q, t in assignments)
        # A repeated parent and rejected hits must not count again or select E.
        text += paf(*assignments[0]) + paf(parents['E'], target_a, 90) + paf(parents['E'], target_a, 100, .5)
    elif fmt in ('uc', 'vsearch'):
        text = ''.join(uc(q, t) for q, t in assignments)
        flags = ['-uparseVer', 'vsearch' if fmt == 'vsearch' else '7']
    elif fmt == 'uparse8':
        text = ''.join(f'{q}\tmatch\t100\t*\t{t}\n' for q, t in assignments)
        text += f'{parents["D"]}\tchimera\t100\t*\t{target_a}(1-100)+{target_b}(101-200)\n'
        flags = ['-uparseVer', '8']
    elif fmt in ('uparse9', 'uparse12'):
        text = ''.join(f'{q}\tmatch\tdqt=0;top={t};size=1;(100%);\n' for q, t in assignments)
        text += f'{parents["E"]}\tperfect_chimera\t*\n'
        flags = ['-uparseVer', '9' if fmt == 'uparse9' else '12']
    elif fmt == 'cdhit':
        text = ''
        for i, (q, t) in enumerate(assignments):
            # CD-HIT represents each cluster with a parent ID. One non-rep line
            # ahead of its rep also exercises the parser's queued-line path.
            text += f'>Cluster {i}\n0\t200nt, >{q}... at +/100%\n1\t200nt, >{q}... *\n'
        flags = ['-uparseVer', 'cdhit']
    elif fmt == 'empty':
        text = ''
    else:
        raise AssertionError(fmt)
    mapping = out/('assign'+suffix)
    mapping.write_text(text)
    # Medium-quality counts come from headers; .rest/reference recovery uses the
    # remaining parent sample vectors. Duplicate REST/primary hits add nothing.
    (out/'add.paf').write_text(paf('S0___medium0', target_a) + paf('S1___medium1', target_b))
    rest = paf(parents['C'], target_a) * 2 + paf(parents['A'], target_a)
    if reference:
        (out/'ref.fna').write_text(f'>reference_F\n{seqs["F"][0][:100]}\n')
        (out/'ref.paf').write_text(paf(parents['F'], 'reference_F') * 2)
        (out/'ref-add.paf').write_text(paf('S1___medium_ref', 'reference_F'))
        (out/'ref-rest.paf').write_text(paf(parents['F'], 'reference_F'))
        flags += ['-OTU_fallback_refclust', out/'ref.fna',
                  '-optimalRead2Cluster_ref', out/'ref.paf',
                  '-ucAdditionalCounts_refclust', out/'ref-add.paf',
                  '-ucAdditionalCounts_refclust1', out/'ref-rest.paf']
    else:
        rest += paf(parents['F'], target_b)
    (out/'rest.paf').write_text(rest)
    outputs = [out/f'seed.{i+1}.fq' for i in range(2 if f['paired'] else 1)]
    args = [*f['common'], *flags, '-i_fastq', ','.join(map(str, f['paths'])),
            '-options', folder/'seed.options', '-seedSubclusters', 1,
            '-optimalRead2Cluster', mapping, '-OTU_fallback', out/'fallback.fna',
            '-derep_map', folder/'derep.map', '-o_fastq', ','.join(map(str, outputs)),
            '-merge_pairs_seed', int(merge), '-otu_matrix', out/'matrix.tsv',
            '-ucAdditionalCounts', out/'add.paf', '-ucAdditionalCounts1', out/'rest.paf',
            '-count_chimeras', 'T' if chimeras else 'F', '-id', f['identity'],
            '-5PR1cut', 70, '-5PR2cut', 70, '-log', out/'seed.log', *extra]
    result = run(binary, folder, name, args, error)
    if error:
        return
    lines = (out/'matrix.tsv').read_text().splitlines()
    samples = lines[0].split('\t')[1:]
    matrices = {line.split('\t')[0]: dict(zip(samples, map(int, line.split('\t')[1:]))) for line in lines[1:]}
    # Associate renamed OTUs with their distinguishable output sequences.
    selected = fq_records(out/'seed.merg.fq') if merge else fq_records(outputs[0])
    if f['paired']:
        selected += fq_records(out/'seed.1.singl.fq')
    by_label = {}
    for otu, sequence, quality in selected:
        matches = [x for x in 'ABE' if sequence.startswith(seqs[x][0][:90])]
        assert len(matches) == 1, (otu, sequence, matches)
        by_label[matches[0]] = (otu.removesuffix('.1'), sequence, quality)
    assert set(by_label) == set('ABE'), (fmt, by_label)
    expected_a = Counter(S0=3, S1=2)
    expected_b = Counter(S0=2, S1=2) if reference else Counter(S0=3, S1=3)
    if fmt == 'uparse8' and chimeras:
        # D has two observations per sample, split equally over A and B.
        expected_a.update(S0=1, S1=1)
        expected_b.update(S0=1, S1=1)
    if fmt == 'empty':
        # A can still be counted through REST; B has only recovered F + medium.
        expected_b = Counter(S0=1, S1=2)
    for label, counts in [('A', expected_a), ('B', expected_b), ('E', Counter(S0=0, S1=0))]:
        assert matrices[by_label[label][0]] == counts, (fmt, label, matrices, counts)
    if reference:
        ref_rows = set(matrices) - {row[0] for row in by_label.values()}
        assert len(ref_rows) == 1
        assert matrices[ref_rows.pop()] == dict(S0=1, S1=2), 'Reference count recovery changed'
    if fmt != 'empty':
        expected_seq = seqs['a'][0]
        if merge:
            # The R2-only mutation lies outside R1; both full mates contribute.
            expected_seq = seqs['a'][0][:150] + reverse_complement(seqs['a'][1])
        assert by_label['A'][1] == expected_seq, (fmt, merge, len(by_label['A'][1]), expected_seq, by_label['A'][1])
        if not merge:
            winner = next(row for rows in f['variants'].values() for row in rows if tuple(r[1] for r in row) == seqs['a'])
            assert by_label['A'][2] == winner[0][2], 'R1 qualities no longer match selected observation'
            if f['paired']:
                mate = {x[0].removesuffix('.2'): x[1:] for x in fq_records(outputs[1])}
                assert mate[by_label['A'][0]] == (winner[1][1], winner[1][2]), 'R2 or its quality vector was split from R1'
        log = (out/'seed.log').read_text()
        if f['paired']:
            assert 'seed candidates had a paired-read overlap' in log
    assert by_label['E'][1] == seqs['E'][0][:100], 'Unmapped fallback changed'
    print(f'PASS {folder.name}/{name}: sample vectors, selected variants and qualities', flush=True)
    return out, args


def legacy_equivalence(binary, baseline, f, seed_args):
    """Default/explicit 0 must retain HQ output; optionally compare a saved binary."""
    out = f['folder'] / 'legacy'
    out.mkdir()
    options = dict(zip(seed_args[::2], seed_args[1::2]))
    options.pop('-seedSubclusters')
    # Quality retention now writes reconstructed variants to the normal HQ
    # paths. Generate ordinary representative HQ input separately to verify
    # that the default/zero seed reader remains compatible with that format.
    stem = out / 'derep'
    run(binary, f['folder'], 'legacy-preprocess', [*f['common'],
        '-i_path', f['folder'], '-map', f['folder']/'samples.tsv',
        '-options', f['folder']/'pre.options', '-o_dereplicate', str(stem)+'.fna',
        '-o_fastq', ','.join(str(out/f'filtered.{m+1}.fq') for m in range(2 if f['paired'] else 1)),
        '-derepStoreQuals', 0, '-derepStoreDiffs', 0, '-derepIdentity', f['identity'],
        '-derepSrchLen', f['cap'], '-min_derep_copies', 2, '-5PR1cut', 5, '-5PR2cut', 5,
        '-derep_reserve', 1, '-log', out/'pre.log'])
    options['-derep_map'] = str(stem)+'.map'
    options['-i_fastq'] = ','.join(str(stem)+f'.{m+1}.hq.fq' for m in range(2 if f['paired'] else 1))
    options['-o_fastq'] = ','.join(str(out/f'seed.{m+1}.fq') for m in range(2 if f['paired'] else 1))
    options['-otu_matrix'] = out/'matrix.tsv'
    options['-log'] = out/'seed.log'
    def execute(executable, label):
        run(executable, f['folder'], label, [v for pair in options.items() for v in pair])
        return {p.name: p.read_bytes() for p in out.iterdir() if p.suffix == '.fq' or p.name == 'matrix.tsv'}
    default = execute(binary, 'legacy-default')
    options['-seedSubclusters'] = 0
    assert execute(binary, 'legacy-explicit-zero') == default
    if baseline:
        del options['-seedSubclusters']
        assert execute(baseline, 'legacy-before-change') == default, 'Legacy HQ output changed'
    print(f'PASS {f["folder"].name}: legacy HQ default/0' + ('/pre-change output identical' if baseline else ''), flush=True)


def reordered_exports(binary, f):
    """Parent map order, mapper order, and FASTQ group order are independent."""
    folder = f['folder']
    lines = (folder/'derep.map').read_text().splitlines()
    shuffled_map = folder/'reordered.map'
    shuffled_map.write_text('\n'.join([lines[0], *reversed(lines[1:])]) + '\n')
    paths = []
    for mate, source in enumerate(f['paths']):
        groups = defaultdict(list)
        for record in fq_records(source):
            groups[parent_id(record[0])].append(record)
        path = folder/f'reordered.{mate+1}.fq'
        write_fq(path, [record for group in reversed(list(groups.values())) for record in group])
        paths.append(path)
    seed(binary, f, tag='independent-orders', reordered=True,
         extra=['-i_fastq', ','.join(map(str, paths)), '-derep_map', shuffled_map])


def invalid_cases(binary, f):
    base = f['folder']/'invalid'
    base.mkdir()
    original = [fq_records(p) for p in f['paths']]
    cases = {
        'missing-last-mate': (original[1][:-1], 'ended at different records'),
        'extra-mate': (original[1]+[original[1][-1]], 'ended at different records'),
        'reordered-mates': ([original[1][1],original[1][0],*original[1][2:]], 'headers do not match'),
        'wrong-count': ([(original[1][0][0].split(';size=')[0] + ';size=99;', *original[1][0][1:]), *original[1][1:]], 'headers do not match'),
        'short-quality': ([(original[1][0][0],original[1][0][1],original[1][0][2][:-1]), *original[1][1:]], 'length mismatch'),
        'empty-mate': ([], 'Could not read complete subcluster'),
    }
    for name, (records, error) in cases.items():
        paths = [base/'r1.fq', base/'r2.fq']
        write_fq(paths[0], original[0]); write_fq(paths[1], records)
        seed(binary, f, tag=name, extra=['-i_fastq', ','.join(map(str, paths))], error=error)
        print(f'PASS rejects {name}', flush=True)
    seed(binary, f, tag='missing-mate-file', extra=['-i_fastq', f'{base}/r1.fq,{base}/absent.fq'], error='Could not find file')
    seed(binary, f, tag='invalid-flag', extra=['-seedSubclusters', 2], error='-seedSubclusters must be 0 or 1')
    seed(binary, f, tag='missing-map', extra=['-derep_map', ''], error='requires the parent -derep_map')
    # Damage both mates equally: pair synchronization alone cannot find this.
    for name, rows, error in [
        ('missing-variant', [[r for r in mate if '.sub2;' not in r[0]] for mate in original], 'Incomplete subclusters'),
        ('duplicate-variant', [[mate[0], *mate] for mate in original], 'Invalid subcluster order'),
        ('unknown-parent', [[('unknown.sub1;size=2;', *mate[0][1:]), *mate[1:]] for mate in original], 'No derep.map parent'),
    ]:
        for path, records in zip(paths, rows): write_fq(path, records)
        seed(binary, f, tag=name, extra=['-i_fastq', ','.join(map(str, paths))], error=error)
        print(f'PASS rejects {name}', flush=True)


def main():
    binary = Path(sys.argv[1] if len(sys.argv)>1 else './sdm').resolve()
    baseline = Path(sys.argv[2]).resolve() if len(sys.argv)>2 else None
    root = Path(tempfile.mkdtemp(prefix='sdm-seed-subclusters-'))
    try:
        help_text = subprocess.run([binary, '-help_flags'], text=True, capture_output=True, check=True).stdout
        assert '-seedSubclusters <0|1>' in help_text and '(default: 0)' in help_text
        for identity in (95, 97.5, 100):
            for paired in (False, True):
                f = fixture(binary, root, paired, identity)
                _, args = seed(binary, f)
                legacy_equivalence(binary, baseline, f, args)
                seed(binary, f, reordered=True)
                if paired:
                    seed(binary, f, merge=True, reordered=True)
                if identity == 97.5:
                    for fmt in ('uc', 'vsearch', 'uparse8', 'uparse9', 'uparse12', 'cdhit', 'empty'):
                        seed(binary, f, fmt=fmt, reordered=True)
                    seed(binary, f, fmt='uparse8', chimeras=False)
                    seed(binary, f, reference=True, tag='reference-count-recovery')
                    reordered_exports(binary, f)
                    if paired:
                        invalid_cases(binary, f)
        for identity in (95, 97.5, 100):
            for paired in (False, True):
                f = fixture(binary, root, paired, identity, cap=-1)
                seed(binary, f, reordered=True)
        # Actual barcodes/primers and fixed 5' cuts are removed once; logical
        # truncation and low-quality tails must not shorten either HQ workflow.
        for identity in (97,100):
            f=fixture(binary,root,True,identity,physical_tags=True)
            _,args=seed(binary,f)
            legacy_equivalence(binary,None,f,args)
            seed(binary,f,merge=True)
    except BaseException:
        print(f'Regression inputs and logs retained in {root}', file=sys.stderr)
        raise
    shutil.rmtree(root)
    print('Native reconstructed seed-extension regressions passed')


if __name__ == '__main__':
    main()
