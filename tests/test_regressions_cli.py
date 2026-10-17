#!/usr/bin/env python3
"""End-to-end regressions for the 2026-09 audit fixes (sdm 3.55-3.60); pass the sdm executable.

Each case builds a small synthetic input with the truth encoded in read names, runs sdm and
checks outputs. Usage: test_regressions_cli.py ./sdm [case-name ...]
Prints PASS/FAIL per case and exits 1 if any case fails. Fixtures of failed cases are kept.
"""
from collections import Counter
import gzip
import os
from pathlib import Path
import random
import re
import shutil
import subprocess
import sys
import tempfile
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_derep_counts import audit  # noqa: E402

SDM = str(Path(sys.argv[1]).resolve()) if len(sys.argv) > 1 else './sdm'
ONLY = set(sys.argv[2:])
COMP = str.maketrans('ACGT', 'TGCA')
FWD, REV = 'GTGCCAGCAGCCGCGGTAA', 'GGACTACCAGGGTATCTAAT'
CASES = []


def case(fn):
    CASES.append(fn)
    return fn


def rc(seq):
    return seq.translate(COMP)[::-1]


def rng(seed):
    return random.Random(seed)


def dna(r, n):
    return ''.join(r.choices('ACGT', k=n))


def fq(records):
    return ''.join(f'@{h}\n{s}\n+\n{q}\n' for h, s, q in records)


def read_fq(path):
    path = Path(path)
    if not path.exists():
        return []
    lines = path.read_text().splitlines()
    return [(lines[i][1:], lines[i + 1], lines[i + 3]) for i in range(0, len(lines), 4)]


def read_fa(path):
    records, name = [], None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:]
            records.append([name, ''])
        elif records:
            records[-1][1] += line
    return [tuple(r) for r in records]


def run(d, *args, expect_rc=0):
    env = dict(os.environ)
    env.setdefault('ASAN_OPTIONS', 'detect_leaks=0:halt_on_error=1')
    env.setdefault('UBSAN_OPTIONS', 'halt_on_error=1:print_stacktrace=1')
    cmd = [SDM, *map(str, args)]
    if '-log' not in cmd:
        cmd += ['-log', str(d / 'run.log')]
    p = subprocess.run(cmd, cwd=d, capture_output=True, text=True, timeout=120, env=env)
    out = p.stdout + p.stderr
    (d / 'command.txt').write_text(' '.join(cmd) + '\n')
    (d / 'console.txt').write_text(out)
    assert 'AddressSanitizer' not in out and 'runtime error:' not in out, 'sanitizer report:\n' + out[-3000:]
    assert 'terminate called' not in out, 'sdm aborted:\n' + out[-2000:]
    if expect_rc is not None:
        assert p.returncode == expect_rc, f'exit code {p.returncode}, expected {expect_rc}\n{out[-2000:]}'
    return p.returncode, out


def opts(d, text, name='opts.txt', crlf=False):
    path = d / name
    body = text.replace('\n', '\r\n') if crlf else text
    path.write_bytes(body.encode())
    return path


def good_reads(n=50, length=150, prefix='read_r', seed=1):
    r = rng(seed)
    return [(f'{prefix}{i}', dna(r, length), ''.join(chr(33 + r.randint(30, 40)) for _ in range(length)))
            for i in range(n)]


# ---------------------------------------------------------------- input handling

@case
def max_seq_length_negative_means_no_limit(d):
    (d / 'in.fq').write_text(fq(good_reads()))
    run(d, '-i', 'in.fq', '-o_fastq', 'o.fq', '-options', opts(d, 'minSeqLength\t30\nmaxSeqLength\t-1\n'))
    assert len(read_fq(d / 'o.fq')) == 50, 'maxSeqLength -1 must not reject reads'


@case
def malformed_fastq_clean_error(d):
    reads = good_reads(10)
    (d / 'in.fq').write_text(fq(reads[:5] + [('bad', 'ACGTACGT', 'IIII')] + reads[5:]))
    _, out = run(d, '-i', 'in.fq', '-o_fastq', 'o.fq', expect_rc=1)
    assert 'Error: FASTQ sequence/quality length mismatch for @bad' in out, out[-500:]


@case
def ignore_io_errors_skips_only_bad_record(d):
    reads = good_reads(50)
    (d / 'in.fq').write_text(fq(reads[:5] + [('bad', 'ACGTACGT', 'IIII')] + reads[5:]))
    _, out = run(d, '-i', 'in.fq', '-o_fastq', 'o.fq', '-ignore_IO_errors', 1)
    assert len(read_fq(d / 'o.fq')) == 50, 'all good reads after the malformed one must be kept'
    assert 'discarding malformed record' in out


@case
def ignore_io_errors_paired_stays_synchronized(d):
    reads = good_reads(50)
    r1 = [(f'{h}/1', s, q) for h, s, q in reads]
    r2 = [(f'{h}/2', rc(s), q) for h, s, q in reads]
    r2[10] = (r2[10][0], r2[10][1], r2[10][2][:40])
    (d / 'in.1.fq').write_text(fq(r1)); (d / 'in.2.fq').write_text(fq(r2))
    run(d, '-paired', 2, '-i', 'in.1.fq,in.2.fq', '-o_fastq', 'o.1.fq,o.2.fq', '-ignore_IO_errors', 1,
        '-options', opts(d, 'minSeqLength\t30\n'))
    a, b = read_fq(d / 'o.1.fq'), read_fq(d / 'o.2.fq')
    assert len(a) == len(b) == 49, (len(a), len(b))
    assert all(x[0].split('/')[0] == y[0].split('/')[0] for x, y in zip(a, b)), 'mates out of sync'
    assert not any(h.startswith('read_r10/') for h, _, _ in a), 'the bad pair must be dropped'


@case
def empty_inputs_report_input_empty(d):
    (d / 'e.fq').write_text(''); (d / 'e2.fq').write_text('')
    with gzip.open(d / 'e.fq.gz', 'wt') as g:
        g.write('')
    for args in (['-i', 'e.fq'], ['-i', 'e.fq.gz'], ['-i_fastq', 'e.fq'], ['-paired', 2, '-i', 'e.fq,e2.fq']):
        out_arg = 'o.1.fq,o.2.fq' if '-paired' in args else 'o.fq'
        _, out = run(d, *args, '-o_fastq', out_arg, expect_rc=1)
        assert 'Input empty' in out, (args, out[-400:])
    (d / 'good.fq').write_text(fq(good_reads(5)))
    (d / 'map.tsv').write_text('#SampleID\tfastqFile\nS1\te.fq\nS2\tgood.fq\n')
    run(d, '-i_path', d, '-map', 'map.tsv', '-o_fastq', 'o.fq')  # one empty sample is only skipped


@case
def empty_fasta_record_keeps_next_record(d):
    r = rng(3)
    s2, s3 = dna(r, 40), dna(r, 40)
    (d / 'in.fa').write_text(f'>r1\n>r2\n{s2}\n>r3\n{s3}\n')
    run(d, '-i_fna', 'in.fa', '-o_fna', 'o.fna')
    got = dict(read_fa(d / 'o.fna'))
    assert got.get('r2') == s2 and got.get('r3') == s3, got


@case
def crlf_options_file_is_parsed(d):
    r = rng(4)
    reads = [(f'p{i}', FWD + dna(r, 150), 'I' * 169) for i in range(20)] + \
            [(f'n{i}', dna(r, 169), 'I' * 169) for i in range(20)]
    (d / 'in.fq').write_text(fq(reads))
    (d / 'map.tsv').write_text(f'#SampleID\tfastqFile\tForwardPrimer\nS1\tin.fq\t{FWD}\n')
    o = opts(d, 'minSeqLength\t100\nRejectSeqWithoutFwdPrim\tT\n', crlf=True)
    run(d, '-i_path', d, '-map', 'map.tsv', '-options', o, '-o_fastq', 'o.fq', '-sample_sep', '___')
    names = [h.split('___')[1].split()[0] for h, _, _ in read_fq(d / 'o.fq')]
    assert len(names) == 20 and all(n.startswith('p') for n in names), '"T\\r" must still mean T'


# ---------------------------------------------------------------- read selection and outputs

@case
def specific_reads_writes_selected_and_excluded(d):
    (d / 'in.fq').write_text(fq(good_reads()))
    (d / 'ids.txt').write_text('read_r1\nread_r5\nread_r7\n')
    run(d, '-i', 'in.fq', '-o_fastq', 'sel.fq', '-specificReads', 'ids.txt', '-excludeFile', 'excl.fq')
    selected = [p for p in d.glob('sel.fq*.fastq')]
    assert selected, 'no selected-reads file'
    assert sorted(h for h, _, _ in read_fq(selected[0])) == ['read_r1', 'read_r5', 'read_r7']
    assert len(read_fq(d / 'excl.fq')) == 47


@case
def mocat_fix_matches_illumina_stems(d):
    r = rng(5)
    reads = [(f'M00123:45:000000000-ABCDE:1:1101:{1000 + i}:{2000 + i} 1:N:0:1', dna(r, 100), 'I' * 100)
             for i in range(20)]
    (d / 'in.fq').write_text(fq(reads))
    (d / 'ids.txt').write_text('M00123:000000000-ABCDE:1:1101:1003\nM00123:000000000-ABCDE:1:1101:1007\n')
    run(d, '-i', 'in.fq', '-o_fastq', 'm.fq', '-specificReads', 'ids.txt', '-mocatFix', 1)
    assert sum(len(read_fq(p)) for p in d.glob('m.fq*.fastq')) == 2
    for p in d.glob('m.fq*'):
        p.unlink()
    run(d, '-i', 'in.fq', '-o_fastq', 'm.fq', '-specificReads', 'ids.txt', '-mocatFix', 0)
    assert sum(len(read_fq(p)) for p in d.glob('m.fq*.fastq')) == 0, '-mocatFix 0 must disable the fix'


@case
def default_log_goes_next_to_output(d):
    (d / 'in.fq').write_text(fq(good_reads(5)))
    (d / 'sub').mkdir()
    subprocess.run([SDM, '-i', 'in.fq', '-o_fastq', './sub/outf'], cwd=d, capture_output=True, timeout=60)
    assert (d / 'sub' / 'outf.log').exists(), sorted(os.listdir(d / 'sub'))
    assert not (d / '.log').exists(), 'log must not be written as ./.log'


@case
def no_barcode_output_single_path(d):
    r = rng(6)
    bc = 'ACGTACGTAC'
    reads = [(f'b{i}', bc + dna(r, 150), 'I' * 160) for i in range(20)] + \
            [(f'n{i}', 'GGCCAATTGG' + dna(r, 150), 'I' * 160) for i in range(10)]
    (d / 'in.fq').write_text(fq(reads))
    (d / 'map.tsv').write_text(f'#SampleID\tBarcodeSequence\nS1\t{bc}\n')
    run(d, '-i_fastq', 'in.fq', '-map', 'map.tsv', '-options', opts(d, 'minSeqLength\t50\n'),
        '-o_fastq', 'o.fq', '-o_fastq_noBC', 'nobc.fq')
    nobc = read_fq(d / 'nobc.fq')
    assert len(nobc) == 10 and all(h.startswith('n') for h, _, _ in nobc), len(nobc)


@case
def paired_single_additional_output_is_an_error(d):
    reads = good_reads(10)
    (d / 'in.1.fq').write_text(fq([(f'{h}/1', s, q) for h, s, q in reads]))
    (d / 'in.2.fq').write_text(fq([(f'{h}/2', s, q) for h, s, q in reads]))
    _, out = run(d, '-paired', 2, '-i', 'in.1.fq,in.2.fq', '-o_fastq', 'o.1.fq,o.2.fq', '-o_fastq2', 'add.fq',
                 '-options', opts(d, 'minSeqLength\t30\n'), expect_rc=57)
    assert 'requires paired additional output file' in out


@case
def paired_input_without_map_or_options_writes_pairs(d):
    reads = good_reads(50)
    (d / 'in.1.fq').write_text(fq([(f'{h}/1', s, q) for h, s, q in reads]))
    (d / 'in.2.fq').write_text(fq([(f'{h}/2', rc(s), q) for h, s, q in reads]))
    for extra in (['-paired', 2], []):
        for p in d.glob('o.*'):
            p.unlink()
        run(d, *extra, '-i', 'in.1.fq,in.2.fq', '-o_fastq', 'o.1.fq,o.2.fq')
        assert len(read_fq(d / 'o.1.fq')) == len(read_fq(d / 'o.2.fq')) == 50, extra
        assert not read_fq(d / 'o.1.singl.fq') and not read_fq(d / 'o.2.singl.fq'), 'pairs went to singletons'


@case
def one_file_per_sample_headers_have_no_junk_tag(d):
    (d / 'a.fq').write_text(fq(good_reads(5, prefix='a')))
    (d / 'b.fq').write_text(fq(good_reads(5, prefix='b', seed=2)))
    (d / 'map.tsv').write_text('#SampleID\tfastqFile\nSA\ta.fq\nSB\tb.fq\n')
    run(d, '-i_path', d, '-map', 'map.tsv', '-options', opts(d, 'minSeqLength\t30\n'), '-o_fastq', 'o.fq',
        '-sample_sep', '___')
    heads = [h for h, _, _ in read_fq(d / 'o.fq')]
    assert len(heads) == 10 and not any('new_bc=FileName' in h for h in heads), heads[:2]
    assert all(h.split('___')[0] == ('SA' if h.split('___')[1].startswith('a') else 'SB') for h in heads)


# ---------------------------------------------------------------- filtering semantics

@case
def barcode_failures_are_not_rescued(d):
    r = rng(7)
    bc = 'ACGTACGTAC'
    cats = {'bchi': (bc, 38), 'bcmid': (bc, 25), 'nobcmid': ('GGCCAATTGG', 25), 'nobchi': ('GGCCAATTGG', 38)}
    reads = [(f'{c}_{i}', p + dna(r, 150), chr(33 + q) * 160) for c, (p, q) in cats.items() for i in range(20)]
    (d / 'in.fq').write_text(fq(reads))
    (d / 'map.tsv').write_text(f'#SampleID\tBarcodeSequence\nS1\t{bc}\n')
    o = opts(d, 'minSeqLength\t100\nminAvgQuality\t30\n*minAvgQuality\t20\nmaxAmbiguousNT\t0\n')
    run(d, '-i_fastq', 'in.fq', '-map', 'map.tsv', '-options', o, '-o_fastq', 'main.fq', '-o_fastq2', 'add.fq',
        '-i_qual_offset', 33)
    count = lambda path: Counter(re.search(r'(bchi|bcmid|nobcmid|nobchi)_', h).group(1) for h, _, _ in read_fq(path))
    assert count(d / 'main.fq') == Counter(bchi=20), count(d / 'main.fq')
    assert count(d / 'add.fq') == Counter(bcmid=20), 'barcode-less reads must not reach the additional output'


@case
def alternative_ambiguity_threshold_is_inclusive(d):
    r = rng(8)
    reads = []
    for n_amb in (0, 1, 2):
        for i in range(10):
            s = list(dna(r, 150))
            for k in range(n_amb):
                s[20 + 40 * k] = 'N'
            reads.append((f'amb{n_amb}_{i}', ''.join(s), 'I' * 150))
    (d / 'in.fq').write_text(fq(reads))
    run(d, '-i', 'in.fq', '-options', opts(d, 'minSeqLength\t30\nmaxAmbiguousNT\t0\n*maxAmbiguousNT\t1\n'),
        '-o_fastq', 'main.fq', '-o_fastq2', 'add.fq')
    main = Counter(h.split('_')[0] for h, _, _ in read_fq(d / 'main.fq'))
    add = Counter(h.split('_')[0] for h, _, _ in read_fq(d / 'add.fq'))
    assert main == Counter(amb0=10) and add == Counter(amb1=10), (main, add)


@case
def min_derep_copies_rules(d):
    r = rng(9)
    counts = {'A': [5, 5, 0, 0, 0], 'B': [9, 0, 0, 0, 0], 'C': [3, 3, 3, 3, 0], 'D': [2, 2, 2, 2, 2], 'E': [1, 1, 1, 0, 0]}
    seqs = {k: dna(r, 150) for k in counts}
    maptext = '#SampleID\tfastqFile\n'
    for s in range(5):
        recs = [(f'S{s}_{k}_{i}', seqs[k], 'I' * 150) for k, c in counts.items() for i in range(c[s])]
        (d / f'S{s}.fq').write_text(fq(recs))
        maptext += f'S{s}\tS{s}.fq\n'
    (d / 'map.tsv').write_text(maptext)
    o = opts(d, 'minSeqLength\t100\nmaxSeqLength\t1000\nminAvgQuality\t20\nmaxAmbiguousNT\t0\n')
    expected = {'8': 'ABCD', '8:1': 'B', '3:4': 'C', '8:1,3:4': 'BC', '8:1;3:4': 'BC'}
    for index, (rule, want) in enumerate(expected.items()):
        out = d / f'rule{index}'
        out.mkdir()
        run(out, '-i_path', d, '-map', d / 'map.tsv', '-options', o, '-o_fastq', 'o.fq', '-o_dereplicate', 'derep.fna',
            '-min_derep_copies', rule, '-sample_sep', '___')
        joined = (out / 'derep.fna').read_text().replace('\n', '')
        got = ''.join(sorted(k for k, s in seqs.items() if s in joined))
        assert got == want, f'rule {rule}: passing {got}, expected {want}'


# ---------------------------------------------------------------- orientation detection

def mixed_orientation_fixture(d, n_per_sample=500):
    r = rng(10)
    templates = [dna(r, 253) for _ in range(12)]
    maptext, swapped = '#SampleID\tfastqFile\tForwardPrimer\tReversePrimer\n', 0
    for s in range(2):
        a, b = [], []
        for n in range(n_per_sample):
            t = r.randrange(len(templates))
            frag = FWD + templates[t] + rc(REV)
            sw = r.random() < 0.5
            swapped += sw
            x, y = (rc(frag)[:250], frag[:250]) if sw else (frag[:250], rc(frag)[:250])
            a.append((f'm{s}_{n}_T{t}_O{int(sw)}/1', x, 'I' * 250))
            b.append((f'm{s}_{n}_T{t}_O{int(sw)}/2', y, 'I' * 250))
        (d / f'S{s}.1.fq').write_text(fq(a)); (d / f'S{s}.2.fq').write_text(fq(b))
        maptext += f'S{s}\tS{s}.1.fq,S{s}.2.fq\t{FWD}\t{REV}\n'
    (d / 'map.tsv').write_text(maptext)
    return templates, swapped


@case
def mixed_orientation_pairs_are_swapped(d):
    templates, swapped = mixed_orientation_fixture(d)
    base = 'minSeqLength\t150\nminAvgQuality\t20\nmaxPrimerErrs\t1\nRejectSeqWithoutFwdPrim\tT\nRejectSeqWithoutRevPrim\tT\n'
    run(d, '-i_path', d, '-map', 'map.tsv', '-options', opts(d, base + 'CheckForMixedPairs\tT\n'), '-paired', 2,
        '-o_fastq', 'o.1.fq,o.2.fq', '-o_dereplicate', 'derep.fna', '-merge_pairs_derep', 1, '-min_derep_copies', 1,
        '-threads', 2, '-log', 'mixed.log', '-sample_sep', '___')
    a, b = read_fq(d / 'o.1.fq'), read_fq(d / 'o.2.fq')
    assert len(a) == len(b) == 1000, (len(a), len(b))
    for (h1, s1, _), (h2, s2, _) in zip(a, b):
        t = templates[int(re.search(r'_T(\d+)_', h1).group(1))]
        fragment = lambda h: re.sub(r'/[12].*', '', h.split('___')[-1])  # R1 carries the sample prefix
        assert fragment(h1) == fragment(h2), (h1, h2)
        assert h1.split()[0].endswith('/1') and h2.split()[0].endswith('/2'), f'mate tags after swap: {h1!r} {h2!r}'
        assert t.startswith(s1) and rc(t).startswith(s2), 'R1 must be normalized to the forward orientation'
    m = re.search(r'Read pairs swapped:\s+([\d,]+)', (d / 'mixed.log').read_text())
    assert m and int(m.group(1).replace(',', '')) == swapped, (m and m.group(1), swapped)
    merged = [s for _, s in read_fa(d / 'derep.merg.fna')]
    assert sorted(merged) == sorted(set(templates)) or set(merged) <= set(templates), 'orientation duplicates in derep'
    # Without the option the reverse-orientation half fails the primer requirement.
    for p in d.glob('o.*'):
        p.unlink()
    run(d, '-i_path', d, '-map', 'map.tsv', '-options', opts(d, base, 'plain.txt'), '-paired', 2,
        '-o_fastq', 'o.1.fq,o.2.fq', '-sample_sep', '___')
    assert len(read_fq(d / 'o.1.fq')) == 1000 - swapped


@case
def reversed_single_reads_are_reverse_complemented(d):
    r = rng(11)
    templates = [dna(r, 253) for _ in range(10)]
    reads, nrev = [], 0
    for n in range(400):
        t = r.randrange(10)
        frag = FWD + templates[t] + rc(REV)
        rev = r.random() < 0.5
        nrev += rev
        s = rc(frag) if rev else frag
        reads.append((f'r{n}_T{t}', s, 'I' * len(s)))
    (d / 'in.fq').write_text(fq(reads))
    (d / 'map.tsv').write_text(f'#SampleID\tfastqFile\tForwardPrimer\tReversePrimer\nS1\tin.fq\t{FWD}\t{REV}\n')
    o = opts(d, 'minSeqLength\t150\nminAvgQuality\t20\nmaxPrimerErrs\t1\nRejectSeqWithoutFwdPrim\tT\n'
                'RejectSeqWithoutRevPrim\tT\nAmpliconShortPE\tT\nCheckForReversedSeqs\tT\n')
    run(d, '-i_path', d, '-map', 'map.tsv', '-options', o, '-paired', 1, '-o_fastq', 'o.fq', '-log', 'rev.log',
        '-sample_sep', '___')
    out = read_fq(d / 'o.fq')
    assert len(out) == 400, len(out)
    assert all(s == templates[int(re.search(r'_T(\d+)', h).group(1))] for h, s, _ in out), 'wrong orientation/trim'
    m = re.search(r'Reads reverse-complemented:\s+([\d,]+)', (d / 'rev.log').read_text())
    assert m and int(m.group(1).replace(',', '')) == nrev, (m and m.group(1), nrev)


@case
def illumina_mate_tags_follow_output_position(d):
    r = rng(13)
    templates = [dna(r, 253) for _ in range(4)]
    a, b = [], []
    for n in range(40):
        frag = FWD + templates[n % 4] + rc(REV)
        x, y = (rc(frag)[:250], frag[:250]) if n % 2 else (frag[:250], rc(frag)[:250])
        a.append((f'M1:1:FC:1:1101:{n}:9 1:N:0:ACGT+TTGA', x, 'I' * 250))
        b.append((f'M1:1:FC:1:1101:{n}:9 2:N:0:ACGT+TTGA', y, 'I' * 250))
    (d / 'in.1.fq').write_text(fq(a)); (d / 'in.2.fq').write_text(fq(b))
    (d / 'map.tsv').write_text(f'#SampleID\tfastqFile\tForwardPrimer\tReversePrimer\nS1\tin.1.fq,in.2.fq\t{FWD}\t{REV}\n')
    base = 'minSeqLength\t100\nRejectSeqWithoutFwdPrim\tT\nRejectSeqWithoutRevPrim\tT\nCheckForMixedPairs\tT\n'
    for fmt in (0, 1):
        o = opts(d, base + f'PEheaderPairFmt\t{fmt}\n', f'o{fmt}.txt')
        for mode, inp in (('nomap', ['-i', 'in.1.fq,in.2.fq']), ('map', ['-i_path', d, '-map', 'map.tsv'])):
            out = f'{mode}{fmt}'
            run(d, *inp, '-options', o, '-paired', 2, '-o_fastq', f'{out}.1.fq,{out}.2.fq', '-sample_sep', '___')
            h1 = [h for h, _, _ in read_fq(d / f'{out}.1.fq')]
            h2 = [h for h, _, _ in read_fq(d / f'{out}.2.fq')]
            assert len(h1) == len(h2) == 40, (out, len(h1), len(h2))
            for x, y in zip(h1, h2):
                if fmt == 1:
                    assert x.split()[0].endswith('/1') and y.split()[0].endswith('/2'), (out, x, y)
                    if mode == 'map':
                        assert x.startswith('S1___'), f'{out}: sample label lost: {x!r}'
                else:
                    assert ' 2:N:0:ACGT+TTGA' in y, (out, y)
                    if mode == 'nomap':
                        assert ' 1:N:0:ACGT+TTGA' in x, (out, x)


# ---------------------------------------------------------------- dereplication ground truth

@case
def paired_derep_ground_truth_and_hq_headers(d):
    r = rng(12)
    templates = [dna(r, 253) for _ in range(15)]
    maptext = '#SampleID\tfastqFile\tForwardPrimer\tReversePrimer\n'
    for s in range(3):
        a, b = [], []
        for n in range(700):
            t = r.randrange(15) if n % 10 else n % 15
            frag = FWD + templates[t] + rc(REV)
            x, y = list(frag[:250]), list(rc(frag)[:250])
            if n % 7 == 0:  # sparse single-base errors create rare variants
                y[150] = 'A' if y[150] != 'A' else 'C'
            a.append((f'p{s}_{n}_T{t}/1', ''.join(x), 'I' * 250))
            b.append((f'p{s}_{n}_T{t}/2', ''.join(y), 'I' * 250))
        (d / f'S{s}.1.fq').write_text(fq(a)); (d / f'S{s}.2.fq').write_text(fq(b))
        maptext += f'S{s}\tS{s}.1.fq,S{s}.2.fq\t{FWD}\t{REV}\n'
    (d / 'map.tsv').write_text(maptext)
    o = opts(d, 'minSeqLength\t150\nminAvgQuality\t20\nmaxPrimerErrs\t1\nRejectSeqWithoutFwdPrim\tT\n'
                'RejectSeqWithoutRevPrim\tT\n')
    results = {}
    for threads in (1, 4):
        out = d / f't{threads}'
        out.mkdir()
        run(out, '-i_path', d, '-map', d / 'map.tsv', '-options', o, '-paired', 2, '-o_fastq', 'o.1.fq,o.2.fq',
            '-o_dereplicate', 'derep.fna', '-merge_pairs_derep', 1, '-min_derep_copies', '8:1,4:2,3:3',
            '-threads', threads, '-sample_sep', '___')
        report = audit(out / 'derep.fna')
        assert report['map_total'] == 2100, report['map_total']
        merged = {s for _, s in read_fa(out / 'derep.merg.fna')}
        assert set(templates) <= merged, 'every template must pass as a merged unique'
        hq1, hq2 = read_fq(out / 'derep.1.hq.fq'), read_fq(out / 'derep.2.hq.fq')
        assert len(hq1) == len(hq2) > 0
        for (h1, _, _), (h2, _, _) in zip(hq1, hq2):
            base1 = re.sub(r';size=\d+;$', '', h1)
            assert h2 == base1[:-2] + '/2', f'reconstructed R2 header {h2!r} vs R1 {h1!r}'
        results[threads] = sorted((n.split(';size=')[1], s) for n, s in read_fa(out / 'derep.merg.fna'))
    assert results[1] == results[4], 'passing merged uniques differ between 1 and 4 threads'


@case
def multithreaded_derep_is_deterministic(d):
    # Failed reads join an exact key only if a passing copy was admitted before
    # them, and the first observation seeds the representative. Many small
    # blocks and failing first copies make both depend on block order, so
    # every thread count must reproduce the single-threaded outputs exactly.
    r = rng(24)
    templates = [dna(r, 150) for _ in range(40)]
    maptext = '#SampleID\tfastqFile\n'
    for s in range(2):
        recs, seen = [], set()
        for n in range(1500):
            t = r.randrange(len(templates))
            seq = list(templates[t])
            if r.random() < 0.2:
                seq[r.randrange(120)] = 'A' if seq[0] != 'A' else 'C'
            seq = ''.join(seq)
            fail = seq not in seen or r.random() < 0.3  # first copy of every key fails
            seen.add(seq)
            qual = ''.join(chr(33 + (r.randint(8, 12) if fail else r.randint(30, 40))) for _ in seq)
            recs.append((f'S{s}_{n}_T{t}_{"F" if fail else "P"}', seq, qual))
        (d / f'S{s}.fq').write_text(fq(recs))
        maptext += f'S{s}\tS{s}.fq\n'
    (d / 'map.tsv').write_text(maptext)
    o = opts(d, 'minSeqLength\t100\nmaxSeqLength\t1000\nminAvgQuality\t25\nmaxAmbiguousNT\t0\n'
                'TruncateSequenceLength\t120\n')
    outputs = {}
    for label, threads in (('t1', 1), ('t4', 4), ('t8', 8), ('t8b', 8)):
        out = d / label
        out.mkdir()
        run(out, '-i_path', d, '-map', d / 'map.tsv', '-options', o, '-o_fastq', 'o.fq', '-o_dereplicate', 'derep.fna',
            '-min_derep_copies', 1, '-threads', threads, '-iniBlockSize', 2, '-sample_sep', '___')
        outputs[label] = {name: (out / name).read_bytes() for name in ('derep.fna', 'derep.fna.rest', 'derep.map', 'derep.1.hq.fq')}
    assert outputs['t1']['derep.map'].count(b'\n') > 40, 'fixture must produce dereplicated keys'
    for label in ('t4', 't8', 't8b'):
        for name, data in outputs['t1'].items():
            assert outputs[label][name] == data, f'{name} differs between 1 and {label} threads'


@case
def documented_ont_example(d):
    # The runnable example in docs/ont.md, with the outputs that page promises.
    fx = Path(__file__).resolve().parent / 'fixtures' / 'ont'
    (d / 'samples').mkdir()
    run(d, '-i_path', fx, '-map', fx / 'map.tsv', '-options', fx / 'options.txt', '-paired', 1, '-ontMode', 1,
        '-barcodeSearchWindow', 200, '-ontBarcodeEnds', 'either', '-threads', 4, '-i_qual_offset', 33,
        '-sample_sep', '___', '-o_fastq', 'full.fq', '-o_demultiplex', f'{d / "samples"}/', '-DemultiBPperSR', 0,
        '-o_dereplicate', 'derep.fas', '-derep_format', 'fa', '-min_derep_copies', 1)
    assert len(read_fq(d / 'full.fq')) == 6
    assert [len(read_fq(d / 'samples' / f'S{s}.fq')) for s in (0, 1)] == [3, 3]
    assert len(read_fa(d / 'derep.fas')) == 2 and len(read_fq(d / 'derep.1.hq.fq')) == 2
    counts = sorted(line.split('\t', 1)[1] for line in (d / 'derep.map').read_text().splitlines()[1:])
    assert counts == ['0:3\t1:2', '1:1'], counts


def main():
    selected = [c for c in CASES if not ONLY or c.__name__ in ONLY]
    if ONLY and len(selected) != len(ONLY):
        print('Unknown case(s):', ', '.join(ONLY - {c.__name__ for c in CASES}))
        return 2
    root = Path(tempfile.mkdtemp(prefix='sdm-cli-regressions-'))
    failed = []
    for fn in selected:
        d = root / fn.__name__
        d.mkdir()
        try:
            fn(d)
            print(f'PASS {fn.__name__}', flush=True)
            shutil.rmtree(d, ignore_errors=True)
        except Exception as e:  # noqa: BLE001 - report every failure, keep going
            failed.append(fn.__name__)
            detail = ''.join(traceback.format_exception_only(type(e), e)).strip()
            print(f'FAIL {fn.__name__}: {detail[:1500]}\n     fixture: {d}', flush=True)
    if failed:
        print(f'{len(failed)} of {len(selected)} CLI regression case(s) failed: {", ".join(failed)}')
        return 1
    shutil.rmtree(root, ignore_errors=True)
    print(f'All {len(selected)} CLI regression cases passed')
    return 0


if __name__ == '__main__':
    sys.exit(main())
