#!/usr/bin/env python3
"""Deterministic single-end ONT integration tests; pass the SDM executable path.

Synthetic adapter/barcode/spacer/primer layouts, not a kit certification.
Uses the existing pipeline test reader/runner; fixtures are retained on failure.
"""
from collections import Counter
import gzip
from pathlib import Path
import random
import re
import shutil
import subprocess
import tempfile
import time

from test_thorough_audit import SDM, records, run, RC

RNG = random.Random(351)
random_dna = lambda n: ''.join(RNG.choices('ACGT', k=n))
BARS = ['ACGTCAGTGCTAGACG', 'TGCATCGACAGTTCGA']
FWD, REV = 'GATCTGACCGTACGTA', 'CTAGGCTACGACTGAC'
BODY = random_dna(420)
SECOND = random_dna(420)


def write_fastq(path, reads):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt') as out:
        for i, (seq, qual) in enumerate(reads):
            out.write(f'@r{i}\n{seq}\n+\n{qual}\n')


def with_quality(seq):
    return seq, ''.join(chr(33 + 25 + i % 15) for i in range(len(seq)))


def reverse_read(read):
    return RC(read[0]), read[1][::-1]


def layout(barcode, body=BODY, rear=None, adapter=37, spacer=45, primers=True):
    left = random_dna(adapter) + barcode + random_dna(spacer)
    right = random_dna(31) + RC(rear) + random_dna(23) if rear else ''
    seq = left + (FWD if primers else '') + body + (RC(REV) if primers else '') + right
    read = with_quality(seq)
    begin = len(left) + (len(FWD) if primers else 0)
    return read, read[1][begin:begin + len(body)]


def options(d, errors=1, primers=True):
    settings = dict(minSeqLength=2, maxSeqLength=20000, minAvgQuality=0,
                    QualWindowThreshhold=0, TrimWindowThreshhold=0,
                    maxAccumulatedError=-1, BinErrorModelAlpha=-1,
                    maxAmbiguousNT=0, maxHomonucleotide=0, trimHomonucleotide=0,
                    TruncateSequenceLength=-1, maxBarcodeErrs=errors, maxPrimerErrs=0,
                    keepBarcodeSeq=0, keepPrimerSeq=0,
                    RejectSeqWithoutFwdPrim='T' if primers else 'F',
                    RejectSeqWithoutRevPrim='T' if primers else 'F')
    (d / 'options.txt').write_text(''.join(f'{k}\t{v}\n' for k, v in settings.items()))


def mapping(d, barcodes=BARS, primers=True, files=None, groups=False):
    columns = ['#SampleID'] + (['BarcodeSequence'] if barcodes is not None else [])
    if primers:
        columns += ['LinkerPrimerSequence', 'ReversePrimer']
    columns += ['fastqFile']
    if groups:
        columns += ['SequencingRun', 'CombineSamples']
    rows = []
    for i in range(len(barcodes) if barcodes is not None else len(files)):
        row = [f'S{i}'] + ([barcodes[i]] if barcodes is not None else [])
        if primers:
            row += [FWD, REV]
        row += [files[i] if files else 'in.fq.gz']
        if groups:
            row += [f'run{i}', f'group{i}']
        rows.append('\t'.join(row))
    (d / 'map.tsv').write_text('\t'.join(columns) + '\n' + '\n'.join(rows) + '\n')


def args(d, extra=()):
    return ['-ontMode', '1', '-sample_sep', '___', '-options', d / 'options.txt', *extra]


def counters(d):
    text = (d / 'filter.log').read_text()
    return {name: int(re.search(r'ONT ' + name + r':\s+([\d,]+)', text).group(1).replace(',', ''))
            for name in ['assigned', 'unassigned', 'ambiguous', 'conflicting ends']}


def matching(d, threads):
    expected = []
    reads = []
    # Exact, internal substitution, insertion and deletion, both read directions.
    bc = BARS[0]
    variants = [bc, bc[:6] + ('A' if bc[6] != 'A' else 'C') + bc[7:],
                bc[:7] + 'T' + bc[7:], bc[:7] + bc[8:]]
    for i, variant in enumerate(variants):
        for reverse in [False, True]:
            read, qual = layout(variant, rear=bc, adapter=17 + i * 11, spacer=32 + i * 9)
            reads.append(reverse_read(read) if reverse else read)
            expected.append(('S0', BODY, qual))
    # Rear-only barcode must also identify a read in either orientation.
    for reverse in [False, True]:
        read, qual = layout('', rear=BARS[1], adapter=0, spacer=18)
        reads.append(reverse_read(read) if reverse else read)
        expected.append(('S1', BODY, qual))
    write_fastq(d / 'in.fq.gz', reads)
    mapping(d); options(d)
    run(d, args(d), threads=threads)
    output = records(d / 'out.fq')
    observed = sorted((name.split('___')[0][1:], seq, ''.join(chr(q + 33) for q in qual))
                      for name, seq, qual in output)
    assert observed == sorted(expected), (observed, expected)
    assert counters(d) == dict(assigned=10, unassigned=0, ambiguous=0, **{'conflicting ends': 0}), counters(d)


def rejection(d, ends):
    bc = BARS[0]
    reads = [layout(bc, rear=bc)[0], layout(bc)[0],
             layout(bc, rear=BARS[1])[0],
             layout('G' * len(bc))[0], layout('G' * len(bc), rear=bc)[0]]
    write_fastq(d / 'in.fq.gz', reads)
    mapping(d); options(d)
    run(d, args(d, ['-ontBarcodeEnds', ends]))
    assigned = 3 if ends == 'either' else 1
    assert len(records(d / 'out.fq')) == assigned
    assert counters(d) == dict(assigned=assigned, unassigned=4-assigned, ambiguous=0,
                               **{'conflicting ends': 1}), counters(d)


def ties(d, exact=False):
    bc = BARS[0]
    if exact:
        # Two different exact sample barcodes in one end window.
        bars, read = BARS, layout(BARS[0] + 'TT' + BARS[1])[0]
    else:
        bars = [bc, bc[:6] + 'C' + bc[7:]]
        read = layout(bc[:6] + 'A' + bc[7:])[0]
    write_fastq(d / 'in.fq.gz', [read, reverse_read(read)])
    mapping(d, bars); options(d)
    run(d, args(d))
    assert records(d / 'out.fq') == []
    assert counters(d)['ambiguous'] == 2, counters(d)


def boundaries(d, offset, window, accepted, short=False):
    body = 'TTAAG' if short else BODY
    seq = 'G' * offset + BARS[0] + body
    write_fastq(d / 'in.fq.gz', [with_quality(seq), reverse_read(with_quality(seq))])
    mapping(d, primers=False); options(d, errors=0, primers=False)
    run(d, args(d, ['-barcodeSearchWindow', str(window)]))
    output = records(d / 'out.fq')
    assert sorted(s for _, s, _ in output) == (sorted([body, RC(body)]) if accepted else []), output
    assert counters(d)['assigned'] == (2 if accepted else 0), counters(d)


def already_trimmed(d, threads, per_file=False, primer_bearing=False):
    if per_file:
        reads = [[with_quality(BODY)] * 3, [with_quality(BODY)] * 2 + [with_quality(SECOND)]]
        files = ['first.fq.gz', 'second.fq']
        for file, rs in zip(files, reads):
            write_fastq(d / file, rs)
        mapping(d, None, primers=False, files=files, groups=True)
    else:
        rs = [with_quality(BARS[0] + BODY)] * 3 + [with_quality(BARS[1] + BODY)] * 2 + [with_quality(BARS[1] + SECOND)]
        if primer_bearing:
            rs = [layout(BARS[0], rear=BARS[0])[0]] * 3 + [reverse_read(layout(BARS[1], rear=BARS[1])[0])] * 2 + [layout(BARS[1], SECOND)[0]]
        write_fastq(d / 'in.fq.gz', rs)
        mapping(d, primers=primer_bearing)
    options(d, errors=0, primers=primer_bearing)
    extra = ['-o_dereplicate', d / 'derep.fas', '-derep_format', 'fa', '-min_derep_copies', '1',
             '-o_demultiplex', str(d) + '/', '-DemultiBPperSR', '0']
    run(d, args(d, extra), threads=threads)
    assert Counter(s for _, s, _ in records(d / 'out.fq')) == Counter({BODY: 5, SECOND: 1})
    prefix = 'S'  # Existing per-sample filenames; headers/map carry CombineSamples.
    for sample, count in [(0, 3), (1, 3)]:
        assert len(records(d / f'{prefix}{sample}.fq')) == count
    derep = (d / 'derep.fas').read_text().splitlines()
    assert len(derep) == 4 and set(derep[1::2]) == {BODY, SECOND}, derep
    hq = records(d / 'derep.1.hq.fq')
    assert len(hq) == 2 and {s for _, s, _ in hq} == {BODY, SECOND}, hq
    table = (d / 'derep.map').read_text().splitlines()
    rows = {line.split('\t')[0].split(';')[0]: line.split('\t')[1:] for line in table[1:]}
    for i in range(0, len(derep), 2):
        key = derep[i][1:].split(';')[0]
        expect = ['0:3', '1:2'] if derep[i+1] == BODY else ['1:1']
        assert sorted(rows[key]) == expect, (rows, expect)
    assert counters(d)['assigned'] == 6, counters(d)
    # Pass the emitted map and HQ representatives through the existing matrix path.
    seed_ids = [name[1:] for name, _, _ in hq]
    (d / 'assign.uc').write_text(''.join(f'{name}\totu{i+1}\t*\n' for i, name in enumerate(seed_ids)))
    (d / 'fallback.fa').write_text(''.join(f'>{name}\n{seq}\n' for (name, seq, _) in hq).replace('>@', '>'))
    command = [SDM, '-i_fastq', str(d / 'derep.1.hq.fq'), '-paired', '1', '-i_qual_offset', '33',
               '-o_fna', str(d / 'seeds.fa'), '-uparseVer', '12', '-optimalRead2Cluster', str(d / 'assign.uc'),
               '-derep_map', str(d / 'derep.map'), '-OTU_fallback', str(d / 'fallback.fa'),
               '-otu_matrix', str(d / 'matrix.tsv'), '-log', str(d / 'matrix.log')]
    result = subprocess.run(command, capture_output=True, text=True, timeout=60)
    (d / 'matrix.console').write_text(result.stdout + result.stderr)
    assert result.returncode == 0, result.stdout + result.stderr
    matrix = (d / 'matrix.tsv').read_text().splitlines()
    assert sorted(line.split('\t')[1:] for line in matrix[1:]) == [['0', '1'], ['3', '2']], matrix


def per_file_primers(d):
    reads = [layout('', adapter=23, spacer=55)[0], reverse_read(layout('', adapter=31, spacer=52)[0])]
    write_fastq(d / 'sample.fq', reads)
    mapping(d, None, files=['sample.fq']); options(d)
    run(d, args(d, ['-ontBarcodeEnds', 'both']))  # No barcode requirements in per-file mode.
    assert [s for _, s, _ in records(d / 'out.fq')] == [BODY, BODY]
    assert counters(d)['assigned'] == 2


def runtime(d):
    read, _ = layout(BARS[0][:7] + 'T' + BARS[0][7:], rear=BARS[0])
    write_fastq(d / 'in.fq.gz', [read, reverse_read(read)] * 2000)
    mapping(d, BARS + [random_dna(24) for _ in range(22)]); options(d)
    start = time.monotonic(); run(d, args(d), threads=4)
    elapsed = time.monotonic() - start
    assert len(records(d / 'out.fq')) == 4000
    print(f'ONT runtime sanity: 4,000 reads with barcode indels, 24-barcode panel, {elapsed:.2f}s', flush=True)



def primer_controls(d):
    read, _ = layout(BARS[0], rear=BARS[0])
    write_fastq(d / 'in.fq.gz', [read, reverse_read(read)])
    mapping(d); options(d)
    path = d / 'options.txt'
    path.write_text(path.read_text().replace('keepPrimerSeq\t0', 'keepPrimerSeq\t1'))
    run(d, args(d))
    assert [seq for _, seq, _ in records(d / 'out.fq')] == [FWD + BODY + RC(REV)] * 2
    # A required but absent primer is an existing primer rejection, after assignment.
    write_fastq(d / 'in.fq.gz', [with_quality(BARS[0] + BODY)])
    options(d)
    run(d, args(d))
    assert not records(d / 'out.fq') and counters(d)['assigned'] == 1
    assert re.search(r'Forward primer detected[^\n]*0 ', (d / 'filter.log').read_text())
    # Reverse-only primer columns use the same trim/orientation machinery.
    read, _ = layout(BARS[0])
    read = with_quality(read[0].replace(FWD, ''))
    write_fastq(d / 'in.fq.gz', [read, reverse_read(read)])
    (d / 'map.tsv').write_text('#SampleID\tBarcodeSequence\tReversePrimer\tfastqFile\n'
                               f'S0\t{BARS[0]}\t{REV}\tin.fq.gz\n')
    run(d, args(d))
    # Without a forward primer the unspecified front spacer remains.
    output = records(d / 'out.fq')
    assert len(output) == 2 and output[0][1] == output[1][1] and output[0][1].endswith(BODY)


def extra_edges(d):
    # Two internal substitutions are beyond budget 1 and accepted at budget 2.
    bc = BARS[0]
    changed = bc[:4] + 'A' + bc[5:10] + 'A' + bc[11:]
    read = layout(changed)[0]
    write_fastq(d / 'in.fq.gz', [read, reverse_read(read)])
    mapping(d); options(d)
    run(d, args(d))
    assert counters(d)['unassigned'] == 2
    options(d, errors=2)
    run(d, args(d))
    assert len(records(d / 'out.fq')) == 2
    # Missing ends and reads shorter than a barcode cannot satisfy both.
    write_fastq(d / 'in.fq.gz', [with_quality('ACGT'), with_quality(BARS[0] + 'TTAAG')])
    mapping(d, primers=False); options(d, errors=0, primers=False)
    run(d, args(d, ['-ontBarcodeEnds', 'both']))
    assert counters(d)['unassigned'] == 2
    # A palindromic barcode seen from both ends is only one physical occurrence.
    write_fastq(d / 'in.fq.gz', [with_quality('ACGTACGTTTAAG')])
    mapping(d, ['ACGTACGT'], primers=False)
    run(d, args(d, ['-ontBarcodeEnds', 'both']))
    assert counters(d)['ambiguous'] == 1 and not records(d / 'out.fq')
    # Inconsistent primer strand evidence is rejected without choosing a rescue.
    write_fastq(d / 'in.fq.gz', [with_quality(BARS[0] + FWD + BODY + RC(FWD))])
    mapping(d); options(d)
    run(d, args(d))
    assert counters(d)['ambiguous'] == 1
    # Mixed primer cells: an empty first sample must not disable another sample.
    write_fastq(d / 'in.fq.gz', [layout(BARS[1])[0]])
    path = d / 'map.tsv'
    path.write_text(path.read_text().replace(FWD + '\t' + REV, '\t', 1))
    run(d, args(d))
    assert [seq for _, seq, _ in records(d / 'out.fq')] == [BODY]



def disabled_mode(d):
    bc = BARS[0]
    substitution = bc[:6] + 'A' + bc[7:]
    insertion = bc[:7] + 'T' + bc[7:]
    reads = [with_quality(x + BODY) for x in [bc, substitution, insertion, 'G' * 60 + bc]]
    write_fastq(d / 'in.fq.gz', reads)
    mapping(d, primers=False); options(d, primers=False)
    run(d, ['-options', d / 'options.txt'])
    default = records(d / 'out.fq')
    run(d, ['-options', d / 'options.txt', '-ontMode', '0'])
    assert records(d / 'out.fq') == default and len(default) == 2
    assert all(seq == BODY for _, seq, _ in default)
    assert 'ONT assigned' not in (d / 'filter.log').read_text()
    run(d, args(d))
    assert len(records(d / 'out.fq')) == 4


def invalid(d):
    write_fastq(d / 'in.fq.gz', [with_quality(BARS[0] + BODY)])
    mapping(d, primers=False); options(d, primers=False)
    cases = [(['-ontMode', 'bad'], '-ontMode'), (['-ontBarcodeEnds', 'bad'], '-ontBarcodeEnds'),
             (['-paired', '2'], 'single-end'), (['-onlyPair', '1'], 'single-end'),
             (['-i_MID_fastq', 'missing.fq'], 'single-end'), (['-GoldenAxe', '1'], 'single-end')]
    cases += [(['-barcodeSearchWindow', x], '-barcodeSearchWindow') for x in ['0', '-1', '1.5', 'x', '10001', '9999999999999999999999999']]
    for extra, message in cases:
        console = run(d, args(d, extra), valid=False)
        assert message in console, console
    for value in ['-1', '4', '1x']:
        options(d, errors=value, primers=False)
        assert 'maxBarcodeErrs' in run(d, args(d), valid=False)
    options(d, primers=False)
    for column in ['Barcode2ndPair', 'MIDfqFile', 'alignmentFile', 'SampleIDinHead']:
        mapping(d, primers=False)
        path = d / 'map.tsv'; lines = path.read_text().splitlines()
        path.write_text(lines[0] + '\t' + column + '\n' + '\n'.join(s + '\t' for s in lines[1:]) + '\n')
        assert 'ONT' in run(d, args(d), valid=False)


def main():
    root = Path(tempfile.mkdtemp(prefix='sdm-ont-pipeline-'))
    cases = [(f'matching_{t}', matching, (t,)) for t in [1, 4]]
    cases += [(f'rejection_{e}', rejection, (e,)) for e in ['either', 'both']]
    cases += [('ties', ties, ()), ('exact_ties', ties, (True,)),
              ('last_start', boundaries, (19, 20, True)), ('outside', boundaries, (20, 20, False)),
              ('one_base_window', boundaries, (0, 1, True)), ('short', boundaries, (0, 200, True, True)),
              ('per_file_primers', per_file_primers, ())]
    cases += [(f'counts_{t}_{f}', already_trimmed, (t, f)) for t in [1, 4] for f in [False, True]]
    cases += [(f'primer_counts_{t}', already_trimmed, (t, False, True)) for t in [1, 4]]
    cases += [('primer_controls', primer_controls, ()), ('extra_edges', extra_edges, ()),
              ('disabled_mode', disabled_mode, ()), ('invalid', invalid, ()), ('runtime', runtime, ())]
    try:
        for name, test, values in cases:
            d = root / name; d.mkdir(); test(d, *values)
            print('PASS', name, flush=True)
    except BaseException:
        print(f'Fixtures retained at {root}', flush=True)
        raise
    shutil.rmtree(root)
    print('ONT pipeline regressions passed')


if __name__ == '__main__':
    main()
