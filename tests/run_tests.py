#!/usr/bin/env python3
"""Run the sdm test suites, continue past failures, and print a summary.

Usage (from the repository root; `make test` does this for you):
  python3 tests/run_tests.py --sdm ./sdm --bin-dir .build/portable-1 [--hts] [--quick] [suite ...]
  python3 tests/run_tests.py --list

Each suite runs in its own process with a timeout. Full output is written to
<log-dir>/<suite>.log; failing suites also print their last lines. Exit status is
0 when every selected suite passed, 1 otherwise.
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys
import time

REPO = Path(__file__).resolve().parent.parent
# Scratch files some C++ tests write below tests/; removed after each run.
GENERATED = ['tests/tmp4.txt', 'tests/tmp_crlf.txt', 'tests/tmp_gz.txt.gz', 'tests/tmp_nonl.txt',
             'tests/tmp_small.txt', 'tests/tmp_small2.txt', 'TestInputHelpers_tmp_crlf_split.txt']
GENERATED_GLOBS = ['tests/reg_fasta_tier.*']


def suites(args):
    """(name, command, slow) for every suite; C++ binaries come from --bin-dir."""
    b, sdm = Path(args.bin_dir), args.sdm
    cpp = ['TestAuditFixes', 'TestInputHelpers', 'TestIfBufStream', 'TestGzipInput', 'TestBgzfOutput',
           'TestRegression', 'TestSearch', 'TestSeedOutputRouting', 'TestCoarseDerep', 'TestPartitionWorkPool']
    result = [(name, [str(b / name)], False) for name in cpp]

    def py(script, slow=False):
        return script.removesuffix('.py'), [sys.executable, f'tests/{script}', sdm], slow

    result += [
        py('test_regressions_cli.py'),
        py('test_audit_pipeline.py'),
        py('test_thorough_audit.py'),
        py('test_merged_pipeline.py'),
        py('test_prefix_derep.py'),
        py('test_ont_pipeline.py'),
        py('test_diff_quality_pipeline.py'),
        py('test_paired_diff_quality_pipeline.py'),
        py('test_seed_extension.py'),
        py('test_seed_subclusters.py', slow=True),
        py('test_storage_only_parity.py', slow=True),
        py('test_derep_behavior_parity.py', slow=True),
        ('test_derep_behavior_parity-fq', [sys.executable, 'tests/test_derep_behavior_parity.py', sdm, '--main-format',
                                           'fq', '--complete-mates', '--consensus-errors', '--require-parity'], True),
    ]
    if args.hts:
        result += [('TestHtsInput', [str(b / 'TestHtsInput')], False),
                   ('test_hts_pipeline', [sys.executable, 'tests/test_hts_pipeline.py', sdm, str(b / 'TestHtsInput')], False)]
    return result


def cleanup():
    for rel in GENERATED:
        (REPO / rel).unlink(missing_ok=True)
    for pattern in GENERATED_GLOBS:
        for path in REPO.glob(pattern):
            path.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('names', nargs='*', help='run only these suites (see --list)')
    parser.add_argument('--sdm', default='./sdm', help='sdm executable under test (default ./sdm)')
    parser.add_argument('--bin-dir', default='.build/portable-1', help='directory with the C++ test binaries')
    parser.add_argument('--hts', action='store_true', help='include HTSlib suites (WITH_HTS=1 builds)')
    parser.add_argument('--quick', action='store_true', help='skip the slow parity/seed-extension suites')
    parser.add_argument('--timeout', type=int, default=1800, help='per-suite timeout in seconds')
    parser.add_argument('--log-dir', default=None, help='where to write per-suite logs (default <bin-dir>/test-logs)')
    parser.add_argument('--list', action='store_true', help='list suites and exit')
    args = parser.parse_args()
    args.sdm = str((REPO / args.sdm).resolve()) if not Path(args.sdm).is_absolute() else args.sdm

    selected = suites(args)
    if args.list:
        for name, command, slow in selected:
            print(f'{name:40} {"(slow)" if slow else "      "} {" ".join(command[1:] if command[0] == sys.executable else command)}')
        return 0
    if args.names:
        unknown = set(args.names) - {name for name, _, _ in selected}
        if unknown:
            print('Unknown suite(s): ' + ', '.join(sorted(unknown)) + ' (see --list)')
            return 2
        selected = [s for s in selected if s[0] in args.names]
    elif args.quick:
        selected = [s for s in selected if not s[2]]

    log_dir = Path(args.log_dir) if args.log_dir else REPO / args.bin_dir / 'test-logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    print(f'Running {len(selected)} test suite(s) against {args.sdm}; logs in {log_dir}', flush=True)
    results = []
    for name, command, _ in selected:
        started = time.time()
        log = log_dir / f'{name}.log'
        if command[0] != sys.executable and not Path(command[0]).exists():
            status, detail = 'FAIL', f'test binary missing: {command[0]} (build it with make)'
            log.write_text(detail + '\n')
        else:
            try:
                p = subprocess.run(command, cwd=REPO, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True, errors='replace', timeout=args.timeout,
                                   env={**os.environ, 'PYTHONDONTWRITEBYTECODE': '1'})  # no tests/__pycache__
                log.write_text(p.stdout)
                status = 'PASS' if p.returncode == 0 else 'FAIL'
                detail = '' if p.returncode == 0 else f'exit code {p.returncode}'
            except subprocess.TimeoutExpired as e:
                output = e.stdout if isinstance(e.stdout, str) else (e.stdout or b'').decode(errors='replace')
                log.write_text(output)
                status, detail = 'FAIL', f'timed out after {args.timeout}s'
        seconds = time.time() - started
        results.append((name, status, seconds, detail))
        print(f'{status} {name:40} {seconds:7.1f}s  {detail}', flush=True)
        if status == 'FAIL':
            tail = log.read_text(errors='replace').splitlines()[-25:]
            print('\n'.join('     | ' + line for line in tail), flush=True)
    cleanup()

    failed = [r for r in results if r[1] == 'FAIL']
    total = sum(r[2] for r in results)
    print(f'\n{len(results) - len(failed)} passed, {len(failed)} failed ({total:.0f}s)')
    for name, _, _, detail in failed:
        print(f'  FAILED: {name} ({detail}); log: {log_dir / (name + ".log")}')
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
