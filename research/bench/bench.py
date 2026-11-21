#!/usr/bin/env python3
"""Benchmark harness.

bench.py --solver "CMD {args}" --set data/sets/X.txt --tag NAME [--jobs 4] [--timeout 600] [--limit K]
  The solver command is a shell-free template; '{args}' is replaced by the numbers (as separate argv).
  For every case: wall time, peak RSS (KB, via wait4), op count, checker verdict are recorded.
  Results -> bench/runs/NAME.csv (+ summary appended to bench/runs/summary.csv).
"""
import argparse, csv, os, shlex, subprocess, sys, time, statistics, json, tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CHECKER = os.path.join(ROOT, 'tools', 'checker')

def run_case(solver_tmpl, nums, timeout):
    argv = []
    for tok in shlex.split(solver_tmpl):
        if tok == '{args}': argv += nums
        else: argv.append(tok)
    with tempfile.TemporaryFile() as out:
        t0 = time.time()
        p = subprocess.Popen(argv, stdout=out, stderr=subprocess.DEVNULL)
        try:
            # wait4 with timeout: poll
            while True:
                pid, status, ru = os.wait4(p.pid, os.WNOHANG)
                if pid: break
                if time.time() - t0 > timeout:
                    p.kill(); pid, status, ru = os.wait4(p.pid, 0); return dict(ok=False, ops=-1, time=timeout, rss=ru.ru_maxrss, verdict='TIMEOUT')
                time.sleep(0.01)
        finally:
            p.returncode = 0
        el = time.time() - t0
        out.seek(0)
        ops = out.read()
    chk = subprocess.run([CHECKER] + nums, input=ops, capture_output=True)
    v = chk.stdout.decode().strip()
    ok = v.startswith('OK')
    cnt = int(v.split()[1]) if ok else -1
    return dict(ok=ok, ops=cnt, time=el, rss=ru.ru_maxrss, verdict=v)

def pct(xs, q):
    xs = sorted(xs); k = (len(xs) - 1) * q; f = int(k); c = min(f + 1, len(xs) - 1)
    return xs[f] + (xs[c] - xs[f]) * (k - f)

def summarize(rows):
    ops = [r['ops'] for r in rows if r['ok']]
    t = [r['time'] for r in rows]
    if not ops: return dict(n_ok=0, n=len(rows))
    return dict(n=len(rows), n_ok=len(ops), min=min(ops), max=max(ops), mean=round(statistics.mean(ops), 2),
                median=statistics.median(ops), std=round(statistics.pstdev(ops), 2), p90=round(pct(ops, .9), 1),
                p95=round(pct(ops, .95), 1), p99=round(pct(ops, .99), 1), time_mean=round(statistics.mean(t), 3),
                time_max=round(max(t), 3), rss_max_kb=max(r['rss'] for r in rows))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--solver', required=True); ap.add_argument('--set', required=True)
    ap.add_argument('--tag', required=True); ap.add_argument('--jobs', type=int, default=4)
    ap.add_argument('--timeout', type=float, default=600); ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--offset', type=int, default=0)
    ap.add_argument('--note', default='')
    a = ap.parse_args()
    cases = [l.split() for l in open(a.set) if l.strip()]
    idx = list(range(len(cases)))[a.offset:]
    if a.limit: idx = idx[:a.limit]
    os.makedirs(os.path.join(HERE, 'runs'), exist_ok=True)
    rows = [None] * len(idx)
    t0 = time.time()
    with ThreadPoolExecutor(a.jobs) as ex:
        futs = {ex.submit(run_case, a.solver, cases[i], a.timeout): k for k, i in enumerate(idx)}
        for f in futs:
            k = futs[f]; r = f.result(); r['case'] = idx[k]; r['n'] = len(cases[idx[k]]); rows[k] = r
    path = os.path.join(HERE, 'runs', a.tag + '.csv')
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['case', 'n', 'ok', 'ops', 'time', 'rss', 'verdict'])
        w.writeheader(); [w.writerow({k: r[k] for k in w.fieldnames}) for r in rows]
    s = summarize(rows)
    s.update(tag=a.tag, set=os.path.basename(a.set), solver=a.solver, jobs=a.jobs, wall=round(time.time() - t0, 1),
             note=a.note, date=time.strftime('%Y-%m-%d %H:%M:%S'))
    sp = os.path.join(HERE, 'runs', 'summary.jsonl')
    with open(sp, 'a') as fh: fh.write(json.dumps(s) + '\n')
    bad = [r for r in rows if not r['ok']]
    print(json.dumps(s))
    if bad: print('FAILURES:', [(r['case'], r['verdict']) for r in bad[:10]], file=sys.stderr); sys.exit(1)

if __name__ == '__main__':
    main()
