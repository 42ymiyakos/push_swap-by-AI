#!/usr/bin/env python3
"""usage: paired.py BASE_TAG TAG [TAG...]  -- paired sign-flip test on bench/runs/<tag>.csv (ops column)"""
import csv, random, statistics as st, sys, os
os.chdir(os.path.join(os.path.dirname(__file__), '..'))
def load(t): return [int(x['ops']) for x in csv.DictReader(open(f'bench/runs/{t}.csv'))]
b = load(sys.argv[1])
print(f"{sys.argv[1]}: mean {st.mean(b):.2f} med {st.median(b)} sd {st.pstdev(b):.1f} max {max(b)}")
for t in sys.argv[2:]:
    x = load(t); d = [p - q for p, q in zip(x, b)]; obs = st.mean(d); rng = random.Random(1)
    cnt = sum(abs(st.mean([v if rng.random() < .5 else -v for v in d])) >= abs(obs) for _ in range(20000))
    print(f"{t}: mean {st.mean(x):.2f} med {st.median(x)} sd {st.pstdev(x):.1f} max {max(x)} | diff {obs:+.2f} "
          f"W {sum(v < 0 for v in d)} T {sum(v == 0 for v in d)} L {sum(v > 0 for v in d)} p {cnt / 20000:.4f}")
