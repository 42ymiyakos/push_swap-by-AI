#!/usr/bin/env python3
"""Deterministic push_swap input generator.

Usage: gen.py --n N --dist DIST --seed S [--count C]  -> prints C lines (one case per line)
Distributions: uniform, sorted, reverse, nearly (k random adjacent swaps), rotated,
  rotated_reverse, two_runs (interleaved sorted runs), blocks_reverse, sawtooth, adversarial_turk
Values: a random set of distinct ints in [-2^31, 2^31-1] mapped by rank (rank pattern is what matters),
  unless --ranks is given (values 0..N-1).
"""
import argparse, random

def perm(n, dist, rng):
    p = list(range(n))
    if dist == 'uniform':
        rng.shuffle(p)
    elif dist == 'sorted':
        pass
    elif dist == 'reverse':
        p.reverse()
    elif dist == 'nearly':
        for _ in range(max(1, n // 20)):
            i = rng.randrange(n - 1) if n > 1 else 0
            if n > 1: p[i], p[i + 1] = p[i + 1], p[i]
    elif dist == 'rotated':
        k = rng.randrange(1, n) if n > 1 else 0
        p = p[k:] + p[:k]
    elif dist == 'rotated_reverse':
        p.reverse(); k = rng.randrange(1, n) if n > 1 else 0
        p = p[k:] + p[:k]
    elif dist == 'two_runs':
        lo, hi = p[: n // 2], p[n // 2:]
        out = []
        while lo or hi:
            src = lo if (lo and (not hi or rng.random() < 0.5)) else hi
            out.append(src.pop(0))
        p = out
    elif dist == 'blocks_reverse':
        b = max(1, int(n ** 0.5)); blocks = [p[i:i + b] for i in range(0, n, b)]
        blocks.reverse(); p = [x for bl in blocks for x in bl]
    elif dist == 'sawtooth':
        b = max(2, n // 10); p = []
        for s in range(0, n, b): p += list(reversed(range(s, min(n, s + b))))
    elif dist == 'local_shuffle':
        w = 5
        for s in range(0, n, w):
            seg = p[s:s + w]; rng.shuffle(seg); p[s:s + w] = seg
    else:
        raise SystemExit('unknown dist ' + dist)
    return p

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, required=True)
    ap.add_argument('--dist', default='uniform')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--count', type=int, default=1)
    ap.add_argument('--ranks', action='store_true')
    a = ap.parse_args()
    for c in range(a.count):
        rng = random.Random(f"{a.n}:{a.dist}:{a.seed + c}")
        p = perm(a.n, a.dist, rng)
        if a.ranks:
            vals = p
        else:
            pool = sorted(rng.sample(range(-2**31, 2**31), a.n))
            vals = [pool[x] for x in p]
        print(' '.join(map(str, vals)))

if __name__ == '__main__':
    main()
