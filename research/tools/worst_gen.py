#!/usr/bin/env python3
"""Structured 'hard' permutations (re-implemented from the description in keusidan/42-push_swap-solver
docs/04_worst_case_analysis.md): coprime stride, bit reversal, alternating bit reversal, golden stride,
Erdos-Szekeres blocks.  usage: worst_gen.py N MODE"""
import math, sys

def coprime(n, s=0):
    if s <= 0:
        if n == 500: s = 49
        else:
            t = max(2, round(2.2 * math.sqrt(n)))
            for d in range(n):
                if t + d < n and math.gcd(t + d, n) == 1: s = t + d; break
                if t - d > 1 and math.gcd(t - d, n) == 1: s = t - d; break
    return [(i * s) % n for i in range(n)]

def bitrev(n):
    bits = 0
    while (1 << bits) < n: bits += 1
    out = []
    for i in range(1 << bits):
        r = int(format(i, f'0{bits}b')[::-1], 2) if bits else 0
        if r < n: out.append(r)
    return out

def altbitrev(n):
    br = bitrev(n); l, r = 0, n - 1; out = []
    for i in range(n):
        if i % 2 == 0: out.append(br[l]); l += 1
        else: out.append(br[r]); r -= 1
    return out

def golden(n):
    s = max(1, round(n * (math.sqrt(5) - 1) / 2))
    while math.gcd(s, n) != 1: s += 1
    return coprime(n, s)

def es(n):
    b = math.ceil(math.sqrt(n)); out = []
    for lo in range(0, n, b): out.extend(range(min(n, lo + b) - 1, lo - 1, -1))
    return out

MODES = {'coprime': coprime, 'bitrev': bitrev, 'alt-bitrev': altbitrev, 'golden': golden, 'es': es}
if __name__ == '__main__':
    n, m = int(sys.argv[1]), sys.argv[2]
    p = MODES[m](n); assert sorted(p) == list(range(n))
    print(*p)
