#!/usr/bin/env python3
"""Analyze an op sequence: phase split (before/after last pb), interleaving, op histogram."""
import sys, collections
ops = [l.strip() for l in open(sys.argv[1]) if l.strip()]
last_pb = max((i for i, o in enumerate(ops) if o == 'pb'), default=-1)
first_pa = min((i for i, o in enumerate(ops) if o == 'pa'), default=len(ops))
p1, p2 = ops[:last_pb + 1], ops[last_pb + 1:]
c1, c2 = collections.Counter(p1), collections.Counter(p2)
print(f"total={len(ops)} lastpb@{last_pb} firstpa@{first_pa}")
print("phase1:", len(p1), dict(c1))
print("phase2:", len(p2), dict(c2))
