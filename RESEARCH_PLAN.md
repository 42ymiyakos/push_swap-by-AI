# RESEARCH PLAN — push_swap autonomous research lab

## Objective
Build a push_swap solver with the lowest possible operation counts (focus N=100, N=500),
always correct, with reasonable runtime. Compare against: human baseline (pending access),
reference implementation (aaax8/push_swap @ 97d38fe), and our own candidates.

## Phases
1. Reconnaissance — done (see RESEARCH_LOG.md, 2026-11-21).
2. Baselines — reference built & measured (normal settings). Human baseline: repository
   `42ymiyakos/push_swap` exists but attaching/using it was blocked by the session permission
   classifier; needs explicit user authorization.
3. Reproduction — independent re-implementation of the reference *ideas*
   (query list + AOrder/BOrder representation, Iterated Greedy with best single-value
   reinsertion, LAHC). Our reinsertion uses an O(Q) DP per (B slot, A slot) instead of the
   reference's Mo's-algorithm/segment-tree machinery.
4. Optimization — candidate reduction, speed, destroy operators, acceptance.
5. Innovation — greedy-reinsertion construction, beam with arc-cover heuristic, two-phase
   constraint, (planned) A2A/sa moves, peephole post-optimizer, adaptive time allocation.
6. Final competition on held-out seeds (data/sets/eval_*).

## Hypotheses (status)
- H1: Single-query SA moves cannot coordinate (B position, B2A time) changes → plateau. CONFIRMED (QSA ~450 @N=100).
- H2: Joint exact reinsertion (i, j, b, t) is the key neighbourhood. CONFIRMED (big jump when all B slots considered).
- H3: Initial solution quality matters under limited time. CONFIRMED (case 2: 571→472 init via K choice → final 431→414).
- H4: N=500 is compute-bound; candidate pruning trades quality for speed. OPEN.
- H5: sa-based A2A moves give ~10-20 ops at N=100. OPEN.

## Budget policy
Level 0 smoke (1–3 cases), Level 1 (3–5 cases, ≤10 s each), Level 2 (≥20 cases dev sets),
Level 4 final on eval sets only. Never more than 4 concurrent solver processes (4 vCPU).
