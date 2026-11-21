# NEXT STEPS (updated 2026-11-21, after v6)

## Current status (resume from here)
- **Best solver: PSIG v6 (C)** in `research/c/` (commit `17c62ae`, results committed in `dcb8528`). Same algorithm as the C++
  research version `research/src/psig.cpp`. With `--iters=N` the two are byte-identical: check with
  `research/tools/cmp_c_cpp.sh SET NCASES ITERS` after **every** change to either version.
- Held-out results (eval3 + eval4, 60 cases per size; reference CSVs already exist: `eval3_ref_*`, `eval4_ref_*`):
  - N=500: 2907 median / 2921 mean vs reference 2967.5 / 2967.5 → mean −46.5, p < 0.0001.
  - N=100: 362.5 / 361.8 vs 363.5 / 363.9 → mean −2.1, p = 0.21 (not significant).
- Default time budget: 15 s for n ≤ 150, 25 s for n ≤ 300, 40 s otherwise. Peephole reserve 0.8 / 1.0 s.
- Machine: 4 vCPU only (cannot be increased). Run at most 4 solver processes at a time, and nothing heavy while a
  time-limited benchmark is running, since timing affects results.

## Key facts learned (do not re-derive)
1. N=500 final ≈ initial − ~410. The initial solution quality drives the result and its variance.
2. IG time at N=500 is worth about 5 ops/s in the 30–40 s regime. Any speed-up of the reinsertion pays off directly.
3. Wall-time split per reinsertion (N=500): sweep ≈ 2.2 ms, exact DP solves ≈ 1.7 ms. Callgrind instruction counts are
   misleading here: measure wall time.
4. Greedy acceptance beats LAHC/racing/restarts. A2A (sa, 1 step) helps at N=100 once the initial solution is good.
5. Failed or neutral (see RESEARCH_LOG 1–41): QSA, LAHC, two-phase constraint, IS+swap init, vectorized/unpruned sweep,
   three-pass sweep, coarse slot sweep, racing, wider init screening, peephole K=8 under a fixed budget.

## Dev/eval protocol
- Develop and compare on `dev40_*` (N=500: `--offset 12 --limit 16`; N=100: `--limit 40`), paired against the previous
  version with a sign-flip test. The python snippet pattern is in RESEARCH_LOG / earlier commits.
- Confirm on eval3 + eval4 (60 cases) only for candidates that win on dev. Do not tune on eval sets. eval/eval2 are
  already contaminated (looked at): do not use them for claims.

## Prioritized next actions
1. **N=500 variance / worst cases**: add more structurally different initial-solution families to the candidate pool
   (e.g. different round orders for the chunk A→B pass, top-only rounds, three chunks per pass, different K choices).
   Measure the mean and the *max* of the initial solution on dev first, then the final result.
2. **Reinsertion speed** (wall time): fewer upper-bound slots (±1 instead of ±2 value neighbours), a cheaper upper bound,
   or skipping the refine solve when the sweep winner equals the upper-bound slot. Every change must keep C/C++ identical.
3. **N=100 significance**: more search diversity per unit time (destroy operators: B-order segments, A-position ranges),
   and A2A with 2 steps (`a2amax`) re-tested on dev.
4. Time budget alignment: the reference takes ~48 s at N=500 and ~17 s at N=100. Decide whether to match it
   (expected −10..−30 at N=500) or keep the current, shorter defaults.
5. Housekeeping: README/ALGORITHM_NOTES sections for v3–v6 (chunk-count search, equal-count chunks, A2A, shared prep).
