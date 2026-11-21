# NEXT STEPS (updated 2026-11-21, after v13)

## Current status (resume from here)
- **Best solver: PSIG v13 (C)** in `research/c/` (commit `b9142b0`). Same algorithm as the C++ research version
  `research/src/psig.cpp`. With `--iters=N` the two are byte-identical: check with
  `research/tools/cmp_c_cpp.sh SET NCASES ITERS` after **every** change to either version. (C-only speed-ups must keep identical output.)
- Held-out results (eval3 + eval4, 60 cases per size; reference CSVs `eval3_ref_*`, `eval4_ref_*` exist):
  - N=100: 356.5 median / 357.7 mean vs 363.5 / 363.9 → −6.18, p = 0.0003 (48 W / 2 T / 10 L).
  - N=500: 2849.5 / 2876.0 vs 2967.5 / 2967.5 → −91.5, p < 0.0001 (49 W / 11 L).
- Defaults: 15 s (n ≤ 150), 25 s (n ≤ 300), 45 s otherwise; irot 16; ubw 1; final beam 1000 / 500; kick 200 (×2); A2A on.
- Machine: 4 vCPU only. At most 4 solver processes; nothing heavy while a time-limited benchmark runs.
- Dev sets: `dev40_u100` (`--limit 40`), `dev40_u500` (`--offset 12 --limit 28`), `dev60_u500` (60 cases) for N=500 decisions.
  N=500 time-limited means vary about ±5 run to run over 28 cases; N=100 results are nearly deterministic.

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

## Prioritized next actions (after v13)
- Remaining weakness: N=500 spread (sd 75 vs 34) driven by the initial solution; worst cases ~3060.
- More C speed (identical output): sweep near-x loop, solve_prep, per-reinsertion candidate building (build_bl/al qsort).
- More diversity at N=100 (kick helped: −1.7 on dev): other perturbation types, adaptive kick strength.

## Older list
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
