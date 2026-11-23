# NEXT STEPS (updated 2026-11-23, after v20)

## Current status (resume from here)
- **Best solver: PSIG v20 (C)** = v19 with kick 100 for N > 150 (neutral at 45 s, −4.7 at 180 s, not significant; log 75).
- **Machine speed changes between containers** (log 75: +36 at 45 s for the same binary). Always rerun the baseline in the
  same container before comparing; never compare against CSVs from an earlier container.
- Previous: PSIG v19 (C): v18 + LAHC with a budget-scaled history (log 74). Held-out N=500 −21.4 vs v16, N=100 −1.8.
  Deep sets vs the reference deep mode: N=500 600 s +35.2 (was +69.7), N=100 300 s −3.0 (same as v18).
- Previous: PSIG v18 (C) (= v17 at default budgets; kick 50 when `--time` > 90 s)
- Previous: PSIG v17 (C) (= v16 on random inputs; IS+swap start group for structured inputs, `--isswap`)
- Previous: PSIG v16 (C) in `research/c/`: v14 + dual start (`--dual`, `--kdir`) + IG alternating between P and P⁻¹ in 8 phases (`--dualig`); C only.
  With `--dual=0 --dualig=1` it is byte-identical to the C++ version (`cmp_c_cpp.sh` passes it). Same algorithm as the C++ research version
  `research/src/psig.cpp`. With `--iters=N` the two are byte-identical: check with
  `research/tools/cmp_c_cpp.sh SET NCASES ITERS` after **every** change to either version. (C-only speed-ups must keep identical output.)
- Held-out results (eval3 + eval4, 60 cases per size; reference CSVs `eval3_ref_*`, `eval4_ref_*` exist):
  - N=100 (v16): 354 median / 353.7 mean vs 363.5 / 363.9 → −10.2, p < 0.0001 (55 W / 5 L).
  - N=500 (v16): 2833 / 2838.9 vs 2967.5 / 2967.5 → −128.6, p < 0.0001 (58 W / 2 L); sd 51.
- Defaults: 15 s (n ≤ 150), 25 s (n ≤ 300), 45 s otherwise; irot 16; ubw 1; final beam 1000 / 500; kick 200 (×2); A2A on.
- Machine: 4 vCPU only. At most 4 solver processes; nothing heavy while a time-limited benchmark runs.
- Dev sets: `dev40_u100` (`--limit 40`), `dev40_u500` (`--offset 12 --limit 28`), `dev60_u500` (60 cases) for N=500 decisions.
  N=500 time-limited means vary about ±5 run to run over 28 cases; N=100 results are nearly deterministic.

- Structured hard inputs: `data/sets/worst_u{100,500}.txt` (tools/worst_gen.py); check them after start-solution changes.

## Key facts learned (do not re-derive)
1. N=500 final ≈ initial − ~410. The initial solution quality drives the result and its variance.
2. IG time at N=500 is worth about 5 ops/s in the 30–40 s regime. Any speed-up of the reinsertion pays off directly.
3. Wall-time split per reinsertion (N=500): sweep ≈ 2.2 ms, exact DP solves ≈ 1.7 ms. Callgrind instruction counts are
   misleading here: measure wall time.
4. LAHC with a short history (≈ iterations/300 at N=500) beats greedy acceptance (log 74); long histories and racing/restarts do not. A2A (sa, 1 step) helps at N=100 once the initial solution is good.
5. Failed or neutral (see RESEARCH_LOG 1–41): QSA, LAHC, two-phase constraint, IS+swap init, vectorized/unpruned sweep,
   three-pass sweep, coarse slot sweep, racing, wider init screening, peephole K=8 under a fixed budget.

## Dev/eval protocol
- Develop and compare on `dev40_*` (N=500: `--offset 12 --limit 16`; N=100: `--limit 40`), paired against the previous
  version with a sign-flip test. The python snippet pattern is in RESEARCH_LOG / earlier commits.
- Confirm on eval3 + eval4 (60 cases) only for candidates that win on dev. Do not tune on eval sets. eval/eval2 are
  already contaminated (looked at): do not use them for claims.

## Prioritized next actions (after v17, agreed with the user)
0. Done (log 68–69): weak-pattern study → losses are mostly reference luck; real outliers come from the start
   choice (probing fixes them but is neutral on average). Rare ops: zero local potential (K=10 exact windows).
1. **Weak-pattern study** (data first): features of held-out cases where we lose to the reference or score badly
   (LIS length, initial-solution quality, local structure); form a hypothesis, generate inputs that amplify it
   (as with es); optionally an adversarial input search with a short-time proxy.
2. **Rarely used ops**: count op usage in our outputs (sb ≈ 0, sa/ss few). Estimate the upper bound of what sb/ss
   could save (e.g. peephole with larger windows / exact local re-solve on sample outputs) before extending the
   representation (B-side swaps, multi-step A2A).
3. **Distance to optimum** at small n (exact BFS for n ≈ 10–12): where is the remaining waste (rotations vs choice
   of pushed values)?
4. Done (log 70): deep comparison. N=100 level with the reference deep mode at 300 s; N=500 +70 at 600 s.
5. **Long-budget N=500 search, continued** (v19 halved the 600 s gap): refine the history rule (adaptive: from the measured
   iteration rate instead of the time), LAHC × kick interplay (kick 50 vs 200 with LAHC at 180 s), destroy-size schedule,
   elite pool / path relinking. The reference keeps more values in A (40–46 vs 25–33) — test a bias towards larger kept sets.
   Earlier: (log 71–72 done: convfix neutral, dualig 32 worse, kick 50 small/unconfirmed): the gain per time doubling drops (25 → 19). Try, on
   dev with 180 s budgets: threshold/late acceptance only after stagnation, growing destroy sizes over time, elite pool
   with restarts/path relinking, more dualig phases for long budgets. Keep the 45 s default unchanged unless it also wins there.

## Older list
1. **N=500 variance / worst cases**: add more structurally different initial-solution families to the candidate pool
   (e.g. different round orders for the chunk A→B pass, top-only rounds, three chunks per pass, different K choices).
   Measure the mean and the *max* of the initial solution on dev first, then the final result.
2. **Reinsertion speed** (wall time): fewer upper-bound slots (±1 instead of ±2 value neighbours), a cheaper upper bound,
   or skipping the refine solve when the sweep winner equals the upper-bound slot. Every change must keep C/C++ identical.
3. **N=100 significance**: more search diversity per unit time (destroy operators: B-order segments, A-position ranges),
   and A2A with 2 steps (`a2amax`) re-tested on dev.
4. Time budget alignment: done in v14 (45 s gave −5.7 at N=500; 16 s gave nothing at N=100).
5. Housekeeping: README/ALGORITHM_NOTES sections for v3–v6 (chunk-count search, equal-count chunks, A2A, shared prep).
