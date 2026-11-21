# EXPERIMENTS

Conventions: sets in `research/data/sets/` (generator `tools/gen.py`, seeds embedded in the file
recipe `tools/make_sets.sh`). train_u100 = seeds 1000.., dev40 = 2000.., eval = 500000.., eval2 = 800000..
Per-case CSVs and summaries: `research/bench/runs/` (`summary.jsonl`). Times are wall-clock per case on
the 4-vCPU container (jobs ≤ 4). "it" = IG iterations.

## E0 Reference reproduction (Level 2)
| Set | Cases | Median | Mean | Min | Max | s/case |
|---|---:|---:|---:|---:|---:|---:|
| dev40_u100 (first 20) | 20 | 362.5 | 361.2 | 343 | 373 | 15.8 |
| dev40_u500 (first 20) | 20 | 2959.5 | 2963.8 | 2922 | 3031 | 43.8 |
Conclusion: README numbers (363 / 2970) reproduced. Here the reference is ~5x faster than reported.

## E1 QSA (single-move SA), N=100 case train#1 (Level 1)
| Config | Result |
|---|---|
| constructor (LIS + 2 chunks + greedy B2A) | 530 |
| 1 s, T0=3→0.2 | 512 |
| 5 s, T0=1.5→0.1 | 456 |
| 20 s, T0=3→0.2 | 458 |
3-case grid at 5 s: T0∈{1.5,0.8,0.5,0.3}: ~445–495. Reference on same case: 359. → abandoned.

## E2 PSIG candidate set for B slots (N=100, cases 1–3, 1500 it)
| B candidates | results | time |
|---|---|---|
| value window ±4 | 465 / 492 / 456 | 5 s |
| value window ±10 | 426 / 484 / 421 | 11 s |
| all slots | 411 / 465 / 411 | 37–40 s |
→ all slots (later made fast by the sweep).

## E3 t-handling (all slots, 1500 it): all t: 411/465/411; provisional t0=first + top-3 refine: 407/469/427;
t0=last: 416/454/416; top-1: 408/468/421 → equivalent; keep t0=first, top-3.

## E4 Acceptance (N=100, 3 cases, 8 s): lahc 1: 412/431/401, 50: 408/458/415, 200: 415/469/429, 1000: 420/474/413.
At 30 s on 8 cases: lahc 1 → 375.0, 30 → 384.1, 300 → 386.4, 3000 → 390.4 (means). → greedy.

## E5 Two-phase constraint (all A2B before B2A), 8 s: 381/441/405 vs free 412/432/401 → mixed; not adopted.

## E6 Initial solution
- K choice (case 2, beam 200): single random LIS 571 → best of 20 starts 472.
- Beam with first heuristic made N=500 *worse* than greedy (2849 vs 2700 phase-2). Arc-cover heuristic:
  N=100 bw 500: 460/517/469 (old bw1000: 468/554/487); N=500 bw 100: 3898 (greedy 4455).
- Greedy-reinsertion construction (from K only): N=100 comparable (453–531), N=500 bad (4705) → off.
- N=500 K screening 3 → 30 starts: outlier case 3826 → 3410.

## E7 Speed
- Separate per-b / per-t precomputation: 2.3x. All-slot sweep + cheap modulo: 2.2x (identical results).
- LB pruning in sweep: +15%. -O3 -march=native: no gain.

## E8 N=500 destroy segment length (30 s, dev case 1): 4: 3097, 8: 3095(40 s), 16: 3086, 32: 3132, 64: 3180.

## E9 N=100, 8 train cases (mean of op counts), 10 s unless noted
| Config | Mean | Median | Max |
|---|---:|---:|---:|
| sweep (baseline) | 377.4 | 379 | 389 |
| + A2A | 380.9 | 383 | 398 |
| + A2A with rotation sync | 378.6 | 382 | 390 |
| + A2A sync, a2amax 1 | 377.4 | 378.5 | 386 |
| + query-destroy 0.5 | 375.9 | 376.5 | 381 |
| + query-destroy 1.0 | 376.0 | 379.5 | 383 |
| qd0.5, topb 10 | 376.0 | 379.5 | 386 |
| qd0.5, t0 last | 376.9 | 379.5 | 386 |
| qd0.5, exact all-t (no sweep) | 379.6 | 383 | 399 |
| 16 s racing 4×, rfrac .5 | 375.8 | 378 | 384 |
| 16 s racing 8×, rfrac .6 | 374.6 | 375.5 | 388 |
| 16 s racing from 4 / 8 different K | 376.6 / 376.4 | | |
| 30 s sweep | 375.0 | 378 | 382 |
| 60 s sweep | 374.3 | 376 | 382 |
| reference (17 s) | 361.8 | 358 | 375 |
Seed variance at 5 s (4 seeds): case2 370–389, case3 372–391.

## E10 Peephole (window BFS, meet in the middle): N=100 saves 3–4 ops (K=4,6,8, ~1.7 s), N=500 saves ~15 (K=4,6, ~2 s).

## E11 Final evaluation v1 (held-out eval sets, commit c6d…; default settings, no peephole)
See BEST_RESULTS.md. N=100: 378 vs ref 361 (wins 2/20). N=500: 3026.5 vs ref 2985 (wins 5/12).

## E12 Final evaluation v2 (fresh `eval2_*`): see BEST_RESULTS.md.

## E13 Initial-solution search over chunk counts (v3)
| Set | Config | Mean init | Init time |
|---|---|---:|---:|
| dev40_u100 (10) | v3 default (m∈{2,3,4}) | 423.3 | 2.6 s |
| dev40_u100 (10) | + m=5 / kper 6 / ktop 6 / bw 3000 / nlis 100 | 422.3–423.3 | up to 5.5 s |
| dev40_u500 (8) | v3 (m∈{4..12 even}, bw 100) | 3398.4 | 4.5 s |
| dev40_u500 (8) | + m 14,16 | 3398.4 | 5.4 s |
| dev40_u500 (8) | kper 5 | 3381.4 | 6.0 s |
| dev40_u500 (8) | final beam 300 | 3344.1 | 8.9 s |

## E14 Final results v3 (dev sets)
| Set | Solver | Median | Mean | Std |
|---|---|---:|---:|---:|
| dev40_u500 (12) | v2 | 3014.0 | 3015.1 | 82.8 |
| dev40_u500 (12) | v3 bw100 | 2927.5 | 2953.6 | 71.1 |
| dev40_u500 (12) | v3 bw300 | 2923.5 | 2940.1 | 63.2 |
| dev40_u500 (12) | v3 bw600 | 2950.5 | 2968.4 | 69.7 |
| dev40_u500 (12) | reference | 2959.5 | 2968.3 | |
| dev40_u100 (20) | v2 10 s | 375.5 | 374.5 | 12.3 |
| dev40_u100 (20) | v3 10 s | 365.0 | 365.6 | 12.0 |
| dev40_u100 (20) | v3 16 s | 362.0 | 363.7 | 11.6 |
| dev40_u100 (20) | reference | 362.5 | 361.2 | 9.0 |
