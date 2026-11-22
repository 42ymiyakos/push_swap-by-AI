# BEST RESULTS (updated 2026-11-22, v17)

## Time–quality curve vs the reference deep mode (v17, deep_u100 / deep_u500 = the reference's own 50+50 deep benchmark inputs)
Reference deep mode: about 8 h per case (its own logs). Ours: `--time=T`, 4 cases in parallel on 4 vCPU.
| N | ours, time per case | mean | median | sd | max | diff vs deep | W/T/L | p |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | reference deep (~8 h) | 345.58 | 342 | 10.7 | 386 | | | |
| 100 | 15 s | 353.76 | 354 | 8.6 | 373 | +8.18 | 8/1/41 | <.0001 |
| 100 | 60 s | 348.98 | 349.5 | 7.1 | 365 | +3.40 | 16/2/32 | .028 |
| 100 | 300 s | **344.34** | 345 | 7.3 | 356 | **−1.24** | 27/2/21 | .40 |
| 500 | reference deep (~8 h) | 2684.26 | 2682.5 | 29.9 | 2765 | | | |
| 500 | 45 s | 2837.30 | 2837.5 | 45.3 | 3035 | +153.0 | 0/0/50 | <.0001 |
| 500 | 180 s | 2786.84 | 2778 | 47.9 | 2948 | +102.6 | 0/0/50 | <.0001 |
| 500 | 600 s | 2753.94 | 2755.5 | 40.8 | 2908 | +69.7 | 1/1/48 | <.0001 |

N=100: at 300 s (≈ 1/100 of the reference's time) we are level with the deep mode (mean slightly better, median 3 worse,
not significant; our max is much lower). N=500: still clearly behind; the gain per time doubling falls from ≈ 25 (45→180 s)
to ≈ 19 (180→600 s). Extrapolating that rate would need roughly 2 h per case to close the gap, so the long-time behaviour
of our IG (greedy acceptance, fixed destroy sizes) is the main thing to fix at N=500.

## v17 = v16 + IS+swap start candidates for structured inputs — structured "hard" inputs (research/tools/worst_gen.py)
Families from keusidan/42-push_swap-solver docs/04 (re-implemented). Random inputs never pass the gate
(IS+swap length ≥ 1.6 × LIS), so v17 is identical to v16 there (the held-out numbers below apply).
| N=500 input | v16 | **v17** | reference | keusidan best |
|---|---:|---:|---:|---:|
| coprime (49i mod 500) | 2980 | 3005 (2964/2967 on reruns; v16 2966/2962) | 3276 | 3477 |
| bitrev | 3074 | 3093 (same algorithm, noise) | 3380 | 3543 |
| alt-bitrev | 2875 | 2905 (same algorithm, noise) | 3079 | 3390 |
| golden | 3246 | 3236 | 3350 | 3529 |
| es (√n descending blocks) | 1968 | **1545** | 1898 | 2179 |
N=100: v17 coprime 388, bitrev 395, alt-bitrev 351, golden 385, es 272 vs reference 418 / 406 / 365 / 418 / 276.


## v16 = v15 + IG alternating between P and P^-1 (8 phases) — held-out eval3 + eval4, 60 cases
| N | Solver | Median | Mean | Std | Max | s/case |
|---:|---|---:|---:|---:|---:|---:|
| 100 | **PSIG v16 (C)** | **354** | **353.7** | 6.7 | 369 | 14.7 |
| 500 | **PSIG v16 (C)** | **2833** | **2838.9** | 50.9 | 3064 | 44.3 |
Paired vs v15: N=100 −1.52 (33 W / 9 T / 18 L, p = 0.032); N=500 −8.2 (40 W / 20 L, p = 0.002).
Paired vs reference: N=100 **−10.2** (55 W / 5 L, p < 0.0001); N=500 **−128.6** (58 W / 2 L, p < 0.0001).
Dev: N=500 phases 2/4/8/16/32 → −5.5/−6.8/−10.0/−7.6/−5.3; N=100 8/16/32 → −2.35/−1.40/−0.72.


## v15 = v14 + dual initial solutions (inverse permutation, 3 full beams per direction) — held-out eval3 + eval4, 60 cases
| N | Solver | Median | Mean | Std | Max | s/case |
|---:|---|---:|---:|---:|---:|---:|
| 100 | **PSIG v15 (C)** | **354.5** | **355.2** | 7.0 | 370 | 14.7 |
| 100 | PSIG v14 | 357.5 | 357.9 | 8.2 | 384 | 15.8 |
| 100 | reference | 363.5 | 363.9 | | | 16.9 |
| 500 | **PSIG v15 (C)** | **2839.5** | **2847.1** | 48.9 | 3052 | 44.3 |
| 500 | PSIG v14 | 2841.5 | 2870.4 | 75.3 | 3060 | 44.4 |
| 500 | reference | 2967.5 | 2967.5 | 34.3 | 3041 | 47.7 |
Paired vs v14: N=100 −2.68 (25 W / 26 T / 9 L, p = 0.004); N=500 −23.3 (41 W / 13 L, p = 0.009).
Paired vs reference: N=100 **−8.7** (51 W / 1 T / 8 L, p < 0.0001); N=500 **−120.4** (57 W / 3 L, p < 0.0001).
Dev (dev60_u500): dual top-3 overall −7.1, top-6 overall −7.4, 2 per direction −10.5, 3 per direction −13.3 (chosen).
The gain is concentrated on bad starts: N=500 sd 75 → 49; median almost unchanged.


## v14 = v13 algorithm, default time 15 s (n ≤ 150) / 45 s (n > 300) — held-out eval3 + eval4, 60 cases
| N | Solver | Median | Mean | Std | Max | s/case |
|---:|---|---:|---:|---:|---:|---:|
| 500 | **PSIG v14 (C, 45 s)** | **2841.5** | **2870.4** | 75.3 | 3060 | 44.4 |
| 500 | reference | 2967.5 | 2967.5 | 34.3 | 3041 | 47.7 |
N=500 paired: mean **−97.2**, 95% CI [−115.4, −77.9], 49 W / 11 L, p < 0.0001 (still faster than the reference).
N=100: 16 s gave the same as 15 s (−6.02 vs −6.18) → 15 s kept; N=100 numbers below (v13) apply.


## v13 (C) vs reference — held-out eval3 + eval4, 60 cases per size, identical inputs
| N | Solver | Median | Mean | Std | Max | s/case |
|---:|---|---:|---:|---:|---:|---:|
| 100 | **PSIG v13 (C, 15 s)** | **356.5** | **357.7** | 8.3 | 384 | 14.8 |
| 100 | reference | 363.5 | 363.9 | 9.3 | 388 | 16.9 |
| 500 | **PSIG v13 (C, 40 s)** | **2849.5** | **2876.0** | 75.1 | 3061 | 39.4 |
| 500 | reference | 2967.5 | 2967.5 | 34.3 | 3041 | 47.7 |

Paired (ours − reference):
- N=100: mean **−6.18**, 95% CI [−9.0, −3.3], 48 W / 2 T / 10 L, **p = 0.0003**.
- N=500: mean **−91.5**, 95% CI [−109.9, −72.2], 49 W / 11 L, **p < 0.0001**.
→ **PSIG v13 (C) is the best known solver for both N=100 and N=500** (faster and with about 1/4–1/5 of the memory).
Commit `b9142b0`.


## v12 (C) vs reference — held-out eval3 + eval4, 60 cases per size
| N | Solver | Median | Mean | Std | Max | s/case |
|---:|---|---:|---:|---:|---:|---:|
| 100 | **PSIG v12 (C, 15 s)** | **361.0** | **361.0** | 9.1 | 388 | 14.7 |
| 100 | reference | 363.5 | 363.9 | 9.3 | 388 | 16.9 |
| 500 | **PSIG v12 (C, 40 s)** | **2849.5** | **2874.6** | 74.7 | 3060 | 39.4 |
| 500 | reference | 2967.5 | 2967.5 | 34.3 | 3041 | 47.7 |

Paired: N=500 mean **−92.9**, 95% CI [−111.4, −73.6], 49 W / 11 L, p < 0.0001.
N=100 mean −2.9, 95% CI [−5.8, +0.1], 37 W / 1 T / 22 L, p = 0.07.
Best solver: PSIG v12 (C), commit `67e937c`.


## v8 (C) vs reference — held-out eval3 + eval4, 60 cases per size
| N | Solver | Median | Mean | Std | Max | s/case |
|---:|---|---:|---:|---:|---:|---:|
| 100 | **PSIG v8 (C, 15 s)** | **361.5** | **361.8** | 9.0 | 389 | 14.8 |
| 100 | reference | 363.5 | 363.9 | 9.3 | 388 | 16.9 |
| 500 | **PSIG v8 (C, 40 s)** | **2884.0** | **2899.2** | 67.4 | 3056 | 39.5 |
| 500 | reference | 2967.5 | 2967.5 | 34.3 | 3041 | 47.7 |

Paired: N=500 mean **−68.3**, 95% CI [−85.4, −51.4], 49 W / 11 L, p < 0.0001.
N=100 mean −2.1, 95% CI [−5.0, +0.9], 36 W / 1 T / 23 L, p = 0.18.
Best solver: PSIG v8 (C), commit `9e638af`.


## v6 (C) vs reference — held-out eval3 + eval4, 60 cases per size, identical inputs
| N | Solver | Median | Mean | Std | Max | s/case |
|---:|---|---:|---:|---:|---:|---:|
| 100 | **PSIG v6 (C, 15 s)** | **362.5** | **361.8** | 9.0 | 389 | 14.8 |
| 100 | reference | 363.5 | 363.9 | 9.3 | 388 | 16.9 |
| 500 | **PSIG v6 (C, 40 s)** | **2907.0** | **2921.1** | 70.7 | 3105 | 39.5 |
| 500 | reference | 2967.5 | 2967.5 | 34.3 | 3041 | 47.7 |

Paired (ours − reference): N=100 mean −2.08, 95% CI [−5.1, +1.1], 37 W / 1 T / 22 L, p = 0.21 (slight edge, not significant);
N=500 mean **−46.5**, 95% CI [−63.9, −28.8], 46 W / 14 L, **p < 0.0001**.
Best solver: PSIG v6 (C), `research/c/`, commit `17c62ae`.


## v5 (C, A2A enabled): N=100 held-out eval3+eval4 (60 cases)
| Solver | Median | Mean | Std | Max |
|---|---:|---:|---:|---:|
| PSIG v5 (C, 15 s) | **362.5** | **362.2** | 9.2 | 389 |
| reference (~17 s) | 363.5 | 363.9 | 9.3 | 388 |
Paired: 37 wins / 1 tie / 22 losses, mean diff −1.73, 95% CI [−4.8, +1.5], p = 0.31 → slight edge, not significant.
v5 vs v4 on dev sets: N=100 −3.15 (p = 0.016, 40 cases); N=500 +4.3 (p = 0.59, 16 cases; C++ gave −7.6) → neutral.


## Final held-out evaluation v4 (C version, fresh seeds eval3 + eval4 = 60 cases per size)
Identical inputs, 4 concurrent processes on a 4-vCPU container. Default settings (15 s for N=100, 40 s for N=500).

| N | Solver | Median | Mean | Std | Max | s/case | Peak RSS |
|---:|---|---:|---:|---:|---:|---:|---:|
| 100 | **PSIG v4 (C)** | 365.0 | 364.5 | 9.6 | 392 | 13.6 | 56 MB |
| 100 | reference | 363.5 | 363.9 | 9.3 | 388 | 16.9 | 355 MB |
| 500 | **PSIG v4 (C)** | **2923.0** | **2941.0** | 66.3 | 3097 | 36.5 | 87 MB |
| 500 | reference | 2967.5 | 2967.5 | 34.3 | 3041 | 47.7 | 415 MB |

Paired comparison (ours − reference):
- N=100: 28 wins / 1 tie / 31 losses, mean diff +0.6, bootstrap 95% CI [−2.5, +3.8], sign-flip p = 0.73 → **parity**.
- N=500: 39 wins / 21 losses, mean diff **−26.5**, 95% CI [−44.1, −9.1], sign-flip p = **0.004** →
  **statistically credible improvement**. Our spread is still about twice the reference's.

### Verdict (v4)
- **N=500: PSIG v4 (C) is the best known solver** (promotion criteria met: correctness verified, identical inputs,
  p = 0.004, faster and smaller than the reference). The source is in `research/c/`, commit `d9d3dc3`.
- **N=100: tie** between PSIG v4 and the reference.

---
Earlier evaluations (kept for history):


## Final held-out evaluation v2 (fresh seeds `eval2_*`, never inspected during development)
Same inputs for both solvers, 4 concurrent processes on a 4-vCPU container.

| N | Solver | Cases | Median | Mean | Std | Min | Max | P90 | s/case | Peak RSS |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | reference (normal) | 20 | **366.5** | 365.3 | 9.0 | 348 | 382 | 379.0 | 16.3 | 347 MB |
| 100 | PSIG v2 (ours, 10 s) | 20 | 376.0 | 377.0 | 11.3 | 358 | 401 | 390.5 | 9.6 | 88 MB |
| 500 | reference (normal) | 12 | **2980.0** | 2986.8 | 55.2 | 2920 | 3157 | 2995.3 | 45.2 | 401 MB |
| 500 | PSIG v2 (ours, 40 s) | 12 | 3076.0 | 3057.8 | 118.5 | **2892** | 3294 | 3204.3 | 38.2 | 52 MB |

| 100 | human baseline (adaptive, default) | 20 | 608 | 807.4 | | 560 | 1084 | 1084 | 0.01 | <2 MB |
| 500 | human baseline (adaptive, default) | 12 | 5349 | 5895.7 | | 5148 | 6784 | 6784 | 0.01 | <2 MB |

Paired: N=100 ours better on 3/20 (1 tie), median diff +9.5; N=500 ours better on 5/12, median diff +101
(very high variance: best case 2892 vs 2954, worst 3294 vs 2996).

Evaluation v1 (`eval_*`, PSIG before K-screening/peephole changes): N=100 378 vs 361 (2/20 wins);
N=500 3026.5 vs 2985 (5/12 wins).

## C port (research/c, default settings, eval2)
| N | Cases | Median | Mean | Min | Max | s/case | Peak RSS |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 100 | 20 | 375.0 | 375.7 | 359 | 395 | 8.7 | 35 MB |
| 500 | 12 | 3029.5 | 3026.8 | 2848 | 3245 | 36.7 | 20 MB |
Byte-identical to the C++ version for fixed iteration budgets; ~13% faster per iteration.

## Verdict
- **Best known solver overall: the reference (aaax8/push_swap @97d38fe), normal settings.**
  Our candidate does not meet the promotion criterion (statistically credible improvement).
- Best own implementation: **PSIG v2** (`research/src/psig.cpp`, default settings), promoted to `research/best/`.
  Strengths: 1.7x faster at N=100, ~7x less memory, exact O(Q) DP reinsertion (much simpler than the
  reference machinery), wins ~40% of N=500 cases. Weakness: variance / basin quality.

## Small n (exact)
n ≤ 8 solved optimally by BFS (n=3: worst 2, n=5: worst 8 over all permutations).
n=10: mean 18.2 / max 21 (train set, 1 s); n=50: mean 150.2 / max 159.

## Correctness
All runs verified by an independent strict checker; special distributions (sorted, reverse, nearly sorted,
rotated, two runs, block-reversed, sawtooth, local shuffle) for n=3,5,10,100,500: all OK.
Valgrind: 0 errors, 0 leaks; ASan/UBSan: clean. Invalid input → "Error" on stderr.
Human baseline (42ymiyakos/push_swap @b9d2c06, unmodified, built in a scratch copy):
200-case `eval_*` sets — adaptive N=100 median 848 (538–1084), N=500 median 6784 (5025–6784);
`--medium` N=100 582.5 (538–626), N=500 5227.5 (4994–5397); `--complex` 1084 / 6784 (fixed, radix);
`--simple` 1459 / 32291.5. Small n (adaptive): n=3 max 2, n=5 max 10, n=10 max 65, n=50 max 467.
All 200/200 correct. Finding: random inputs have disorder ≈ 0.5, exactly the DISORDER_HIGH threshold,
so ~50% of random inputs go to radix (100/200 at N=100, 101/200 at N=500), which is why adaptive is
bimodal and worse than `--medium` on half the cases.
Vs human adaptive on eval2: PSIG v2 −38% (N=100, 376 vs 608) and −42% (N=500, 3076 vs 5349).
