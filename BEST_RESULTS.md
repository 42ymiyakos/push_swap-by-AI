# BEST RESULTS (as of 2026-11-21, session 1)

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
