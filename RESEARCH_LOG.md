# RESEARCH LOG (chronological)

## 2026-11-21 — session 1
1. Recon: empty repo (branch claude/practical-hamilton-31vohm). 4 vCPU Xeon 2.1 GHz, 15 GB RAM,
   gcc 13.3, clang 18, cmake/ninja, valgrind, python 3.13. No perf.
2. Human baseline: repository `42ymiyakos/push_swap` is listed for this account, but attaching it
   (add_repo) and building a read-only clone were both blocked by the session's permission classifier.
   **Not used.** Needs explicit authorization from the user.
3. Reference aaax8/push_swap @97d38fe built (Release, -O2). Normal settings reproduced the README:
   N=100 median 362.5 (20 dev cases, 16 s/case here), N=500 median 2959.5 (20 cases, 44 s/case).
   (README: 363 / 2970, author reported ~90 s / 140 s per case on Windows.)
4. Read the reference report (docs/push_swap_report.qmd): query list + AOrder/BOrder, IG with best
   (A2B, B2A, BOrder) reinsertion, LAHC, chunk+beam initial solution, LCS sync, bidirectional BFS peephole.
5. QSA (single-move SA over the representation, incremental evaluation with checkpoints): ~450 @N=100.
   Diagnosis: reference spends more in phase 1 to make phase 2 cheap; single moves cannot coordinate.
6. PSIG: exact O(Q) DP reinsertion for fixed (B slot, A slot), decomposition delta = A_b(i) + B_bt(j).
   Bug found via delta-check mode: A-side candidate keys ignored absent values' init keys → key
   collisions → respace changed the order. Fixed; delta now exact (verified on thousands of reinsertions).
7. All B slots instead of value-neighbour window: large gain (N=100 @1500 it: 465→411 on case 1).
8. LAHC history>1 hurt at every budget tried (10–30 s); greedy (≤ current) best.
9. Multi-start LIS (distinct LIS from many start points, screened by small beam) — init 571→472 on case 2.
10. Beam heuristic: arc-cover of A positions that still need an insertion (+ reaching min) → better init.
11. All-slot time-major sweep (slot arrays updated per query) + cheap modulo: 2.2x faster, identical results.
12. N=500 with all slots: 3191 @30 s (windowed: 3343 @30 s, 3316 @120 s).
13. A2A queries (sa-based moves inside A) + sync of A2A rotations with the next push: correct but no
    measurable gain at 10 s on N=100 (377.4 → 377–380). Kept as option (off by default).
14. Query-destroy (destroyed values stay in A while others are reinserted): small gain (377.4 → 375.9).
15. Racing restarts / restarts from different K: no significant gain.
16. Lower-bound pruning in the sweep with an upper bound from neighbour slots: ~15% faster at N=500.
17. Research review: N=100 plateaus at ~375 (ref 362) regardless of time (10→60 s: 377→374).
    Per-case gap 0–25; ref uses ~10 swaps per case; our solutions use fewer swaps and more pushes.
18. Final eval v1 (eval sets): N=100 378 vs 361, N=500 3026.5 vs 2985. Outlier N=500 case traced to a bad
    initial K (3826); K screening 3→30 starts fixed it (3410). Since this used an eval case, a fresh
    eval2 set was generated for an uncontaminated evaluation.
19. Peephole optimizer (meet-in-the-middle BFS on local stack windows): −4 @N=100, −15 @N=500.
20. Final eval v2 (eval2, fresh): N=100 376 vs 366.5; N=500 3076 vs 2980 (std 118 vs 55; 5/12 wins).
    Conclusion: competitive but not better; main weakness at N=500 is run-to-run variance.
21. Human baseline measured (user authorized clone): adaptive N=100 median 848 / N=500 6784 on 200-case
    eval sets; --medium 582.5 / 5227.5. Adaptive sends ~50% of random inputs to radix (disorder ≈ 0.5 =
    DISORDER_HIGH).
22. C port of PSIG (research/c, ~1300 lines, default configuration only). C++ sorts made deterministic
    (tie-breaks in beam/IG/peephole) so both versions are byte-identical with fixed --iters (verified on
    n=5,10,50,100,500 + special sets). ASan/UBSan/Valgrind clean. eval2: N=100 375, N=500 3029.5.

## 2026-11-21 — session 1 (continued, after repository hand-over)
23. Diagnosis of N=500 variance (eval2): final result tracks the initial solution almost monotonically
    (init 3305–3330 → final 2848–2905; init 3626–3705 → final 3038–3245).
24. Initial solution: chunk count matters per instance (case: m=4..10 → 3683/3946/3705/3650/3495/3528).
    New construction = greedy screening of all LIS candidates per chunk count, small beam on the top 3 per
    count, full beam on the overall top 3. N=500 (counts 4,6,8,10,12, small beam 20): bad cases −200..−250.
    N=100 (counts 2,3,4, small beam 100): mean init 445 → 423 (= reference initial median).
25. dev40_u500 (12 cases, 40 s): v2 3014 → v3 2927.5 median (mean 3015 → 2954); ref 2959.5 / 2968; 9/12 wins.
    Final beam 300 instead of 100: mean 2953.6 → 2940.1; 600 worse (less IG time).
26. dev40_u100 (20 cases): v2 375.5 → v3 365.0 (10 s); 16 s: 362.0 (ref 362.5).
27. Bug: passing any large-n option disabled all size-dependent defaults; fixed (defaults first, explicit
    options override).
28. Fresh held-out eval3 (20+20 cases, C port v3 vs reference): N=100 median 364.0 vs 364.5, mean 363.2 vs
    365.8, 11/20 wins (13.7 s vs 16.1 s); N=500 median 2915 vs 2963.5, mean 2944 vs 2960, 12/20 wins
    (36.7 s vs 44.4 s). Paired sign-flip p ≈ 0.39 / 0.43 → parity, not a statistically credible win.
    Weakness: N=500 spread (sd 89 vs 40).
29. Worst dev case (init 3575): 120 s still 3049 (ref 3028); wider chunk-count search no help.
    Equal-count chunk boundaries (cmode 1) give 3285 on it. v4 = both chunk modes in the candidate pool:
    dev40_u500 median 2927 / mean 2938 / sd 49 / max 3075 (v3: 2923.5 / 2940 / 63 / 3109);
    dev40_u100 (15 s) 362.0 / 363.3; N=100 mean init 422.6 → 417.7.
30. C beam heuristic: removed qsort (points already ascending), incremental B-value bitset, and B-side marking when
    |B| < |A| → N=500 init 9.8 → 8.5 s, identical output. -O3 -march=native: about 11% faster.
31. Large held-out eval (eval3 + eval4, 60 cases per size, C v4): N=100 parity (mean +0.6, p = 0.73);
    N=500 mean −26.5 (CI [−44, −9], p = 0.004), 39/60 wins → v4 promoted as best for N=500.
32. FAILED: branch-free, unpruned "vectorizable" sweep: 300 iterations 18.7 s → 24.7 s (the lost pruning costs more than
    SIMD gains). Reverted (sources restored to d9d3dc3). The lower-bound pruning is worth keeping.
33. Racing from top-2/3 initial solutions at N=500: worse (mean 2939 → 2959/2974).
34. Sweep: c1 tabulated by B-rank (O(sb) per step): about 5% faster, identical output. Tabulating all six quantities was slower.
35. A2A (sa moves, a2amax 1) re-tested with the better init: N=100 −3.3 (C++, 40 cases, p = 0.012), −3.15 (C, p = 0.016);
    N=500 neutral. Ported A2A + rotation sync to C (byte-identical to C++). v5 = A2A on by default.
    Held-out N=100 (60 cases): 362.5 vs ref 363.5 median, mean −1.73, p = 0.31.
    Valgrind reports a 'definite' leak only for -march=native builds (LSan and generic -O3 are clean): treated as a false positive.
36. Init screening accuracy at N=500 (bwm 50, ktop 5, ...): at most −5 mean init for +4–11 s → not adopted.
37. FAILED: IS+swap initial solutions (LIS with swappable consecutive pairs; DP + Fenwick; sa applied when a pair reaches
    the top; converter maps sa to an A2A query). Kept values 43 → 53 at N=500, but mean init only −15 (8 cases) and
    final +7.4 (16 cases, p = 0.41). Option --isswap=1 kept in C++ (default off), not ported to C.
