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
38. Coarse slot sweep (every k-th slot + exact refine of neighbours): IG at most 16% faster, quality equal → option only.
39. Peephole K=8 at N=500 saves 13 more ops but costs ~10 s; at equal total time it was +45 worse (16/16 losses):
    IG time at N=500 is worth about 5 ops/s in the 30–40 s regime.
40. Wall-time profile (C): two exact DP solves per reinsertion ≈ 50% of IG time (callgrind Ir suggested 20%).
    Shared prep (query pass, per-target arrays, keep/A2A options computed once; slot ranks from the query list in O(Q))
    → identical output, ~10% faster. Peephole reserve 2/4 s → 0.8/1.0 s (measured 0.4–0.5 s).
    v6 vs v5 (dev, time-limited): N=500 −24.7 (16/16 wins, p < 0.001), N=100 −1.1 (11 W / 1 L, p = 0.002).
41. FAILED: three-pass (vectorizable) sweep with snapshot pruning: not faster (19.5 → 20.1 s / 300 iterations).
42. New chunk families (half-shifted ranges, mirrored ranges/counts): individually worse (3416–3450 vs 3339 mean init);
    adding them to the pool gives −7..−8 for +3..5 s of init → rejected (worse than 5 ops/s).
43. Upper-bound slots ±1 instead of ±2 value neighbours (ubw=1): about 13% faster IG; dev N=500 −1.7 (11 W / 4 L), N=100 −0.1 →
    adopted as default.
44. Final solutions: A rotation before a pb is small (mean 1.4, p99 9, max 21 at N=500) → option irot: in the sweep,
    A2B(v) is not considered at times where A's top is farther than irot from v. At fixed iterations irot=32 gives the same
    result 21% faster, irot=16 is 30% faster with similar quality, and irot=8 is worse. Dev test pending.
45. irot dev results (time-limited, vs irot unlimited): N=100 no effect; N=500 irot=16 −10.75 (16 cases) and −5.4 (12 more),
    combined −8.5 (17 W / 11 L, p = 0.10); irot=12 −9.4 on the second slice. Adopted irot=16 as default (v7).
46. v8 (C only, identical output to C++): far A2B times (82% of times at N=500) are skipped as a whole when a lower bound
    on min_s(prefix-min + middle-effect sum) proves that no slot passes the B2A test; the middle-effect and rank updates are
    deferred as range adds (difference arrays) and materialized only when needed. 300 iterations: 15.5 → 13.5 s.
    Dev N=500 (28 cases): −13.6, 28 W / 0 L.
47. v8 held-out (60 cases): N=500 2884 median, mean −68.3 vs reference (49 W / 11 L); N=100 361.5, −2.1 (p = 0.18).
48. FAILED: BOrder-segment destroy (dbo 2/4): N=100 +0.6 / +2.55.
    Parameter re-tune at N=100 (pqd 0.25/0.75, seglen 5/12, a2amax 2): all worse (+1.1..+3.6). The time-limited N=100 result
    is nearly deterministic (base rerun differs in 1 of 40 cases), so these differences are real: the defaults are a local optimum.
    6 s vs 15 s at N=100: 361.9 vs 356.4, and the best of 2 seeds × 6 s is 359.3 → multi-start rejected; time (speed) matters at N=100 too.
49. v9 (C, identical output): per-slot prep arrays cached (no recomputation for the top-3 in tmode 1); exact fast modulo; comb_cost
    fast path when there are no pending A2A rotations. About 12% faster at N=100.
50. v9 vs v8 (dev): N=500 −4.3 (22 W / 0 L), N=100 −0.35 (6 W / 1 L) → main.
    FAILED: skipping the exact re-decode when no A2A query is present: solutions almost always contain A2A queries → no gain; reverted.
51. Cost-biased destroy window (dbias = best of k random windows by value cost): N=100 +1.85 (k=2) / +2.58 (k=3), N=500 −3.8 (p = 0.52)
    → not adopted (option kept, default 1).
52. v10 (C, identical output): with A2A on, the full_cost calls after the destroy and after each query-destroy removal are unused
    (the cost is recomputed after every reinsertion) → removed; a few % faster.
53. N=500 parameter re-tune vs v10 base (seglen 10/24, pqd 0.3/0.7): −4..−8, none significant → defaults unchanged.
    N=500 time-limited runs vary run to run by about ±5 in mean over 28 cases.
54. Init stage timing (C v10): N=500 stage1 0.43 s, stage2 (30 beams width 20) 3.0 s, stage3 (3 beams width 300) 4.7 s;
    N=100 total 3.6 s. Beam profile: qsort ≈ 25% → v11: stamp dedup, top-k selection, inline merge sort (identical output):
    init 9.6 → 5.0 s (N=500), 4.4 → 2.9 s (N=100). Dev: N=500 −3.2 (20 W / 5 L, p = 0.011), N=100 unchanged.
55. With the faster beam: final beam 500 −9.6 (p = 0.43), stage-2 beam 40 −4.4 (p = 0.66) on 28 cases → 60-case test running.
56. 60-case dev (dev60_u500): final beam 500 vs 300: −8.8 (34 W / 25 L, p = 0.21), consistent with −9.6 on 28 cases → v12: bw=500 for N>300.
57. N=100 greedy stage-1 screening (bws 0): init 423 → 429, final +1.7 (p = 0.21) → keep bws 20 at N ≤ 300.
58. v12 held-out (60 cases): N=500 2849.5 median, −92.9 vs reference (49 W / 11 L); N=100 361.0, −2.9 (p = 0.07).
59. Kick: after K iterations without improving the current solution, one destroy with seglen/vwin/rmax × M is accepted
    unconditionally (the best is kept separately). Dev N=100 vs v12: K=200 M=3 −1.40 (p = 0.07), K=500 −1.18, K=100 −1.10,
    K=200 M=2 −1.73 (24 W / 10 L, p = 0.033), M=5 +0.05. N=500 (K=200, M=2): −1.0 (11 W / 5 L, rarely triggers).
    v13: kick 200 / ×2 default.
60. Kick from the best solution instead of the current one (kickbest): −0.17 (p = 0.85) → option only. Kick thresholds 100/300 (×2) ≈ 200.
