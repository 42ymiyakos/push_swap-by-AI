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
