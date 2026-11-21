# NEXT STEPS (prioritized)

1. **Human baseline**: get user authorization to clone/build `42ymiyakos/push_swap` (read-only) and benchmark it.
2. ~~Peephole post-optimizer~~ (done: −4 / −15).
3. **N=500 variance**: restarts with K diversity at N=500, adaptive time allocation (stop early when stagnant, restart).
   Original item 2:**Peephole post-optimizer** (bidirectional BFS over windows of the final command list, as in the
   reference ch.31; reference gains ~2 @N=100, ~15 @N=500). Cheap, independent of the search.
3. **Why A2A does not pay off**: reference uses ~10 sa/ss per N=100 case. Hypotheses: (a) our A2A option
   is evaluated with stale anchors (DP inexact around A2A: ~1–2% of reinsertions mismatch), (b) the
   IS+swap structure needs to be present in the *initial* solution (reference builds IS+swap first),
   (c) B-side swaps (sb) are not representable at all.
4. **N=500 speed**: sweep is O(Q·m) comb evaluations per value; next: incremental per-slot evaluation
   only where effects change (piecewise-constant structure), or reference-style candidate reduction.
5. **Multi-trip values / B2B queries** (representation limits shared with the reference).
6. Stronger initial solution for N=500 (reference 3253 vs ours ~3500).
7. Larger final evaluation (≥100 cases per size) and adversarial inputs once the solver stabilizes.
