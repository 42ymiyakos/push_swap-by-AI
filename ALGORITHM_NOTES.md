# ALGORITHM NOTES

## Problem facts
- 11 ops; between two pushes only rotations/swaps happen, so with fixed "which element / which slot"
  the rotation cost of a push is exactly 1 + min(max(ra,rb), max(rra,rrb), ra+rrb, rra+rb)  (proven: the
  four direction combos are the only non-dominated ways; state after the push is identical).
- Info-theoretic lower bound: log(n!)/log(11): N=100 ≥ 152, N=500 ≥ 1089 (weak).
- Exact optimum for n ≤ 7 by BFS over (A,B) states (src/exact.hpp).

## Representation (src/rep.hpp) — independent reproduction of the reference idea
- K: values that never move (circularly increasing subsequence of the input).
- Query list q: A2B(v), B2A(v) for v ∉ K (single trip each).
- AOrder keys: ik[v] (initial position), tk[v] (target); BOrder keys bk[v]. Each stack is always
  cyclically sorted by key; an insertion goes where the cyclic order is kept.
- Validity: any query order with A2B(v) before B2A(v), any BOrder, and targets cyclically sorted by value
  → always decodes to a correct sort. Decoder: sorted arrays + top index, O(log n) per query.
- Converter: any swap-free single-trip command sequence → representation via ghost lists.

## QSA (src/qsa.cpp) — single-move simulated annealing (candidate, superseded)
Moves: shift query, move in BOrder, move target in AOrder, K add/remove. Incremental evaluation with
checkpoints + early stop when tops re-converge. ~200k moves/s at N=100 but plateaus ~450 (N=100).

## PSIG (src/psig.cpp + src/ig.hpp) — Iterated Greedy with exact DP reinsertion
- Destroy: values of a random query segment / value range / random set; values removed completely
  ("stack destroy").
- Reinsert each removed value (highest previous cost first): for candidate B slot b and A target slot t,
  every other query's cost with v present is recomputed exactly from stored per-query info
  ("anchors": the top after a query is either the pushed element or the successor of the popped key).
  Total delta(i,j) = A_b(i) + B_{b,t}(j) separates → best (i<j) by prefix minimum in O(Q).
  Also j==i (immediate round trip) and "keep in A" (K) options. Verified exact (check mode).
- t handling: provisional t0, then full t×j DP for the top-3 b.
- Acceptance: LAHC (history 1 = greedy ≤ current was best at short budgets).
- Construction: LIS from many start points (distinct sets), chunk A2B (middle-out top/bottom chunks with
  delayed rb merged into rr), B2A by beam search with an arc-cover heuristic (minimal arc of A that must
  still be visited + reaching the minimum).

## Additions v3–v8 (current best solver, `research/c/`)

### Initial solution (v3, v4)
- Candidate K sets: LIS read cyclically from 50 (N ≤ 300) / 30 (N > 300) start positions, deduplicated.
- Chunk A→B pass: chunk count m is searched per instance (N ≤ 150: {2,3,4}; N ≤ 300: {3,4,5,6}; N > 300: {4,6,8,10,12})
  and two boundary families: equal value ranges (cmode 0) and equal counts of moved values (cmode 1).
- Three stages: greedy B→A for every (K, m, family) → small beam (width 100 / 20) on the 3 best K per (m, family)
  → full beam (width 1000 / 300) on the 3 best overall.
- Empirical: at N=500 the final result is close to "initial − 410", so initial quality drives the result and its spread.

### A2A moves (v5)
- Query type 2: value u moves inside A by d steps with sa: forward = rotate u to top, sa,(ra,sa)^(d−1); backward =
  rotate pred(u) to top, sa,(rra,sa)^(d−1). Default: only d = 1 (a2amax 1).
- The rotations of an A2A run before a push can be paired with that push's B rotations (rr/rrr): comb_sync.
- The reinsertion DP offers "v as A2A at time i with target t" next to keep / A2B+B2A. The DP is inexact around A2A
  queries (anchors are approximations there), so after each reinsertion the total cost is recomputed exactly.

### Reinsertion speed (v6–v8; all changes keep C output byte-identical to C++ for fixed iterations)
- v6: the slot-independent part (query pass, prefix effect P, keep/A2A options, per-target arrays) is computed once per
  reinsertion and shared by the upper-bound DP, the sweep and the refine DP. Slot ranks over time come from the query
  list in O(Q) per slot. Peephole time reserve reduced to the measured ~0.5 s.
- v7: ubw = 1 (upper-bound slots next to the ±1 value neighbours). irot = 16: in the sweep, A2B(v) is only considered
  at times where A's top is within 16 rotations of v (in final solutions the A rotation before a pb is ≤ 9 for 99%).
- v8: at "far" times (A2B not allowed) a slot can only matter through the B2A test
  `pm[s] + G[s] + lbC2 + lbTail < bestSoFar`. A lower bound lbMin ≤ min_s(pm[s] + G[s]) is maintained (exact after each
  full step, then increased by the minimum middle effect). When it proves that no slot passes, the step is skipped:
  the middle-effect update of G (piecewise constant over slots, 3 intervals) and the rank update are deferred into
  difference arrays and materialized only before the next full step.
- Value of IG time at N=500: about 5 ops per second in the 30–40 s regime (measured), so every speed-up pays off directly.
