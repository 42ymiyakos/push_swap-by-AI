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
