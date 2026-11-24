# push_swap solver research — consolidated report (v20, 2026-11-24)

This document summarizes where the project stands: how the solver works, how good it is, how far it is from
the optimum (as far as that can be measured), what did not work, and which directions are still open.
Chronological details: [RESEARCH_LOG.md](RESEARCH_LOG.md) (entries 1–78). All numbers: [BEST_RESULTS.md](BEST_RESULTS.md).

## 1. Summary

| | N=100 | N=500 |
|---|---:|---:|
| This solver, default budget (15 s / 45 s), held-out eval3+eval4, mean (v19 = v20 on random inputs at these budgets) | **351.9** | **2817.5** |
| Reference aaax8/push_swap, normal mode, same inputs | 363.9 | 2967.5 |
| Human implementation (42 validated), eval2 median | 608 | 5349 |
| Reference deep mode (~8 h per case), its own 50+50 deep inputs, mean | 345.58 | 2684.26 |
| This solver on the same deep inputs, longest run (N=100: 300 s, N=500: 600 s) | **342.58** | 2719.42 |
| Best known per case on the deep inputs (min over all our runs and the reference deep mode), mean | 338.74 | 2680.96 |
| Proven lower bound (information-theoretic, log(N!)/log 11) | 152 | 1089 |

- Random inputs, n ≤ 8: the solver returns the exact optimum on **every** permutation (checked over all n! inputs for n = 3..8).
- n = 9, 10, 11 (200 random inputs each): mean gap to the exact optimum +0.105 / +0.145 / +0.205 operations
  (+0.75% / +0.89% / +1.10%). Every one of the gap cases needs a swap of B (sb/ss) in all of its optimal solutions (section 4.2).
- N=100: better than the reference deep mode at 300 s (−3.0, p .033). N=500: still +35 behind the deep mode at 600 s and
  +38 behind the best known values; the gap shrinks by about 25 per doubling of time.
- **Optimality is not claimed anywhere above n = 8.** For N ≥ 12 the only proven lower bound is the very weak
  information-theoretic one; the true distance to the optimum at N=100/500 is unknown.

## 2. Algorithm design (research/c, ~2400 lines of C)

### 2.1 Representation (independent reproduction of the reference's idea)
A solution is not a command list. It is:
- **K**: values that never leave A (a cyclically increasing subsequence of the input);
- a **query list**: for each other value v, A2B(v) and later B2A(v) (each value makes exactly one round trip), plus
  optional **A2A** queries (move a value one step inside A with `sa`);
- **keys**: AOrder (initial and target positions in A) and BOrder (position in B). Each stack is always cyclically
  sorted by key, so the place where a value is inserted is implied.

Any query order that respects "A2B before B2A", with any BOrder, decodes to a correct sort. The decoder emits the
rotations between consecutive pushes with the cheapest of the four direction combinations
(1 + min(max(ra,rb), max(rra,rrb), ra+rrb, rra+rb)). This choice is provably optimal for a fixed pair of push targets.
Rotations of A2A runs are merged with the next push (rr/rrr). A converter maps any swap-free single-trip command
sequence back into this representation.

### 2.2 Construction (≈ 28% of the time at N=500)
1. Candidate K sets: LIS read cyclically from 50 (N ≤ 300) / 30 (N > 300) start positions, deduplicated. For
   structured inputs, "increasing after swapping adjacent pairs" (IS+swap) candidates are added as a separate group.
2. A→B: chunk passes, middle-out, with delayed rb merged into rr. The chunk count is searched per instance
   ({2,3,4} / {3..6} / {4,6,8,10,12}), with two boundary families (equal value ranges, equal moved counts).
3. B→A: beam search with an arc-cover heuristic (the minimal arc of A that must still be visited, plus reaching the minimum).
   The stages are a greedy screen of every (K, m, family), then small beams on the top 3 per (m, family), then full beams
   (1000 / 300 / 500 wide) on the overall top 3.
4. **Dual start**: everything above also runs on P⁻¹. Its op sequence, reversed and with every op inverted, sorts P.
   The shorter result is kept. This removes most bad starts (N=500 worst case 3510 → 3408).

### 2.3 Improvement: Iterated Greedy with exact DP reinsertion (≈ 70% of the time)
- **Destroy**: remove a few values. The removal is a random query segment, value range or random set (segment length
  ~16 at N=500). Removed values disappear completely from the solution.
- **Repair**: values are reinserted one by one, the most expensive first. For one value v, the choice is joint:
  - when it goes to B (time i), when it comes back (time j), where it sits in B (slot b) and where it lands in A (t);
  - or "keep in A" (add to K);
  - or "A2A at time i".

  For fixed (b, t) the cost change of every other query is computed exactly from stored per-query anchors. The total
  delta splits as A_b(i) + B_{b,t}(j), so the best i < j is a prefix minimum: O(Q) per (b, t). All B slots are swept
  at once (time-major sweep). The sweep has lower-bound pruning, an upper bound from the ±1 value-neighbour slots,
  and a cap on A rotation (irot 16). Far times whose B2A test provably fails are skipped, with deferred range updates.
  The exact t×j DP runs only for the top-3 slots. The total cost is recomputed exactly after each reinsertion
  (the DP is approximate around A2A queries).
- **Acceptance**: late acceptance (LAHC): accept if cost ≤ the cost h iterations ago. The history grows with the
  budget: h = 0.25·T (N > 300), 0.5·T (N ≤ 300), 2·T (N ≤ 150, T ≤ 90 s), plain descent for N ≤ 150 above 90 s.
- **Kick**: after K iterations without improvement, one destroy twice as large is accepted unconditionally (K = 100 for
  N > 150; 200, or 50 above 90 s, for N ≤ 150). The best solution is kept separately.
- **dualig**: the budget is split into 8 phases. Between phases the best solution is converted to P⁻¹ (reverse + invert +
  convert) and the search continues there, so the search alternates between the two equivalent problems.

### 2.4 Post-processing and small n
- **Peephole**: windows of the op sequence are replaced by shorter equivalent sequences found by meet-in-the-middle BFS
  over all 11 ops. This is the only place where sb/ss/rb-only sequences can appear.
- n ≤ 8: exact BFS. n ≤ 5: hard bounds (n=3 ≤ 2 ops, n=5 ≤ 8 ops).

### 2.5 Where the time goes (gprof, N=500, 30 s run)
| Function | Share |
|---|---:|
| beam_solve | 27% |
| solve_slots (reinsertion sweep) | 27% |
| chunk_construct | 22% |
| prep_b | 6% |
| lb_ll | 5% |

Measured worth of IG time at N=500: about 5 ops per second in the 30–45 s regime. A container change that made the
machine slower cost +36 ops at 45 s.

## 3. Results

### 3.1 Held-out evaluation (eval3+eval4, 60 cases per size, seeds never used in development)
| Version | Main idea | N=100 mean | N=500 mean |
|---|---|---:|---:|
| v4 | construction search | parity with reference (+0.6) | 2944 |
| v8 | fast reinsertion | 361.5 (median) | 2884 (median) |
| v13 | kick | 356.5 (median) | 2849.5 (median) |
| v16 | dual start + dualig | 353.7 | 2838.9 |
| v19 / v20 | late acceptance (+ kick 100) | **351.9** | **2817.5** |

Paired against the reference: N=100 −12.0 (53 W / 7 L), N=500 −150.0 (58 W / 2 L), both p < 0.0001.

### 3.2 Time–quality curve (the reference's own deep inputs; reference deep ≈ 8 h per case)
| N | Budget | Ours (mean) | vs reference deep |
|---:|---:|---:|---:|
| 100 | 15 s | 351.96 | +6.4 |
| 100 | 60 s | 347.94 | +2.4 |
| 100 | 300 s | **342.58** | **−3.0** (p .033) |
| 500 | 45 s | 2815.0 | +130.7 |
| 500 | 180 s | 2760.4 (v20 on a slower machine: 2761.8) | +76.2 |
| 500 | 600 s | 2719.4 | +35.2 |

## 4. Distance to the optimum

### 4.1 Exact optimum for small n (research/tools/optbfs.c)
A BFS from the sorted state gives the exact distance of every state (the op set is closed under inversion).
n=11 has 479 M states (one byte each) and builds in about 17 min on 4 vCPU. Statistics over all n! inputs:

| n | Mean optimum | Max optimum (diameter over inputs) | Solver vs optimum |
|---:|---:|---:|---|
| 3 | 1.1667 | 2 | identical on all inputs |
| 4 | 2.8750 | 5 | identical on all inputs |
| 5 | 4.8583 | 8 | identical on all inputs |
| 6 | 7.0597 | 10 | identical on all inputs |
| 7 | 9.3062 | 13 | identical on all inputs |
| 8 | 11.5904 | 16 | identical on all 40 320 inputs |
| 9 | 13.9110 | 19 | +0.105 (+0.75%), 200 random inputs |
| 10 | 16.2756 | 22 | +0.145 (+0.89%), 200 random inputs |
| 11 | 18.6816 | 25 | +0.205 (+1.10%), 200 random inputs |

For n = 9..11, the solver at 1 s and at the default 15 s gives the same mean: more time does not close the gap.
Every output was verified to end sorted.

### 4.2 Where the small-n gap comes from (optbfs `model` mode, IDA* with the exact table as heuristic)
Shortest solutions under two restrictions, on the same 200 inputs per n:
- **M1**: each value is pushed to B at most once. The optimum is unchanged in all 600 cases. Multiple round trips are
  never needed at n ≤ 11.
- **M2**: M1 and no sb/ss, which is the solver's representation before peephole. The optimum is lost in 33 / 50 / 62 of 200
  cases (mean +0.17 / +0.26 / +0.33).

| n | Cases where ours > optimum | …of which M2 > optimum | Cases where ours < M2 (peephole found a B swap) |
|---:|---:|---:|---:|
| 9 | 20 | **20** | 13 |
| 10 | 28 | **28** | 22 |
| 11 | 38 | **38** | 24 |

**Every remaining gap case at n = 9..11 is one where every optimal solution swaps B.** Inside the solver's own model
(M2) the solver is optimal in all cases but one, and often beats M2 thanks to the peephole.
This contradicts the earlier N=100 conclusion (log 69: "no evidence for B-side swaps"). That test only looked for local,
≤ 10-op windows inside a finished solution. The small-n evidence says B-side swaps matter, but they must be planned
inside the search, not patched in afterwards. The share of such cases grows with n (17% → 25% → 31%).

### 4.3 N=100 and N=500
- **Best known values** (deep inputs, mean of the per-case minimum over all our runs and the reference deep mode):
  N=100 **338.74** (the reference deep mode is strictly better than all our runs on only 13/50 cases), N=500 **2680.96** (the reference deep mode is best on
  39/50 cases).
- **Proven lower bounds**: only 152 / 1089. That is far too weak to say anything, so the true gap is unknown.
- **Estimate (not a measurement)**: at N=500 the gain per doubling of time fell from ≈ 27.5 (45 → 180 s) to ≈ 23.6
  (180 → 600 s). If that ratio (≈ 0.86 per doubling) continued, this search would level off around 2720 − 145 ≈ **2575**.
  That is about 4% below today's best known value, and it is only a guess about where this search would saturate, not a
  bound on the optimum. The small-n trend (+0.75% → +1.10% gap from n = 9 to 11, all from B swaps) suggests the
  optimum at N=500 is lower still by an unknown amount.

## 5. What did not work (all measured; log entry in brackets)
| Idea | Result |
|---|---|
| Single-move simulated annealing over the representation (QSA) | plateau ≈ 450 at N=100 (5) |
| LAHC with long histories at short budgets (first attempt) | worse at 10–30 s (8). Later rehabilitated with budget-scaled histories (74) |
| Racing / restarts from different K / top-2-3 starts | no gain or worse (15, 33) |
| Branch-free "vectorizable" sweeps | slower than the pruned sweep (32, 41) |
| IS+swap starts on random inputs | +7.4 (37). Useful only on structured inputs (67) |
| Larger peephole windows (K=8) at N=500 | the time is worth more in IG: +45 (39) |
| New chunk families (shifted / mirrored) | not worth their time (42) |
| BOrder-segment destroy, cost-biased destroy windows | worse (48, 51) |
| Parameter re-tunes (seglen, pqd, a2amax, dualig 4/16/32, kdir, beam width) | defaults are a local optimum (48, 53, 64, 71, 76, 77) |
| Multi-start with equal total time (N=100) | worse than one long run (48) |
| Probe-based start selection | neutral (68) |
| Conversion fix for dualig | no effect (71) |
| Automatic LAHC history (rate-based) | no gain over the fixed rule (75) |
| Keep-in-A bonus | much worse (75) |
| Moving time between construction and IG | no gain (77) |

## 6. Open problems / where to look next
1. **B-side swaps inside the search** (strongest evidence, section 4.2). The representation has A2A (sa) but no B2B.
   A query type "swap two adjacent values in B" (sb, merged with sa into ss when both happen), with a DP-reinsertion
   option analogous to A2A, would cover every remaining small-n gap case. Cheap first check: rerun the 4.2 table with the
   new query type and confirm that M2 then equals the optimum.
2. **Search stagnation at N=500**: the gain per time doubling is falling. Population or elite-based restarts, path relinking
   between elites, or larger structured destroys (whole chunks, whole BOrder blocks) are untested at long budgets.
3. **Speed**: construction plus beam still take about half of the profile. IG iterations translate almost directly into ops
   (≈ 5 ops/s). Candidates are incremental prep, SIMD on the slot sweep (the pruned version must stay), and parallel
   phases on 4 cores (the solver is single-threaded).
4. **Lower bounds**: nothing useful exists above n = 11. A relaxation (e.g. sum of per-value rotation lower bounds over a
   fixed K, or an LP over the query model) would give the first real optimality-gap numbers for N=100/500.
5. **Exact data beyond n = 11**: n = 12 needs 6.2 G states (≈ 6 GB). That is feasible with a bit-packed table on a larger
   machine and would extend the curve in 4.1.
6. **Structured inputs**: the reference still wins on one worst-case family (es: 1968 vs 1898). C++ IS+swap gives 1596 there,
   but in C it is not selected among the start candidates (67).

## 7. Reproducing
- Solver: `cd research/c && make && ./push_swap --time=15 <numbers>`.
- Benchmarks: `research/bench/` (CSV per run in `research/bench/runs/`), paired tests: `research/tools/paired.py`.
- Exact tables: `gcc -O2 -o optbfs research/tools/optbfs.c`, then `./optbfs 10 t10.bin build`, then
  `./optbfs 10 t10.bin stats | analyze | model`.
- Small-n data for 4.1/4.2: `research/bench/runs/optimum/` (200-case inputs per n as ranks 0..n−1, v20 outputs at
  1 s and 15 s as "input | ops" lines for `optbfs analyze`, `model` results, and the n=10/11 `stats` histograms).
