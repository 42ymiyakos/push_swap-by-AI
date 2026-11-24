# push_swap solver research — consolidated report (v20, 2026-11-24)

This document summarizes where the project stands: why it was started, how the research was run (an AI-driven
research loop with a human making the strategic decisions), how the solver works, how good it is, how far it is from
the optimum (as far as that can be measured), what did not work, and which directions are still open.
Chronological details: [RESEARCH_LOG.md](RESEARCH_LOG.md) (entries 1–78). All numbers: [BEST_RESULTS.md](BEST_RESULTS.md).

## 1. Background and purpose

### 1.1 The problem
push_swap is a sorting exercise from the 42 school curriculum. The input is a list of distinct integers on **stack A**;
**stack B** starts empty. The program must print a sequence of operations that leaves A sorted in ascending order (smallest
on top) and B empty. Only 11 operations exist:

| Operation | Effect |
|---|---|
| `sa`, `sb`, `ss` | swap the top two elements of A, of B, or of both |
| `pa`, `pb` | move the top of B onto A / the top of A onto B |
| `ra`, `rb`, `rr` | rotate A, B or both upwards (the top element goes to the bottom) |
| `rra`, `rrb`, `rrr` | rotate A, B or both downwards (the bottom element goes to the top) |

The score is the **number of operations**; running time is not part of the score. In the curriculum, a project passes
the standard benchmark with fewer than 700 operations for 100 random numbers and fewer than 5500 for 500. Usual student
solutions are hand-designed greedy, chunk or radix strategies.

Mathematically it is a shortest-path problem: each state is a pair of stacks, and each operation is an edge.
The state space grows like N!·(N+1), so exact search is only possible for very small N (here: N ≤ 11, section 6.1). For
N = 100 and 500 the problem becomes a large combinatorial optimization problem, with no known way to compute or certify
the optimum.

### 1.2 Starting point
| | N=100 | N=500 |
|---|---:|---:|
| Human implementation (the repository owner's own push_swap, validated at 42; adaptive mode, 200 random inputs, median) | 848 | 6784 |
| Same implementation, best fixed strategy (`--medium`) | 582.5 | 5227.5 |
| Strongest public reference found, aaax8/push_swap (iterated greedy + LAHC, normal mode, median, reproduced here) | 362.5 | 2959.5 |

The human implementation meets the curriculum thresholds. The public reference shows that search-based methods do far
better than hand-designed strategies.

### 1.3 Why minimize the operation count further
- The pass thresholds measure whether an exercise is done, not how good an answer is. Once you ask "how many operations
  are actually necessary?", push_swap becomes a clean, well-defined optimization problem with an easy-to-verify objective.
  It has a small, closed set of moves, and an exact optimum is computable for small N. That makes it a good testbed for
  search algorithms.
- The project set out to build the strongest solver it could. That meant beating the best public reference under
  comparable time, measuring progress with statistically sound held-out tests, and documenting what works and what does
  not, so that later push_swap researchers can start from here instead of from scratch.
- Constraints fixed at the start: always correct output (checked by an independent strict checker); the human baseline is
  never modified; no optimality claims without proof; development and evaluation inputs are kept strictly separate.

## 2. How the research was done: an AI-run research loop

### 2.1 Roles
The work was done by Claude (an AI model, running as Claude Code in a cloud container with 4 vCPU) acting as an
autonomous research lab, with a human supervisor (the repository owner) setting the direction.

| Claude | Human |
|---|---|
| Surveyed prior work (built the reference and read its report) and wrote the research plan | Set the goal ("the strongest push_swap solver") and the ground rules (no changes to the baseline, no force-push, no unproven optimality claims, separate train/eval data) |
| Proposed hypotheses (RESEARCH_PLAN.md H1–H5 and every log entry) | Granted access when needed (cloning the human baseline, log 21; handing over the repository, after log 22) |
| Implemented everything: C++ research version, C port, tools (checker, generators, exact BFS, statistics) | Pointed to a friend's solver as a source of hints only, with the instruction not to follow it too closely. Two ideas came from its documentation (dual start, log 63; worst-case generator, log 67); its code was never run |
| Designed and ran the experiments, including scheduling 4 parallel jobs and multi-hour background runs | Encouraged continuing when progress looked finished ("まだいけるぞ" — there is still more to gain) |
| Did the statistics (paired tests), made the adopt/reject decisions per version, and kept the docs and commit history | Made the strategic call at the plateau: stop parameter tuning, measure the distance to the optimum, and write this report so that others can continue |

In numbers: 4 days (2026-11-21 to 11-24), 105 commits, 78 log entries, 20 solver versions, about 300 benchmark result
files. The human wrote no solver code and chose no parameters. The human's input was a handful of short messages at
decision points.

### 2.2 How Claude worked: the project in phases
**Day 1: from an empty repository to parity with the reference.**
Claude started with an empty repository. It surveyed the machine, built the public reference and reproduced its published
numbers before trusting them, and read the reference's design report.
- The first own design, simulated annealing with single moves (QSA), stalled at about 450 operations for N=100. Claude
  traced this to single moves being unable to change "when a value goes to B" and "where it sits in B" together. It
  replaced them with a reinsertion that chooses all four decisions jointly.
- To make that affordable, Claude derived that the cost change separates into two terms, A_b(i) + B_{b,t}(j). That turns
  the choice into an O(Q) dynamic program.
- It then built a check mode that compares every predicted delta with a full re-evaluation. That check found a real bug
  (key collisions), which was fixed until the deltas matched exactly on thousands of reinsertions.
- That evening came the first honest verdict: "competitive but not better" than the reference. The main weakness was
  diagnosed as run-to-run variance at N=500.

**Days 1–2: diagnosis-driven improvement.**
Instead of guessing parameters, Claude measured where the losses came from.
- Plotting initial versus final cost showed an almost monotone relation at N=500. That led to the construction search
  (chunk counts, two boundary families, staged beams), which brought a significant held-out win (log 31).
- Profiling showed that two DP solves took half the search time. A series of exactness-preserving speed-ups followed:
  shared preparation, an upper bound from neighbour slots, a cap on A rotations, and lower-bound skipping of 82% of the
  sweep steps. Each one was verified to produce byte-identical output to the C++ version. Over v6–v12 these came to more
  than 90 operations ahead of the reference at N=500.
- Two clever-looking ideas were slower in practice: SIMD-friendly unpruned sweeps (log 32, 41). They were measured,
  reverted and documented.

**Days 2–3: diversification and long budgets.**
- Kick (v13) made the N=100 gain statistically significant for the first time.
- The dual start and alternating search on P⁻¹ (v15/v16) came from one paragraph in a friend's documentation. They were
  re-implemented from the idea alone.
- Claude then compared against the reference's 8-hour deep mode on its own inputs. That exposed stagnation at long
  budgets, which led back to late acceptance, previously rejected, now with a budget-scaled history (v19).
- When a container restart made the machine slower, Claude noticed that the same binary lost 36 operations. It
  discarded every cross-container comparison and reran the baselines.

**Day 4: plateau and the optimum study.**
- When every remaining parameter tested neutral, the human proposed measuring the distance to the optimum.
- Claude built exact distance tables for every state up to n = 11 (479 million states).
- To find why the remaining gap exists, it wrote a restricted-model IDA* search. The search showed that every gap case
  needs B-side swaps. This is a concrete, testable lead for the next researcher, and it corrected one of Claude's own
  earlier conclusions.

**Throughout:**
- Hours-long benchmark jobs ran in the background, detached so they would survive tool time limits.
- Never more than four solver processes ran at a time, so measurements were not disturbed.
- Every negative result went into the log. About 25 ideas failed, and section 7 lists them.
- Claims were downgraded whenever the held-out data did not support them.

### 2.3 The loop
Each idea went through the same cycle, recorded in RESEARCH_LOG.md:
1. **Hypothesis** from a diagnosis. Example: "the N=500 final result tracks the initial solution, so better starts give
   better finals" (log 23 → v3/v4).
2. **Implementation** behind a command-line option, default off. An option keeps the old behaviour reproducible.
3. **Quick check** on a few cases, then a **development test** on 20–60 fixed dev cases. The comparison is always paired:
   both versions run the same inputs in the same container.
4. **Statistics**: mean paired difference, wins/ties/losses, and a paired sign-flip permutation test
   (`research/tools/paired.py`).
5. **Adoption** only after a **held-out test** on evaluation inputs never used during development (eval3+eval4, 60 cases
   per size). Adopted versions were promoted to `main`.
6. **Documentation** of every result, including failures (section 7), so later work does not repeat them.

### 2.4 Safeguards against fooling ourselves
- **Data contamination was caught and fixed.** At log 18 an evaluation case was used to diagnose a bug. That set was
  then treated as contaminated, and fresh evaluation seeds were generated (eval2, later eval3/eval4).
- **Noise was measured, not assumed.** The same configuration was rerun to size the run-to-run noise (about ±8 at
  N=500, 45 s, log 71). Some dev wins did not reproduce on held-out data and were reported as such (log 72: v18 kept
  only as "neutral or better").
- **Machine effects**: after a container restart the same binary was 36 operations worse (log 75). From then on only
  comparisons within one container were used, and the baselines were rerun.
- **Exactness checks**: a delta-check mode verified the DP reinsertion against full re-evaluation (it found a real bug,
  log 6). The C and C++ versions are byte-identical for a fixed number of iterations (checked after every change). Every
  output is replayed by an independent strict checker. Sanitizers and Valgrind are clean.
- **No unproven claims**: optimality is stated only where exhaustive search proves it (n ≤ 8). Estimates are labelled
  as such.

### 2.5 Self-correction
Several conclusions were later overturned by better evidence. The log keeps both the old and the new result:
- Late acceptance was rejected early (log 8, short budgets, long histories). Revisited with budget-scaled histories, it
  became the largest late-stage gain (log 74: −21 at N=500).
- "B-side swaps are not needed" (log 69, based on windowed re-solves at N=100) was overturned by the exact small-n study
  (log 78): every remaining small-n gap needs them.
- A dev-set gain for kick 50 at long budgets did not transfer to held-out data (log 72); the claim was downgraded.

### 2.6 What this setup was good and bad at
- **Good at**: running many careful, controlled experiments quickly; keeping negative results; implementing intricate exact
  machinery (the O(Q) reinsertion DP, the exact BFS tables) and checking it automatically.
- **Limited by**: compute (one 4-vCPU container; N=500 runs take 45 s to 10 min per case, so 60-case tests take hours) and
  session constraints (some actions needed explicit human permission; long jobs had to be run detached).
- **Where the human mattered most**: deciding what was worth doing. That covered granting access, pointing to outside
  ideas while limiting how far to follow them, refusing to stop early, and finally recognizing the plateau and turning
  the work into a write-up for future researchers.

## 3. Summary

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
  (+0.75% / +0.89% / +1.10%). Every one of the gap cases needs a swap of B (sb/ss) in all of its optimal solutions (section 6.2).
- N=100: better than the reference deep mode at 300 s (−3.0, p .033). N=500: still +35 behind the deep mode at 600 s and
  +38 behind the best known values; the gap shrinks by about 25 per doubling of time.
- **Optimality is not claimed anywhere above n = 8.** For N ≥ 12 the only proven lower bound is the very weak
  information-theoretic one; the true distance to the optimum at N=100/500 is unknown.

## 4. Algorithm design (research/c, ~2400 lines of C)

### 4.1 Representation (independent reproduction of the reference's idea)
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

### 4.2 Construction (≈ 28% of the time at N=500)
1. Candidate K sets: LIS read cyclically from 50 (N ≤ 300) / 30 (N > 300) start positions, deduplicated. For
   structured inputs, "increasing after swapping adjacent pairs" (IS+swap) candidates are added as a separate group.
2. A→B: chunk passes, middle-out, with delayed rb merged into rr. The chunk count is searched per instance
   ({2,3,4} / {3..6} / {4,6,8,10,12}), with two boundary families (equal value ranges, equal moved counts).
3. B→A: beam search with an arc-cover heuristic (the minimal arc of A that must still be visited, plus reaching the minimum).
   The stages are a greedy screen of every (K, m, family), then small beams on the top 3 per (m, family), then full beams
   (1000 / 300 / 500 wide) on the overall top 3.
4. **Dual start**: everything above also runs on P⁻¹. Its op sequence, reversed and with every op inverted, sorts P.
   The shorter result is kept. This removes most bad starts (N=500 worst case 3510 → 3408).

### 4.3 Improvement: Iterated Greedy with exact DP reinsertion (≈ 70% of the time)
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

### 4.4 Post-processing and small n
- **Peephole**: windows of the op sequence are replaced by shorter equivalent sequences found by meet-in-the-middle BFS
  over all 11 ops. This is the only place where sb/ss/rb-only sequences can appear.
- n ≤ 8: exact BFS. n ≤ 5: hard bounds (n=3 ≤ 2 ops, n=5 ≤ 8 ops).

### 4.5 Where the time goes (gprof, N=500, 30 s run)
| Function | Share |
|---|---:|
| beam_solve | 27% |
| solve_slots (reinsertion sweep) | 27% |
| chunk_construct | 22% |
| prep_b | 6% |
| lb_ll | 5% |

Measured worth of IG time at N=500: about 5 ops per second in the 30–45 s regime. A container change that made the
machine slower cost +36 ops at 45 s.

## 5. Results

### 5.1 Held-out evaluation (eval3+eval4, 60 cases per size, seeds never used in development)
| Version | Main idea | N=100 mean | N=500 mean |
|---|---|---:|---:|
| v4 | construction search | parity with reference (+0.6) | 2944 |
| v8 | fast reinsertion | 361.5 (median) | 2884 (median) |
| v13 | kick | 356.5 (median) | 2849.5 (median) |
| v16 | dual start + dualig | 353.7 | 2838.9 |
| v19 / v20 | late acceptance (+ kick 100) | **351.9** | **2817.5** |

Paired against the reference: N=100 −12.0 (53 W / 7 L), N=500 −150.0 (58 W / 2 L), both p < 0.0001.

### 5.2 Time–quality curve (the reference's own deep inputs; reference deep ≈ 8 h per case)
| N | Budget | Ours (mean) | vs reference deep |
|---:|---:|---:|---:|
| 100 | 15 s | 351.96 | +6.4 |
| 100 | 60 s | 347.94 | +2.4 |
| 100 | 300 s | **342.58** | **−3.0** (p .033) |
| 500 | 45 s | 2815.0 | +130.7 |
| 500 | 180 s | 2760.4 (v20 on a slower machine: 2761.8) | +76.2 |
| 500 | 600 s | 2719.4 | +35.2 |

## 6. Distance to the optimum

### 6.1 Exact optimum for small n (research/tools/optbfs.c)
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

### 6.2 Where the small-n gap comes from (optbfs `model` mode, IDA* with the exact table as heuristic)
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

### 6.3 N=100 and N=500
- **Best known values** (deep inputs, mean of the per-case minimum over all our runs and the reference deep mode):
  N=100 **338.74** (the reference deep mode is strictly better than all our runs on only 13/50 cases), N=500 **2680.96** (the reference deep mode is best on
  39/50 cases).
- **Proven lower bounds**: only 152 / 1089. That is far too weak to say anything, so the true gap is unknown.
- **Estimate (not a measurement)**: at N=500 the gain per doubling of time fell from ≈ 27.5 (45 → 180 s) to ≈ 23.6
  (180 → 600 s). If that ratio (≈ 0.86 per doubling) continued, this search would level off around 2720 − 145 ≈ **2575**.
  That is about 4% below today's best known value, and it is only a guess about where this search would saturate, not a
  bound on the optimum. The small-n trend (+0.75% → +1.10% gap from n = 9 to 11, all from B swaps) suggests the
  optimum at N=500 is lower still by an unknown amount.

## 7. What did not work (all measured; log entry in brackets)
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

## 8. Open problems / where to look next
1. **B-side swaps inside the search** (strongest evidence, section 6.2). The representation has A2A (sa) but no B2B.
   A query type "swap two adjacent values in B" (sb, merged with sa into ss when both happen), with a DP-reinsertion
   option analogous to A2A, would cover every remaining small-n gap case. Cheap first check: rerun the 6.2 table with the
   new query type and confirm that M2 then equals the optimum.
2. **Search stagnation at N=500**: the gain per time doubling is falling. Population or elite-based restarts, path relinking
   between elites, or larger structured destroys (whole chunks, whole BOrder blocks) are untested at long budgets.
3. **Speed**: construction plus beam still take about half of the profile. IG iterations translate almost directly into ops
   (≈ 5 ops/s). Candidates are incremental prep, SIMD on the slot sweep (the pruned version must stay), and parallel
   phases on 4 cores (the solver is single-threaded).
4. **Lower bounds**: nothing useful exists above n = 11. A relaxation (e.g. sum of per-value rotation lower bounds over a
   fixed K, or an LP over the query model) would give the first real optimality-gap numbers for N=100/500.
5. **Exact data beyond n = 11**: n = 12 needs 6.2 G states (≈ 6 GB). That is feasible with a bit-packed table on a larger
   machine and would extend the curve in 6.1.
6. **Structured inputs**: the reference still wins on one worst-case family (es: 1968 vs 1898). C++ IS+swap gives 1596 there,
   but in C it is not selected among the start candidates (67).

## 9. Reproducing
- Solver: `cd research/c && make && ./push_swap --time=15 <numbers>`.
- Benchmarks: `research/bench/` (CSV per run in `research/bench/runs/`), paired tests: `research/tools/paired.py`.
- Exact tables: `gcc -O2 -o optbfs research/tools/optbfs.c`, then `./optbfs 10 t10.bin build`, then
  `./optbfs 10 t10.bin stats | analyze | model`.
- Small-n data for 6.1/6.2: `research/bench/runs/optimum/` (200-case inputs per n as ranks 0..n−1, v20 outputs at
  1 s and 15 s as "input | ops" lines for `optbfs analyze`, `model` results, and the n=10/11 `stats` histograms).
