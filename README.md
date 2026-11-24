# push_swap-by-AI

An autonomous research project to build the strongest possible push_swap solver: as few
operations as possible, always correct.

## Results (held-out inputs eval3+eval4, 60 cases per size, mean operation count; human implementation: eval2 median)

| N | This repository (C version, v19; v20 = current default) | Reference [aaax8/push_swap](https://github.com/aaax8/push_swap) | Human implementation (42 validated) |
|---:|---:|---:|---:|
| 100 | **351.9** (15 s) | 363.9 (about 17 s) | 608 |
| 500 | **2817.5** (45 s) | 2967.5 (about 48 s) | 5349 |

- Both sizes are significantly better than the reference normal mode: N=100 −12.0 (53 W / 7 L), N=500 −150.0 (58 W / 2 L), p < 0.0001.
- Against the reference deep mode (~8 h per case): N=100 better at 300 s (−3.0, p .033); N=500 still +35 at 600 s.
- v20 (current default) differs from v19 only for N > 150. It was measured on a slower machine, where it is −0.73 vs v19 at N=500 (not significant).
- Exact optimum on every input for n ≤ 8 (exhaustive check); on 200-input samples for n = 9–11, +0.75–1.1% above the exact optimum, every gap case needing B-side swaps.
- Structured "es" input (N=500): 1545 vs the reference's 1898.
- **Overview of the design, the distance to the optimum, failed ideas and open problems: [REPORT.md](REPORT.md).**
- Full figures and per-case data: [BEST_RESULTS.md](BEST_RESULTS.md) and `research/bench/runs/`.

## Quick start (C version)

```sh
cd research/c && make
./push_swap 3 2 1
./push_swap --time=5 $(python3 ../tools/gen.py --n 100 --seed 1)
```

Default time limit: 15 s for n ≤ 150, 25 s for n ≤ 300, 45 s above that. `--iters=N --seed=S` gives a
deterministic run. Invalid input prints `Error`. This is not a norminette-compliant 42 submission.

## Algorithm summary

1. **Representation**: a solution is not stored as the command list itself. It is stored as which values stay
   in A (K), the order in which each value moves (query list), and its cyclic position inside each stack
   (AOrder/BOrder keys). Any such combination always decodes to a correct sort.
2. **Initial solution**: LIS candidates from many start positions, chunked A→B, then beam search for B→A
   (with an arc-cover estimate of the remaining rotations). The same construction is also run on the inverse
   permutation P⁻¹; reversing that op sequence and inverting each op (pa↔pb, ra↔rra, …) also sorts P, and the shorter
   one is kept ("dual" start, idea from keusidan/42-push_swap-solver). This mostly removes bad starting solutions.
   For structured inputs where "increasing after swapping adjacent pairs" keeps far more values in A than the LIS
   (e.g. descending blocks), IS+swap starts are added as an extra candidate group.
3. **Iterated Greedy** (the time is split into 8 phases; between phases the best solution is converted to the
   inverse permutation and back, so the search alternates between the two equivalent problems): remove a few values, then reinsert each one at the jointly optimal choice of
   "when it goes to B, when it comes back, and where it sits in B and in A". The cost difference of each
   candidate is computed exactly, which turns the choice into an O(Q) DP. All B slots are evaluated at once
   with a time-major sweep plus lower-bound pruning.
4. **Post-processing**: replace short windows by shorter equivalent sequences found with meet-in-the-middle
   BFS (peephole).

Details: [REPORT.md](REPORT.md) (current design), [ALGORITHM_NOTES.md](ALGORITHM_NOTES.md) (derivations)

## Layout

| Path | Contents |
|---|---|
| `research/c/` | C implementation (main deliverable) |
| `research/src/` | C++ research version (psig = same algorithm plus experimental options, qsa = simulated annealing prototype) |
| `research/tools/` | Strict checker, input generator, op-sequence analyzer, C/C++ equivalence check |
| `research/bench/` | Benchmark harness and all measurement results |
| `research/data/sets/` | Benchmark inputs (development / held-out sets kept separate) |
| `*.md` (root) | Research documents: plan, log, experiments, algorithm notes, results, next steps, resource usage |
