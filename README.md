# push_swap-by-AI

An autonomous research project to build the strongest possible push_swap solver: as few
operations as possible, always correct.

## Results (held-out inputs eval3+eval4, 60 cases per size, median operation count; human implementation measured on eval2)

| N | This repository (C version) | Reference [aaax8/push_swap](https://github.com/aaax8/push_swap) | Human implementation (42 validated) |
|---:|---:|---:|---:|
| 100 | **354.5** (about 15 s) | 363.5 (about 17 s) | 608 |
| 500 | **2839.5** (about 44 s) | 2967.5 (about 48 s) | 5349 |

- n ≤ 8 is solved optimally by exhaustive search (n=3: at most 2 ops, n=5: at most 8 ops).
- Both sizes are significantly better than the reference: N=500 by 120.4 operations on average (p < 0.0001, 57 wins out of 60, worst case 3052), N=100 by 8.7 operations (p < 0.0001, 51 wins, 1 tie, 8 losses).
- Full figures, conditions and per-case data: [BEST_RESULTS.md](BEST_RESULTS.md) and `research/bench/runs/`.

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
3. **Iterated Greedy**: remove a few values, then reinsert each one at the jointly optimal choice of
   "when it goes to B, when it comes back, and where it sits in B and in A". The cost difference of each
   candidate is computed exactly, which turns the choice into an O(Q) DP. All B slots are evaluated at once
   with a time-major sweep plus lower-bound pruning.
4. **Post-processing**: replace short windows by shorter equivalent sequences found with meet-in-the-middle
   BFS (peephole).

Details: [ALGORITHM_NOTES.md](ALGORITHM_NOTES.md)

## Layout

| Path | Contents |
|---|---|
| `research/c/` | C implementation (main deliverable) |
| `research/src/` | C++ research version (psig = same algorithm plus experimental options, qsa = simulated annealing prototype) |
| `research/tools/` | Strict checker, input generator, op-sequence analyzer, C/C++ equivalence check |
| `research/bench/` | Benchmark harness and all measurement results |
| `research/data/sets/` | Benchmark inputs (development / held-out sets kept separate) |
| `*.md` (root) | Research documents: plan, log, experiments, algorithm notes, results, next steps, resource usage |
