# push_swap research lab

Autonomous research project: build the strongest push_swap solver.
Top-level research documents live in the repository root (RESEARCH_PLAN.md, RESEARCH_LOG.md,
EXPERIMENTS.md, ALGORITHM_NOTES.md, NEXT_STEPS.md, BEST_RESULTS.md, RESOURCE_USAGE.md).

## Layout
- `c/` — **C implementation** of the solver (see c/README.md).
- `src/` — solver sources (C++17, header-based): `psig.cpp` (best solver), `qsa.cpp` (SA prototype),
  `rep.hpp` (representation/decoder), `ig.hpp` (exact DP reinsertion), `construct.hpp`/`beam.hpp`
  (initial solutions), `exact.hpp` (BFS optimum for n<=8), `peephole.hpp` (window BFS post-optimizer).
- `tools/` — strict checker (C), deterministic generator, set builder, op-sequence analyzer.
- `bench/` — benchmark harness (`bench.py`), `runs/` (per-case CSV + `summary.jsonl`), `final_eval.sh`.
- `data/sets/` — benchmark inputs. `train_*`/`dev40_*` used for development, `eval_*` held out.
- `baselines/` — notes on reference / human baselines. `best/` — pointer to the promoted solver.

## Build & run
```
make                   # bin/psig, bin/qsa, tools/checker
bin/psig 3 2 1         # prints operations
bin/psig --time=5 $(python3 tools/gen.py --n 100 --seed 1) | tools/checker $(python3 tools/gen.py --n 100 --seed 1)
```
Options: `--time=S` (default 10 s for n<=150, 40 s for n>300), `--seed=`, `--peep=4,6` (peephole windows),
`--a2a=1` (sa moves), see `src/psig.cpp` Params for the full list.

## Benchmark
```
python3 bench/bench.py --solver "bin/psig --time=10 {args}" --set data/sets/dev40_u100.txt --tag mytag --jobs 4
```
