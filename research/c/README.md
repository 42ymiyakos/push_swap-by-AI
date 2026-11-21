# push_swap — C implementation (port of PSIG)

Plain C (gnu11), not norminette-compliant. Same algorithm and default settings as `research/src/psig.cpp`
(query representation, Iterated Greedy with exact DP reinsertion over all B slots, query-destroy,
chunk+beam initial solution with multi-start LIS, peephole post-optimization, exact BFS for n ≤ 8).
Optional experimental features of the C++ version (A2A/sa moves, restarts, two-phase, LAHC>1, ...) were not ported.

```
make            # -> ./push_swap
./push_swap 3 2 1
./push_swap --time=5 <numbers>      # default 10 s for n<=150, 20 s for n<=300, 40 s otherwise
./push_swap --iters=N --seed=S ...  # fixed iteration budget (deterministic)
make asan       # ASan/UBSan build
```

Equivalence: with `--iters` fixed, the output is byte-identical to the C++ version
(`tools/cmp_c_cpp.sh`, checked on n=5/10/50/100/500 and special distributions).
Files: `main.c` (CLI, flow), `rep.c` (representation/decoder/converter), `construct.c` + `beam.c`
(initial solution), `ig.c` (reinsertion DP, sweep, IG loop), `peep.c` (peephole), `exact.c` (BFS), `util.c`.
