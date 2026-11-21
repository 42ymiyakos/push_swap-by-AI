#!/bin/bash
# Level-4 evaluation on held-out seeds (eval_* sets never used for tuning).
cd "$(dirname "$0")/.."
python3 bench/bench.py --solver "bin/psig_final_eval {args}" --set data/sets/eval_u100.txt --limit 20 --jobs 4 --tag final_psig_eval_u100 --note "default time 10s"
python3 bench/bench.py --solver "bin/psig_final_eval {args}" --set data/sets/eval_u500.txt --limit 12 --jobs 4 --tag final_psig_eval_u500 --note "default time 40s"
python3 bench/bench.py --solver "bin/ref_push_swap {args}" --set data/sets/eval_u100.txt --limit 20 --jobs 4 --tag final_ref_eval_u100 --timeout 1800
python3 bench/bench.py --solver "bin/ref_push_swap {args}" --set data/sets/eval_u500.txt --limit 12 --jobs 4 --tag final_ref_eval_u500 --timeout 1800
