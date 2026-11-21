#!/bin/bash
# Level-4 evaluation v2 on fresh held-out seeds (eval2_*), never looked at during development.
cd "$(dirname "$0")/.."
python3 bench/bench.py --solver "bin/psig_v2 {args}" --set data/sets/eval2_u100.txt --jobs 4 --tag final2_psig_u100 --note "v2 default 10s"
python3 bench/bench.py --solver "bin/psig_v2 {args}" --set data/sets/eval2_u500.txt --jobs 4 --tag final2_psig_u500 --note "v2 default 40s"
python3 bench/bench.py --solver "bin/ref_push_swap {args}" --set data/sets/eval2_u100.txt --jobs 4 --tag final2_ref_u100 --timeout 1800
python3 bench/bench.py --solver "bin/ref_push_swap {args}" --set data/sets/eval2_u500.txt --jobs 4 --tag final2_ref_u500 --timeout 1800
