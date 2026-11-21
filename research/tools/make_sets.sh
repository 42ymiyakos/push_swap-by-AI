#!/bin/bash
# Build the canonical benchmark sets (deterministic). Train and eval seeds are disjoint.
set -e
D=$(dirname "$0")/../data/sets; mkdir -p "$D"; G="python3 $(dirname "$0")/gen.py"
for n in 3 5 10 50 100 500; do
  $G --n $n --dist uniform --seed 1000 --count 200 > $D/train_u$n.txt   # development
  $G --n $n --dist uniform --seed 500000 --count 200 > $D/eval_u$n.txt  # held-out final eval
done
for n in 100 500; do
  $G --n $n --dist uniform --seed 2000 --count 40 > $D/dev40_u$n.txt   # quick dev subset
  for d in sorted reverse nearly rotated rotated_reverse two_runs blocks_reverse sawtooth local_shuffle; do
    $G --n $n --dist $d --seed 7000 --count 3; done > $D/special_$n.txt
done
for n in 3 5 10; do for d in sorted reverse nearly rotated rotated_reverse; do $G --n $n --dist $d --seed 7000 --count 2; done > $D/special_$n.txt; done
