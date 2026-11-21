#!/bin/bash
# mean initial-solution length over the first N cases of a set: init_mean.sh SET N "opts"
cd "$(dirname "$0")/.."
head -n $2 $1 | while read -r L; do bin/psig --iters=0 --verbose=1 $3 $L 2>&1 >/dev/null | grep -o 'init=[0-9]* valid=1 t=[0-9.]*'; done | awk '{split($1,a,"=");split($3,b,"="); s+=a[2]; t+=b[2]; n++} END {printf "mean init %.1f  mean t %.2f\n", s/n, t/n}'
