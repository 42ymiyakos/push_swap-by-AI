#!/bin/bash
# parallel mean init: init_mean_par.sh SET OFFSET N "opts"
cd "$(dirname "$0")/.."
sed -n "$(( $2 + 1 )),$(( $2 + $3 ))p" $1 | xargs -P 4 -I{} sh -c "bin/psig --iters=0 --verbose=1 $4 {} 2>&1 >/dev/null | grep -o 'init=[0-9]* valid=1 t=[0-9.]*'" | awk '{split($1,a,"=");split($3,b,"="); s+=a[2]; t+=b[2]; n++} END {printf "mean init %.1f  mean t %.2f (n=%d)\n", s/n, t/n, n}'
