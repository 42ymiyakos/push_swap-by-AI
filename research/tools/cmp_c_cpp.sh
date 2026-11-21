#!/bin/bash
# Compare C port vs C++ PSIG outputs byte-for-byte with a fixed iteration budget.
# usage: cmp_c_cpp.sh SETFILE NCASES ITERS
set -u
cd "$(dirname "$0")/.."
same=0; diff=0
while read -r L; do
  a=$(bin/psig --iters=$3 $L | md5sum); b=$(c/push_swap --iters=$3 $L | md5sum)
  if [ "$a" == "$b" ]; then same=$((same+1)); else diff=$((diff+1)); fi
done < <(head -n $2 $1)
echo "$1 iters=$3: identical=$same different=$diff"
