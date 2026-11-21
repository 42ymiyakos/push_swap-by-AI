# Reference baseline: aaax8/push_swap @ 97d38feaec41522272b4558581cb65f322a5ed3f
Not vendored. Build (as done here):
```
git clone https://github.com/aaax8/push_swap && cd push_swap && git checkout 97d38fe
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS=-O2
cmake --build build --target push_swap checker
```
Normal settings (no options). Results: see BEST_RESULTS.md / bench/runs/ref_*.csv.
