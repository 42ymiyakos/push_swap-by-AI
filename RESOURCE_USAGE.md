# RESOURCE USAGE

Machine: 4 vCPU (Xeon @2.1GHz), 15 GB RAM. Peak RSS of our solver < 50 MB (N=500); reference ~350–415 MB.

| Item | CPU time (approx) |
|---|---|
| Reference build | ~3 min |
| Reference N=100, 20+10 cases (normal) | ~8 CPU-min |
| Reference N=500, 20 cases (normal) | ~15 CPU-min |
| QSA experiments | ~5 CPU-min |
| PSIG development experiments (Level 0–2) | ~3 CPU-hours |
| Final evaluation (held-out) | ~12 CPU-min ours + ~14 CPU-min reference |

Policy followed: ≤4 concurrent solver processes; no unbounded runs (all runs time- or iteration-capped).
