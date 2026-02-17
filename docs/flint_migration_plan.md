# FLINT Migration Plan (P=200 to P=500)

## Why migrate

Measured in the current Wolfram implementation:
- The Newton linear solve is negligible.
- The expensive part is repeated objective/gradient/Hessian assembly.
- Conditioning worsens rapidly as `P` grows, which amplifies precision loss.

The paper's implementation also identifies derivative assembly as the bottleneck and uses precomputation plus parallelism.

## What is now implemented

A first C++/FLINT executable:
- `cpp/src/constants_objective.cpp`
- Binary: `cpp/build/constants_objective`

It computes:
- finite objective: `1/2 + sum_{k=1..Nmax} Sk(k)^4`
- mod-4 asymptotic tail objective
- finite + tail gradient
- finite + tail Hessian
- total objective and `sqrt(objective)`
- constrained Newton optimization with damping
- staged continuation via `tools/run_flint_bootstrap.py`
  - includes retry-based escalation (precision/iterations/damping floor)

Utility scripts:
- `tools/build_flint_cpp.sh`
- `tools/run_flint_objective.sh`
- `tools/test_flint_objective.sh`
- `tools/test_flint_derivatives.sh`
- `tools/test_flint_optimize_smoke.sh`
- `tools/run_flint_bootstrap.py`
- `tools/test_flint_bootstrap_retry.sh`
- `tools/run_flint_calibrate_stage.py`
- `tools/test_flint_calibrate_stage.sh`

Current parity status:
- Matches known objective values for:
  - `P=1, a={1}`
  - `P=2, a={0.986,0.014}`
  to the printed precision.
- Matches WL objective/gradient/Hessian at a fixed `P=4` coefficient vector
  with max absolute errors around `1e-70` (objective/gradient) and `1e-69` (Hessian).
- Newton continuation reproduces the expected `P=1,2,4,8` trajectory in FLINT.

## Next implementation steps

1. Add persistent precompute tables (disk cache):
   - `B(p,k)` for `k <= Nmax`
   - asymptotic basis tensors by residue class mod 4
   - Hurwitz-zeta weights for `(s, n0 + r/4)`

2. Optimize derivative kernels for large `P`:
   - reduce temporary convolutions in tail Hessian
   - thread parallelism over `(i,l)` blocks
   - benchmark memory/compute tradeoffs for pairwise basis products

3. Harden optimization for large `P`:
   - improve conditioning strategy (projected coordinates or regularized KKT)
   - smarter precision scheduling across iterations/stages
   - failure recovery (retry with higher precision / smaller step floor)
   - objective stability checks across `(N,K,prec)` sweeps

4. Add benchmarking harness:
   - per-stage timing for finite/tail grad/hess
   - cache hit/miss counts
   - iteration history and residual norms

5. Run progression:
   - `P=8 -> 12 -> 16 -> ... -> 64` to calibrate precision schedule
   - then push to `P=101`, `P=160`, and `P=200`
   - reassess before attempting `P=500`

## Practical expectation

- `P=200` should be feasible with FLINT + caching + optimized derivative kernels.
- `P=500` is likely feasible only with careful precompute reuse and parallelized kernels; brute-force recomputation each iteration will be too slow.
