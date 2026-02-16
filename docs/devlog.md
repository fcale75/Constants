# Devlog (Constants project)

## Current architecture
- `wl/Constants.wl` — core setup + precision control + heartbeat logging; cached `BesselBlock[p,k]` implementing \mathcal{J}(p,k); mod-4 Kummer tail for Objective only (`Objective[a,Nmax,K]`) and residue-class series coefficients
- `wl/Newton.wl` — finite-sum objective/derivatives (`NewtonObjective`, `NewtonDerivatives`); constrained KKT Newton step `NewtonStep`
- `wl/TailDerivatives.wl` — mod-4 Kummer tail contributions for gradient/Hessian
- `wl/Driver.wl` — iteration harness `NewtonOptimize[P,Nmax; MaxIterations,Tolerance,Log]`

## Tests
Run:
```wl
math -file tests/runAll.wl
```
This loads all harnesses: `runTests`, `newtonTests`, `driverTests`, `tailTests`.

## Near-term todos
- Integrate tail derivative terms into the main driver path (fully tail-aware Newton loop).
- Improve precision strategy & stability diagnostics (condition numbers, damping, stopping criteria).
- Parameterize K, P, N schedule and add reproducible seed/logging for debugging runs.

## Note (2026-02-15)
We had a multi-hour period with no commits and no WL harness updates. Future workflow should rely on small, verifiable commits and CI/log artifacts rather than assumed background work.

## Update (2026-02-15)
Progress:
- WL harness is now executable end-to-end for a minimal N=8 run via `tests/runAll.wl` and `tools/iterate_fix.py`.
- `NewtonStep` now returns an association and `NewtonOptimize` returns an association with keys (including `aFinal`), enabling the harness to read `aFinal`.
- The harness summary now validates `Total[aFinal]≈1` and `Length[aFinal]=8`, and the Python controller includes a timeout to avoid idle hangs.
- Package loading is stabilized by explicit `Get["wl/*.wl"]` in the harness path.

Next goal:
- Verify the WL implementation matches the Rechnitzer paper exactly (objective functional, constraint(s), residuals, gradient, Hessian, and KKT system).
- Add a “paper fidelity” test that checks the scaled sequence `8^n a_n (n+1)(-1)^n` remains near constant (paper heuristic) and flags mismatches.
- Add instrumentation to log per-iteration diagnostics for debugging (norms, condition numbers/damping decisions, and fidelity metrics).
