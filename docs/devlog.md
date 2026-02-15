# Devlog (Constants project)

## Current architecture
- `wl/Constants.wl`
  - core setup + precision control + heartbeat logging
  - cached `BesselBlock[p,k]` implementing \mathcal{J}(p,k)
  - mod-4 Kummer tail for Objective only (`Objective[a,Nmax,K]`) and residue-class series coefficients
- `wl/Newton.wl`
  - finite-sum objective/derivatives (`NewtonObjective`, `NewtonDerivatives`)
  - constrained KKT Newton step `NewtonStep`
- `wl/TailDerivatives.wl`
  - mod-4 Kummer tail contributions for gradient/Hessian
- `wl/Driver.wl`
  - iteration harness `NewtonOptimize[P,Nmax; MaxIterations,Tolerance,Log]`

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

