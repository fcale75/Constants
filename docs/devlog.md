# Devlog (Constants project)

## Current architecture (as of 2026-02-16)

### Wolfram Language code (`wl/`)

- `wl/Constants.wl`
  - Core definitions for the Rechnitzer ansatz and accelerated objective.
  - Implements `BesselBlock[p,k] = J_p(pi k/2) p! (4/(pi k))^p` with **precision-safe caching keyed by (p,k,prec)**.
  - Implements the mod-4 asymptotic series `SkSeriesCoefficients[a,r,K]` and the accelerated tail `TailObjective[a,N,K]`.
  - `Objective[a,N,K] = 1/2 + Sum_{k=1..N} Sk[a,k]^4 + TailObjective[a,N,K]`.

- `wl/Newton.wl`
  - Finite-sum (k <= N) objective + derivatives (`NewtonObjective`, `NewtonDerivatives`).
  - Provides a single constrained KKT Newton step (`NewtonStep`) for the truncated system.
  - This module is intentionally "no tail" (tail terms are added in the Driver).

- `wl/TailDerivatives.wl`
  - Tail contributions for gradient and Hessian:
    - `TailGradient[a,N,K]` for `Sum_{k>N} 4 J(l,k) Sk^3`
    - `TailHessian[a,N,K]` for `Sum_{k>N} 12 J(i,k) J(l,k) Sk^2`
  - Uses the same residue-class asymptotic series machinery as `TailObjective`.

- `wl/Driver.wl`
  - Main Newton optimization harness:
    - Uses the **KKT residual** `g(a) - lambda*1` (not `||g||`) as the stationarity measure.
    - Default is **paper-faithful**: includes tail objective + tail derivatives (`UseTail -> True`).
    - Includes a simple backtracking damping step for stability.
  - Entry point: `Constants`Driver`NewtonOptimize[P, N, opts]` (with option `K -> ...`).

### Test harness (`tests/`)

- `tests/runAll.wl` calls `tests/runHarness.wl`.
- `tests/runHarness.wl` runs:
  - `runTests.wl` (core constants),
  - `newtonTests.wl` (finite-sum derivatives),
  - `driverTests.wl` (driver integration),
  - `tailTests.wl` (tail derivative + asymptotic sanity).
- The harness writes a brief key=value summary to `logs/latest_summary.txt`.

### Python controller (`tools/iterate_fix.py`)

Runs `wolframscript -file tests/runAll.wl`, prints stdout/stderr, and prints the parsed summary.

## Key fixes made (2026-02-16)

- Fixed a major asymptotic scaling bug in `SkSeriesCoefficients`: for `z = pi k/2`, the Bessel asymptotic series uses powers of `(pi k)^{-n}` (not `(2/pi)^n`).
- Fixed precision poisoning: `BesselBlock[p,k]` caching is now keyed by working precision.
- Driver now uses the correct stationarity test for the constrained KKT system: `||g(a) - lambda||` plus the constraint residual.
- Test scripts no longer call `Exit[]` mid-run; the orchestrator (`runHarness.wl`) controls the process exit code.
- Fixed `tools/iterate_fix.py` syntax/structure issues.

## Near-term next steps

1. Add a "paper fidelity" regression:
   - check `u_n := 8^n (n+1) (-1)^n a_n` is approximately constant for the early coefficients.
2. Add objective/gradient/hessian finite-difference checks at a parameter set closer to the paper (e.g. P=8, N>=128, K>=16).
3. Add robust logging of line-search decisions, residual norms, and basic conditioning diagnostics in `Driver.wl`.
4. Start implementing the lower-bound dual construction (Section 4 of Rechnitzer).

