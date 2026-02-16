# FLINT Backend

This directory contains the C++/FLINT backend for the Constants project.

Current executable:
- `constants_objective`

What it does today:
- Parses a coefficient vector `a_0,...,a_{P-1}`.
- Computes
  - finite objective `1/2 + sum_{k=1..Nmax} Sk(k)^4`
  - optional mod-4 asymptotic tail objective (Kummer-style split)
  - optional gradient and Hessian (finite + tail).
- Reports objective and `sqrt(objective)` as a `nu2` candidate.

Build:
```bash
/Users/fcale/Dropbox/ChatGPT/Constants/tools/build_flint_cpp.sh
```

Run:
```bash
/Users/fcale/Dropbox/ChatGPT/Constants/tools/run_flint_objective.sh \
  --coeffs 1 \
  --nmax 256 \
  --k 32 \
  --prec-dps 120 \
  --print-dps 60 \
  --tail true \
  --derivatives true \
  --hessian true
```

Notes:
- FLINT (with Arb) is required (`brew install flint`).
- Tail derivatives reuse precomputed basis/zeta tables inside a run.
- Regressions:
  - objective parity: `/Users/fcale/Dropbox/ChatGPT/Constants/tools/test_flint_objective.sh`
  - derivative parity: `/Users/fcale/Dropbox/ChatGPT/Constants/tools/test_flint_derivatives.sh`
- Next major step is Newton/continuation and persistent on-disk precompute tables.
