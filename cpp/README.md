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
- Supports a constrained Newton optimize mode (`sum a_j = 1`) with KKT solve and damping.

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

Run Newton optimization:
```bash
/Users/fcale/Dropbox/ChatGPT/Constants/tools/run_flint_objective.sh \
  --coeffs 1,0 \
  --nmax 256 \
  --k 32 \
  --prec-dps 120 \
  --print-dps 60 \
  --tail true \
  --optimize true \
  --max-it 20 \
  --tol 1e-24 \
  --constraint-tol 1e-24 \
  --damping true \
  --step-min-exp 20 \
  --log false
```

Run staged continuation:
```bash
python3 /Users/fcale/Dropbox/ChatGPT/Constants/tools/run_flint_bootstrap.py \
  --stage-list "1:256:32:120:20:1e-24:1e-24,2:256:32:120:20:1e-24:1e-24,4:256:32:120:20:1e-24:1e-24,8:256:32:120:20:1e-24:1e-24" \
  --max-retries 2 \
  --retry-prec-step 20 \
  --retry-max-it-step 8 \
  --records /Users/fcale/Dropbox/ChatGPT/Constants/checkpoints/flint_bootstrap_records.json
```

Notes:
- FLINT (with Arb) is required (`brew install flint`).
- Tail derivatives reuse precomputed basis/zeta tables inside a run.
- Regressions:
  - objective parity: `/Users/fcale/Dropbox/ChatGPT/Constants/tools/test_flint_objective.sh`
  - derivative parity: `/Users/fcale/Dropbox/ChatGPT/Constants/tools/test_flint_derivatives.sh`
  - optimize smoke: `/Users/fcale/Dropbox/ChatGPT/Constants/tools/test_flint_optimize_smoke.sh`
- runner reliability:
  - retry escalation: `/Users/fcale/Dropbox/ChatGPT/Constants/tools/test_flint_bootstrap_retry.sh`
- Next major step is persistent on-disk precompute tables and parallelized Hessian kernels for large `P`.
