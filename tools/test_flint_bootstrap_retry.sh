#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/fcale/Dropbox/ChatGPT/Constants"
BOOT="${ROOT}/tools/run_flint_bootstrap.py"
TMP_OUT="$(mktemp)"
trap 'rm -f "$TMP_OUT"' EXIT

python3 "${BOOT}" \
  --p-list 2 \
  --nmax 256 \
  --k 32 \
  --prec-dps 100 \
  --print-dps 45 \
  --max-it 1 \
  --tol 1e-30 \
  --constraint-tol 1e-30 \
  --max-retries 2 \
  --retry-prec-step 10 \
  --retry-max-it-step 5 \
  --retry-step-min-exp-step 2 \
  --log false \
  --records /tmp/flint_bootstrap_retry_test.json \
  --no-save \
  > "${TMP_OUT}"

python3 - "${TMP_OUT}" <<'PY'
import re
import sys
from decimal import Decimal, getcontext

getcontext().prec = 120
text = open(sys.argv[1], "r", encoding="utf-8").read()

if "retry 1/2" not in text:
    raise SystemExit("Expected retry message not found")

m = re.search(
    r"P=2\s+N=\d+\s+K=\d+\s+prec=(\d+)\s+init=\S+\s+converged=(\w+)\s+kktSolveOk=(\w+)\s+iters=(\S+)\s+objective=([0-9.]+)",
    text,
)
if not m:
    raise SystemExit("Final stage summary line not found")

prec = int(m.group(1))
converged = m.group(2) == "true"
kkt_ok = m.group(3) == "true"
iters = m.group(4)
obj = Decimal(m.group(5))

if prec <= 100:
    raise SystemExit(f"Expected escalated precision, got {prec}")
if not (converged and kkt_ok):
    raise SystemExit("Stage did not converge after retry escalation")
if iters == "1":
    raise SystemExit("Expected more than one iteration after retries")

ref = Decimal("0.57463961900793418221607541986723768774586354806362979837544737291407986934717195")
if abs(obj - ref) > Decimal("1e-30"):
    raise SystemExit("Final objective is outside expected tolerance")

print("FLINT bootstrap retry test passed.")
print(f"final_prec={prec} final_iters={iters} objective={obj}")
PY
