#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/fcale/Dropbox/ChatGPT/Constants"
CAL="${ROOT}/tools/run_flint_calibrate_stage.py"
OUT="$(mktemp /tmp/flint_calibrate_XXXXXX.json)"
TMP_STDOUT="$(mktemp /tmp/flint_calibrate_stdout_XXXXXX.txt)"
trap 'rm -f "$OUT" "$TMP_STDOUT"' EXIT

python3 "${CAL}" \
  --p 8 \
  --n-list 256,384 \
  --k-list 32 \
  --prec-list 120 \
  --max-it 12 \
  --tol 1e-24 \
  --constraint-tol 1e-24 \
  --stability-threshold 1e-6 \
  --require-consecutive 1 \
  --out "$OUT" \
  > "$TMP_STDOUT"

python3 - "$OUT" <<'PY'
import json
import sys
from decimal import Decimal

path = sys.argv[1]
data = json.load(open(path, "r", encoding="utf-8"))
rows = data.get("rows", [])
if len(rows) < 2:
    raise SystemExit("Expected at least 2 calibration rows")

for row in rows:
    if row.get("P") != 8:
        raise SystemExit("Unexpected P in calibration rows")
    if not row.get("converged", False):
        raise SystemExit("Non-converged row in calibration test")
    if not row.get("kktSolveOk", False):
        raise SystemExit("kktSolveOk false in calibration test")

if data.get("recommended") is None:
    raise SystemExit("Expected a recommended setting")

obj0 = Decimal(rows[0]["objective"])
obj1 = Decimal(rows[1]["objective"])
if abs(obj1 - obj0) > Decimal("1e-2"):
    raise SystemExit("Calibration produced implausibly large objective jump")

print("FLINT calibration stage test passed.")
print(f"rows={len(rows)} recommended={data['recommended']}")
PY
