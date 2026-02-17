#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/fcale/Dropbox/ChatGPT/Constants"
BOOT="${ROOT}/tools/run_flint_bootstrap.py"
TMP_OUT="$(mktemp)"
trap 'rm -f "$TMP_OUT"' EXIT

python3 "${BOOT}" \
  --p-list 1,2,4,8 \
  --nmax 256 \
  --k 32 \
  --prec-dps 120 \
  --print-dps 70 \
  --max-it 20 \
  --tol 1e-24 \
  --constraint-tol 1e-24 \
  --log false \
  --records /tmp/flint_bootstrap_smoke.json \
  --no-save \
  > "${TMP_OUT}"

python3 - "${TMP_OUT}" <<'PY'
import re
import sys
from decimal import Decimal, getcontext

getcontext().prec = 120
text = open(sys.argv[1], "r", encoding="utf-8").read()

line_re = re.compile(
    r"P=(\d+)\s+N=\d+\s+K=\d+\s+prec=\d+\s+init=\S+\s+converged=(\w+)\s+kktSolveOk=(\w+)\s+iters=\S+\s+objective=([0-9.]+)"
)
rows = []
for m in line_re.finditer(text):
    p = int(m.group(1))
    conv = m.group(2) == "true"
    kkt_ok = m.group(3) == "true"
    obj = Decimal(m.group(4))
    rows.append((p, conv, kkt_ok, obj))

if [p for p, _, _, _ in rows] != [1, 2, 4, 8]:
    raise SystemExit(f"Unexpected stage rows: {rows}")

for p, conv, kkt_ok, _ in rows:
    if not (conv and kkt_ok):
        raise SystemExit(f"Stage P={p} did not converge cleanly")

obj = {p: v for p, _, _, v in rows}
if not (obj[1] > obj[2] > obj[4] >= obj[8]):
    raise SystemExit(f"Objective monotonicity failed: {obj}")

ref_p1 = Decimal("0.574694861894825706792169188192308052244615365442183425654224244188404730797633575601903795")
ref_p2 = Decimal("0.574639619007934182216075419867237687745863548063629798375447372914079869347171954272272976578")

if abs(obj[1] - ref_p1) > Decimal("1e-40"):
    raise SystemExit("P=1 objective mismatch")
if abs(obj[2] - ref_p2) > Decimal("1e-35"):
    raise SystemExit("P=2 objective mismatch")

print("FLINT optimize smoke passed.")
for p in [1, 2, 4, 8]:
    print(f"P={p} objective={obj[p]}")
PY
