#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/fcale/Dropbox/ChatGPT/Constants"
RUN="${ROOT}/tools/run_flint_objective.sh"

extract_objective() {
  awk -F= '/^objective=/{print $2}'
}

obj_p1="$(
  "${RUN}" --coeffs 1 --nmax 256 --k 32 --prec-dps 120 --print-dps 90 --tail true \
    | extract_objective
)"

obj_p2_fixed="$(
  "${RUN}" --coeffs 0.986,0.014 --nmax 256 --k 32 --prec-dps 120 --print-dps 90 --tail true \
    | extract_objective
)"

python3 - <<PY
from decimal import Decimal, getcontext
getcontext().prec = 120

obj_p1 = Decimal("${obj_p1}")
obj_p2 = Decimal("${obj_p2_fixed}")

ref_p1 = Decimal("0.574694861894825706792169188192308052244615365442183425654224244188404730797633575601903795")
ref_p2 = Decimal("0.574639647705652190514028806318248367678818613688424419427493565264744056487304757774002157")

err1 = abs(obj_p1 - ref_p1)
err2 = abs(obj_p2 - ref_p2)

print(f"P=1  objective={obj_p1}")
print(f"P=1  absErr={err1}")
print(f"P=2  objective={obj_p2}")
print(f"P=2  absErr={err2}")

tol = Decimal("1e-70")
if err1 > tol or err2 > tol:
    raise SystemExit("FLINT objective regression failed")
PY

echo "FLINT objective regression passed."
