#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/fcale/Dropbox/ChatGPT/Constants"
RUN_FLINT="${ROOT}/tools/run_flint_objective.sh"

NMAX=256
K=32
PREC=110
PRINT=70

COEFFS_CSV="0.986020809644789355689012836072932717254029034983701632178039766988506187064702897092183357128581,0.014894771298517023283538955113605724901665282561444506716994265687921326857733437091781803016145,-0.000989950749235760146946391506037425262624327911249788818325577261366934867249833622042701015273,0.000074369805929381174394600319498983106930010366103649923291544584939420944813499438077540870546"

TMP_CPP="$(mktemp)"
TMP_WL="$(mktemp)"
TMP_WL_SCRIPT="$(mktemp /tmp/flint_deriv_XXXXXX.wl)"
trap 'rm -f "$TMP_CPP" "$TMP_WL" "$TMP_WL_SCRIPT"' EXIT

"${RUN_FLINT}" \
  --coeffs "${COEFFS_CSV}" \
  --nmax "${NMAX}" \
  --k "${K}" \
  --prec-dps "${PREC}" \
  --print-dps "${PRINT}" \
  --tail true \
  --derivatives true \
  --hessian true \
  > "${TMP_CPP}"

cat > "${TMP_WL_SCRIPT}" <<'EOF'
repoRoot = "/Users/fcale/Dropbox/ChatGPT/Constants";
AppendTo[$Path, FileNameJoin[{repoRoot, "wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "Constants.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "Newton.wl"}]];
Get[FileNameJoin[{repoRoot, "wl", "TailDerivatives.wl"}]];

Off[N::meprec];

prec = 110;
ConstantsSetPrecision[prec];

a = {
  0.986020809644789355689012836072932717254029034983701632178039766988506187064702897092183357128581,
  0.014894771298517023283538955113605724901665282561444506716994265687921326857733437091781803016145,
  -0.000989950749235760146946391506037425262624327911249788818325577261366934867249833622042701015273,
  0.000074369805929381174394600319498983106930010366103649923291544584939420944813499438077540870546
};
P = Length[a];
Nmax = 256;
K = 32;

repFin = Constants`Newton`NewtonDerivatives[a, P, Nmax, 0];
obj = repFin["objective"] + Constants`TailObjective[a, Nmax, K];
grad = repFin["gpart"] + Constants`TailDerivatives`TailGradient[a, Nmax, K];
hess = repFin["hess"] + Constants`TailDerivatives`TailHessian[a, Nmax, K];

toS[x_] := StringReplace[ToString[N[x, 70], InputForm], RegularExpression["`[0-9.]+"] -> ""];

json = <|
  "objective" -> toS[obj],
  "gradient" -> (toS /@ grad),
  "hessian" -> ((toS /@ #) & /@ hess)
|>;

Print[ExportString[json, "JSON"]];
Exit[0];
EOF

wolframscript -file "${TMP_WL_SCRIPT}" > "${TMP_WL}"

python3 - "${TMP_CPP}" "${TMP_WL}" <<'PY'
import json
import sys
from decimal import Decimal, getcontext

getcontext().prec = 120

cpp_path, wl_path = sys.argv[1], sys.argv[2]

cpp = {}
with open(cpp_path, "r", encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if "=" in line:
            k, v = line.split("=", 1)
            cpp[k] = v

def to_dec(s: str) -> Decimal:
    return Decimal(s.replace("*^", "e"))

def parse_cpp_vec(s: str):
    return [to_dec(x) for x in s.split(",") if x]

def parse_cpp_mat(s: str):
    return [parse_cpp_vec(row) for row in s.split(";") if row]

if "objective" not in cpp or "gradientCsv" not in cpp or "hessianRows" not in cpp:
    raise SystemExit("Missing expected FLINT outputs")

wl_raw = open(wl_path, "r", encoding="utf-8").read().strip()
wl = json.loads(wl_raw)

obj_cpp = to_dec(cpp["objective"])
obj_wl = to_dec(wl["objective"])
grad_cpp = parse_cpp_vec(cpp["gradientCsv"])
grad_wl = [to_dec(x) for x in wl["gradient"]]
hess_cpp = parse_cpp_mat(cpp["hessianRows"])
hess_wl = [[to_dec(x) for x in row] for row in wl["hessian"]]

if len(grad_cpp) != len(grad_wl):
    raise SystemExit("Gradient length mismatch")
if len(hess_cpp) != len(hess_wl):
    raise SystemExit("Hessian row count mismatch")

max_grad_err = Decimal(0)
for a, b in zip(grad_cpp, grad_wl):
    e = abs(a - b)
    if e > max_grad_err:
        max_grad_err = e

max_hess_err = Decimal(0)
for row_c, row_w in zip(hess_cpp, hess_wl):
    if len(row_c) != len(row_w):
        raise SystemExit("Hessian column count mismatch")
    for a, b in zip(row_c, row_w):
        e = abs(a - b)
        if e > max_hess_err:
            max_hess_err = e

obj_err = abs(obj_cpp - obj_wl)

print(f"objective absErr={obj_err}")
print(f"gradient maxAbsErr={max_grad_err}")
print(f"hessian maxAbsErr={max_hess_err}")

tol = Decimal("1e-42")
if obj_err > tol or max_grad_err > tol or max_hess_err > tol:
    raise SystemExit("FLINT derivative regression failed")
PY

echo "FLINT derivative regression passed."
