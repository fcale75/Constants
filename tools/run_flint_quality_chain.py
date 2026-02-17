#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation, getcontext
from pathlib import Path
from typing import Dict, List, Optional, Tuple

ROOT = Path("/Users/fcale/Dropbox/ChatGPT/Constants")
RUNNER = ROOT / "tools" / "run_flint_objective.sh"


@dataclass
class StageSpec:
    p: int
    nmax: int
    k: int
    prec_dps: int
    max_it: int
    tol: str
    constraint_tol: str
    step_min_exp: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run FLINT continuation in P with strict quality gates. "
            "Each P-stage runs a base solve and a stricter check solve; "
            "the stage only passes when achieved digits meet expectation."
        )
    )
    parser.add_argument(
        "--stage-list",
        default="32:8192:128:340:24:1e-28:1e-36:40,50:12288:128:360:28:1e-30:1e-40:42,100:16384:128:384:32:1e-34:1e-44:45",
        help=(
            "Comma-separated stage specs: "
            "P:N:K:PREC:MAXIT:TOL:CTOL:STEPMINEXP"
        ),
    )
    parser.add_argument(
        "--seed-coeffs",
        default="",
        help="CSV coefficients for first stage; when omitted, seed records or asymptotic seed is used.",
    )
    parser.add_argument(
        "--seed-records",
        default=str(ROOT / "checkpoints" / "flint_chain_p50_p100_p150_reg1e22_tol1e13.json"),
        help="JSON checkpoint file used to warm-start when --seed-coeffs is omitted.",
    )
    parser.add_argument(
        "--seed-p",
        type=int,
        default=32,
        help="Preferred P from --seed-records (fallback: highest P <= first stage P).",
    )
    parser.add_argument("--print-dps", type=int, default=120)
    parser.add_argument("--min-it", type=int, default=1)
    parser.add_argument("--damping", choices=["true", "false"], default="true")
    parser.add_argument("--kkt-scaling", choices=["asymptotic", "none"], default="asymptotic")
    parser.add_argument("--kkt-reg", default="1e-24")
    parser.add_argument("--max-abs-a", type=float, default=0.0)
    parser.add_argument("--step-max-u", type=float, default=0.0)
    parser.add_argument("--nonnegative-grid", type=int, default=400)
    parser.add_argument("--nonnegative-tol", type=float, default=1e-28)
    parser.add_argument("--log", choices=["true", "false"], default="false")

    parser.add_argument("--check-n-delta", type=int, default=4096)
    parser.add_argument("--check-k-delta", type=int, default=32)
    parser.add_argument("--check-prec-delta", type=int, default=40)
    parser.add_argument("--check-max-it-delta", type=int, default=8)
    parser.add_argument(
        "--check-tol-decades",
        type=int,
        default=4,
        help="Check solve uses tolerance scaled by 10^(-decades).",
    )
    parser.add_argument(
        "--check-constraint-tol-decades",
        type=int,
        default=4,
        help="Check solve uses constraint tolerance scaled by 10^(-decades).",
    )

    parser.add_argument("--max-attempts", type=int, default=4)
    parser.add_argument("--escalate-n-mult", type=float, default=1.5)
    parser.add_argument("--escalate-k-add", type=int, default=16)
    parser.add_argument("--escalate-prec-add", type=int, default=20)
    parser.add_argument("--escalate-max-it-add", type=int, default=6)
    parser.add_argument("--escalate-step-min-exp-add", type=int, default=2)
    parser.add_argument("--escalate-tol-decades", type=int, default=1)
    parser.add_argument("--escalate-constraint-tol-decades", type=int, default=1)
    parser.add_argument("--reg-decay", default="0.1")
    parser.add_argument("--reg-min", default="1e-30")

    parser.add_argument(
        "--expected-digits-map",
        default="32:40,50:62,100:128,150:190",
        help="Optional per-stage nu2 digits target as P:digits list.",
    )
    parser.add_argument("--expected-digits-slope", type=float, default=1.25)
    parser.add_argument("--expected-digits-intercept", type=float, default=0.0)
    parser.add_argument("--expected-digits-floor", type=float, default=20.0)
    parser.add_argument(
        "--coeff-digits-offset",
        type=float,
        default=8.0,
        help="Required coeff digits = max(min-coeff-digits, expected-nu2-digits - offset).",
    )
    parser.add_argument("--min-coeff-digits", type=float, default=24.0)

    parser.add_argument(
        "--out",
        default=str(ROOT / "checkpoints" / "flint_quality_chain.json"),
        help="JSON report path (updated after each attempt).",
    )
    return parser.parse_args()


def parse_stage_list(raw: str) -> List[StageSpec]:
    stages: List[StageSpec] = []
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        parts = token.split(":")
        if len(parts) != 8:
            raise ValueError(f"Bad stage spec: {token}")
        stages.append(
            StageSpec(
                p=int(parts[0]),
                nmax=int(parts[1]),
                k=int(parts[2]),
                prec_dps=int(parts[3]),
                max_it=int(parts[4]),
                tol=parts[5],
                constraint_tol=parts[6],
                step_min_exp=int(parts[7]),
            )
        )
    if not stages:
        raise ValueError("Empty --stage-list")
    return stages


def parse_expected_digits_map(raw: str) -> Dict[int, float]:
    out: Dict[int, float] = {}
    if not raw.strip():
        return out
    for token in raw.split(","):
        token = token.strip()
        if not token:
            continue
        parts = token.split(":")
        if len(parts) != 2:
            raise ValueError(f"Bad expected digits map token: {token}")
        out[int(parts[0])] = float(parts[1])
    return out


def parse_kv_output(stdout: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for line in stdout.splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def parse_csv_decimals(csv: str) -> List[Decimal]:
    if not csv:
        return []
    return [Decimal(x) for x in csv.split(",")]


def format_coeffs(vals: List[Decimal], digits: int) -> str:
    out: List[str] = []
    for v in vals:
        s = format(v, f".{digits}f")
        if "." in s:
            s = s.rstrip("0").rstrip(".")
        out.append(s if s else "0")
    return ",".join(out)


def project_sum_to_one(vals: List[Decimal]) -> List[Decimal]:
    diff = (Decimal(1) - sum(vals)) / Decimal(len(vals))
    return [v + diff for v in vals]


def warm_start(prev: List[Decimal], p: int) -> List[Decimal]:
    if len(prev) >= p:
        vals = prev[:p]
    else:
        vals = prev + [Decimal(0)] * (p - len(prev))
    return project_sum_to_one(vals)


def asymptotic_seed(p: int) -> List[Decimal]:
    vals: List[Decimal] = []
    for n in range(p):
        term = Decimal((-1) ** n) / (Decimal(8) ** n * Decimal(n + 1))
        vals.append(term)
    total = sum(vals)
    vals = [v / total for v in vals]
    return project_sum_to_one(vals)


def collect_seed_candidates(node, out: List[Tuple[int, str, str]]) -> None:
    if isinstance(node, dict):
        p = node.get("P")
        a_csv = node.get("aCsv", "")
        if not a_csv:
            result = node.get("result")
            if isinstance(result, dict):
                a_csv = result.get("aCsv", "")
        if p is not None and isinstance(a_csv, str) and a_csv.strip():
            ts = str(node.get("timestamp", ""))
            out.append((int(p), ts, a_csv))
        for value in node.values():
            collect_seed_candidates(value, out)
        return
    if isinstance(node, list):
        for item in node:
            collect_seed_candidates(item, out)


def choose_seed(
    seed_coeffs: str,
    seed_records: Path,
    seed_p: int,
    first_stage_p: int,
) -> Tuple[List[Decimal], str]:
    if seed_coeffs.strip():
        vals = project_sum_to_one(parse_csv_decimals(seed_coeffs.strip()))
        return vals, "explicit-seed"

    if seed_records.exists():
        payload = json.loads(seed_records.read_text())
        candidates: List[Tuple[int, str, str]] = []
        collect_seed_candidates(payload, candidates)
        if candidates:
            candidates = sorted(candidates, key=lambda x: (x[0], x[1]))
            exact = [c for c in candidates if c[0] == seed_p]
            if exact:
                chosen = exact[-1]
            else:
                below = [c for c in candidates if c[0] <= first_stage_p]
                chosen = below[-1] if below else candidates[-1]
            vals = project_sum_to_one(parse_csv_decimals(chosen[2]))
            return vals, f"records:{seed_records.name}:P{chosen[0]}"

    return asymptotic_seed(first_stage_p), "asymptotic-seed"


def run_flint(
    coeffs: List[Decimal],
    stage: StageSpec,
    args: argparse.Namespace,
    kkt_reg: str,
) -> Dict[str, str]:
    coeffs_csv = format_coeffs(coeffs, max(stage.prec_dps, args.print_dps))
    cmd = [
        str(RUNNER),
        "--coeffs", coeffs_csv,
        "--nmax", str(stage.nmax),
        "--k", str(stage.k),
        "--prec-dps", str(stage.prec_dps),
        "--print-dps", str(args.print_dps),
        "--tail", "true",
        "--optimize", "true",
        "--max-it", str(stage.max_it),
        "--min-it", str(args.min_it),
        "--tol", str(stage.tol),
        "--constraint-tol", str(stage.constraint_tol),
        "--damping", str(args.damping),
        "--step-min-exp", str(stage.step_min_exp),
        "--kkt-scaling", str(args.kkt_scaling),
        "--kkt-reg", str(kkt_reg),
        "--max-abs-a", str(args.max_abs_a),
        "--step-max-u", str(args.step_max_u),
        "--nonnegative-grid", str(args.nonnegative_grid),
        "--nonnegative-tol", str(args.nonnegative_tol),
        "--log", str(args.log),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"FLINT run failed for P={stage.p}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    out = parse_kv_output(proc.stdout)
    out["_stdout"] = proc.stdout
    return out


def to_decimal(s: str) -> Optional[Decimal]:
    if s is None:
        return None
    s = str(s).strip()
    if not s:
        return None
    try:
        return Decimal(s)
    except InvalidOperation:
        return None


def digits_from_diff(diff: Optional[Decimal]) -> Optional[float]:
    if diff is None:
        return None
    if diff == 0:
        return float("inf")
    if diff < 0:
        diff = -diff
    try:
        return float(-diff.log10())
    except (InvalidOperation, ValueError):
        return None


def max_abs_diff(a: List[Decimal], b: List[Decimal]) -> Optional[Decimal]:
    if len(a) != len(b) or not a:
        return None
    out = Decimal(0)
    for x, y in zip(a, b):
        cur = abs(x - y)
        if cur > out:
            out = cur
    return out


def tighten_tol(raw: str, decades: int) -> str:
    val = Decimal(raw)
    out = val.scaleb(-decades)
    return f"{out:.6E}".replace("E+", "E")


def expected_digits_for_p(
    p: int,
    expected_map: Dict[int, float],
    slope: float,
    intercept: float,
    floor_val: float,
) -> float:
    if p in expected_map:
        return expected_map[p]
    return max(floor_val, slope * p + intercept)


def maybe_inf(v: Optional[float]):
    if v is None:
        return None
    if v == float("inf"):
        return "inf"
    return v


def write_report(path: Path, meta: Dict, stages: List[Dict]) -> None:
    payload = {
        "meta": meta,
        "stages": stages,
        "overallStatus": "passed" if stages and all(s.get("status") == "passed" for s in stages) else "failed",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2))


def main() -> int:
    args = parse_args()
    stages = parse_stage_list(args.stage_list)
    expected_map = parse_expected_digits_map(args.expected_digits_map)

    max_prec = max(s.prec_dps for s in stages) + args.check_prec_delta + args.escalate_prec_add * max(args.max_attempts, 1)
    getcontext().prec = max(200, max_prec + 80)

    coeffs, seed_source = choose_seed(
        args.seed_coeffs,
        Path(args.seed_records),
        args.seed_p,
        stages[0].p,
    )

    meta = {
        "timestampStart": datetime.now(timezone.utc).isoformat(),
        "seedSource": seed_source,
        "stageList": [s.__dict__ for s in stages],
        "checkDeltas": {
            "n": args.check_n_delta,
            "k": args.check_k_delta,
            "prec": args.check_prec_delta,
            "maxIt": args.check_max_it_delta,
            "tolDecades": args.check_tol_decades,
            "constraintTolDecades": args.check_constraint_tol_decades,
        },
        "escalation": {
            "maxAttempts": args.max_attempts,
            "nMult": args.escalate_n_mult,
            "kAdd": args.escalate_k_add,
            "precAdd": args.escalate_prec_add,
            "maxItAdd": args.escalate_max_it_add,
            "stepMinExpAdd": args.escalate_step_min_exp_add,
            "tolDecades": args.escalate_tol_decades,
            "constraintTolDecades": args.escalate_constraint_tol_decades,
            "regDecay": args.reg_decay,
            "regMin": args.reg_min,
        },
        "expectedDigitsMap": expected_map,
        "expectedDigitsSlope": args.expected_digits_slope,
        "expectedDigitsIntercept": args.expected_digits_intercept,
        "expectedDigitsFloor": args.expected_digits_floor,
        "coeffDigitsOffset": args.coeff_digits_offset,
        "minCoeffDigits": args.min_coeff_digits,
        "global": {
            "printDps": args.print_dps,
            "minIt": args.min_it,
            "damping": args.damping,
            "kktScaling": args.kkt_scaling,
            "kktRegInitial": args.kkt_reg,
            "maxAbsA": args.max_abs_a,
            "stepMaxU": args.step_max_u,
            "nonnegativeGrid": args.nonnegative_grid,
            "nonnegativeTol": args.nonnegative_tol,
            "log": args.log,
        },
    }

    stage_reports: List[Dict] = []
    out_path = Path(args.out)

    reg_decay = Decimal(str(args.reg_decay))
    reg_min = Decimal(str(args.reg_min))

    for base_stage in stages:
        expected_nu2 = expected_digits_for_p(
            base_stage.p,
            expected_map,
            args.expected_digits_slope,
            args.expected_digits_intercept,
            args.expected_digits_floor,
        )
        required_coeff_digits = max(args.min_coeff_digits, expected_nu2 - args.coeff_digits_offset)

        stage_report: Dict = {
            "P": base_stage.p,
            "expectedNu2Digits": expected_nu2,
            "requiredCoeffDigits": required_coeff_digits,
            "attempts": [],
            "status": "failed",
        }

        # Guarantee the current seed matches this stage P.
        coeffs = warm_start(coeffs, base_stage.p)

        stage = StageSpec(**base_stage.__dict__)
        current_reg = Decimal(str(args.kkt_reg))
        stage_passed = False
        final_coeffs: Optional[List[Decimal]] = None

        for attempt in range(args.max_attempts + 1):
            check_stage = StageSpec(
                p=stage.p,
                nmax=stage.nmax + args.check_n_delta,
                k=stage.k + args.check_k_delta,
                prec_dps=stage.prec_dps + args.check_prec_delta,
                max_it=stage.max_it + args.check_max_it_delta,
                tol=tighten_tol(stage.tol, args.check_tol_decades),
                constraint_tol=tighten_tol(stage.constraint_tol, args.check_constraint_tol_decades),
                step_min_exp=min(60, stage.step_min_exp + 2),
            )

            base_out = run_flint(coeffs, stage, args, f"{current_reg:.3E}".replace("E+", "E"))
            base_csv = base_out.get("aCsv", "")
            if not base_csv:
                raise RuntimeError(f"Missing aCsv in base stage output for P={stage.p}")
            base_coeffs = parse_csv_decimals(base_csv)

            check_out = run_flint(base_coeffs, check_stage, args, f"{current_reg:.3E}".replace("E+", "E"))
            check_csv = check_out.get("aCsv", "")
            if not check_csv:
                raise RuntimeError(f"Missing aCsv in check stage output for P={stage.p}")
            check_coeffs = parse_csv_decimals(check_csv)

            base_nu2 = to_decimal(base_out.get("nu2Candidate", ""))
            check_nu2 = to_decimal(check_out.get("nu2Candidate", ""))
            base_obj = to_decimal(base_out.get("objective", ""))
            check_obj = to_decimal(check_out.get("objective", ""))
            nu2_diff = None if (base_nu2 is None or check_nu2 is None) else abs(check_nu2 - base_nu2)
            obj_diff = None if (base_obj is None or check_obj is None) else abs(check_obj - base_obj)
            coeff_diff = max_abs_diff(base_coeffs, check_coeffs)

            nu2_digits = digits_from_diff(nu2_diff)
            obj_digits = digits_from_diff(obj_diff)
            coeff_digits = digits_from_diff(coeff_diff)

            base_conv = base_out.get("converged", "false") == "true"
            base_kkt = base_out.get("kktSolveOk", "false") == "true"
            check_conv = check_out.get("converged", "false") == "true"
            check_kkt = check_out.get("kktSolveOk", "false") == "true"

            pass_nu2 = nu2_digits is not None and nu2_digits >= expected_nu2
            pass_coeff = coeff_digits is not None and coeff_digits >= required_coeff_digits
            passed = base_conv and base_kkt and check_conv and check_kkt and pass_nu2 and pass_coeff

            attempt_row = {
                "attempt": attempt,
                "kktReg": f"{current_reg:.3E}".replace("E+", "E"),
                "baseStage": stage.__dict__,
                "checkStage": check_stage.__dict__,
                "base": {
                    "converged": base_conv,
                    "kktSolveOk": base_kkt,
                    "iterations": base_out.get("iterations", ""),
                    "resInfUpper": base_out.get("resInfUpper", ""),
                    "constraintAbsUpper": base_out.get("constraintAbsUpper", ""),
                    "timeSec": base_out.get("timeSec", ""),
                    "objective": base_out.get("objective", ""),
                    "nu2Candidate": base_out.get("nu2Candidate", ""),
                },
                "check": {
                    "converged": check_conv,
                    "kktSolveOk": check_kkt,
                    "iterations": check_out.get("iterations", ""),
                    "resInfUpper": check_out.get("resInfUpper", ""),
                    "constraintAbsUpper": check_out.get("constraintAbsUpper", ""),
                    "timeSec": check_out.get("timeSec", ""),
                    "objective": check_out.get("objective", ""),
                    "nu2Candidate": check_out.get("nu2Candidate", ""),
                },
                "drift": {
                    "nu2AbsDiff": str(nu2_diff) if nu2_diff is not None else "",
                    "objectiveAbsDiff": str(obj_diff) if obj_diff is not None else "",
                    "coeffInfAbsDiff": str(coeff_diff) if coeff_diff is not None else "",
                    "nu2Digits": maybe_inf(nu2_digits),
                    "objectiveDigits": maybe_inf(obj_digits),
                    "coeffDigits": maybe_inf(coeff_digits),
                },
                "targets": {
                    "expectedNu2Digits": expected_nu2,
                    "requiredCoeffDigits": required_coeff_digits,
                },
                "passFlags": {
                    "baseConverged": base_conv,
                    "baseKkt": base_kkt,
                    "checkConverged": check_conv,
                    "checkKkt": check_kkt,
                    "nu2DigitsOk": pass_nu2,
                    "coeffDigitsOk": pass_coeff,
                },
                "passed": passed,
            }
            stage_report["attempts"].append(attempt_row)
            write_report(out_path, meta, stage_reports + [stage_report])

            print(
                f"P={stage.p} attempt={attempt} reg={attempt_row['kktReg']} "
                f"base(conv={base_conv},kkt={base_kkt},it={base_out.get('iterations','?')},res={base_out.get('resInfUpper','?')}) "
                f"check(conv={check_conv},kkt={check_kkt},it={check_out.get('iterations','?')},res={check_out.get('resInfUpper','?')}) "
                f"nu2Digits={attempt_row['drift']['nu2Digits']} coeffDigits={attempt_row['drift']['coeffDigits']} "
                f"targetNu2={expected_nu2:.2f} targetCoeff={required_coeff_digits:.2f} pass={passed}",
                flush=True,
            )

            if passed:
                stage_passed = True
                final_coeffs = check_coeffs
                stage_report["status"] = "passed"
                stage_report["final"] = {
                    "attempt": attempt,
                    "aCsv": check_csv,
                    "nu2Candidate": check_out.get("nu2Candidate", ""),
                    "objective": check_out.get("objective", ""),
                    "resInfUpper": check_out.get("resInfUpper", ""),
                    "constraintAbsUpper": check_out.get("constraintAbsUpper", ""),
                }
                break

            if attempt >= args.max_attempts:
                break

            stage.nmax = int(round(stage.nmax * args.escalate_n_mult))
            stage.k = stage.k + args.escalate_k_add
            stage.prec_dps = stage.prec_dps + args.escalate_prec_add
            stage.max_it = stage.max_it + args.escalate_max_it_add
            stage.step_min_exp = min(60, stage.step_min_exp + args.escalate_step_min_exp_add)
            stage.tol = tighten_tol(stage.tol, args.escalate_tol_decades)
            stage.constraint_tol = tighten_tol(stage.constraint_tol, args.escalate_constraint_tol_decades)
            current_reg = max(reg_min, current_reg * reg_decay)

            # Keep warm start from the stricter run even on failed gate.
            coeffs = check_coeffs

        stage_reports.append(stage_report)
        write_report(out_path, meta, stage_reports)

        if not stage_passed or final_coeffs is None:
            print(f"stage_failed P={base_stage.p} wrote={out_path}", flush=True)
            return 2

        coeffs = final_coeffs

    print(f"all_stages_passed wrote={out_path}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
