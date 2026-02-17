#!/usr/bin/env python3
import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from decimal import Decimal, getcontext
from pathlib import Path

ROOT = Path("/Users/fcale/Dropbox/ChatGPT/Constants")
RUNNER = ROOT / "tools" / "run_flint_objective.sh"
DEFAULT_RECORDS = ROOT / "checkpoints" / "flint_bootstrap_records.json"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Run FLINT Newton continuation in P with zero-padding warm starts."
    )
    p.add_argument("--p-list", default="1,2,4,8", help="Comma-separated list of P values.")
    p.add_argument(
        "--stage-list",
        default="",
        help=(
            "Optional comma-separated stage specs. "
            "Each stage is P[:N[:K[:PREC[:MAXIT[:TOL[:CTOL]]]]]]. "
            "If omitted, global options are used."
        ),
    )
    p.add_argument("--nmax", type=int, default=256)
    p.add_argument("--k", type=int, default=32)
    p.add_argument("--prec-dps", type=int, default=120)
    p.add_argument("--print-dps", type=int, default=70)
    p.add_argument("--max-it", type=int, default=30)
    p.add_argument("--tol", default="1e-24")
    p.add_argument("--constraint-tol", default="1e-24")
    p.add_argument("--damping", choices=["true", "false"], default="true")
    p.add_argument("--step-min-exp", type=int, default=20)
    p.add_argument("--log", choices=["true", "false"], default="false")
    p.add_argument("--records", default=str(DEFAULT_RECORDS))
    p.add_argument("--no-save", action="store_true")
    p.add_argument("--resume-from-records", choices=["true", "false"], default="true")
    p.add_argument("--start-coeffs", default="", help="Optional CSV coefficients for the first stage.")
    p.add_argument("--max-retries", type=int, default=2, help="Retry attempts per stage on failure.")
    p.add_argument("--retry-prec-step", type=int, default=20, help="Extra precision digits per retry.")
    p.add_argument("--retry-max-it-step", type=int, default=8, help="Extra max iterations per retry.")
    p.add_argument(
        "--retry-step-min-exp-step",
        type=int,
        default=2,
        help="Extra damping floor exponent per retry.",
    )
    return p.parse_args()


def parse_kv_output(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_csv_numbers(csv: str) -> list[Decimal]:
    if not csv:
        return []
    return [Decimal(x) for x in csv.split(",")]


def format_coeffs(vals: list[Decimal], digits: int) -> str:
    getcontext().prec = max(getcontext().prec, digits + 20)
    out = []
    for v in vals:
        # Fixed-point string (no scientific notation) for robust parser input.
        s = format(v, f".{digits}f")
        # Trim trailing zeros but keep at least one digit after decimal for stability.
        if "." in s:
            s = s.rstrip("0").rstrip(".")
        out.append(s if s else "0")
    return ",".join(out)


def project_sum_to_one(vals: list[Decimal]) -> list[Decimal]:
    n = Decimal(len(vals))
    diff = (Decimal(1) - sum(vals)) / n
    return [v + diff for v in vals]


def asymptotic_seed(p: int) -> list[Decimal]:
    vals = []
    for n in range(p):
        term = Decimal((-1) ** n) / (Decimal(8) ** n * Decimal(n + 1))
        vals.append(term)
    s = sum(vals)
    vals = [v / s for v in vals]
    return project_sum_to_one(vals)


def warm_start(prev: list[Decimal], p: int) -> list[Decimal]:
    if len(prev) >= p:
        vals = prev[:p]
    else:
        vals = prev + [Decimal(0)] * (p - len(prev))
    return project_sum_to_one(vals)


def parse_stage_list(spec: str, args: argparse.Namespace) -> list[dict]:
    if not spec.strip():
        p_list = [int(x) for x in args.p_list.split(",") if x.strip()]
        if not p_list:
            raise ValueError("Empty --p-list")
        return [
            {
                "P": p,
                "Nmax": args.nmax,
                "K": args.k,
                "precDps": args.prec_dps,
                "maxIt": args.max_it,
                "tol": str(args.tol),
                "constraintTol": str(args.constraint_tol),
                "stepMinExp": args.step_min_exp,
            }
            for p in p_list
        ]

    stages: list[dict] = []
    for raw_stage in spec.split(","):
        raw_stage = raw_stage.strip()
        if not raw_stage:
            continue
        parts = raw_stage.split(":")
        if len(parts) > 7:
            raise ValueError(f"Bad stage spec: {raw_stage}")
        stage = {
            "P": int(parts[0]),
            "Nmax": args.nmax,
            "K": args.k,
            "precDps": args.prec_dps,
            "maxIt": args.max_it,
            "tol": str(args.tol),
            "constraintTol": str(args.constraint_tol),
            "stepMinExp": args.step_min_exp,
        }
        if len(parts) >= 2 and parts[1]:
            stage["Nmax"] = int(parts[1])
        if len(parts) >= 3 and parts[2]:
            stage["K"] = int(parts[2])
        if len(parts) >= 4 and parts[3]:
            stage["precDps"] = int(parts[3])
        if len(parts) >= 5 and parts[4]:
            stage["maxIt"] = int(parts[4])
        if len(parts) >= 6 and parts[5]:
            stage["tol"] = parts[5]
        if len(parts) >= 7 and parts[6]:
            stage["constraintTol"] = parts[6]
        stages.append(stage)

    if not stages:
        raise ValueError("Empty --stage-list")
    return stages


def run_stage(coeffs: list[Decimal], stage: dict, args: argparse.Namespace) -> dict[str, str]:
    coeffs_csv = format_coeffs(coeffs, max(stage["precDps"], args.print_dps))
    cmd = [
        str(RUNNER),
        "--coeffs", coeffs_csv,
        "--nmax", str(stage["Nmax"]),
        "--k", str(stage["K"]),
        "--prec-dps", str(stage["precDps"]),
        "--print-dps", str(args.print_dps),
        "--tail", "true",
        "--optimize", "true",
        "--max-it", str(stage["maxIt"]),
        "--tol", str(stage["tol"]),
        "--constraint-tol", str(stage["constraintTol"]),
        "--damping", args.damping,
        "--step-min-exp", str(stage["stepMinExp"]),
        "--log", args.log,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"FLINT stage failed for P={stage['P']}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    out = parse_kv_output(proc.stdout)
    out["_stdout"] = proc.stdout
    return out


def run_stage_with_retries(
    coeffs: list[Decimal], stage: dict, args: argparse.Namespace
) -> tuple[dict[str, str], dict, list[str]]:
    stage_run = dict(stage)
    attempt_logs: list[dict] = []
    messages: list[str] = []
    last_out: dict[str, str] | None = None

    for attempt in range(args.max_retries + 1):
        try:
            out = run_stage(coeffs, stage_run, args)
            converged = out.get("converged", "false") == "true"
            kkt_ok = out.get("kktSolveOk", "true") == "true"
            last_out = out

            attempt_logs.append(
                {
                    "attempt": attempt,
                    "params": dict(stage_run),
                    "status": "ok",
                    "converged": converged,
                    "kktSolveOk": kkt_ok,
                    "iterations": out.get("iterations", ""),
                    "objective": out.get("objective", ""),
                    "resInfUpper": out.get("resInfUpper", ""),
                    "constraintAbsUpper": out.get("constraintAbsUpper", ""),
                    "timeSec": out.get("timeSec", ""),
                }
            )

            if converged and kkt_ok:
                out["_attemptLogs"] = json.dumps(attempt_logs)
                return out, stage_run, messages

            if attempt < args.max_retries:
                msg = (
                    f"retry {attempt + 1}/{args.max_retries}: "
                    f"P={stage_run['P']} did not fully converge "
                    f"(converged={converged}, kktSolveOk={kkt_ok}); "
                    "escalating precision/iterations."
                )
                messages.append(msg)
                stage_run["precDps"] = int(stage_run["precDps"]) + int(args.retry_prec_step)
                stage_run["maxIt"] = int(stage_run["maxIt"]) + int(args.retry_max_it_step)
                stage_run["stepMinExp"] = min(
                    60, int(stage_run["stepMinExp"]) + int(args.retry_step_min_exp_step)
                )
            continue
        except Exception as exc:
            attempt_logs.append(
                {
                    "attempt": attempt,
                    "params": dict(stage_run),
                    "status": "error",
                    "error": str(exc),
                }
            )
            if attempt < args.max_retries:
                msg = (
                    f"retry {attempt + 1}/{args.max_retries}: "
                    f"P={stage_run['P']} stage error ({exc}); escalating and retrying."
                )
                messages.append(msg)
                stage_run["precDps"] = int(stage_run["precDps"]) + int(args.retry_prec_step)
                stage_run["maxIt"] = int(stage_run["maxIt"]) + int(args.retry_max_it_step)
                stage_run["stepMinExp"] = min(
                    60, int(stage_run["stepMinExp"]) + int(args.retry_step_min_exp_step)
                )
                continue
            raise RuntimeError(
                f"FLINT stage failed after retries for P={stage_run['P']}: {exc}"
            ) from exc

    if last_out is None:
        raise RuntimeError(f"FLINT stage produced no output for P={stage['P']}")

    last_out["_attemptLogs"] = json.dumps(attempt_logs)
    return last_out, stage_run, messages


def load_records(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return json.loads(path.read_text())


def save_records(path: Path, records: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(records, indent=2))


def find_latest_record(records: list[dict], stage: dict) -> dict | None:
    candidates = []
    for r in records:
        if r.get("P") != stage["P"]:
            continue
        params = r.get("params", {})
        if int(params.get("Nmax", -1)) != int(stage["Nmax"]):
            continue
        if int(params.get("K", -1)) != int(stage["K"]):
            continue
        if int(params.get("precDps", -1)) > int(stage["precDps"]):
            continue
        candidates.append(r)
    if not candidates:
        return None
    candidates.sort(key=lambda x: x.get("timestamp", ""))
    return candidates[-1]


def main() -> int:
    args = parse_args()
    stages = parse_stage_list(args.stage_list, args)
    max_prec = max(int(s["precDps"]) for s in stages) + int(args.max_retries) * int(
        args.retry_prec_step
    )
    getcontext().prec = max(max_prec + 20, 180)

    records_path = Path(args.records)
    records = load_records(records_path)

    prev_coeffs: list[Decimal] | None = None
    start_p_list = [int(s["P"]) for s in stages]
    print(f"FLINT bootstrap start: stages={start_p_list}")

    if args.start_coeffs.strip():
        prev_coeffs = parse_csv_numbers(args.start_coeffs.strip())

    for stage in stages:
        p = int(stage["P"])
        if prev_coeffs is None:
            if args.resume_from_records == "true":
                rec = find_latest_record(records, stage)
                if rec is not None:
                    init = parse_csv_numbers(rec["result"]["aCsv"])
                    init_source = "saved-record"
                else:
                    init = asymptotic_seed(p)
                    init_source = "asymptotic-seed"
            else:
                init = asymptotic_seed(p)
                init_source = "asymptotic-seed"
        else:
            init = warm_start(prev_coeffs, p)
            init_source = "previous-stage"

        out, stage_used, retry_messages = run_stage_with_retries(init, stage, args)
        for msg in retry_messages:
            print(msg)

        converged = out.get("converged", "false")
        iters = out.get("iterations", "?")
        obj = out.get("objective", "")
        res = out.get("resInfUpper", "")
        cabs = out.get("constraintAbsUpper", "")
        tsec = out.get("timeSec", "")
        kkt_ok = out.get("kktSolveOk", "true")

        print(
            f"P={p} N={stage_used['Nmax']} K={stage_used['K']} prec={stage_used['precDps']} "
            f"init={init_source} converged={converged} kktSolveOk={kkt_ok} iters={iters} "
            f"objective={obj} resInf={res} constraintAbs={cabs} timeSec={tsec}"
        )

        a_csv = out.get("aCsv", "")
        if not a_csv:
            raise RuntimeError(f"Missing aCsv in stage output for P={p}")
        prev_coeffs = parse_csv_numbers(a_csv)

        rec = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "P": p,
            "initSource": init_source,
            "params": {
                "Nmax": stage_used["Nmax"],
                "K": stage_used["K"],
                "precDps": stage_used["precDps"],
                "maxIt": stage_used["maxIt"],
                "tol": str(stage_used["tol"]),
                "constraintTol": str(stage_used["constraintTol"]),
                "damping": args.damping,
                "stepMinExp": stage_used["stepMinExp"],
            },
            "result": {
                "converged": converged == "true",
                "kktSolveOk": kkt_ok == "true",
                "iterations": int(iters) if iters.isdigit() else iters,
                "objective": obj,
                "nu2Candidate": out.get("nu2Candidate", ""),
                "lambda": out.get("lambda", ""),
                "resInfUpper": res,
                "constraintAbsUpper": cabs,
                "timeSec": tsec,
                "aCsv": a_csv,
                "retryMessages": retry_messages,
                "attemptLogs": json.loads(out.get("_attemptLogs", "[]"))
                if out.get("_attemptLogs")
                else [],
            },
        }
        records.append(rec)

    if not args.no_save:
        save_records(records_path, records)
        print(f"Saved records: {records_path}")

    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
