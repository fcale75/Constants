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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Run a fixed chain of FLINT optimize stages with explicit warm starts "
            "and checkpoint writes after every stage."
        )
    )
    p.add_argument(
        "--stage-list",
        default="50:16384:64:260:8,100:32768:64:260:8,150:32768:64:260:8",
        help="Comma list of P:N:K:PREC:MAXIT stage specs.",
    )
    p.add_argument(
        "--seed-coeffs",
        default="",
        help="Optional CSV coefficients for the first stage.",
    )
    p.add_argument(
        "--seed-records",
        default=str(ROOT / "checkpoints" / "flint_bootstrap_records_p50_candidate.json"),
        help="Records file used when --seed-coeffs is omitted.",
    )
    p.add_argument(
        "--seed-p",
        type=int,
        default=50,
        help="Preferred P to pull from --seed-records (fallback: latest <= first stage P).",
    )
    p.add_argument("--print-dps", type=int, default=70)
    p.add_argument("--min-it", type=int, default=3)
    p.add_argument("--tol", default="1e-13")
    p.add_argument("--constraint-tol", default="1e-24")
    p.add_argument("--damping", choices=["true", "false"], default="true")
    p.add_argument("--step-min-exp", type=int, default=35)
    p.add_argument("--kkt-scaling", choices=["asymptotic", "none"], default="asymptotic")
    p.add_argument("--kkt-reg", default="1e-22")
    p.add_argument("--max-abs-a", type=float, default=0.0)
    p.add_argument("--step-max-u", type=float, default=0.0)
    p.add_argument("--nonnegative-grid", type=int, default=200)
    p.add_argument("--nonnegative-tol", type=float, default=1e-18)
    p.add_argument("--log", choices=["true", "false"], default="false")
    p.add_argument(
        "--out",
        default=str(ROOT / "checkpoints" / "flint_chain_records.json"),
        help="JSON output path (rewritten after each stage).",
    )
    return p.parse_args()


def parse_stage_list(spec: str) -> list[dict]:
    out: list[dict] = []
    for raw in spec.split(","):
        raw = raw.strip()
        if not raw:
            continue
        parts = raw.split(":")
        if len(parts) != 5:
            raise ValueError(f"Bad stage spec: {raw}")
        out.append(
            {
                "P": int(parts[0]),
                "Nmax": int(parts[1]),
                "K": int(parts[2]),
                "precDps": int(parts[3]),
                "maxIt": int(parts[4]),
            }
        )
    if not out:
        raise ValueError("Empty --stage-list")
    return out


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
    out = []
    for v in vals:
        s = format(v, f".{digits}f")
        if "." in s:
            s = s.rstrip("0").rstrip(".")
        out.append(s if s else "0")
    return ",".join(out)


def project_sum_to_one(vals: list[Decimal]) -> tuple[list[Decimal], Decimal]:
    n = Decimal(len(vals))
    diff = (Decimal(1) - sum(vals)) / n
    return [v + diff for v in vals], diff


def warm_start(prev: list[Decimal], p: int) -> tuple[list[Decimal], Decimal]:
    if len(prev) >= p:
        vals = prev[:p]
    else:
        vals = prev + [Decimal(0)] * (p - len(prev))
    return project_sum_to_one(vals)


def choose_seed_from_records(records_path: Path, seed_p: int, first_p: int) -> tuple[list[Decimal], str]:
    if not records_path.exists():
        raise FileNotFoundError(f"Seed records not found: {records_path}")
    records = json.loads(records_path.read_text())

    candidates: list[tuple[int, str, str]] = []
    for r in records:
        p = r.get("P")
        ts = r.get("timestamp", "")
        a_csv = r.get("result", {}).get("aCsv", "")
        if p is None or not a_csv:
            continue
        candidates.append((int(p), ts, a_csv))
    if not candidates:
        raise ValueError(f"No usable rows in seed records: {records_path}")

    candidates.sort(key=lambda x: (x[0], x[1]))

    exact = [c for c in candidates if c[0] == seed_p]
    if exact:
        chosen = exact[-1]
    else:
        below = [c for c in candidates if c[0] <= first_p]
        chosen = below[-1] if below else candidates[-1]

    vals = parse_csv_numbers(chosen[2])
    vals, _ = project_sum_to_one(vals)
    return vals, f"records:{records_path.name}:P{chosen[0]}"


def run_stage(init: list[Decimal], stage: dict, args: argparse.Namespace) -> dict[str, str]:
    coeffs_csv = format_coeffs(init, max(stage["precDps"], args.print_dps))
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
        "--min-it", str(args.min_it),
        "--tol", str(args.tol),
        "--constraint-tol", str(args.constraint_tol),
        "--damping", args.damping,
        "--step-min-exp", str(args.step_min_exp),
        "--kkt-scaling", str(args.kkt_scaling),
        "--kkt-reg", str(args.kkt_reg),
        "--max-abs-a", str(args.max_abs_a),
        "--step-max-u", str(args.step_max_u),
        "--nonnegative-grid", str(args.nonnegative_grid),
        "--nonnegative-tol", str(args.nonnegative_tol),
        "--log", args.log,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"FLINT stage failed (P={stage['P']}).\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    out = parse_kv_output(proc.stdout)
    out["_stdout"] = proc.stdout
    return out


def write_rows(path: Path, meta: dict, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "meta": meta,
        "rows": rows,
    }
    path.write_text(json.dumps(payload, indent=2))


def main() -> int:
    args = parse_args()
    stages = parse_stage_list(args.stage_list)
    getcontext().prec = max(int(s["precDps"]) for s in stages) + 60

    if args.seed_coeffs.strip():
        seed, seed_source = project_sum_to_one(parse_csv_numbers(args.seed_coeffs.strip()))
        seed_source = "explicit-seed"
    else:
        seed, seed_source = choose_seed_from_records(
            Path(args.seed_records), args.seed_p, stages[0]["P"]
        )

    rows: list[dict] = []
    prev = seed

    meta = {
        "timestampStart": datetime.now(timezone.utc).isoformat(),
        "stageList": stages,
        "seedSource": seed_source,
        "tol": str(args.tol),
        "constraintTol": str(args.constraint_tol),
        "minIt": args.min_it,
        "kktReg": str(args.kkt_reg),
        "kktScaling": args.kkt_scaling,
        "stepMinExp": args.step_min_exp,
        "nonnegativeGrid": args.nonnegative_grid,
        "nonnegativeTol": args.nonnegative_tol,
        "damping": args.damping,
        "log": args.log,
    }

    print(
        f"FLINT chain start: stages={[s['P'] for s in stages]} seed={seed_source}",
        flush=True,
    )

    out_path = Path(args.out)

    for idx, stage in enumerate(stages):
        if idx == 0:
            init = prev
            proj_diff = Decimal(0)
        else:
            init, proj_diff = warm_start(prev, stage["P"])

        init_sum = sum(init)
        print(
            f"start P={stage['P']} initLen={len(init)} initSum={init_sum} projDiff={proj_diff}",
            flush=True,
        )

        out = run_stage(init, stage, args)
        a_csv = out.get("aCsv", "")
        if not a_csv:
            raise RuntimeError(f"Missing aCsv for P={stage['P']}")

        coeffs = parse_csv_numbers(a_csv)
        coeff_sum = sum(coeffs)
        max_abs = max(abs(v) for v in coeffs)

        row = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "P": stage["P"],
            "Nmax": stage["Nmax"],
            "K": stage["K"],
            "precDps": stage["precDps"],
            "maxIt": stage["maxIt"],
            "initLen": len(init),
            "initSum": str(init_sum),
            "initProjDiff": str(proj_diff),
            "converged": out.get("converged", ""),
            "kktSolveOk": out.get("kktSolveOk", ""),
            "iterations": out.get("iterations", ""),
            "objective": out.get("objective", ""),
            "nu2Candidate": out.get("nu2Candidate", ""),
            "lambda": out.get("lambda", ""),
            "resInfUpper": out.get("resInfUpper", ""),
            "constraintAbsUpper": out.get("constraintAbsUpper", ""),
            "timeSec": out.get("timeSec", ""),
            "timeBessel": out.get("timeBessel", ""),
            "timeTailPrecompute": out.get("timeTailPrecompute", ""),
            "timeOptimize": out.get("timeOptimize", ""),
            "sumA": str(coeff_sum),
            "maxAbsA": format(max_abs, "E"),
            "aCsv": a_csv,
        }
        rows.append(row)

        write_rows(out_path, meta, rows)
        print(
            f"done P={stage['P']} conv={row['converged']} kkt={row['kktSolveOk']} "
            f"iters={row['iterations']} res={row['resInfUpper']} "
            f"obj={row['objective'][:34]} time={row['timeSec']}",
            flush=True,
        )

        prev = coeffs

    print(f"wrote {out_path}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
