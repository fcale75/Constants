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
DEFAULT_BOOTSTRAP_RECORDS = ROOT / "checkpoints" / "flint_bootstrap_records.json"
DEFAULT_OUT = ROOT / "checkpoints" / "flint_precision_ladder.json"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Run fixed-(P,N,K) FLINT optimization across increasing precision/tolerances "
            "with warm starts."
        )
    )
    p.add_argument("--p", type=int, required=True, help="Coefficient count P.")
    p.add_argument("--nmax", type=int, required=True, help="Finite cutoff N.")
    p.add_argument("--k", type=int, required=True, help="Asymptotic terms K.")
    p.add_argument("--print-dps", type=int, default=90)
    p.add_argument("--min-it", type=int, default=1)
    p.add_argument(
        "--stages",
        default="220:1e-18:1e-24:40,280:1e-24:1e-30:60,340:1e-30:1e-36:80",
        help=(
            "Comma-separated stage specs PREC:TOL:CTOL:MAXIT[:STEPMINEXP]. "
            "Example: 240:1e-20:1e-28:50,300:1e-26:1e-34:70:40"
        ),
    )
    p.add_argument("--damping", choices=["true", "false"], default="true")
    p.add_argument("--step-min-exp", type=int, default=30)
    p.add_argument("--kkt-scaling", choices=["asymptotic", "none"], default="asymptotic")
    p.add_argument("--max-abs-a", type=float, default=0.0)
    p.add_argument("--step-max-u", type=float, default=0.0)
    p.add_argument("--nonnegative-grid", type=int, default=0)
    p.add_argument("--nonnegative-tol", type=float, default=0.0)
    p.add_argument("--kkt-reg", type=float, default=0.0)
    p.add_argument("--log", choices=["true", "false"], default="false")
    p.add_argument("--bootstrap-records", default=str(DEFAULT_BOOTSTRAP_RECORDS))
    p.add_argument("--start-coeffs", default="")
    p.add_argument("--out", default=str(DEFAULT_OUT))
    p.add_argument("--append", action="store_true")
    return p.parse_args()


def parse_stage_spec(spec: str, default_step_min_exp: int) -> list[dict]:
    out: list[dict] = []
    for raw in spec.split(","):
        raw = raw.strip()
        if not raw:
            continue
        parts = raw.split(":")
        if len(parts) < 4 or len(parts) > 5:
            raise ValueError(f"Invalid stage spec: {raw}")
        stage = {
            "precDps": int(parts[0]),
            "tol": parts[1],
            "constraintTol": parts[2],
            "maxIt": int(parts[3]),
            "stepMinExp": default_step_min_exp,
        }
        if len(parts) == 5 and parts[4]:
            stage["stepMinExp"] = int(parts[4])
        out.append(stage)
    if not out:
        raise ValueError("Empty --stages")
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


def project_sum_to_one(vals: list[Decimal]) -> list[Decimal]:
    n = Decimal(len(vals))
    diff = (Decimal(1) - sum(vals)) / n
    return [v + diff for v in vals]


def warm_start(prev: list[Decimal], p: int) -> list[Decimal]:
    if len(prev) >= p:
        vals = prev[:p]
    else:
        vals = prev + [Decimal(0)] * (p - len(prev))
    return project_sum_to_one(vals)


def asymptotic_seed(p: int) -> list[Decimal]:
    vals = []
    for n in range(p):
        term = Decimal((-1) ** n) / (Decimal(8) ** n * Decimal(n + 1))
        vals.append(term)
    s = sum(vals)
    vals = [v / s for v in vals]
    return project_sum_to_one(vals)


def load_bootstrap_records(path: Path) -> list[dict]:
    if not path.exists():
        return []
    return json.loads(path.read_text())


def pick_seed_from_records(records: list[dict], p: int) -> tuple[list[Decimal], str] | None:
    candidates: list[tuple[int, str, str]] = []
    for r in records:
        p_rec = r.get("P")
        ts = r.get("timestamp", "")
        a_csv = r.get("result", {}).get("aCsv", "")
        if p_rec is None or not a_csv:
            continue
        candidates.append((int(p_rec), ts, a_csv))
    if not candidates:
        return None
    candidates.sort(key=lambda x: (x[0], x[1]))
    below = [x for x in candidates if x[0] <= p]
    chosen = below[-1] if below else candidates[-1]
    coeffs = parse_csv_numbers(chosen[2])
    return warm_start(coeffs, p), f"bootstrap-record-P{chosen[0]}"


def run_stage(
    coeffs: list[Decimal],
    stage: dict,
    args: argparse.Namespace,
) -> dict[str, str]:
    coeffs_csv = format_coeffs(coeffs, max(stage["precDps"], args.print_dps))
    cmd = [
        str(RUNNER),
        "--coeffs", coeffs_csv,
        "--nmax", str(args.nmax),
        "--k", str(args.k),
        "--prec-dps", str(stage["precDps"]),
        "--print-dps", str(args.print_dps),
        "--tail", "true",
        "--optimize", "true",
        "--max-it", str(stage["maxIt"]),
        "--min-it", str(args.min_it),
        "--tol", str(stage["tol"]),
        "--constraint-tol", str(stage["constraintTol"]),
        "--damping", args.damping,
        "--step-min-exp", str(stage["stepMinExp"]),
        "--kkt-scaling", str(args.kkt_scaling),
        "--max-abs-a", str(args.max_abs_a),
        "--step-max-u", str(args.step_max_u),
        "--nonnegative-grid", str(args.nonnegative_grid),
        "--nonnegative-tol", str(args.nonnegative_tol),
        "--kkt-reg", str(args.kkt_reg),
        "--log", args.log,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"FLINT stage failed.\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    out = parse_kv_output(proc.stdout)
    out["_stdout"] = proc.stdout
    return out


def main() -> int:
    args = parse_args()
    if args.p <= 0:
        raise ValueError("--p must be positive")
    if args.nmax <= 0 or args.k <= 0:
        raise ValueError("--nmax and --k must be positive")

    stages = parse_stage_spec(args.stages, args.step_min_exp)
    max_prec = max(int(s["precDps"]) for s in stages)
    getcontext().prec = max(200, max_prec + 40)

    if args.start_coeffs.strip():
        current = project_sum_to_one(parse_csv_numbers(args.start_coeffs.strip()))
        init_source = "explicit-start"
    else:
        records = load_bootstrap_records(Path(args.bootstrap_records))
        pick = pick_seed_from_records(records, args.p)
        if pick is None:
            current = asymptotic_seed(args.p)
            init_source = "asymptotic-seed"
        else:
            current, init_source = pick

    print(
        f"FLINT precision ladder start: "
        f"P={args.p} N={args.nmax} K={args.k} seed={init_source} "
        f"stages={[s['precDps'] for s in stages]}"
    )

    rows = []
    for idx, stage in enumerate(stages):
        out = run_stage(current, stage, args)
        converged = out.get("converged", "false") == "true"
        kkt_ok = out.get("kktSolveOk", "true") == "true"

        row = {
            "idx": idx,
            "P": args.p,
            "Nmax": args.nmax,
            "K": args.k,
            "precDps": stage["precDps"],
            "tol": stage["tol"],
            "constraintTol": stage["constraintTol"],
            "maxIt": stage["maxIt"],
            "stepMinExp": stage["stepMinExp"],
            "converged": converged,
            "kktSolveOk": kkt_ok,
            "iterations": out.get("iterations", ""),
            "objective": out.get("objective", ""),
            "nu2Candidate": out.get("nu2Candidate", ""),
            "resInfUpper": out.get("resInfUpper", ""),
            "constraintAbsUpper": out.get("constraintAbsUpper", ""),
            "timeSec": out.get("timeSec", ""),
            "aCsv": out.get("aCsv", ""),
        }
        rows.append(row)

        print(
            f"[{idx:02d}] prec={stage['precDps']} tol={stage['tol']} ctol={stage['constraintTol']} "
            f"conv={converged} kkt={kkt_ok} iters={row['iterations']} "
            f"resInf={row['resInfUpper']} obj={row['objective']} time={row['timeSec']}"
        )

        if not row["aCsv"]:
            raise RuntimeError("Missing aCsv in stage output")
        current = parse_csv_numbers(row["aCsv"])

    result = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "config": {
            "P": args.p,
            "Nmax": args.nmax,
            "K": args.k,
            "stages": stages,
            "damping": args.damping,
            "minIt": args.min_it,
            "kktScaling": args.kkt_scaling,
            "maxAbsA": args.max_abs_a,
            "stepMaxU": args.step_max_u,
            "nonnegativeGrid": args.nonnegative_grid,
            "nonnegativeTol": args.nonnegative_tol,
            "kktReg": args.kkt_reg,
            "initSource": init_source,
        },
        "rows": rows,
        "final": rows[-1] if rows else None,
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if args.append and out_path.exists():
        prior = json.loads(out_path.read_text())
        if isinstance(prior, list):
            prior.append(result)
            out_path.write_text(json.dumps(prior, indent=2))
        else:
            out_path.write_text(json.dumps([prior, result], indent=2))
    else:
        out_path.write_text(json.dumps(result, indent=2))

    print(f"Wrote precision ladder report: {out_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
