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
DEFAULT_OUT = ROOT / "checkpoints" / "flint_stage_calibration.json"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Calibrate FLINT stage settings by sweeping (N,K,prec) for fixed P "
            "and tracking objective stability."
        )
    )
    p.add_argument("--p", type=int, required=True, help="Target coefficient count P.")
    p.add_argument("--n-list", default="256,384,512", help="Comma list of N values.")
    p.add_argument("--k-list", default="32,48,64", help="Comma list of K values.")
    p.add_argument("--prec-list", default="120,140,170", help="Comma list of precision digits.")
    p.add_argument("--print-dps", type=int, default=70)
    p.add_argument("--max-it", type=int, default=30)
    p.add_argument("--tol", default="1e-26")
    p.add_argument("--constraint-tol", default="1e-26")
    p.add_argument("--damping", choices=["true", "false"], default="true")
    p.add_argument("--step-min-exp", type=int, default=20)
    p.add_argument("--log", choices=["true", "false"], default="false")
    p.add_argument(
        "--stability-threshold",
        default="1e-12",
        help="Absolute objective drift threshold for declaring stability.",
    )
    p.add_argument(
        "--require-consecutive",
        type=int,
        default=2,
        help="Consecutive stable drifts required to mark recommended setting.",
    )
    p.add_argument(
        "--stability-mode",
        choices=["by-n", "sequential"],
        default="by-n",
        help=(
            "by-n: compare objective drift only across increasing N for fixed (K,prec). "
            "sequential: compare against immediately previous grid point."
        ),
    )
    p.add_argument("--bootstrap-records", default=str(DEFAULT_BOOTSTRAP_RECORDS))
    p.add_argument("--start-coeffs", default="", help="Optional CSV to override initial coefficients.")
    p.add_argument("--out", default=str(DEFAULT_OUT))
    p.add_argument("--append", action="store_true")
    return p.parse_args()


def parse_csv_ints(s: str) -> list[int]:
    out = [int(x.strip()) for x in s.split(",") if x.strip()]
    if not out:
        raise ValueError("Empty list")
    return out


def parse_kv_output(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in text.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def parse_csv_decimals(csv: str) -> list[Decimal]:
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


def choose_seed_from_records(records: list[dict], p: int) -> tuple[list[Decimal], str] | None:
    candidates = []
    for r in records:
        rr = r.get("result", {})
        params = r.get("params", {})
        p_rec = r.get("P")
        a_csv = rr.get("aCsv", "")
        if p_rec is None or not a_csv:
            continue
        candidates.append((int(p_rec), r))
    if not candidates:
        return None
    candidates.sort(key=lambda x: (x[0], x[1].get("timestamp", "")))

    # prefer highest P <= target, else highest available P.
    below = [x for x in candidates if x[0] <= p]
    chosen = below[-1] if below else candidates[-1]
    p_rec, rec = chosen
    coeffs = parse_csv_decimals(rec["result"]["aCsv"])
    seed = warm_start(coeffs, p)
    return seed, f"bootstrap-record-P{p_rec}"


def run_stage(
    coeffs: list[Decimal],
    p: int,
    n: int,
    k: int,
    prec: int,
    args: argparse.Namespace,
) -> dict[str, str]:
    coeffs_csv = format_coeffs(coeffs, max(prec, args.print_dps))
    cmd = [
        str(RUNNER),
        "--coeffs", coeffs_csv,
        "--nmax", str(n),
        "--k", str(k),
        "--prec-dps", str(prec),
        "--print-dps", str(args.print_dps),
        "--tail", "true",
        "--optimize", "true",
        "--max-it", str(args.max_it),
        "--tol", str(args.tol),
        "--constraint-tol", str(args.constraint_tol),
        "--damping", args.damping,
        "--step-min-exp", str(args.step_min_exp),
        "--log", args.log,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"FLINT stage failed for P={p} N={n} K={k} prec={prec}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )
    out = parse_kv_output(proc.stdout)
    out["_stdout"] = proc.stdout
    return out


def combos(
    n_list: list[int],
    k_list: list[int],
    prec_list: list[int],
    stability_mode: str,
) -> list[tuple[int, int, int]]:
    if stability_mode == "by-n":
        # Keep (K, prec) fixed while N increases so stability checks are meaningful.
        return [
            (n, k, prec)
            for prec in sorted(prec_list)
            for k in sorted(k_list)
            for n in sorted(n_list)
        ]

    grid = [(n, k, prec) for prec in sorted(prec_list) for n in sorted(n_list) for k in sorted(k_list)]
    # Monotone by increasing computational cost.
    grid.sort(key=lambda t: (t[2], t[0], t[1], t[0] * t[1]))
    return grid


def main() -> int:
    args = parse_args()
    if args.p <= 0:
        raise ValueError("--p must be positive")

    n_list = parse_csv_ints(args.n_list)
    k_list = parse_csv_ints(args.k_list)
    prec_list = parse_csv_ints(args.prec_list)
    stability_threshold = Decimal(args.stability_threshold)
    require_consecutive = int(args.require_consecutive)
    if require_consecutive <= 0:
        raise ValueError("--require-consecutive must be positive")

    max_prec = max(prec_list)
    getcontext().prec = max(200, max_prec + 40)

    if args.start_coeffs.strip():
        current = project_sum_to_one(parse_csv_decimals(args.start_coeffs.strip()))
        init_source = "explicit-start"
    else:
        recs = load_bootstrap_records(Path(args.bootstrap_records))
        pick = choose_seed_from_records(recs, args.p)
        if pick is None:
            current = asymptotic_seed(args.p)
            init_source = "asymptotic-seed"
        else:
            current, init_source = pick

    print(
        f"FLINT calibration start: P={args.p}, seed={init_source}, "
        f"N={sorted(n_list)}, K={sorted(k_list)}, prec={sorted(prec_list)}"
    )

    rows = []
    seq_last_obj: Decimal | None = None
    seq_stable_streak = 0
    by_n_last_obj: dict[tuple[int, int], Decimal] = {}
    by_n_stable_streak: dict[tuple[int, int], int] = {}
    recommended_idx: int | None = None

    for idx, (n, k, prec) in enumerate(combos(n_list, k_list, prec_list, args.stability_mode)):
        out = run_stage(current, args.p, n, k, prec, args)
        converged = out.get("converged", "false") == "true"
        kkt_ok = out.get("kktSolveOk", "true") == "true"
        obj = Decimal(out.get("objective", "nan"))

        if args.stability_mode == "by-n":
            key = (k, prec)
            prev_obj = by_n_last_obj.get(key)
            drift = None if prev_obj is None else abs(obj - prev_obj)
            stable_streak = by_n_stable_streak.get(key, 0)
            if drift is not None and drift <= stability_threshold and converged and kkt_ok:
                stable_streak += 1
            else:
                stable_streak = 0
            by_n_stable_streak[key] = stable_streak
            by_n_last_obj[key] = obj
            drift_scope = f"byN(K={k},prec={prec})"
        else:
            drift = None if seq_last_obj is None else abs(obj - seq_last_obj)
            if drift is not None and drift <= stability_threshold and converged and kkt_ok:
                seq_stable_streak += 1
            else:
                seq_stable_streak = 0
            stable_streak = seq_stable_streak
            seq_last_obj = obj
            drift_scope = "sequential"

        if recommended_idx is None and stable_streak >= require_consecutive:
            recommended_idx = idx

        row = {
            "idx": idx,
            "P": args.p,
            "Nmax": n,
            "K": k,
            "precDps": prec,
            "converged": converged,
            "kktSolveOk": kkt_ok,
            "iterations": out.get("iterations", ""),
            "objective": str(obj),
            "nu2Candidate": out.get("nu2Candidate", ""),
            "resInfUpper": out.get("resInfUpper", ""),
            "constraintAbsUpper": out.get("constraintAbsUpper", ""),
            "timeSec": out.get("timeSec", ""),
            "stabilityMode": args.stability_mode,
            "driftScope": drift_scope,
            "driftFromPrev": str(drift) if drift is not None else "",
            "stableStreak": stable_streak,
        }
        rows.append(row)

        print(
            f"[{idx:02d}] P={args.p} N={n} K={k} prec={prec} "
            f"conv={converged} kkt={kkt_ok} obj={row['objective']} "
            f"drift={row['driftFromPrev'] or 'NA'} ({row['driftScope']}) "
            f"time={row['timeSec']}"
        )

        a_csv = out.get("aCsv", "")
        if a_csv:
            current = parse_csv_decimals(a_csv)

    result = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "config": {
            "P": args.p,
            "nList": sorted(n_list),
            "kList": sorted(k_list),
            "precList": sorted(prec_list),
            "stabilityThreshold": str(stability_threshold),
            "requireConsecutive": require_consecutive,
            "stabilityMode": args.stability_mode,
            "initSource": init_source,
        },
        "rows": rows,
        "recommendedIndex": recommended_idx,
        "recommended": rows[recommended_idx] if recommended_idx is not None else None,
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

    if recommended_idx is not None:
        rec = rows[recommended_idx]
        print(
            "Recommended settings: "
            f"N={rec['Nmax']} K={rec['K']} prec={rec['precDps']} "
            f"objective={rec['objective']}"
        )
    else:
        print("No stable recommendation found; expand N/K/prec grid.")

    print(f"Wrote calibration report: {out_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(1)
