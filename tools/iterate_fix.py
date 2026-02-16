#!/usr/bin/env python3
"""
tools/iterate_fix.py

Small controller script to run the Wolfram Language test harness:

  wolframscript -file tests/runAll.wl

It prints:
  - the wolframscript combined stdout/stderr
  - the parsed key=value summary from logs/latest_summary.txt (if present)
"""

from __future__ import annotations

import pathlib
import subprocess
import sys
from typing import Dict, Tuple


ROOT = pathlib.Path(__file__).resolve().parents[1]
TEST_FILE = ROOT / "tests" / "runAll.wl"
SUMMARY_FILE = ROOT / "logs" / "latest_summary.txt"


def parse_summary(text: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        k, v = line.split("=", 1)
        out[k.strip()] = v.strip()
    return out


def run_once(timeout_s: int = 600) -> Tuple[int, Dict[str, str], str]:
    cmd = ["wolframscript", "-file", str(TEST_FILE)]
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(ROOT),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout_s,
            check=False,
        )
    except FileNotFoundError:
        msg = (
            "ERROR: wolframscript not found on PATH.\n"
            "Install WolframScript (or add it to PATH) and re-run.\n"
        )
        return 127, {}, msg
    except subprocess.TimeoutExpired as e:
        out = e.stdout or ""
        out += f"\n[timed out after {timeout_s}s]"
        return 124, {}, out

    summary: Dict[str, str] = {}
    if SUMMARY_FILE.exists():
        try:
            summary = parse_summary(SUMMARY_FILE.read_text(errors="ignore"))
        except Exception:
            summary = {}

    return proc.returncode, summary, proc.stdout


def main() -> int:
    status, summary, stdout = run_once()

    print("=== wolframscript stdout/stderr ===")
    print(stdout)

    print("=== parsed summary (logs/latest_summary.txt) ===")
    if summary:
        for k in sorted(summary):
            print(f"{k}={summary[k]}")
    else:
        print("(no summary file)")

    return status


if __name__ == "__main__":
    raise SystemExit(main())

