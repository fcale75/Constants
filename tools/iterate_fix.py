#!/usr/bin/env python3
"""
tools/iterate_fix.py

Runs:
  wolframscript -file tests/runAll.wl

Features:
  - Streams WL output live (no buffering)
  - Enforces a real timeout
  - Force-kills WL on timeout
  - Prints parsed logs/latest_summary.txt if present
"""

from __future__ import annotations

import argparse
import pathlib
import subprocess
import sys
import time
from typing import Dict


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


def run_with_timeout(timeout_s: int) -> int:
    cmd = ["wolframscript", "-file", str(TEST_FILE)]

    try:
        proc = subprocess.Popen(
            cmd,
            cwd=str(ROOT),
            stdout=sys.stdout,   # live streaming
            stderr=sys.stderr,
            text=True,
        )
    except FileNotFoundError:
        print("ERROR: wolframscript not found on PATH.")
        return 127

    deadline = time.time() + timeout_s

    while True:
        rc = proc.poll()
        if rc is not None:
            return rc

        if time.time() >= deadline:
            print(f"\nERROR: timed out after {timeout_s} seconds. Terminating...")
            try:
                proc.terminate()
                time.sleep(1.0)
                if proc.poll() is None:
                    proc.kill()
            except Exception:
                pass
            return 124

        time.sleep(0.2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--timeout",
        type=int,
        default=120,
        help="Timeout in seconds (default 120)"
    )
    args = parser.parse_args()

    rc = run_with_timeout(args.timeout)

    print("\n=== parsed summary (logs/latest_summary.txt) ===")

    if SUMMARY_FILE.exists():
        try:
            summary = parse_summary(SUMMARY_FILE.read_text(errors="ignore"))
        except Exception:
            summary = {}
        if summary:
            for k in sorted(summary):
                print(f"{k}={summary[k]}")
        else:
            print("(summary file exists but could not parse)")
    else:
        print("(no summary file)")

    return rc


if __name__ == "__main__":
    raise SystemExit(main())

