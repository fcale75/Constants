#!/usr/bin/env python3
import subprocess
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
TEST_FILE = ROOT / "tests" / "runAll.wl"
SUMMARY_FILE = ROOT / "logs" / "latest_summary.txt"


def run_once():
    proc = subprocess.run(
        ["wolframscript", "-file", str(TEST_FILE)],
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )

    summary = {}
    if SUMMARY_FILE.exists():
        for line in SUMMARY_FILE.read_text(errors="ignore").splitlines():
            if "=" in line:
                k, v = line.split("=", 1)
                summary[k.strip()] = v.strip()

    return proc.returncode, summary, proc.stdout


def main():
    status, summary, stdout = run_once()
    print("=== stdout ===")
    print(stdout)
    print("=== summary ===")
    if summary:
        for k, v in summary.items():
            print(f"{k}={v}")
    else:
        print("(no summary file)")

    sys.exit(status)


if __name__ == "__main__":
    main()
