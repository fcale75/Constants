# Session Handoff (FLINT Quality-Gated Runs)

Last updated: 2026-02-17

## Canonical repo state

- Commit: `177b164`
- Branch: `main`
- Remote: `origin/main` (pushed)

Essential files committed in this handoff commit:
- `/Users/fcale/Dropbox/ChatGPT/Constants/tools/run_flint_quality_chain.py`
- `/Users/fcale/Dropbox/ChatGPT/Constants/checkpoints/flint_quality_p32_gate.json`

## What was done

1. Added a new gated runner:
   - `tools/run_flint_quality_chain.py`
   - For each stage P, it runs:
     - a base solve
     - a stricter check solve
   - It computes achieved digits from drift in `nu2Candidate` and coefficient `L_inf` drift, then only passes a stage if thresholds are met.

2. Ran gated `P=32` attempts (then stopped manually on user request):
   - Checkpoint file:
     - `checkpoints/flint_quality_p32_gate.json`

## Current numerical status (`P=32` gate)

From `checkpoints/flint_quality_p32_gate.json`:

- Attempt 0:
  - Base: iterations `24`, `resInfUpper=2.53925e-13`
  - Check: iterations `34`, `resInfUpper=2.37495e-14`
  - Achieved digits:
    - `nu2Digits=14.565823783224188`
    - `coeffDigits=6.850665444089551`
  - Target was `nu2 >= 40`, coeff digits `>= 32` -> failed

- Attempt 1:
  - Base: iterations `30`, `resInfUpper=4.05241e-14`
  - Check: iterations `40`, `resInfUpper=1.80856e-14`
  - Achieved digits:
    - `nu2Digits=15.1417631995922`
    - `coeffDigits=7.366841963042819`
  - Target still failed

Stage status in checkpoint: `failed`  
Overall status in checkpoint: `failed`

## Why runtimes increased vs earlier short runs

- Earlier exploratory runs often stopped in very few iterations (e.g. `3/3/5`), which was suspicious.
- The gated runs intentionally enforce harder settings:
  - larger `N`, `K`, precision
  - larger `max-it`
  - stricter check pass per attempt
- This causes many more Newton iterations and much longer runtime per attempt.

## Process state at handoff

- All active optimization processes were stopped before handoff commit.
- No active `constants_objective` runs remain.

## Safe multi-machine workflow (important)

- Use one active writer at a time (especially in Dropbox-backed folders).
- Before switching machines:
  1. stop runs
  2. wait for Dropbox sync
  3. commit/push
  4. pull on the other machine

## Resume commands

From repo root:

```bash
cd /Users/fcale/Dropbox/ChatGPT/Constants
git fetch origin
git checkout main
git pull --ff-only
```

Re-run the same gated `P=32` workflow:

```bash
python3 tools/run_flint_quality_chain.py \
  --stage-list "32:8192:128:340:24:1e-20:1e-28:40" \
  --seed-records checkpoints/flint_bootstrap_records_p32_N16384_guard2.json \
  --seed-p 32 \
  --expected-digits-map "32:40" \
  --check-n-delta 8192 \
  --check-k-delta 64 \
  --check-prec-delta 60 \
  --check-max-it-delta 10 \
  --check-tol-decades 4 \
  --check-constraint-tol-decades 4 \
  --max-attempts 4 \
  --escalate-n-mult 1.5 \
  --escalate-k-add 16 \
  --escalate-prec-add 20 \
  --escalate-max-it-add 6 \
  --escalate-step-min-exp-add 2 \
  --escalate-tol-decades 1 \
  --escalate-constraint-tol-decades 1 \
  --kkt-reg 1e-22 \
  --reg-decay 0.1 \
  --reg-min 1e-30 \
  --out checkpoints/flint_quality_p32_gate.json
```

If you want persistent console logs on rerun:

```bash
python3 tools/run_flint_quality_chain.py ... | tee logs/flint_quality_p32_gate.log
```

## Known repo noise (not part of this commit)

This repo has many pre-existing untracked checkpoint files and unrelated local changes.  
Handoff commit intentionally included only the essential runner + checkpoint snapshot.
