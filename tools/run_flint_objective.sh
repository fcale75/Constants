#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/fcale/Dropbox/ChatGPT/Constants"
BIN="${ROOT}/cpp/build/constants_objective"

if [[ ! -x "${BIN}" ]]; then
  "${ROOT}/tools/build_flint_cpp.sh"
fi

"${BIN}" "$@"
