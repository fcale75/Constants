#!/usr/bin/env bash
set -euo pipefail

ROOT="/Users/fcale/Dropbox/ChatGPT/Constants"
BUILD_DIR="${ROOT}/cpp/build"
SRC="${ROOT}/cpp/src/constants_objective.cpp"
BIN="${BUILD_DIR}/constants_objective"

mkdir -p "${BUILD_DIR}"

if command -v cmake >/dev/null 2>&1; then
  cmake -S "${ROOT}/cpp" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
  cmake --build "${BUILD_DIR}" --config Release -j
else
  CXX="${CXX:-clang++}"
  "${CXX}" \
    -std=c++17 -O3 -DNDEBUG \
    -I/opt/homebrew/opt/flint/include \
    -I/opt/homebrew/include \
    "${SRC}" \
    -L/opt/homebrew/opt/flint/lib \
    -L/opt/homebrew/lib \
    -lflint -lmpfr -lgmp -lm \
    -Wl,-rpath,/opt/homebrew/opt/flint/lib \
    -Wl,-rpath,/opt/homebrew/lib \
    -o "${BIN}"
fi

echo "Built: ${BIN}"
