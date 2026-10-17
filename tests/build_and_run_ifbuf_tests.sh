#!/usr/bin/env bash
set -euo pipefail
CXX=${CXX:-g++}
SRC_DIR="$(dirname "$0")/.."
TEST_DIR="$(dirname "$0")"
OUT=test_ifbuf_runner

mkdir -p "$TEST_DIR"

echo "Compiling ifbufstream tests..."
SOURCES=()
for src in "$SRC_DIR"/*.cpp; do
  [[ "$(basename "$src")" == "Demultipl.cpp" ]] || SOURCES+=("$src")
done
$CXX -std=c++20 -O2 -I"$SRC_DIR" \
  "$TEST_DIR/TestIfBufStream.cpp" "$TEST_DIR/TestSupport.cpp" \
  "${SOURCES[@]}" -pthread -lz -o "$TEST_DIR/$OUT"

echo "Running tests..."
"$TEST_DIR/$OUT"

echo "Done." 
