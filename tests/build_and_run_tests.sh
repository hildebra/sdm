#!/usr/bin/env bash
# Simple script to compile and run the tests locally using g++.
# Usage: cd sdm && ./tests/build_and_run_tests.sh

set -euo pipefail

# Attempt to find g++
CXX=${CXX:-g++}

SRC_DIR="$(dirname "$0")/.."
TEST_DIR="$(dirname "$0")"
OUT=test_input_helpers_runner

echo "Compiling tests..."
SOURCES=()
for src in "$SRC_DIR"/*.cpp; do
  [[ "$(basename "$src")" == "Demultipl.cpp" ]] || SOURCES+=("$src")
done
$CXX -std=c++20 -I"$SRC_DIR" -O2 \
  "$TEST_DIR/TestInputHelpers.cpp" "$TEST_DIR/TestSupport.cpp" \
  "${SOURCES[@]}" -pthread -lz -o "$TEST_DIR/$OUT"

echo "Running tests..."
"$TEST_DIR/$OUT"

echo "Done."
