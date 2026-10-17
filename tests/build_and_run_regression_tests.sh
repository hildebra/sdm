#!/usr/bin/env bash
set -euo pipefail

CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:--std=c++20 -O1 -g -Wall -Wextra}

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)

sources=()
for source in "$REPO_DIR"/*.cpp; do
    [[ "$source" == "$REPO_DIR/Demultipl.cpp" ]] && continue
    sources+=("$source")
done

"$CXX" $CXXFLAGS -I"$REPO_DIR" \
    "$SCRIPT_DIR/TestRegression.cpp" "$SCRIPT_DIR/TestSupport.cpp" "${sources[@]}" \
    -pthread -lz -o "$SCRIPT_DIR/TestRegression"
cd "$REPO_DIR"
"$SCRIPT_DIR/TestRegression" "$@"
