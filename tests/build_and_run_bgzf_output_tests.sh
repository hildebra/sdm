#!/usr/bin/env bash
set -euo pipefail

CXX=${CXX:-g++}
CXXFLAGS=${CXXFLAGS:--std=c++20 -O1 -g -Wall -Wextra}

SCRIPT_DIR=$(cd -- "$(dirname -- "$0")" && pwd)
REPO_DIR=$(cd -- "$SCRIPT_DIR/.." && pwd)

"$CXX" $CXXFLAGS -I"$REPO_DIR" \
    "$SCRIPT_DIR/TestBgzfOutput.cpp" \
    "$REPO_DIR/OutputFileStream.cpp" \
    "$REPO_DIR/InputFileStream.cpp" \
    "$REPO_DIR/Common.cpp" \
    -pthread -lz -o "$SCRIPT_DIR/TestBgzfOutput"

cd "$REPO_DIR"
"$SCRIPT_DIR/TestBgzfOutput"
