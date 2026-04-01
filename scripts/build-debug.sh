#!/usr/bin/env bash
set -e
mkdir -p build_debug
cd build_debug
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j"$(nproc)"
