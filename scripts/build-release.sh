#!/usr/bin/env bash
set -e
mkdir -p build_release
cd build_release
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j"$(nproc)"
