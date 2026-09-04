#!/usr/bin/env bash
set -e
mkdir -p build-release
cd build-release
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15
make -j"$(nproc)"
