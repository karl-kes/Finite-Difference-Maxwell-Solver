#!/usr/bin/env bash
set -e
mkdir -p build-debug
cd build-debug
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=g++-15
make -j"$(nproc)"
