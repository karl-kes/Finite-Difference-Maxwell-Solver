#!/usr/bin/env bash
set -e
mkdir -p build_test
cd build_test
cmake .. -DCMAKE_BUILD_TYPE=Test
make -j"$(nproc)" run_tests
ctest --output-on-failure
