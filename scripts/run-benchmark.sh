#!/usr/bin/env bash
set -e
mkdir -p build_release
cd build_release
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j"$(nproc)" benchmark
cd ..

echo ""
echo "Running benchmark..."
echo ""

# Default: PML on, 100 steps, max grid 200^3
# Pass args through: ./scripts/run-benchmark.sh --steps 200 --max-n 256 --no-pml
./build_release/benchmark "$@"
