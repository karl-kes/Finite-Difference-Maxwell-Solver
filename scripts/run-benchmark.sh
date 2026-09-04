#!/usr/bin/env bash
set -e

# Build both benchmark binaries:
#   - `benchmark`        : OpenMP-enabled, used for the multi-thread run and
#                          the "1 thread OMP" column.
#   - `benchmark-serial` : core rebuilt without -fopenmp; true serial baseline.
mkdir -p build-release
cd build-release
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15
make -j"$(nproc)" benchmark benchmark-serial
cd ..

mkdir -p benchmark/results
STAMP="$(date +%Y%m%d-%H%M%S)"
HOST="$(hostname | tr -d '[:space:]')"
TAG="${HOST}-${STAMP}"

# Pin OMP threads to physical cores to cut variance. Override by exporting
# OMP_PLACES / OMP_PROC_BIND before invoking this script if you need a
# different policy.
export OMP_PLACES="${OMP_PLACES:-cores}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"

echo ""
echo "=== OpenMP benchmark (--threads = max) ==="
./build-release/benchmark --output "benchmark/results/omp-${TAG}.csv" "$@"

echo ""
echo "=== Serial benchmark (true single-thread, no -fopenmp) ==="
./build-release/benchmark-serial --output "benchmark/results/serial-${TAG}.csv" "$@"

echo ""
echo "Results written to benchmark/results/"
echo "  benchmark/results/omp-${TAG}.csv"
echo "  benchmark/results/serial-${TAG}.csv"
