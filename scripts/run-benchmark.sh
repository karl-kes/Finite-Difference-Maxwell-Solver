#!/usr/bin/env bash
set -e

# Build both benchmark binaries:
#   - `benchmark`        : OpenMP-enabled, used for the multi-thread run and
#                          the "1 thread OMP" column.
#   - `benchmark_serial` : core rebuilt without -fopenmp; true serial baseline.
mkdir -p build_release
cd build_release
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j"$(nproc)" benchmark benchmark_serial
cd ..

mkdir -p benchmark/results
STAMP="$(date +%Y%m%d_%H%M%S)"
HOST="$(hostname | tr -d '[:space:]')"
TAG="${HOST}_${STAMP}"

# Pin OMP threads to physical cores to cut variance. Override by exporting
# OMP_PLACES / OMP_PROC_BIND before invoking this script if you need a
# different policy.
export OMP_PLACES="${OMP_PLACES:-cores}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"

echo ""
echo "=== OpenMP benchmark (--threads = max) ==="
./build_release/benchmark --output "benchmark/results/omp_${TAG}.csv" "$@"

echo ""
echo "=== Serial benchmark (true single-thread, no -fopenmp) ==="
./build_release/benchmark_serial --output "benchmark/results/serial_${TAG}.csv" "$@"

echo ""
echo "Results written to benchmark/results/"
echo "  benchmark/results/omp_${TAG}.csv"
echo "  benchmark/results/serial_${TAG}.csv"