#!/usr/bin/env bash
set -e

# Build both benchmark binaries:
#   - `kestrel-maxwell-benchmark`        : OpenMP-enabled.
#   - `kestrel-maxwell-benchmark-serial` : core rebuilt without -fopenmp.
mkdir -p build-release
cd build-release
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15
make -j"$(nproc)" kestrel-maxwell-benchmark kestrel-maxwell-benchmark-serial
cd ..

mkdir -p examples/maxwell/benchmark/results
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
./build-release/examples/maxwell/kestrel-maxwell-benchmark \
    --output "examples/maxwell/benchmark/results/omp-${TAG}.csv" "$@"

echo ""
echo "=== Serial benchmark (true single-thread, no -fopenmp) ==="
./build-release/examples/maxwell/kestrel-maxwell-benchmark-serial \
    --output "examples/maxwell/benchmark/results/serial-${TAG}.csv" "$@"

echo ""
echo "Results written to examples/maxwell/benchmark/results/"
echo "  examples/maxwell/benchmark/results/omp-${TAG}.csv"
echo "  examples/maxwell/benchmark/results/serial-${TAG}.csv"
