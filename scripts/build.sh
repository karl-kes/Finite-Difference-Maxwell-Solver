#!/usr/bin/env bash
set -e

PROJECT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_DIR"

JOBS=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

for MODE in Debug Test Release; do
    DIR="build_$(echo "$MODE" | tr '[:upper:]' '[:lower:]')"

    echo ""
    echo "================================"
    echo "  Building: $MODE -> $DIR"
    echo "================================"

    rm -rf "$DIR"
    cmake -B "$DIR" -DCMAKE_BUILD_TYPE="$MODE"
    cmake --build "$DIR" -j"$JOBS"

    if [ "$MODE" != "Release" ]; then
        echo ""
        echo "--- Running tests ($MODE) ---"
        ctest --test-dir "$DIR" --output-on-failure
    fi
done

echo ""
echo "================================"
echo "  All builds complete."
echo "  build_debug/    -> Debug"
echo "  build_test/     -> Test"
echo "  build_release/  -> Release"
echo "================================"