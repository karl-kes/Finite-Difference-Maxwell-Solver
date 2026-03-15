#!/usr/bin/env bash
set -euo pipefail

MODE="${1:-release}"
BUILD_DIR="build"
JOBS="$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)"

case "${MODE,,}" in
    release|r)
        BUILD_TYPE="Release"
        ;;
    debug|d)
        BUILD_TYPE="Debug"
        ;;
    clean)
        echo "Cleaning build directory..."
        rm -rf "$BUILD_DIR"
        echo "Done."
        exit 0
        ;;
    *)
        echo "Usage: $0 [release|debug|clean]"
        exit 1
        ;;
esac

echo "=== $BUILD_TYPE Build ==="

# Configure (only if needed or build type changed):
if [ ! -f "$BUILD_DIR/CMakeCache.txt" ] || \
   ! grep -q "CMAKE_BUILD_TYPE:STRING=$BUILD_TYPE" "$BUILD_DIR/CMakeCache.txt" 2>/dev/null; then
    echo "Configuring ($BUILD_TYPE)..."
    rm -rf "$BUILD_DIR"
    cmake -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="$BUILD_TYPE"
fi

# Build:
echo "Building with $JOBS threads..."
cmake --build "$BUILD_DIR" -j "$JOBS"

echo ""
echo "=== Build Complete ==="
echo "Binary: $BUILD_DIR/main"

# Run:
echo ""
"./$BUILD_DIR/main"