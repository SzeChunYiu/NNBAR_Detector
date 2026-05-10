#!/bin/bash
# ============================================================================
# NNBAR Detector Simulation - Build Script
# Automatically configures and builds with bundled dependencies
# ============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
EXTERNAL_DIR="$PROJECT_DIR/external"
BUILD_DIR="$PROJECT_DIR/build"
OS_NAME="$(uname -s)"

echo "============================================="
echo "NNBAR Build Script"
echo "============================================="

# Load local package paths when available.
# shellcheck disable=SC1091
source "$SCRIPT_DIR/setup_environment.sh"

# Check if dependencies are set up. On macOS, do not use checked-in Linux ELF
# bundles; prefer active/system package configs such as conda's nnbar_env.
if [ -z "${Arrow_DIR:-}" ] || [ -z "${Parquet_DIR:-}" ] || [ -z "${spdlog_DIR:-}" ]; then
    if [ "$OS_NAME" != "Darwin" ] && [ -d "$EXTERNAL_DIR/arrow-install-linux" ]; then
        :
    elif [ ! -d "$EXTERNAL_DIR/arrow-install" ]; then
        echo "Dependencies not found. Running setup script first..."
        bash "$SCRIPT_DIR/setup_dependencies.sh"
    fi
fi

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# If this checkout moved, an old CMake cache points at the previous source tree
# and CMake refuses to reconfigure. Keep the stale files for inspection and let
# this build directory regenerate cleanly.
if [ -f CMakeCache.txt ]; then
    CACHED_SOURCE="$(grep '^CMAKE_HOME_DIRECTORY:INTERNAL=' CMakeCache.txt | cut -d= -f2- || true)"
    if [ -n "$CACHED_SOURCE" ] && [ "$CACHED_SOURCE" != "$PROJECT_DIR" ]; then
        STALE_CACHE_DIR="$BUILD_DIR/stale-cmake-cache-$(date +%Y%m%d%H%M%S)"
        mkdir -p "$STALE_CACHE_DIR"
        mv CMakeCache.txt "$STALE_CACHE_DIR/"
        if [ -d CMakeFiles ]; then
            mv CMakeFiles "$STALE_CACHE_DIR/"
        fi
        echo ">>> Moved stale CMake cache from $CACHED_SOURCE to $STALE_CACHE_DIR"
    fi
fi

# Configure with bundled dependencies
echo ""
echo ">>> Configuring CMake..."

CMAKE_ARGS=(
    -DCMAKE_BUILD_TYPE=Release
    -DMCPL_BUILD=OFF
    -DTARGET_BUILD=ON
    -DWITH_SCINTILLATION=OFF
    -DWITH_GARFIELD=OFF
    -DWITH_GARFIELD_GPU=OFF
    -DWITH_CELERITAS=OFF
    -DWITH_OPTICKS=OFF
    -DWITH_DASHBOARD=OFF
)

# Add bundled dependency paths if they exist
if [ -d "$EXTERNAL_DIR/arrow-install" ]; then
    CMAKE_ARGS+=(-DArrow_DIR="$EXTERNAL_DIR/arrow-install/lib/cmake/Arrow")
    CMAKE_ARGS+=(-DParquet_DIR="$EXTERNAL_DIR/arrow-install/lib/cmake/Parquet")
elif [ "$OS_NAME" != "Darwin" ] && [ -d "$EXTERNAL_DIR/arrow-install-linux" ]; then
    CMAKE_ARGS+=(-DArrow_DIR="$EXTERNAL_DIR/arrow-install-linux/lib/cmake/Arrow")
    CMAKE_ARGS+=(-DParquet_DIR="$EXTERNAL_DIR/arrow-install-linux/lib/cmake/Parquet")
elif [ -n "${Arrow_DIR:-}" ] && [ -n "${Parquet_DIR:-}" ]; then
    CMAKE_ARGS+=(-DArrow_DIR="$Arrow_DIR")
    CMAKE_ARGS+=(-DParquet_DIR="$Parquet_DIR")
fi

if [ -d "$EXTERNAL_DIR/spdlog-install" ]; then
    CMAKE_ARGS+=(-Dspdlog_DIR="$EXTERNAL_DIR/spdlog-install/lib/cmake/spdlog")
elif [ "$OS_NAME" != "Darwin" ] && [ -d "$EXTERNAL_DIR/spdlog-install-linux" ]; then
    CMAKE_ARGS+=(-Dspdlog_DIR="$EXTERNAL_DIR/spdlog-install-linux/lib/cmake/spdlog")
elif [ -n "${spdlog_DIR:-}" ]; then
    CMAKE_ARGS+=(-Dspdlog_DIR="$spdlog_DIR")
fi

if [ -d "$EXTERNAL_DIR/json" ]; then
    CMAKE_ARGS+=(-Dnlohmann_json_DIR="$EXTERNAL_DIR/json")
fi

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --scintillation)
            CMAKE_ARGS+=(-DWITH_SCINTILLATION=ON)
            echo ">>> Enabling scintillation (optical photon generation)"
            shift
            ;;
        --garfield)
            CMAKE_ARGS+=(-DWITH_GARFIELD=ON)
            echo ">>> Enabling Garfield++ TPC simulation"
            shift
            ;;
        --opticks)
            CMAKE_ARGS+=(-DWITH_OPTICKS=ON)
            echo ">>> Enabling Opticks GPU optical photon propagation"
            shift
            ;;
        --debug)
            CMAKE_ARGS+=(-DDEBUG_VERBOSE=ON)
            echo ">>> Enabling debug output"
            shift
            ;;
        --mcpl)
            CMAKE_ARGS+=(-DMCPL_BUILD=ON)
            echo ">>> Enabling MCPL particle source"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            echo "Usage: $0 [--scintillation] [--garfield] [--opticks] [--debug] [--mcpl]"
            exit 1
            ;;
    esac
done

# Run CMake
cmake "${CMAKE_ARGS[@]}" ..

# Build
echo ""
echo ">>> Building..."
if command -v nproc >/dev/null 2>&1; then
    JOBS="$(nproc)"
else
    JOBS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
fi
make -j"$JOBS"

echo ""
echo "============================================="
echo "Build complete!"
echo "Executable: $BUILD_DIR/nnbar-detector-simulation"
echo "============================================="
echo ""
echo "To run a smoke test: ./nnbar-detector-simulation -m macro/quick_test.mac"
echo "To use an MCPL source: ./scripts/build.sh --mcpl"
echo ""
