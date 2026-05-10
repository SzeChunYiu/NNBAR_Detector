#!/bin/bash
# ============================================================================
# setup_dependencies.sh - Download and setup all external dependencies
# ============================================================================
# This script downloads and installs all required dependencies for the
# NNBAR detector simulation to the external/ directory.
#
# Dependencies installed:
#   - Apache Arrow/Parquet (columnar data format)
#   - nlohmann_json (header-only JSON library)
#   - spdlog (fast logging library)
#
# MCPL and parquet-writer are already bundled in the repository.
#
# Usage: ./scripts/setup_dependencies.sh [options]
# Options:
#   --clean       Remove all dependencies and start fresh
#   --parallel N  Use N parallel jobs for compilation (default: nproc)
#   --skip-arrow  Skip Arrow build (if already installed system-wide)
# ============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# Versions
ARROW_VERSION="17.0.0"
JSON_VERSION="3.11.3"
SPDLOG_VERSION="1.12.0"

# Directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
EXTERNAL_DIR="${PROJECT_DIR}/external"

# Parse arguments
CLEAN=false
PARALLEL=$(nproc 2>/dev/null || echo 4)
SKIP_ARROW=false

while [[ $# -gt 0 ]]; do
    case $1 in
        --clean)
            CLEAN=true
            shift
            ;;
        --parallel)
            PARALLEL="$2"
            shift 2
            ;;
        --skip-arrow)
            SKIP_ARROW=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --clean       Remove all dependencies and start fresh"
            echo "  --parallel N  Use N parallel jobs (default: $(nproc 2>/dev/null || echo 4))"
            echo "  --skip-arrow  Skip Arrow build (if already installed system-wide)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Banner
echo ""
echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║${NC}${BOLD}         NNBAR Detector Simulation - Dependency Setup                      ${NC}${CYAN}║${NC}"
echo -e "${CYAN}║${NC}         Bundling all required libraries for portable builds               ${CYAN}║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${BLUE}Dependencies to install:${NC}"
echo -e "  • Apache Arrow/Parquet ${GREEN}v${ARROW_VERSION}${NC} - Columnar data format"
echo -e "  • nlohmann_json        ${GREEN}v${JSON_VERSION}${NC}  - JSON parsing (header-only)"
echo -e "  • spdlog               ${GREEN}v${SPDLOG_VERSION}${NC} - Fast logging"
echo ""
echo -e "${BLUE}Already bundled:${NC}"
echo -e "  • MCPL                 ${GREEN}bundled${NC}    - Monte Carlo Particle List"
echo -e "  • parquet-writer       ${GREEN}bundled${NC}    - Parquet output helper"
echo ""

mkdir -p "${EXTERNAL_DIR}"
cd "${EXTERNAL_DIR}"

# ============================================================================
# nlohmann_json (Header-only - just download)
# ============================================================================
echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}[1/3] nlohmann_json (Header-only JSON library)${NC}"
echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

JSON_DIR="${EXTERNAL_DIR}/json"

if [[ "$CLEAN" == "true" ]] && [[ -d "${JSON_DIR}" ]]; then
    echo -e "  Cleaning previous installation..."
    rm -rf "${JSON_DIR}"
fi

if [[ -f "${JSON_DIR}/include/nlohmann/json.hpp" ]]; then
    echo -e "  ${GREEN}✓${NC} nlohmann_json already installed"
else
    echo -e "  Downloading nlohmann_json v${JSON_VERSION}..."

    mkdir -p "${JSON_DIR}/include/nlohmann"

    if command -v wget &> /dev/null; then
        wget -q --show-progress -O "${JSON_DIR}/include/nlohmann/json.hpp" \
            "https://github.com/nlohmann/json/releases/download/v${JSON_VERSION}/json.hpp"
    else
        curl -L -# -o "${JSON_DIR}/include/nlohmann/json.hpp" \
            "https://github.com/nlohmann/json/releases/download/v${JSON_VERSION}/json.hpp"
    fi

    echo -e "  ${GREEN}✓${NC} nlohmann_json installed to ${JSON_DIR}"
fi
echo ""

# ============================================================================
# spdlog (Fast logging library)
# ============================================================================
echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}[2/3] spdlog (Fast logging library)${NC}"
echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

SPDLOG_INSTALL="${EXTERNAL_DIR}/spdlog-install"
SPDLOG_SRC="${EXTERNAL_DIR}/spdlog-src"
SPDLOG_BUILD="${EXTERNAL_DIR}/spdlog-build"

if [[ "$CLEAN" == "true" ]]; then
    rm -rf "${SPDLOG_INSTALL}" "${SPDLOG_SRC}" "${SPDLOG_BUILD}"
fi

if [[ -f "${SPDLOG_INSTALL}/lib/cmake/spdlog/spdlogConfig.cmake" ]]; then
    echo -e "  ${GREEN}✓${NC} spdlog already installed"
else
    echo -e "  Downloading spdlog v${SPDLOG_VERSION}..."

    SPDLOG_URL="https://github.com/gabime/spdlog/archive/refs/tags/v${SPDLOG_VERSION}.tar.gz"

    if command -v wget &> /dev/null; then
        wget -q --show-progress -O "spdlog-${SPDLOG_VERSION}.tar.gz" "${SPDLOG_URL}"
    else
        curl -L -# -o "spdlog-${SPDLOG_VERSION}.tar.gz" "${SPDLOG_URL}"
    fi

    echo -e "  Extracting..."
    tar -xzf "spdlog-${SPDLOG_VERSION}.tar.gz"
    mv "spdlog-${SPDLOG_VERSION}" "${SPDLOG_SRC}"
    rm "spdlog-${SPDLOG_VERSION}.tar.gz"

    echo -e "  Configuring..."
    mkdir -p "${SPDLOG_BUILD}"
    cd "${SPDLOG_BUILD}"

    cmake "${SPDLOG_SRC}" \
        -DCMAKE_INSTALL_PREFIX="${SPDLOG_INSTALL}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DSPDLOG_BUILD_EXAMPLE=OFF \
        -DSPDLOG_BUILD_TESTS=OFF \
        -DSPDLOG_BUILD_SHARED=ON \
        > cmake_output.log 2>&1

    echo -e "  Building..."
    make -j${PARALLEL} > build_output.log 2>&1

    echo -e "  Installing..."
    make install > install_output.log 2>&1

    cd "${EXTERNAL_DIR}"
    echo -e "  ${GREEN}✓${NC} spdlog installed to ${SPDLOG_INSTALL}"

    # Cleanup source and build
    rm -rf "${SPDLOG_SRC}" "${SPDLOG_BUILD}"
fi
echo ""

# ============================================================================
# Apache Arrow/Parquet
# ============================================================================
echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}[3/3] Apache Arrow/Parquet (Columnar data format)${NC}"
echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

ARROW_INSTALL="${EXTERNAL_DIR}/arrow-install"
ARROW_SRC="${EXTERNAL_DIR}/arrow-src"
ARROW_BUILD="${EXTERNAL_DIR}/arrow-build"

if [[ "$SKIP_ARROW" == "true" ]]; then
    echo -e "  ${YELLOW}○${NC} Skipped (--skip-arrow flag)"
elif [[ -f "${ARROW_INSTALL}/lib/cmake/Arrow/ArrowConfig.cmake" ]] && [[ "$CLEAN" == "false" ]]; then
    echo -e "  ${GREEN}✓${NC} Arrow already installed"
else
    if [[ "$CLEAN" == "true" ]]; then
        rm -rf "${ARROW_INSTALL}" "${ARROW_SRC}" "${ARROW_BUILD}"
    fi

    echo -e "  ${YELLOW}Note:${NC} Arrow build takes 5-15 minutes"
    echo ""

    # Download
    echo -e "  Downloading Apache Arrow v${ARROW_VERSION}..."
    ARROW_URL="https://github.com/apache/arrow/archive/refs/tags/apache-arrow-${ARROW_VERSION}.tar.gz"

    if command -v wget &> /dev/null; then
        wget -q --show-progress -O "apache-arrow-${ARROW_VERSION}.tar.gz" "${ARROW_URL}"
    else
        curl -L -# -o "apache-arrow-${ARROW_VERSION}.tar.gz" "${ARROW_URL}"
    fi

    echo -e "  Extracting..."
    tar -xzf "apache-arrow-${ARROW_VERSION}.tar.gz"
    mv "arrow-apache-arrow-${ARROW_VERSION}" "${ARROW_SRC}"
    rm "apache-arrow-${ARROW_VERSION}.tar.gz"

    # Check for compression libraries
    SNAPPY_OPT="OFF"
    ZLIB_OPT="OFF"
    if pkg-config --exists snappy 2>/dev/null; then
        SNAPPY_OPT="ON"
    fi
    if pkg-config --exists zlib 2>/dev/null || [[ -f /usr/include/zlib.h ]]; then
        ZLIB_OPT="ON"
    fi

    echo -e "  Configuring..."
    mkdir -p "${ARROW_BUILD}"
    cd "${ARROW_BUILD}"

    cmake "${ARROW_SRC}/cpp" \
        -DCMAKE_INSTALL_PREFIX="${ARROW_INSTALL}" \
        -DCMAKE_BUILD_TYPE=Release \
        -DARROW_BUILD_STATIC=OFF \
        -DARROW_BUILD_SHARED=ON \
        -DARROW_PARQUET=ON \
        -DARROW_WITH_SNAPPY=${SNAPPY_OPT} \
        -DARROW_WITH_ZLIB=${ZLIB_OPT} \
        -DARROW_WITH_ZSTD=OFF \
        -DARROW_WITH_LZ4=OFF \
        -DARROW_WITH_BROTLI=OFF \
        -DARROW_COMPUTE=OFF \
        -DARROW_DATASET=OFF \
        -DARROW_FILESYSTEM=ON \
        -DARROW_JSON=OFF \
        -DARROW_CSV=OFF \
        -DARROW_FLIGHT=OFF \
        -DARROW_GANDIVA=OFF \
        -DARROW_ORC=OFF \
        -DARROW_PLASMA=OFF \
        -DARROW_BUILD_UTILITIES=OFF \
        -DARROW_BUILD_TESTS=OFF \
        -DARROW_BUILD_EXAMPLES=OFF \
        -DARROW_BUILD_BENCHMARKS=OFF \
        -DARROW_BUILD_INTEGRATION=OFF \
        -DPARQUET_BUILD_EXECUTABLES=OFF \
        -DPARQUET_BUILD_EXAMPLES=OFF \
        -DPARQUET_REQUIRE_ENCRYPTION=OFF \
        > cmake_output.log 2>&1

    if [[ $? -ne 0 ]]; then
        echo -e "${RED}Error: CMake configuration failed${NC}"
        echo "Check ${ARROW_BUILD}/cmake_output.log for details"
        exit 1
    fi

    echo -e "  Building (using ${PARALLEL} parallel jobs)..."
    make -j${PARALLEL} 2>&1 | tee build_output.log | grep -E "^\[|Built target"

    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
        echo -e "${RED}Error: Build failed${NC}"
        echo "Check ${ARROW_BUILD}/build_output.log for details"
        exit 1
    fi

    echo -e "  Installing..."
    make install > install_output.log 2>&1

    cd "${EXTERNAL_DIR}"
    echo -e "  ${GREEN}✓${NC} Arrow installed to ${ARROW_INSTALL}"

    # Cleanup source and build (keep for debugging if needed)
    # rm -rf "${ARROW_SRC}" "${ARROW_BUILD}"
fi
echo ""

# ============================================================================
# Summary
# ============================================================================
echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║${NC}${BOLD}                    Dependency Setup Complete!                             ${NC}${CYAN}║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${BOLD}Installed to external/:${NC}"
echo -e "  ${GREEN}✓${NC} json/           - nlohmann_json headers"
echo -e "  ${GREEN}✓${NC} spdlog-install/ - spdlog library"
if [[ "$SKIP_ARROW" != "true" ]]; then
    echo -e "  ${GREEN}✓${NC} arrow-install/  - Arrow/Parquet libraries"
fi
echo -e "  ${GREEN}✓${NC} mcpl/           - MCPL (already bundled)"
echo -e "  ${GREEN}✓${NC} parquet-writer/ - Parquet writer (already bundled)"
echo ""

echo -e "${BOLD}Library locations for manual linking (if needed):${NC}"
echo -e "  Arrow:    ${CYAN}${ARROW_INSTALL}/lib/cmake/Arrow${NC}"
echo -e "  Parquet:  ${CYAN}${ARROW_INSTALL}/lib/cmake/Parquet${NC}"
echo -e "  spdlog:   ${CYAN}${SPDLOG_INSTALL}/lib/cmake/spdlog${NC}"
echo -e "  json:     ${CYAN}${JSON_DIR}/include${NC}"
echo ""

echo -e "${BOLD}Next steps:${NC}"
echo -e "  1. Configure: ${CYAN}mkdir -p build && cd build && cmake ..${NC}"
echo -e "  2. Build:     ${CYAN}make -j\$(nproc)${NC}"
echo -e "  3. Run:       ${CYAN}./nnbar-calo-sim${NC}"
echo ""

# Check if everything is ready
ALL_READY=true
[[ ! -f "${JSON_DIR}/include/nlohmann/json.hpp" ]] && ALL_READY=false
[[ ! -f "${SPDLOG_INSTALL}/lib/cmake/spdlog/spdlogConfig.cmake" ]] && ALL_READY=false
if [[ "$SKIP_ARROW" != "true" ]]; then
    [[ ! -f "${ARROW_INSTALL}/lib/cmake/Arrow/ArrowConfig.cmake" ]] && ALL_READY=false
fi

if [[ "$ALL_READY" == "true" ]]; then
    echo -e "${GREEN}All dependencies are ready! You can now build the simulation.${NC}"
else
    echo -e "${YELLOW}Some dependencies may not have installed correctly. Check the logs above.${NC}"
    exit 1
fi
