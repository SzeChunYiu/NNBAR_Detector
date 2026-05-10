#!/bin/bash
# ============================================================================
# build_arrow.sh - Download and build Apache Arrow with Parquet support
# ============================================================================
# This script downloads and builds Apache Arrow as a bundled dependency
# for the NNBAR detector simulation.
#
# Usage: ./scripts/build_arrow.sh [options]
# Options:
#   --clean       Remove existing Arrow build and start fresh
#   --parallel N  Use N parallel jobs (default: nproc)
# ============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# Configuration
ARROW_VERSION="17.0.0"
ARROW_URL="https://github.com/apache/arrow/archive/refs/tags/apache-arrow-${ARROW_VERSION}.tar.gz"

# Directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
EXTERNAL_DIR="${PROJECT_DIR}/external"
ARROW_SOURCE_DIR="${EXTERNAL_DIR}/arrow-src"
ARROW_BUILD_DIR="${EXTERNAL_DIR}/arrow-build"
ARROW_INSTALL_DIR="${EXTERNAL_DIR}/arrow-install"

# Parse arguments
CLEAN=false
PARALLEL=$(nproc 2>/dev/null || echo 4)

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
        -h|--help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --clean       Remove existing Arrow build and start fresh"
            echo "  --parallel N  Use N parallel jobs (default: $(nproc 2>/dev/null || echo 4))"
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
echo -e "${CYAN}║${NC}${BOLD}           Apache Arrow/Parquet - Dependency Builder                       ${NC}${CYAN}║${NC}"
echo -e "${CYAN}║${NC}           Version: ${GREEN}${ARROW_VERSION}${NC}                                               ${CYAN}║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Check if already installed
if [[ -f "${ARROW_INSTALL_DIR}/lib/cmake/Arrow/ArrowConfig.cmake" ]] && [[ "$CLEAN" == "false" ]]; then
    echo -e "${GREEN}✓${NC} Arrow is already installed at: ${ARROW_INSTALL_DIR}"
    echo -e "  Use ${CYAN}--clean${NC} to rebuild from scratch."
    echo ""
    echo -e "${GREEN}Arrow installation complete!${NC}"
    exit 0
fi

# Clean if requested
if [[ "$CLEAN" == "true" ]]; then
    echo -e "${YELLOW}[1/5]${NC} Cleaning previous build..."
    rm -rf "${ARROW_SOURCE_DIR}" "${ARROW_BUILD_DIR}" "${ARROW_INSTALL_DIR}"
    echo -e "  ${GREEN}✓${NC} Cleaned"
fi

# Check dependencies
echo -e "${BLUE}[1/5]${NC} Checking build dependencies..."

check_command() {
    if command -v "$1" &> /dev/null; then
        echo -e "  ${GREEN}✓${NC} $1"
        return 0
    else
        echo -e "  ${RED}✗${NC} $1 - ${RED}not found${NC}"
        return 1
    fi
}

DEPS_OK=true
check_command cmake || DEPS_OK=false
check_command g++ || check_command clang++ || DEPS_OK=false
check_command make || DEPS_OK=false
check_command wget || check_command curl || DEPS_OK=false

if [[ "$DEPS_OK" == "false" ]]; then
    echo ""
    echo -e "${RED}Error: Missing build dependencies${NC}"
    echo "Please install: cmake, g++ or clang++, make, wget or curl"
    exit 1
fi

# Check for optional compression libraries
echo ""
echo -e "  ${CYAN}Optional compression libraries:${NC}"
if pkg-config --exists snappy 2>/dev/null; then
    echo -e "  ${GREEN}✓${NC} snappy (compression)"
    SNAPPY_OPT="ON"
else
    echo -e "  ${YELLOW}○${NC} snappy - not found (will build without)"
    SNAPPY_OPT="OFF"
fi

if pkg-config --exists zlib 2>/dev/null || [[ -f /usr/include/zlib.h ]]; then
    echo -e "  ${GREEN}✓${NC} zlib (compression)"
    ZLIB_OPT="ON"
else
    echo -e "  ${YELLOW}○${NC} zlib - not found (will build without)"
    ZLIB_OPT="OFF"
fi

# Download Arrow source
echo ""
echo -e "${BLUE}[2/5]${NC} Downloading Apache Arrow ${GREEN}v${ARROW_VERSION}${NC}..."

mkdir -p "${EXTERNAL_DIR}"
cd "${EXTERNAL_DIR}"

if [[ ! -d "${ARROW_SOURCE_DIR}" ]]; then
    TARBALL="apache-arrow-${ARROW_VERSION}.tar.gz"

    if command -v wget &> /dev/null; then
        wget -q --show-progress -O "${TARBALL}" "${ARROW_URL}"
    else
        curl -L -# -o "${TARBALL}" "${ARROW_URL}"
    fi

    echo -e "  Extracting..."
    tar -xzf "${TARBALL}"
    mv "arrow-apache-arrow-${ARROW_VERSION}" "${ARROW_SOURCE_DIR}"
    rm "${TARBALL}"
    echo -e "  ${GREEN}✓${NC} Download complete"
else
    echo -e "  ${GREEN}✓${NC} Source already downloaded"
fi

# Configure Arrow
echo ""
echo -e "${BLUE}[3/5]${NC} Configuring Arrow build..."

mkdir -p "${ARROW_BUILD_DIR}"
cd "${ARROW_BUILD_DIR}"

cmake "${ARROW_SOURCE_DIR}/cpp" \
    -DCMAKE_INSTALL_PREFIX="${ARROW_INSTALL_DIR}" \
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
    2>&1 | tee cmake_output.log

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo -e "${RED}Error: CMake configuration failed${NC}"
    echo "Check ${ARROW_BUILD_DIR}/cmake_output.log for details"
    exit 1
fi

echo -e "  ${GREEN}✓${NC} Configuration complete"

# Build Arrow
echo ""
echo -e "${BLUE}[4/5]${NC} Building Arrow (using ${PARALLEL} parallel jobs)..."
echo -e "  ${YELLOW}This may take 5-15 minutes...${NC}"

make -j${PARALLEL} 2>&1 | tee build_output.log

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo -e "${RED}Error: Build failed${NC}"
    echo "Check ${ARROW_BUILD_DIR}/build_output.log for details"
    exit 1
fi

echo -e "  ${GREEN}✓${NC} Build complete"

# Install Arrow
echo ""
echo -e "${BLUE}[5/5]${NC} Installing Arrow to external/arrow-install/..."

make install 2>&1 | tee install_output.log

if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo -e "${RED}Error: Installation failed${NC}"
    exit 1
fi

echo -e "  ${GREEN}✓${NC} Installation complete"

# Verify installation
echo ""
echo -e "${CYAN}═══════════════════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}${BOLD}Apache Arrow/Parquet successfully installed!${NC}"
echo -e "${CYAN}═══════════════════════════════════════════════════════════════════════════${NC}"
echo ""
echo -e "  ${BOLD}Installation directory:${NC}"
echo -e "    ${ARROW_INSTALL_DIR}"
echo ""
echo -e "  ${BOLD}CMake will automatically detect Arrow at:${NC}"
echo -e "    ${ARROW_INSTALL_DIR}/lib/cmake/Arrow/ArrowConfig.cmake"
echo -e "    ${ARROW_INSTALL_DIR}/lib/cmake/Parquet/ParquetConfig.cmake"
echo ""
echo -e "  ${BOLD}Library files:${NC}"
ls -la "${ARROW_INSTALL_DIR}/lib/"*.so* 2>/dev/null | head -6 || true
echo ""
echo -e "  ${BOLD}Next step:${NC} Run CMake to build the simulation"
echo -e "    ${CYAN}cd build && cmake .. && make -j\$(nproc)${NC}"
echo ""

# Clean up source and build directories to save space (optional)
echo -e "  ${YELLOW}Note:${NC} Source and build dirs kept for debugging."
echo -e "        Delete them with: ${CYAN}rm -rf external/arrow-{src,build}${NC}"
echo ""
