#!/bin/bash
# ============================================================================
# NNBAR Detector Simulation - Interactive Configuration Script
# ============================================================================
# This script helps configure the build with user-friendly prompts
# Run this BEFORE cmake to set up your build options
# ============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
EXTERNAL_DIR="$PROJECT_DIR/external"
BUILD_DIR="$PROJECT_DIR/build"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color
BOLD='\033[1m'

# Default options
OPT_SCINTILLATION="OFF"
OPT_GARFIELD="OFF"
OPT_OPTICKS="OFF"
OPT_DASHBOARD="OFF"
OPT_MCPL="OFF"
OPT_DEBUG="OFF"

echo ""
echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║${NC}${BOLD}        NNBAR/HIBEAM Detector Simulation Configuration            ${NC}${CYAN}║${NC}"
echo -e "${CYAN}║${NC}        Geant4-based Monte Carlo for n-nbar oscillation           ${CYAN}║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# ============================================================================
# Check Dependencies
# ============================================================================
echo -e "${BOLD}Checking dependencies...${NC}"
echo ""

# Function to check if a dependency exists
check_dep() {
    local name=$1
    local bundled_path=$2
    local system_check=$3

    if [ -d "$bundled_path" ]; then
        echo -e "  ${GREEN}✓${NC} $name: ${GREEN}Found (bundled)${NC}"
        return 0
    elif eval "$system_check" 2>/dev/null; then
        echo -e "  ${GREEN}✓${NC} $name: ${GREEN}Found (system)${NC}"
        return 0
    else
        echo -e "  ${YELLOW}○${NC} $name: ${YELLOW}Not found (will be downloaded)${NC}"
        return 1
    fi
}

echo -e "${BLUE}Required Dependencies:${NC}"

# Check Geant4
if command -v geant4-config &> /dev/null; then
    G4_VERSION=$(geant4-config --version)
    echo -e "  ${GREEN}✓${NC} Geant4: ${GREEN}Found v${G4_VERSION}${NC}"
else
    echo -e "  ${RED}✗${NC} Geant4: ${RED}NOT FOUND - Required!${NC}"
    echo -e "    ${YELLOW}Please install Geant4 first: https://geant4.web.cern.ch/${NC}"
    exit 1
fi

# Check Arrow/Parquet
check_dep "Arrow/Parquet" "$EXTERNAL_DIR/arrow-install" "pkg-config --exists arrow-glib"
ARROW_STATUS=$?

# Check nlohmann_json
check_dep "nlohmann_json" "$EXTERNAL_DIR/json" "pkg-config --exists nlohmann_json"
JSON_STATUS=$?

# Check spdlog
check_dep "spdlog" "$EXTERNAL_DIR/spdlog-install" "pkg-config --exists spdlog"
SPDLOG_STATUS=$?

# Check MCPL (always bundled)
if [ -f "$EXTERNAL_DIR/mcpl/mcpl.c" ]; then
    echo -e "  ${GREEN}✓${NC} MCPL: ${GREEN}Found (bundled)${NC}"
else
    echo -e "  ${RED}✗${NC} MCPL: ${RED}Missing from external/mcpl${NC}"
fi

echo ""
echo -e "${BLUE}Optional Dependencies:${NC}"

# Check Garfield++
if [ -n "$GARFIELD_HOME" ] && [ -d "$GARFIELD_HOME" ]; then
    echo -e "  ${GREEN}✓${NC} Garfield++: ${GREEN}Found at $GARFIELD_HOME${NC}"
    GARFIELD_AVAILABLE=1
else
    echo -e "  ${YELLOW}○${NC} Garfield++: ${YELLOW}Not found (TPC drift simulation disabled)${NC}"
    GARFIELD_AVAILABLE=0
fi

# Check Opticks/CUDA
if command -v nvcc &> /dev/null; then
    CUDA_VERSION=$(nvcc --version | grep release | awk '{print $5}' | cut -d',' -f1)
    echo -e "  ${GREEN}✓${NC} CUDA: ${GREEN}Found v${CUDA_VERSION}${NC}"
    if [ -n "$OPTICKS_HOME" ]; then
        echo -e "  ${GREEN}✓${NC} Opticks: ${GREEN}Found at $OPTICKS_HOME${NC}"
        OPTICKS_AVAILABLE=1
    else
        echo -e "  ${YELLOW}○${NC} Opticks: ${YELLOW}Not found (GPU optical disabled)${NC}"
        OPTICKS_AVAILABLE=0
    fi
else
    echo -e "  ${YELLOW}○${NC} CUDA/Opticks: ${YELLOW}Not found (GPU acceleration disabled)${NC}"
    OPTICKS_AVAILABLE=0
fi

# Check Qt
if command -v qmake &> /dev/null || command -v qmake6 &> /dev/null; then
    echo -e "  ${GREEN}✓${NC} Qt: ${GREEN}Found${NC}"
    QT_AVAILABLE=1
else
    echo -e "  ${YELLOW}○${NC} Qt: ${YELLOW}Not found (Dashboard disabled)${NC}"
    QT_AVAILABLE=0
fi

# ============================================================================
# Install Missing Dependencies
# ============================================================================
if [ $ARROW_STATUS -ne 0 ] || [ $JSON_STATUS -ne 0 ] || [ $SPDLOG_STATUS -ne 0 ]; then
    echo ""
    echo -e "${YELLOW}Some dependencies are missing.${NC}"
    read -p "Would you like to download and install them now? [Y/n] " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]] || [[ -z $REPLY ]]; then
        echo -e "${BLUE}Running dependency setup script...${NC}"
        bash "$SCRIPT_DIR/setup_dependencies.sh"
    fi
fi

# ============================================================================
# Interactive Configuration
# ============================================================================
echo ""
echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║${NC}${BOLD}                    Build Configuration                            ${NC}${CYAN}║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Scintillation
echo -e "${BOLD}1. Optical Photon Simulation (Scintillation/Cerenkov)${NC}"
echo -e "   Generates optical photons in scintillators and lead glass."
echo -e "   ${YELLOW}Warning: Significantly slower simulation (~10-100x)${NC}"
read -p "   Enable optical photons? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    OPT_SCINTILLATION="ON"
    echo -e "   ${GREEN}→ Optical photons: ENABLED${NC}"
else
    echo -e "   ${BLUE}→ Optical photons: DISABLED (fast mode)${NC}"
fi
echo ""

# Garfield++
if [ $GARFIELD_AVAILABLE -eq 1 ]; then
    echo -e "${BOLD}2. Garfield++ TPC Simulation${NC}"
    echo -e "   Uses Garfield++ for realistic electron drift in TPC gas."
    echo -e "   Provides accurate ionization and diffusion modeling."
    read -p "   Enable Garfield++ TPC? [y/N] " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        OPT_GARFIELD="ON"
        echo -e "   ${GREEN}→ Garfield++ TPC: ENABLED${NC}"
    else
        echo -e "   ${BLUE}→ Garfield++ TPC: DISABLED (using simple TPC model)${NC}"
    fi
else
    echo -e "${BOLD}2. Garfield++ TPC Simulation${NC}"
    echo -e "   ${YELLOW}→ Not available (set GARFIELD_HOME to enable)${NC}"
fi
echo ""

# Opticks GPU
if [ $OPTICKS_AVAILABLE -eq 1 ]; then
    echo -e "${BOLD}3. Opticks GPU Optical Photon Acceleration${NC}"
    echo -e "   Uses NVIDIA GPU for optical photon propagation."
    echo -e "   ${GREEN}Provides 50-200x speedup for optical simulations.${NC}"
    read -p "   Enable Opticks GPU? [y/N] " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        OPT_OPTICKS="ON"
        echo -e "   ${GREEN}→ Opticks GPU: ENABLED${NC}"
    else
        echo -e "   ${BLUE}→ Opticks GPU: DISABLED (using CPU optical physics)${NC}"
    fi
else
    echo -e "${BOLD}3. GPU Optical Photon Acceleration${NC}"
    echo -e "   ${YELLOW}→ Not available (requires CUDA + Opticks)${NC}"
fi
echo ""

# Qt Dashboard
if [ $QT_AVAILABLE -eq 1 ]; then
    echo -e "${BOLD}4. Qt Monitoring Dashboard${NC}"
    echo -e "   Real-time visualization of particles, energy deposits, and statistics."
    read -p "   Enable Qt Dashboard? [y/N] " -n 1 -r
    echo ""
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        OPT_DASHBOARD="ON"
        echo -e "   ${GREEN}→ Qt Dashboard: ENABLED${NC}"
    else
        echo -e "   ${BLUE}→ Qt Dashboard: DISABLED${NC}"
    fi
else
    echo -e "${BOLD}4. Qt Monitoring Dashboard${NC}"
    echo -e "   ${YELLOW}→ Not available (install Qt5 or Qt6 to enable)${NC}"
fi
echo ""

# MCPL Input
echo -e "${BOLD}5. Particle Source${NC}"
echo -e "   Choose between built-in particle gun or MCPL file input."
echo -e "   MCPL files contain pre-generated particle lists from external sources."
read -p "   Use MCPL file input? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    OPT_MCPL="ON"
    echo -e "   ${GREEN}→ Particle source: MCPL file${NC}"
else
    echo -e "   ${BLUE}→ Particle source: Built-in particle gun${NC}"
fi
echo ""

# Debug mode
echo -e "${BOLD}6. Debug Output${NC}"
read -p "   Enable verbose debug output? [y/N] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    OPT_DEBUG="ON"
    echo -e "   ${GREEN}→ Debug output: ENABLED${NC}"
else
    echo -e "   ${BLUE}→ Debug output: DISABLED${NC}"
fi

# ============================================================================
# Configuration Summary
# ============================================================================
echo ""
echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║${NC}${BOLD}                    Configuration Summary                          ${NC}${CYAN}║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "  ${BOLD}Option                    Value${NC}"
echo -e "  ─────────────────────────────────────"
echo -e "  WITH_SCINTILLATION      ${OPT_SCINTILLATION}"
echo -e "  WITH_GARFIELD           ${OPT_GARFIELD}"
echo -e "  WITH_OPTICKS            ${OPT_OPTICKS}"
echo -e "  WITH_DASHBOARD          ${OPT_DASHBOARD}"
echo -e "  MCPL_BUILD              ${OPT_MCPL}"
echo -e "  DEBUG_VERBOSE           ${OPT_DEBUG}"
echo ""

# ============================================================================
# Generate CMake Command
# ============================================================================
CMAKE_CMD="cmake"
CMAKE_CMD+=" -DCMAKE_BUILD_TYPE=Release"
CMAKE_CMD+=" -DWITH_SCINTILLATION=${OPT_SCINTILLATION}"
CMAKE_CMD+=" -DWITH_GARFIELD=${OPT_GARFIELD}"
CMAKE_CMD+=" -DWITH_OPTICKS=${OPT_OPTICKS}"
CMAKE_CMD+=" -DWITH_DASHBOARD=${OPT_DASHBOARD}"
CMAKE_CMD+=" -DMCPL_BUILD=${OPT_MCPL}"
CMAKE_CMD+=" -DDEBUG_VERBOSE=${OPT_DEBUG}"
CMAKE_CMD+=" -DTARGET_BUILD=ON"

# Add bundled dependency paths if they exist
if [ -d "$EXTERNAL_DIR/arrow-install" ]; then
    CMAKE_CMD+=" -DArrow_DIR=$EXTERNAL_DIR/arrow-install/lib/cmake/Arrow"
    CMAKE_CMD+=" -DParquet_DIR=$EXTERNAL_DIR/arrow-install/lib/cmake/Parquet"
fi
if [ -d "$EXTERNAL_DIR/spdlog-install" ]; then
    CMAKE_CMD+=" -Dspdlog_DIR=$EXTERNAL_DIR/spdlog-install/lib/cmake/spdlog"
fi

CMAKE_CMD+=" .."

# ============================================================================
# Run CMake
# ============================================================================
read -p "Proceed with configuration? [Y/n] " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Nn]$ ]]; then
    echo -e "${YELLOW}Configuration cancelled.${NC}"
    echo ""
    echo "To configure manually, run:"
    echo "  mkdir -p build && cd build"
    echo "  $CMAKE_CMD"
    exit 0
fi

# Create and enter build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

echo ""
echo -e "${BLUE}Running CMake...${NC}"
echo ""
eval $CMAKE_CMD

echo ""
echo -e "${GREEN}╔══════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║${NC}${BOLD}                Configuration Complete!                            ${NC}${GREEN}║${NC}"
echo -e "${GREEN}╚══════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "To build the simulation:"
echo -e "  ${CYAN}cd build && make -j\$(nproc)${NC}"
echo ""
echo -e "To run:"
echo -e "  ${CYAN}./nnbar-calo-sim macro/signal/run_signal.mac${NC}"
echo ""
