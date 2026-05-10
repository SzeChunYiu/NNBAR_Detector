#!/bin/bash
# ============================================================================
# Opticks Rebuild Guide for Geant4 11.2.2
# ============================================================================
#
# IMPORTANT: Rebuilding Opticks for a new Geant4 version is complex!
#
# The current Opticks installation was built against Geant4 10.7.4.
# To use Opticks with Geant4 11.2.2, the ENTIRE Opticks stack needs
# to be rebuilt from scratch, not just U4 and G4CX.
#
# Required components in order:
#   1. CLHEP (standalone, not Geant4's internal CLHEP)
#   2. Geant4 11.2.2 (linked against standalone CLHEP)
#   3. Opticks core: OKConf, SysRap, NPY, AnalyticMesh
#   4. Opticks geometry: CSG, CSGOptiX, GDXML
#   5. Opticks Geant4 integration: U4, G4CX
#
# The Opticks build system (opticks-full) handles this complexity.
#
# Quick Rebuild Option:
#   If you have the opticks build system set up, run:
#   $ source /path/to/opticks/bin/opticks-env.sh
#   $ opticks-full
#
# Alternative: Download pre-built Opticks for Geant4 11.2.x from:
#   https://github.com/simoncblyth/opticks/releases
#
# For now, Celeritas GPU acceleration is available and working!
# Opticks can be enabled later once rebuilt.
#
# ============================================================================

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  Opticks Status Check for NNBAR${NC}"
echo -e "${CYAN}========================================${NC}"

# Environment setup
export PATH=/usr/local/bin:/usr/local/cuda-12.4/bin:/usr/bin:/bin:$PATH
export CUDA_HOME=/usr/local/cuda-12.4
unset LD_LIBRARY_PATH

# Directories
GEANT4_INSTALL=/home/billy/nnbar/simulation/GEANT4_Packages/install/geant4-11.2.2
OPTICKS_INSTALL=/home/billy/nnbar/simulation/GEANT4_Packages/install/opticks
OPTIX_INSTALL=/home/billy/nnbar/simulation/GEANT4_Packages/install/optix

echo ""
echo -e "${GREEN}Checking environment...${NC}"
echo ""

# Check Geant4
if [ -f "${GEANT4_INSTALL}/bin/geant4.sh" ]; then
    source ${GEANT4_INSTALL}/bin/geant4.sh
    echo -e "  Geant4:   ${GREEN}$(geant4-config --version)${NC} at ${GEANT4_INSTALL}"
else
    echo -e "  Geant4:   ${RED}NOT FOUND${NC}"
fi

# Check CUDA
if command -v nvcc &> /dev/null; then
    echo -e "  CUDA:     ${GREEN}$(nvcc --version 2>&1 | grep release | sed 's/.*release //' | sed 's/,.*//')${NC}"
else
    echo -e "  CUDA:     ${RED}NOT FOUND${NC}"
fi

# Check OptiX
if [ -d "${OPTIX_INSTALL}/include" ]; then
    OPTIX_VER=$(grep "OPTIX_VERSION" ${OPTIX_INSTALL}/include/optix.h 2>/dev/null | head -1 | grep -oP '\d+' || echo "unknown")
    echo -e "  OptiX:    ${GREEN}${OPTIX_VER}${NC} at ${OPTIX_INSTALL}"
else
    echo -e "  OptiX:    ${RED}NOT FOUND${NC}"
fi

# Check Opticks
if [ -f "${OPTICKS_INSTALL}/lib/libG4CX.so" ]; then
    echo -e "  Opticks:  ${GREEN}INSTALLED${NC} at ${OPTICKS_INSTALL}"

    # Check what Geant4 version Opticks was built against
    G4_LINK=$(ldd "${OPTICKS_INSTALL}/lib/libG4CX.so" 2>/dev/null | grep libG4global | head -1)
    if echo "$G4_LINK" | grep -q "geant4-11.2.2"; then
        echo -e "            ${GREEN}Linked to Geant4 11.2.2${NC}"
    elif echo "$G4_LINK" | grep -q "geant4-10"; then
        G4_OLD=$(echo "$G4_LINK" | grep -oP 'geant4-\d+\.\d+\.\d+')
        echo -e "            ${YELLOW}WARNING: Linked to ${G4_OLD}${NC}"
        echo -e "            ${YELLOW}Needs rebuild for Geant4 11.2.2${NC}"
    else
        echo -e "            ${YELLOW}Cannot determine linked Geant4 version${NC}"
    fi
else
    echo -e "  Opticks:  ${RED}NOT FOUND${NC}"
fi

echo ""
echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  Current GPU Acceleration Status${NC}"
echo -e "${CYAN}========================================${NC}"
echo ""
echo -e "  ${GREEN}[WORKING]${NC}  Celeritas - GPU EM physics offloading"
echo -e "  ${YELLOW}[DISABLED]${NC} Opticks   - GPU optical photon propagation"
echo ""
echo -e "${YELLOW}To enable Opticks:${NC}"
echo "  1. Rebuild entire Opticks stack against Geant4 11.2.2"
echo "  2. See https://bitbucket.org/simoncblyth/opticks for build instructions"
echo "  3. After rebuild, rerun: cmake .. -DWITH_OPTICKS=ON"
echo ""

# ============================================================================
# Partial rebuild attempt (for reference - usually fails due to CLHEP issue)
# ============================================================================

if [ "$1" == "--attempt-rebuild" ]; then
    echo -e "${YELLOW}Attempting partial rebuild (experimental)...${NC}"
    echo ""

    OPTICKS_SRC=/home/billy/nnbar/simulation/GEANT4_Packages/src/opticks
    BUILD_DIR=/home/billy/nnbar/simulation/GEANT4_Packages/build/opticks_rebuild

    export OPTICKS_OPTIX_PREFIX=${OPTIX_INSTALL}
    export CMAKE_PREFIX_PATH=${OPTICKS_INSTALL}:${OPTICKS_INSTALL}/lib/cmake:${OPTICKS_INSTALL}/externals:${GEANT4_INSTALL}
    export CMAKE_MODULE_PATH=${OPTICKS_SRC}/cmake/Modules
    export OPTICKS_PREFIX=${OPTICKS_INSTALL}

    mkdir -p ${BUILD_DIR}/sysrap
    cd ${BUILD_DIR}/sysrap

    echo "Building SysRap..."
    cmake ${OPTICKS_SRC}/sysrap \
        -DCMAKE_INSTALL_PREFIX=${OPTICKS_INSTALL} \
        -DCMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH}" \
        -DCMAKE_MODULE_PATH="${CMAKE_MODULE_PATH}" \
        -DOPTICKS_PREFIX="${OPTICKS_PREFIX}" \
        -DCMAKE_CXX_FLAGS="-Wno-error" \
        -DCMAKE_BUILD_TYPE=Release

    make -j$(nproc) SysRap
    make install

    echo ""
    echo -e "${GREEN}SysRap rebuilt.${NC}"
    echo -e "${YELLOW}U4 and G4CX require full Opticks rebuild due to CLHEP dependency.${NC}"
fi
