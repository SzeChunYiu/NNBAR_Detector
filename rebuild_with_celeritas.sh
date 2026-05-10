#!/bin/bash
# ============================================================================
# NNBAR Detector Simulation - Full GPU Rebuild Script
# This rebuilds with Geant4 11.2.2, Celeritas, and GarfieldGPU CUDA
# ============================================================================

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  NNBAR Full GPU Rebuild${NC}"
echo -e "${CYAN}========================================${NC}"

# ============================================================================
# Setup paths - CUDA must be first!
# ============================================================================
export PATH=/usr/local/cuda-12.4/bin:/usr/local/cuda/bin:/usr/local/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda-12.4/lib64:/usr/local/cuda/lib64:$LD_LIBRARY_PATH

# Verify CUDA
if ! command -v nvcc &> /dev/null; then
    echo -e "${RED}ERROR: nvcc not found. Please install CUDA.${NC}"
    exit 1
fi
echo -e "${GREEN}CUDA:${NC} $(nvcc --version | grep release)"

# ============================================================================
# Geant4 11.2.2
# ============================================================================
GEANT4_DIR="/home/billy/nnbar/simulation/GEANT4_Packages/install/geant4-11.2.2"
if [ -f "${GEANT4_DIR}/bin/geant4.sh" ]; then
    source ${GEANT4_DIR}/bin/geant4.sh
    export Geant4_DIR=${GEANT4_DIR}/lib/cmake/Geant4
    export LD_LIBRARY_PATH=${GEANT4_DIR}/lib:$LD_LIBRARY_PATH
    echo -e "${GREEN}Geant4:${NC} $(geant4-config --version)"
else
    echo -e "${RED}ERROR: Geant4 11.2.2 not found at ${GEANT4_DIR}${NC}"
    exit 1
fi

# ============================================================================
# Celeritas
# ============================================================================
CELERITAS_DIR="/home/billy/nnbar/simulation/GEANT4_Packages/install/celeritas"
if [ -d "${CELERITAS_DIR}/lib/cmake/Celeritas" ]; then
    export Celeritas_DIR=${CELERITAS_DIR}/lib/cmake/Celeritas
    export LD_LIBRARY_PATH=${CELERITAS_DIR}/lib:$LD_LIBRARY_PATH
    CELERITAS_FLAG="-DCeleritas_DIR=${CELERITAS_DIR}/lib/cmake/Celeritas"
    echo -e "${GREEN}Celeritas:${NC} Available"
else
    echo -e "${YELLOW}WARNING: Celeritas not found, building without it${NC}"
    CELERITAS_FLAG=""
fi

# ============================================================================
# Opticks (GPU Optical Photons)
# ============================================================================
OPTICKS_DIR="/home/billy/nnbar/simulation/GEANT4_Packages/install/opticks"
OPTIX_DIR="/home/billy/nnbar/simulation/GEANT4_Packages/install/optix"
if [ -f "${OPTICKS_DIR}/lib/libG4CX.so" ] && [ -f "${OPTICKS_DIR}/lib/libU4.so" ]; then
    export LD_LIBRARY_PATH=${OPTICKS_DIR}/lib:$LD_LIBRARY_PATH
    export OPTICKS_HOME=${OPTICKS_DIR}
    echo -e "${GREEN}Opticks:${NC} G4CX/U4 Available"
else
    echo -e "${YELLOW}WARNING: Opticks G4CX/U4 not found${NC}"
fi

# ============================================================================
# Clean and enter build directory
# ============================================================================
BUILD_DIR="/home/billy/nnbar/simulation/NNBAR_Detector/build"
echo ""
echo -e "${YELLOW}Cleaning build directory...${NC}"
rm -rf ${BUILD_DIR}/CMakeCache.txt ${BUILD_DIR}/CMakeFiles
cd ${BUILD_DIR}

# ============================================================================
# Configure CMake
# ============================================================================
echo ""
echo -e "${CYAN}Configuring CMake with GPU acceleration...${NC}"
echo ""

cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DWITH_DASHBOARD=ON \
    -DWITH_GARFIELD_GPU=ON \
    -DWITH_CELERITAS=ON \
    -DWITH_OPTICKS=OFF \
    -DCMAKE_CUDA_COMPILER=/usr/local/cuda-12.4/bin/nvcc \
    -DCMAKE_CUDA_ARCHITECTURES=86 \
    -DGeant4_DIR=${Geant4_DIR} \
    -DThrust_DIR=/usr/local/cuda-12.4/targets/x86_64-linux/lib/cmake/thrust \
    ${CELERITAS_FLAG} \
    -DMCPL_BUILD=ON

# ============================================================================
# Build
# ============================================================================
echo ""
echo -e "${CYAN}Building (this may take a few minutes)...${NC}"
make -j$(nproc)

# ============================================================================
# Summary
# ============================================================================
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Build Complete!${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "To run the simulation, first source the environment:"
echo ""
echo "  export PATH=/usr/local/cuda-12.4/bin:\$PATH"
echo "  export LD_LIBRARY_PATH=/usr/local/cuda-12.4/lib64:\$LD_LIBRARY_PATH"
echo "  source /home/billy/nnbar/simulation/GEANT4_Packages/install/geant4-11.2.2/bin/geant4.sh"
echo "  export LD_LIBRARY_PATH=/home/billy/nnbar/simulation/GEANT4_Packages/install/celeritas/lib:\$LD_LIBRARY_PATH"
echo ""
echo "Then run:"
echo "  ./nnbar-detector-simulation"
echo ""
