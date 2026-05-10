#!/bin/bash
# Step 5: Install Celeritas (GPU EM Physics)
set -e

INSTALL_DIR="/home/billy/nnbar/simulation/GEANT4_Packages"
NPROC=$(nproc)

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║          Installing Celeritas (GPU EM Physics)               ║"
echo "╚══════════════════════════════════════════════════════════════╝"

# Check CUDA
if ! command -v nvcc &> /dev/null; then
    echo "ERROR: CUDA not found. Please run install_step1_cuda.sh first."
    exit 1
fi

# Source Geant4
if [ -f "$INSTALL_DIR/install/geant4-11.2.2/bin/geant4.sh" ]; then
    source "$INSTALL_DIR/install/geant4-11.2.2/bin/geant4.sh"
    export Geant4_DIR="$INSTALL_DIR/install/geant4-11.2.2/lib/cmake/Geant4"
fi

cd "$INSTALL_DIR/src"

# Clone Celeritas
if [ ! -d "celeritas" ]; then
    echo "Cloning Celeritas..."
    git clone --depth 1 https://github.com/celeritas-project/celeritas.git
fi

# Build
mkdir -p "$INSTALL_DIR/build/celeritas"
cd "$INSTALL_DIR/build/celeritas"

echo "Configuring Celeritas..."
cmake "$INSTALL_DIR/src/celeritas" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/install/celeritas" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCELERITAS_USE_CUDA=ON \
    -DCELERITAS_USE_Geant4=ON \
    -DCELERITAS_USE_ROOT=OFF \
    -DCELERITAS_USE_MPI=OFF \
    -DCELERITAS_USE_OpenMP=ON \
    -DCELERITAS_BUILD_TESTS=OFF \
    -DCELERITAS_BUILD_DEMOS=OFF

echo "Building Celeritas with $NPROC cores..."
make -j$NPROC

echo "Installing Celeritas..."
make install

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║          Celeritas Installation Complete!                    ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Installed to: $INSTALL_DIR/install/celeritas"
