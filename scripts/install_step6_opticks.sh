#!/bin/bash
# Step 6: Install Opticks (GPU Optical Photons)
set -e

INSTALL_DIR="/home/billy/nnbar/simulation/GEANT4_Packages"
NPROC=$(nproc)

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║          Installing Opticks (GPU Optical Photons)            ║"
echo "╚══════════════════════════════════════════════════════════════╝"

# Check CUDA
if ! command -v nvcc &> /dev/null; then
    echo "ERROR: CUDA not found. Please run install_step1_cuda.sh first."
    exit 1
fi

# Check for OptiX
OPTIX_DIR=""
for dir in /opt/optix* /usr/local/optix* ~/NVIDIA-OptiX* /home/billy/NVIDIA-OptiX*; do
    if [ -d "$dir" ] && [ -f "$dir/include/optix.h" ]; then
        OPTIX_DIR="$dir"
        break
    fi
done

if [ -z "$OPTIX_DIR" ]; then
    echo ""
    echo "═══════════════════════════════════════════════════════════════"
    echo "  OptiX SDK Not Found!"
    echo "═══════════════════════════════════════════════════════════════"
    echo ""
    echo "Opticks requires NVIDIA OptiX 7.x SDK."
    echo ""
    echo "To install OptiX:"
    echo "  1. Go to: https://developer.nvidia.com/designworks/optix/download"
    echo "  2. Create/login to NVIDIA Developer account"
    echo "  3. Download OptiX 7.7 SDK for Linux"
    echo "  4. Run: sh NVIDIA-OptiX-SDK-7.7.0-linux64-x86_64.sh --prefix=/opt/optix"
    echo ""
    echo "After installing OptiX, run this script again."
    echo ""
    exit 1
fi

echo "Found OptiX at: $OPTIX_DIR"

# Source Geant4
if [ -f "$INSTALL_DIR/install/geant4-11.2.2/bin/geant4.sh" ]; then
    source "$INSTALL_DIR/install/geant4-11.2.2/bin/geant4.sh"
    export Geant4_DIR="$INSTALL_DIR/install/geant4-11.2.2/lib/cmake/Geant4"
fi

cd "$INSTALL_DIR/src"

# Clone Opticks
if [ ! -d "opticks" ]; then
    echo "Cloning Opticks..."
    git clone --depth 1 https://bitbucket.org/simoncblyth/opticks.git
fi

# Build
mkdir -p "$INSTALL_DIR/build/opticks"
cd "$INSTALL_DIR/build/opticks"

echo "Configuring Opticks..."
cmake "$INSTALL_DIR/src/opticks" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/install/opticks" \
    -DCMAKE_BUILD_TYPE=Release \
    -DOptiX_INSTALL_DIR="$OPTIX_DIR"

echo "Building Opticks with $NPROC cores..."
make -j$NPROC

echo "Installing Opticks..."
make install

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║          Opticks Installation Complete!                      ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Installed to: $INSTALL_DIR/install/opticks"
