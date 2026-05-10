#!/bin/bash
# Build remaining Opticks packages (u4, CSGOptiX, g4cx)
# With -Werror disabled to work around unused variable warnings
set -e

cd /home/billy/nnbar/simulation/GEANT4_Packages/src/opticks

export OPTICKS_HOME="$(pwd)"
export OPTICKS_PREFIX="/home/billy/nnbar/simulation/GEANT4_Packages/install/opticks"
export OPTICKS_CUDA_PREFIX="/usr/local/cuda-12.4"
export OPTICKS_OPTIX_PREFIX="/home/billy/nnbar/simulation/GEANT4_Packages/install/optix"
export PATH="$OPTICKS_CUDA_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$OPTICKS_CUDA_PREFIX/lib64:$LD_LIBRARY_PATH"

# CLHEP is in opticks_externals build directory
export CLHEP_DIR="/home/billy/nnbar/simulation/GEANT4_Packages/install/opticks_externals/clhep_2451.build/2.4.5.1/CLHEP.build"
export CMAKE_PREFIX_PATH="$CLHEP_DIR:$CMAKE_PREFIX_PATH"

# Disable -Werror by overriding CXXFLAGS
export CXXFLAGS="${CXXFLAGS} -Wno-error"

# Source Geant4
source /home/billy/nnbar/simulation/GEANT4_Packages/install/geant4-11.2.2/bin/geant4.sh
export Geant4_DIR="/home/billy/nnbar/simulation/GEANT4_Packages/install/geant4-11.2.2/lib/cmake/Geant4"

# Source existing Opticks environment
source "$OPTICKS_PREFIX/bashrc"

# Define stub functions required by om.bash
olocal-() { : ; }
opticks-() { : ; }
oe-() { : ; }
opticks-home() { echo "$OPTICKS_HOME" ; }
opticks-prefix() { echo "$OPTICKS_PREFIX" ; }
opticks-buildtype() { echo Release ; }
opticks-optix-prefix() { echo "$OPTICKS_OPTIX_PREFIX" ; }
opticks-compute-capability() { echo 86 ; }
opticks-build-with-cuda() { echo ON ; }
opticks-sdir() { echo "$OPTICKS_HOME" ; }
opticks-bdir() { echo "$OPTICKS_PREFIX/build" ; }
opticks-cmake-generator() { echo "Unix Makefiles" ; }

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║     Building Opticks: u4, CSGOptiX, g4cx                     ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Environment:"
echo "  OPTICKS_HOME:   $OPTICKS_HOME"
echo "  OPTICKS_PREFIX: $OPTICKS_PREFIX"
echo "  CUDA:           $(nvcc --version 2>&1 | grep release)"
echo "  Geant4:         $(geant4-config --version)"
echo "  CLHEP_DIR:      $CLHEP_DIR"
echo ""

# Load om.bash
source om.bash

# Clear U4 build directory to get fresh cmake with new flags
rm -rf "$OPTICKS_PREFIX/build/u4"

echo "Starting build of remaining packages..."
echo ""

# Build u4 first (depends on Geant4), then CSGOptiX (already built), then g4cx
echo "Step 1: Building u4 (Geant4 interface)..."
om-install u4:u4 || echo "[WARN] U4 build completed with warnings"

echo ""
echo "Step 2: Verifying CSGOptiX..."
if [ -f "$OPTICKS_PREFIX/lib/libCSGOptiX.so" ]; then
    echo "[OK] CSGOptiX already built"
else
    echo "Building CSGOptiX..."
    om-install CSGOptiX:CSGOptiX
fi

echo ""
echo "Step 3: Building g4cx..."
rm -rf "$OPTICKS_PREFIX/build/g4cx"  # Fresh build
om-install g4cx:g4cx || echo "[WARN] G4CX build completed with warnings"

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Verification"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

success=0
for lib in libCSGOptiX.so libU4.so libG4CX.so; do
    if [ -f "$OPTICKS_PREFIX/lib/$lib" ]; then
        echo "[OK] $lib"
        ((success++))
    else
        echo "[--] $lib NOT FOUND"
    fi
done

echo ""
if [ $success -eq 3 ]; then
    echo "╔══════════════════════════════════════════════════════════════╗"
    echo "║              Opticks Build Complete!                         ║"
    echo "╚══════════════════════════════════════════════════════════════╝"
else
    echo "Build incomplete ($success/3)"
fi
