#!/bin/bash
# Complete Opticks Build using Native System with all required functions
# This script provides all stub functions and builds remaining packages

INSTALL_DIR="/home/billy/nnbar/simulation/GEANT4_Packages"
OPTICKS_SRC="$INSTALL_DIR/src/opticks"
OPTICKS_PREFIX="$INSTALL_DIR/install/opticks"

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║     Completing Opticks Build (Native System)                 ║"
echo "╚══════════════════════════════════════════════════════════════╝"

# Setup CUDA
export CUDA_HOME="/usr/local/cuda-12.4"
export OPTICKS_CUDA_PREFIX="$CUDA_HOME"
export PATH="$CUDA_HOME/bin:$PATH"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:$LD_LIBRARY_PATH"

# Setup OptiX
export OPTICKS_OPTIX_PREFIX="$INSTALL_DIR/install/optix"
export OptiX_INSTALL_DIR="$OPTICKS_OPTIX_PREFIX"

# Setup Opticks environment
export OPTICKS_HOME="$OPTICKS_SRC"
export OPTICKS_PREFIX="$OPTICKS_PREFIX"
export OPTICKS_BUILD_WITH_CUDA="ON"

# GPU Compute Capability for RTX A3000 (Ampere)
export OPTICKS_COMPUTE_CAPABILITY=86

# Source Geant4 first
if [ -f "$INSTALL_DIR/install/geant4-11.2.2/bin/geant4.sh" ]; then
    source "$INSTALL_DIR/install/geant4-11.2.2/bin/geant4.sh"
    export Geant4_DIR="$INSTALL_DIR/install/geant4-11.2.2/lib/cmake/Geant4"
    echo "[OK] Geant4 environment loaded ($(geant4-config --version 2>/dev/null || echo 'version'))"
else
    echo "[ERROR] Geant4 not found"
    exit 1
fi

# === Define ALL required opticks- functions ===

# Core functions
olocal-() { : ; }
opticks-() { : ; }
oe-() { : ; }

opticks-source() { echo "$OPTICKS_SRC/opticks.bash" ; }
opticks-home() { echo "$OPTICKS_HOME" ; }
opticks-prefix() { echo "$OPTICKS_PREFIX" ; }
opticks-sdir() { echo "$OPTICKS_HOME" ; }
opticks-bdir() { echo "$OPTICKS_PREFIX/build" ; }

opticks-buildtype() { echo Release ; }
opticks-optix-prefix() { echo "$OPTICKS_OPTIX_PREFIX" ; }
opticks-compute-capability() { echo "$OPTICKS_COMPUTE_CAPABILITY" ; }
opticks-build-with-cuda() { echo "${OPTICKS_BUILD_WITH_CUDA:-ON}" ; }

opticks-git-clone() { git clone "$@" ; }
opticks-curl() { curl -L -O "$@" ; }

opticks-cmake-generator() { echo "Unix Makefiles" ; }

# Geocache functions
opticks-user-home() { echo "${HOME}" ; }
opticks-geocache-prefix() { echo "${OPTICKS_GEOCACHE_PREFIX:-$HOME/.opticks}" ; }
opticks-rngcache-prefix() { echo "${OPTICKS_RNGCACHE_PREFIX:-$HOME/.opticks}" ; }
opticks-geocachedir() { echo "$(opticks-geocache-prefix)/geocache" ; }
opticks-rngcachedir() { echo "$(opticks-rngcache-prefix)/rngcache" ; }

# Export all functions
export -f olocal-
export -f opticks-
export -f oe-
export -f opticks-source
export -f opticks-home
export -f opticks-prefix
export -f opticks-sdir
export -f opticks-bdir
export -f opticks-buildtype
export -f opticks-optix-prefix
export -f opticks-compute-capability
export -f opticks-build-with-cuda
export -f opticks-git-clone
export -f opticks-curl
export -f opticks-cmake-generator
export -f opticks-user-home
export -f opticks-geocache-prefix
export -f opticks-rngcache-prefix
export -f opticks-geocachedir
export -f opticks-rngcachedir

# Source the existing Opticks bashrc if available
if [ -f "$OPTICKS_PREFIX/bashrc" ]; then
    echo "[INFO] Sourcing $OPTICKS_PREFIX/bashrc"
    source "$OPTICKS_PREFIX/bashrc"
fi

# Navigate to Opticks source
cd "$OPTICKS_SRC"

# Display environment
echo ""
echo "Environment:"
echo "  OPTICKS_HOME:         $OPTICKS_HOME"
echo "  OPTICKS_PREFIX:       $OPTICKS_PREFIX"
echo "  OPTICKS_CUDA_PREFIX:  $OPTICKS_CUDA_PREFIX"
echo "  OPTICKS_OPTIX_PREFIX: $OPTICKS_OPTIX_PREFIX"
echo "  OPTICKS_COMPUTE_CAPABILITY: $OPTICKS_COMPUTE_CAPABILITY"
echo "  Geant4_DIR:           $Geant4_DIR"
echo ""

# Source om.bash
echo "Sourcing om.bash..."
source "$OPTICKS_SRC/om.bash" 2>&1 | grep -v "^$"

# Check if om functions exist
if ! type om-subs &>/dev/null; then
    echo "[ERROR] om.bash functions not loaded properly"
    exit 1
fi
echo "[OK] om.bash functions loaded"

# List packages
echo ""
echo "Checking which packages need building..."
om-subs 2>/dev/null | tail -20

echo ""
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
echo "  Building remaining Opticks packages"
echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

# Build function
build_with_om() {
    local pkg=$1
    echo ""
    echo "Building: $pkg"
    echo "---"
    cd "$OPTICKS_SRC/$pkg"
    om-conf 2>&1 | tail -20
    if [ $? -eq 0 ]; then
        om-make 2>&1 | tail -10
        return $?
    fi
    return 1
}

# Check and build CSGOptiX
echo ""
echo "Step 1/3: CSGOptiX"
if [ ! -f "$OPTICKS_PREFIX/lib/libCSGOptiX.so" ]; then
    cd "$OPTICKS_SRC/CSGOptiX"
    echo "Directory: $(pwd)"
    build_with_om "CSGOptiX" || echo "[WARN] CSGOptiX build issues"
else
    echo "[OK] CSGOptiX already exists"
fi

# Check and build U4
echo ""
echo "Step 2/3: U4"
if [ ! -f "$OPTICKS_PREFIX/lib/libU4.so" ]; then
    cd "$OPTICKS_SRC/u4"
    echo "Directory: $(pwd)"
    build_with_om "u4" || echo "[WARN] U4 build issues"
else
    echo "[OK] U4 already exists"
fi

# Check and build G4CX
echo ""
echo "Step 3/3: G4CX"
if [ ! -f "$OPTICKS_PREFIX/lib/libG4CX.so" ]; then
    cd "$OPTICKS_SRC/g4cx"
    echo "Directory: $(pwd)"
    build_with_om "g4cx" || echo "[WARN] G4CX build issues"
else
    echo "[OK] G4CX already exists"
fi

# Verify
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
    echo ""
    echo "To use Opticks in your simulation:"
    echo "  export OPTICKS_HOME=$OPTICKS_PREFIX"
    echo "  cd /home/billy/nnbar/simulation/NNBAR_Detector/build"
    echo "  cmake .. -DWITH_OPTICKS=ON"
    echo "  make -j\$(nproc)"
else
    echo "[INFO] Build incomplete ($success/3 libraries)"
    echo ""
    echo "Check the build logs above for specific errors."
    echo "Common issues:"
    echo "  - Missing CMake modules (try adding CMAKE_MODULE_PATH)"
    echo "  - Missing dependencies (PLog, Custom4, etc.)"
    echo ""
    echo "You may need to build interactively:"
    echo "  cd $OPTICKS_SRC"
    echo "  source om.bash"
    echo "  cd CSGOptiX && om"
    echo "  cd ../u4 && om"
    echo "  cd ../g4cx && om"
fi
