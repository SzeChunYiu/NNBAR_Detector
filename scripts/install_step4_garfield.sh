#!/bin/bash
# Step 4: Install Garfield++ (requires ROOT - already installed)
set -e

INSTALL_DIR="/home/billy/nnbar/simulation/GEANT4_Packages"
NPROC=$(nproc)

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║          Installing Garfield++ (TPC Simulation)              ║"
echo "╚══════════════════════════════════════════════════════════════╝"

# Source Geant4
if [ -f "$INSTALL_DIR/install/geant4-11.2.2/bin/geant4.sh" ]; then
    source "$INSTALL_DIR/install/geant4-11.2.2/bin/geant4.sh"
    export Geant4_DIR="$INSTALL_DIR/install/geant4-11.2.2/lib/cmake/Geant4"
fi

cd "$INSTALL_DIR/src"

# Clone Garfield++
if [ ! -d "garfieldpp" ]; then
    echo "Cloning Garfield++..."
    git clone --depth 1 https://gitlab.cern.ch/garfield/garfieldpp.git
fi

# Build
mkdir -p "$INSTALL_DIR/build/garfield"
cd "$INSTALL_DIR/build/garfield"

echo "Configuring Garfield++..."
cmake "$INSTALL_DIR/src/garfieldpp" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/install/garfield" \
    -DCMAKE_BUILD_TYPE=Release \
    -DWITH_GEANT4=ON \
    -DWITH_EXAMPLES=OFF

echo "Building Garfield++ with $NPROC cores..."
make -j$NPROC

echo "Installing Garfield++..."
make install

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║          Garfield++ Installation Complete!                   ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Installed to: $INSTALL_DIR/install/garfield"
