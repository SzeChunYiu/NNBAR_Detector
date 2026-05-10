#!/bin/bash
# Step 3: Download and build Geant4 11.2.2
set -e

INSTALL_DIR="/home/billy/nnbar/simulation/GEANT4_Packages"
GEANT4_VERSION="11.2.2"
NPROC=$(nproc)

echo "╔══════════════════════════════════════════════════════════════╗"
echo "║          Installing Geant4 ${GEANT4_VERSION}                            ║"
echo "╚══════════════════════════════════════════════════════════════╝"

mkdir -p "$INSTALL_DIR"/{src,build,install}
cd "$INSTALL_DIR/src"

# Download Geant4
if [ ! -f "geant4-v${GEANT4_VERSION}.tar.gz" ]; then
    echo "Downloading Geant4 ${GEANT4_VERSION}..."
    wget --progress=bar:force https://gitlab.cern.ch/geant4/geant4/-/archive/v${GEANT4_VERSION}/geant4-v${GEANT4_VERSION}.tar.gz
fi

# Extract
if [ ! -d "geant4-v${GEANT4_VERSION}" ]; then
    echo "Extracting Geant4..."
    tar -xzf geant4-v${GEANT4_VERSION}.tar.gz
fi

# Build
mkdir -p "$INSTALL_DIR/build/geant4"
cd "$INSTALL_DIR/build/geant4"

echo "Configuring Geant4..."
cmake "$INSTALL_DIR/src/geant4-v${GEANT4_VERSION}" \
    -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR/install/geant4-${GEANT4_VERSION}" \
    -DCMAKE_BUILD_TYPE=Release \
    -DGEANT4_INSTALL_DATA=ON \
    -DGEANT4_USE_GDML=ON \
    -DGEANT4_USE_QT=ON \
    -DGEANT4_USE_OPENGL_X11=ON \
    -DGEANT4_USE_RAYTRACER_X11=ON \
    -DGEANT4_BUILD_MULTITHREADED=ON \
    -DGEANT4_USE_SYSTEM_EXPAT=ON \
    -DGEANT4_USE_SYSTEM_ZLIB=ON

echo ""
echo "Building Geant4 with $NPROC cores..."
echo "This will take 20-40 minutes..."
echo ""
make -j$NPROC

echo "Installing Geant4..."
make install

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║          Geant4 ${GEANT4_VERSION} Installation Complete!              ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""
echo "Installed to: $INSTALL_DIR/install/geant4-${GEANT4_VERSION}"
