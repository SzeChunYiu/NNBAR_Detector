#!/bin/bash
# ============================================================================
# NNBAR Detector Simulation - Complete GPU Environment Installation
# ============================================================================
# This script installs:
#   1. CUDA Toolkit 12.x
#   2. Geant4 11.2.2 (with Qt, OpenGL, GDML)
#   3. Celeritas (GPU EM physics)
#   4. Opticks + OptiX 7.x (GPU optical photons)
#   5. Garfield++ (realistic TPC simulation)
# ============================================================================

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Installation directory
INSTALL_DIR="${1:-/home/billy/nnbar/simulation/GEANT4_Packages}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
NPROC=$(nproc)

echo -e "${CYAN}"
echo "╔══════════════════════════════════════════════════════════════════════════╗"
echo "║     NNBAR Detector - Complete GPU Environment Installation               ║"
echo "╚══════════════════════════════════════════════════════════════════════════╝"
echo -e "${NC}"

echo -e "${BLUE}Installation directory: ${INSTALL_DIR}${NC}"
echo -e "${BLUE}Build parallelism: ${NPROC} cores${NC}"
echo ""

# Create directories
mkdir -p "$INSTALL_DIR"/{src,build,install}
cd "$INSTALL_DIR"

# ============================================================================
# Helper functions
# ============================================================================

log_step() {
    echo -e "\n${GREEN}════════════════════════════════════════════════════════════════${NC}"
    echo -e "${GREEN}  $1${NC}"
    echo -e "${GREEN}════════════════════════════════════════════════════════════════${NC}\n"
}

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

check_command() {
    if command -v "$1" &> /dev/null; then
        return 0
    else
        return 1
    fi
}

# ============================================================================
# Step 1: Install CUDA Toolkit
# ============================================================================

install_cuda() {
    log_step "Step 1: Installing CUDA Toolkit 12.x"

    if check_command nvcc; then
        CUDA_VERSION=$(nvcc --version | grep release | awk '{print $5}' | cut -d',' -f1)
        log_info "CUDA already installed: version $CUDA_VERSION"
        return 0
    fi

    log_info "Installing CUDA Toolkit..."

    # Add NVIDIA package repository
    wget -q https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-keyring_1.1-1_all.deb
    sudo dpkg -i cuda-keyring_1.1-1_all.deb
    rm cuda-keyring_1.1-1_all.deb

    sudo apt-get update
    sudo apt-get install -y cuda-toolkit-12-4

    # Add to PATH
    echo 'export PATH=/usr/local/cuda/bin:$PATH' >> ~/.bashrc
    echo 'export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
    export PATH=/usr/local/cuda/bin:$PATH
    export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

    log_info "CUDA Toolkit installed successfully"
}

# ============================================================================
# Step 2: Install/Upgrade Geant4 11.2.2
# ============================================================================

install_geant4() {
    log_step "Step 2: Installing Geant4 11.2.2"

    GEANT4_VERSION="11.2.2"
    GEANT4_DIR="$INSTALL_DIR/install/geant4-${GEANT4_VERSION}"

    if [ -f "$GEANT4_DIR/bin/geant4-config" ]; then
        log_info "Geant4 ${GEANT4_VERSION} already installed at $GEANT4_DIR"
        return 0
    fi

    cd "$INSTALL_DIR/src"

    # Download Geant4
    if [ ! -f "geant4-v${GEANT4_VERSION}.tar.gz" ]; then
        log_info "Downloading Geant4 ${GEANT4_VERSION}..."
        wget -q --show-progress https://gitlab.cern.ch/geant4/geant4/-/archive/v${GEANT4_VERSION}/geant4-v${GEANT4_VERSION}.tar.gz
    fi

    # Extract
    if [ ! -d "geant4-v${GEANT4_VERSION}" ]; then
        log_info "Extracting Geant4..."
        tar -xzf geant4-v${GEANT4_VERSION}.tar.gz
    fi

    # Build
    mkdir -p "$INSTALL_DIR/build/geant4"
    cd "$INSTALL_DIR/build/geant4"

    log_info "Configuring Geant4 (this may take a moment)..."
    cmake "$INSTALL_DIR/src/geant4-v${GEANT4_VERSION}" \
        -DCMAKE_INSTALL_PREFIX="$GEANT4_DIR" \
        -DCMAKE_BUILD_TYPE=Release \
        -DGEANT4_INSTALL_DATA=ON \
        -DGEANT4_USE_GDML=ON \
        -DGEANT4_USE_QT=ON \
        -DGEANT4_USE_OPENGL_X11=ON \
        -DGEANT4_USE_RAYTRACER_X11=ON \
        -DGEANT4_BUILD_MULTITHREADED=ON \
        -DGEANT4_USE_SYSTEM_EXPAT=ON \
        -DGEANT4_USE_SYSTEM_ZLIB=ON

    log_info "Building Geant4 (using $NPROC cores - this will take 20-40 minutes)..."
    make -j$NPROC

    log_info "Installing Geant4..."
    make install

    log_info "Geant4 ${GEANT4_VERSION} installed successfully"

    # Create environment script
    cat > "$INSTALL_DIR/setup_geant4.sh" << EOF
# Geant4 ${GEANT4_VERSION} Environment
source ${GEANT4_DIR}/bin/geant4.sh
export Geant4_DIR=${GEANT4_DIR}/lib/cmake/Geant4
EOF
}

# ============================================================================
# Step 3: Install Celeritas
# ============================================================================

install_celeritas() {
    log_step "Step 3: Installing Celeritas (GPU EM Physics)"

    CELERITAS_DIR="$INSTALL_DIR/install/celeritas"

    if [ -f "$CELERITAS_DIR/lib/cmake/Celeritas/CeleritasConfig.cmake" ]; then
        log_info "Celeritas already installed at $CELERITAS_DIR"
        return 0
    fi

    # Check CUDA
    if ! check_command nvcc; then
        log_warn "CUDA not found - skipping Celeritas installation"
        log_warn "Run this script again after installing CUDA"
        return 1
    fi

    cd "$INSTALL_DIR/src"

    # Clone Celeritas
    if [ ! -d "celeritas" ]; then
        log_info "Cloning Celeritas repository..."
        git clone --depth 1 https://github.com/celeritas-project/celeritas.git
    fi

    # Build
    mkdir -p "$INSTALL_DIR/build/celeritas"
    cd "$INSTALL_DIR/build/celeritas"

    # Source Geant4 environment
    source "$INSTALL_DIR/setup_geant4.sh" 2>/dev/null || true

    log_info "Configuring Celeritas..."
    cmake "$INSTALL_DIR/src/celeritas" \
        -DCMAKE_INSTALL_PREFIX="$CELERITAS_DIR" \
        -DCMAKE_BUILD_TYPE=Release \
        -DCELERITAS_USE_CUDA=ON \
        -DCELERITAS_USE_Geant4=ON \
        -DCELERITAS_USE_ROOT=OFF \
        -DCELERITAS_USE_MPI=OFF \
        -DCELERITAS_USE_OpenMP=ON \
        -DCELERITAS_BUILD_TESTS=OFF \
        -DCELERITAS_BUILD_DEMOS=OFF

    log_info "Building Celeritas (using $NPROC cores)..."
    make -j$NPROC

    log_info "Installing Celeritas..."
    make install

    log_info "Celeritas installed successfully"

    # Add to environment script
    cat >> "$INSTALL_DIR/setup_gpu.sh" << EOF
# Celeritas
export Celeritas_DIR=${CELERITAS_DIR}/lib/cmake/Celeritas
export LD_LIBRARY_PATH=${CELERITAS_DIR}/lib:\$LD_LIBRARY_PATH
EOF
}

# ============================================================================
# Step 4: Install OptiX 7.x and Opticks
# ============================================================================

install_opticks() {
    log_step "Step 4: Installing Opticks (GPU Optical Photons)"

    OPTICKS_DIR="$INSTALL_DIR/install/opticks"

    if [ -f "$OPTICKS_DIR/lib/libOpticks.so" ] || [ -f "$OPTICKS_DIR/lib64/libOpticks.so" ]; then
        log_info "Opticks already installed at $OPTICKS_DIR"
        return 0
    fi

    # Check CUDA
    if ! check_command nvcc; then
        log_warn "CUDA not found - skipping Opticks installation"
        return 1
    fi

    cd "$INSTALL_DIR/src"

    # Check for OptiX
    OPTIX_DIR=""
    for dir in /opt/optix* /usr/local/optix* ~/NVIDIA-OptiX*; do
        if [ -d "$dir" ] && [ -f "$dir/include/optix.h" ]; then
            OPTIX_DIR="$dir"
            break
        fi
    done

    if [ -z "$OPTIX_DIR" ]; then
        log_warn "OptiX SDK not found!"
        log_warn "Please download OptiX 7.x from: https://developer.nvidia.com/designworks/optix/download"
        log_warn "Install it to /opt/optix or ~/NVIDIA-OptiX-SDK-7.x"
        log_warn "Then run this script again."

        # Create placeholder script
        cat > "$INSTALL_DIR/install_optix_manually.sh" << 'EOF'
#!/bin/bash
echo "To install OptiX:"
echo "1. Go to https://developer.nvidia.com/designworks/optix/download"
echo "2. Download OptiX 7.7 or later (requires NVIDIA developer account)"
echo "3. Run: sudo sh NVIDIA-OptiX-SDK-7.7.0-linux64-x86_64.sh --prefix=/opt/optix"
echo "4. Then run install_all_gpu.sh again"
EOF
        chmod +x "$INSTALL_DIR/install_optix_manually.sh"
        return 1
    fi

    log_info "Found OptiX at: $OPTIX_DIR"

    # Clone Opticks
    if [ ! -d "opticks" ]; then
        log_info "Cloning Opticks repository..."
        git clone --depth 1 https://bitbucket.org/simoncblyth/opticks.git
    fi

    # Build Opticks
    mkdir -p "$INSTALL_DIR/build/opticks"
    cd "$INSTALL_DIR/build/opticks"

    # Source environments
    source "$INSTALL_DIR/setup_geant4.sh" 2>/dev/null || true

    log_info "Configuring Opticks..."
    cmake "$INSTALL_DIR/src/opticks" \
        -DCMAKE_INSTALL_PREFIX="$OPTICKS_DIR" \
        -DCMAKE_BUILD_TYPE=Release \
        -DOptiX_INSTALL_DIR="$OPTIX_DIR" \
        -DOPTICKS_GEANT4_HOME="$INSTALL_DIR/install/geant4-11.2.2"

    log_info "Building Opticks (using $NPROC cores)..."
    make -j$NPROC

    log_info "Installing Opticks..."
    make install

    log_info "Opticks installed successfully"

    # Add to environment script
    cat >> "$INSTALL_DIR/setup_gpu.sh" << EOF
# Opticks
export OPTICKS_HOME=${OPTICKS_DIR}
export Opticks_DIR=${OPTICKS_DIR}/lib/cmake/Opticks
export LD_LIBRARY_PATH=${OPTICKS_DIR}/lib:\$LD_LIBRARY_PATH
EOF
}

# ============================================================================
# Step 5: Install Garfield++
# ============================================================================

install_garfield() {
    log_step "Step 5: Installing Garfield++ (Realistic TPC Simulation)"

    GARFIELD_DIR="$INSTALL_DIR/install/garfield"

    if [ -f "$GARFIELD_DIR/lib/cmake/Garfield/GarfieldConfig.cmake" ]; then
        log_info "Garfield++ already installed at $GARFIELD_DIR"
        return 0
    fi

    # Check ROOT
    if ! check_command root-config; then
        log_warn "ROOT not found - skipping Garfield++ installation"
        log_warn "Install ROOT first: sudo apt install root-system"
        return 1
    fi

    cd "$INSTALL_DIR/src"

    # Clone Garfield++
    if [ ! -d "garfieldpp" ]; then
        log_info "Cloning Garfield++ repository..."
        git clone --depth 1 https://gitlab.cern.ch/garfield/garfieldpp.git
    fi

    # Build
    mkdir -p "$INSTALL_DIR/build/garfield"
    cd "$INSTALL_DIR/build/garfield"

    # Source Geant4 environment
    source "$INSTALL_DIR/setup_geant4.sh" 2>/dev/null || true

    log_info "Configuring Garfield++..."
    cmake "$INSTALL_DIR/src/garfieldpp" \
        -DCMAKE_INSTALL_PREFIX="$GARFIELD_DIR" \
        -DCMAKE_BUILD_TYPE=Release \
        -DWITH_GEANT4=ON \
        -DWITH_EXAMPLES=OFF

    log_info "Building Garfield++ (using $NPROC cores)..."
    make -j$NPROC

    log_info "Installing Garfield++..."
    make install

    log_info "Garfield++ installed successfully"

    # Add to environment script
    cat >> "$INSTALL_DIR/setup_gpu.sh" << EOF
# Garfield++
export GARFIELD_INSTALL=${GARFIELD_DIR}
export Garfield_DIR=${GARFIELD_DIR}/lib/cmake/Garfield
export HEED_DATABASE=${GARFIELD_DIR}/share/Heed/database
export LD_LIBRARY_PATH=${GARFIELD_DIR}/lib:\$LD_LIBRARY_PATH
EOF
}

# ============================================================================
# Step 6: Create Master Environment Script
# ============================================================================

create_environment() {
    log_step "Step 6: Creating Environment Setup"

    cat > "$INSTALL_DIR/setup_env.sh" << EOF
#!/bin/bash
# ============================================================================
# NNBAR Detector Simulation - GPU Environment
# ============================================================================
# Source this file before building/running the simulation:
#   source ${INSTALL_DIR}/setup_env.sh
# ============================================================================

# CUDA
if [ -d "/usr/local/cuda" ]; then
    export PATH=/usr/local/cuda/bin:\$PATH
    export LD_LIBRARY_PATH=/usr/local/cuda/lib64:\$LD_LIBRARY_PATH
fi

# Geant4 11.2.2
if [ -f "${INSTALL_DIR}/install/geant4-11.2.2/bin/geant4.sh" ]; then
    source ${INSTALL_DIR}/install/geant4-11.2.2/bin/geant4.sh
    export Geant4_DIR=${INSTALL_DIR}/install/geant4-11.2.2/lib/cmake/Geant4
fi

# Celeritas
if [ -d "${INSTALL_DIR}/install/celeritas" ]; then
    export Celeritas_DIR=${INSTALL_DIR}/install/celeritas/lib/cmake/Celeritas
    export LD_LIBRARY_PATH=${INSTALL_DIR}/install/celeritas/lib:\$LD_LIBRARY_PATH
fi

# Opticks
if [ -d "${INSTALL_DIR}/install/opticks" ]; then
    export OPTICKS_HOME=${INSTALL_DIR}/install/opticks
    export Opticks_DIR=${INSTALL_DIR}/install/opticks/lib/cmake/Opticks
    export LD_LIBRARY_PATH=${INSTALL_DIR}/install/opticks/lib:\$LD_LIBRARY_PATH
fi

# Garfield++
if [ -d "${INSTALL_DIR}/install/garfield" ]; then
    export GARFIELD_INSTALL=${INSTALL_DIR}/install/garfield
    export Garfield_DIR=${INSTALL_DIR}/install/garfield/lib/cmake/Garfield
    export HEED_DATABASE=${INSTALL_DIR}/install/garfield/share/Heed/database
    export LD_LIBRARY_PATH=${INSTALL_DIR}/install/garfield/lib:\$LD_LIBRARY_PATH
fi

echo "NNBAR GPU Environment loaded:"
echo "  Geant4:    \$(geant4-config --version 2>/dev/null || echo 'not found')"
echo "  CUDA:      \$(nvcc --version 2>/dev/null | grep release | awk '{print \$5}' | cut -d',' -f1 || echo 'not found')"
echo "  Celeritas: \$([ -d \"\$Celeritas_DIR\" ] && echo 'available' || echo 'not found')"
echo "  Opticks:   \$([ -d \"\$Opticks_DIR\" ] && echo 'available' || echo 'not found')"
echo "  Garfield:  \$([ -d \"\$Garfield_DIR\" ] && echo 'available' || echo 'not found')"
EOF

    chmod +x "$INSTALL_DIR/setup_env.sh"

    log_info "Environment script created: $INSTALL_DIR/setup_env.sh"
}

# ============================================================================
# Step 7: Rebuild NNBAR Simulation
# ============================================================================

rebuild_simulation() {
    log_step "Step 7: Rebuilding NNBAR Simulation with GPU Support"

    cd "$PROJECT_DIR/build"

    # Source environment
    source "$INSTALL_DIR/setup_env.sh"

    log_info "Cleaning previous build..."
    rm -rf CMakeCache.txt CMakeFiles/

    log_info "Configuring with GPU support..."
    cmake .. \
        -DWITH_CELERITAS=ON \
        -DWITH_OPTICKS=ON \
        -DWITH_GARFIELD=ON \
        -DWITH_DASHBOARD=ON \
        -DWITH_SCINTILLATION=ON

    log_info "Building simulation..."
    make -j$NPROC

    log_info "Build complete!"
}

# ============================================================================
# Main Installation Flow
# ============================================================================

main() {
    echo -e "${YELLOW}This will install:${NC}"
    echo "  1. CUDA Toolkit 12.4"
    echo "  2. Geant4 11.2.2 (upgrade from 10.7.4)"
    echo "  3. Celeritas (GPU EM physics - 10-100x speedup)"
    echo "  4. Opticks (GPU optical photons - 50-200x speedup)"
    echo "  5. Garfield++ (realistic TPC simulation)"
    echo ""
    echo -e "${YELLOW}Estimated time: 1-2 hours${NC}"
    echo -e "${YELLOW}Disk space required: ~15 GB${NC}"
    echo ""

    read -p "Continue with installation? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Installation cancelled."
        exit 0
    fi

    # Install system dependencies
    log_step "Installing System Dependencies"
    sudo apt-get update
    sudo apt-get install -y \
        build-essential \
        cmake \
        git \
        wget \
        libxerces-c-dev \
        libexpat1-dev \
        zlib1g-dev \
        libgl1-mesa-dev \
        libglu1-mesa-dev \
        libxt-dev \
        libxmu-dev \
        libxi-dev \
        qtbase5-dev \
        libqt5opengl5-dev \
        libgsl-dev \
        libboost-dev

    # Run installation steps
    install_cuda
    install_geant4
    install_garfield
    install_celeritas
    install_opticks
    create_environment

    echo ""
    log_step "Installation Summary"

    source "$INSTALL_DIR/setup_env.sh" 2>/dev/null || true

    echo -e "${GREEN}╔══════════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║                    Installation Complete!                                ║${NC}"
    echo -e "${GREEN}╚══════════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "${CYAN}To use the GPU-accelerated simulation:${NC}"
    echo ""
    echo "  1. Source the environment:"
    echo -e "     ${YELLOW}source $INSTALL_DIR/setup_env.sh${NC}"
    echo ""
    echo "  2. Rebuild the simulation:"
    echo -e "     ${YELLOW}cd $PROJECT_DIR/build${NC}"
    echo -e "     ${YELLOW}rm -rf CMakeCache.txt CMakeFiles/${NC}"
    echo -e "     ${YELLOW}cmake .. -DWITH_CELERITAS=ON -DWITH_OPTICKS=ON -DWITH_GARFIELD=ON${NC}"
    echo -e "     ${YELLOW}make -j\$(nproc)${NC}"
    echo ""
    echo "  3. Run the simulation:"
    echo -e "     ${YELLOW}./nnbar-detector-simulation macro/signal/run_signal.mac${NC}"
    echo ""
}

# Run main
main "$@"
