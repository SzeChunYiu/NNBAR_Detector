#!/bin/bash
# ============================================================================
# install_gpu_packages.sh - Master script for GPU acceleration packages
# ============================================================================
# Installs CUDA, Celeritas, Opticks, and Garfield++ for the NNBAR simulation.
#
# Installation directory: ../GEANT4_Packages/
#
# Usage: ./scripts/install_gpu_packages.sh [options]
# Options:
#   --cuda        Install CUDA Toolkit only
#   --celeritas   Install Celeritas only
#   --opticks     Install Opticks only
#   --garfield    Install Garfield++ only
#   --all         Install all packages (default)
#   --geant4      Upgrade Geant4 to 11.x
#   --help        Show this help message
# ============================================================================

set -e

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
NC='\033[0m'
BOLD='\033[1m'

# Directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PACKAGES_DIR="$(dirname "$PROJECT_DIR")/GEANT4_Packages"
PARALLEL=$(nproc 2>/dev/null || echo 4)

# Package versions
GEANT4_VERSION="11.2.1"
CELERITAS_VERSION="0.5.0"
GARFIELD_VERSION="master"

# Parse arguments
INSTALL_CUDA=false
INSTALL_CELERITAS=false
INSTALL_OPTICKS=false
INSTALL_GARFIELD=false
INSTALL_GEANT4=false
INSTALL_ALL=true

while [[ $# -gt 0 ]]; do
    case $1 in
        --cuda)      INSTALL_CUDA=true; INSTALL_ALL=false; shift ;;
        --celeritas) INSTALL_CELERITAS=true; INSTALL_ALL=false; shift ;;
        --opticks)   INSTALL_OPTICKS=true; INSTALL_ALL=false; shift ;;
        --garfield)  INSTALL_GARFIELD=true; INSTALL_ALL=false; shift ;;
        --geant4)    INSTALL_GEANT4=true; INSTALL_ALL=false; shift ;;
        --all)       INSTALL_ALL=true; shift ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo "Options:"
            echo "  --cuda        Install CUDA Toolkit"
            echo "  --celeritas   Install Celeritas"
            echo "  --opticks     Install Opticks"
            echo "  --garfield    Install Garfield++"
            echo "  --geant4      Upgrade Geant4 to 11.x"
            echo "  --all         Install all packages (default)"
            exit 0
            ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [[ "$INSTALL_ALL" == "true" ]]; then
    INSTALL_CUDA=true
    INSTALL_CELERITAS=true
    INSTALL_OPTICKS=true
    INSTALL_GARFIELD=true
fi

# Banner
echo ""
echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║${NC}${BOLD}           NNBAR Simulation - GPU Package Installation                     ${NC}${CYAN}║${NC}"
echo -e "${CYAN}║${NC}           Installing advanced physics packages for GPU acceleration         ${CYAN}║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${BLUE}Installation directory:${NC} ${PACKAGES_DIR}"
echo ""

# Create packages directory
mkdir -p "${PACKAGES_DIR}"
cd "${PACKAGES_DIR}"

# ============================================================================
# Check Prerequisites
# ============================================================================
echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}[0/5] Checking Prerequisites${NC}"
echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

# Check for NVIDIA GPU
if nvidia-smi &>/dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader | head -1)
    echo -e "  ${GREEN}✓${NC} NVIDIA GPU detected: ${GPU_NAME}"
else
    echo -e "  ${RED}✗${NC} No NVIDIA GPU detected"
    echo -e "  ${YELLOW}Note:${NC} GPU packages require an NVIDIA GPU"
    exit 1
fi

# Check for CUDA
if which nvcc &>/dev/null; then
    CUDA_VER=$(nvcc --version | grep release | awk '{print $5}' | tr -d ',')
    echo -e "  ${GREEN}✓${NC} CUDA installed: v${CUDA_VER}"
    CUDA_INSTALLED=true
else
    echo -e "  ${YELLOW}○${NC} CUDA not installed"
    CUDA_INSTALLED=false
fi

# Check Geant4 version
G4_VER=$(geant4-config --version 2>/dev/null || echo "not found")
echo -e "  ${GREEN}✓${NC} Geant4 version: ${G4_VER}"

if [[ "${G4_VER}" < "11.0" ]] && [[ "$INSTALL_CELERITAS" == "true" ]]; then
    echo -e "  ${YELLOW}!${NC} Celeritas requires Geant4 >= 11.0"
    echo -e "     Run with ${CYAN}--geant4${NC} to upgrade Geant4 first"
fi

echo ""

# ============================================================================
# Install CUDA Toolkit
# ============================================================================
if [[ "$INSTALL_CUDA" == "true" ]] && [[ "$CUDA_INSTALLED" == "false" ]]; then
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BOLD}[1/5] Installing CUDA Toolkit${NC}"
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    echo -e "  ${YELLOW}Note:${NC} CUDA installation requires sudo"
    echo ""
    echo -e "  To install CUDA manually:"
    echo -e "    ${CYAN}sudo apt install nvidia-cuda-toolkit${NC}"
    echo -e "  Or download from: ${CYAN}https://developer.nvidia.com/cuda-downloads${NC}"
    echo ""

    read -p "  Attempt automatic CUDA install? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        sudo apt update
        sudo apt install -y nvidia-cuda-toolkit
        echo -e "  ${GREEN}✓${NC} CUDA installed"
    else
        echo -e "  ${YELLOW}○${NC} Skipping CUDA installation"
    fi
    echo ""
fi

# ============================================================================
# Install Geant4 11.x
# ============================================================================
if [[ "$INSTALL_GEANT4" == "true" ]]; then
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BOLD}[2/5] Installing Geant4 ${GEANT4_VERSION}${NC}"
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    GEANT4_INSTALL="${PACKAGES_DIR}/geant4-install"

    if [[ -f "${GEANT4_INSTALL}/bin/geant4-config" ]]; then
        echo -e "  ${GREEN}✓${NC} Geant4 already installed"
    else
        echo -e "  Downloading Geant4 ${GEANT4_VERSION}..."

        wget -q --show-progress -O "geant4-v${GEANT4_VERSION}.tar.gz" \
            "https://gitlab.cern.ch/geant4/geant4/-/archive/v${GEANT4_VERSION}/geant4-v${GEANT4_VERSION}.tar.gz"

        tar -xzf "geant4-v${GEANT4_VERSION}.tar.gz"
        rm "geant4-v${GEANT4_VERSION}.tar.gz"

        mkdir -p geant4-build && cd geant4-build

        echo -e "  Configuring..."
        cmake "../geant4-v${GEANT4_VERSION}" \
            -DCMAKE_INSTALL_PREFIX="${GEANT4_INSTALL}" \
            -DCMAKE_BUILD_TYPE=Release \
            -DGEANT4_USE_GDML=ON \
            -DGEANT4_USE_QT=ON \
            -DGEANT4_USE_OPENGL_X11=ON \
            -DGEANT4_INSTALL_DATA=ON \
            > cmake.log 2>&1

        echo -e "  Building (this takes 30-60 minutes)..."
        make -j${PARALLEL} > build.log 2>&1

        echo -e "  Installing..."
        make install > install.log 2>&1

        cd "${PACKAGES_DIR}"
        rm -rf geant4-build "geant4-v${GEANT4_VERSION}"

        echo -e "  ${GREEN}✓${NC} Geant4 ${GEANT4_VERSION} installed to ${GEANT4_INSTALL}"
    fi
    echo ""
fi

# ============================================================================
# Install Celeritas
# ============================================================================
if [[ "$INSTALL_CELERITAS" == "true" ]]; then
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BOLD}[3/5] Installing Celeritas${NC}"
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    CELERITAS_INSTALL="${PACKAGES_DIR}/celeritas-install"

    if [[ -f "${CELERITAS_INSTALL}/lib/cmake/Celeritas/CeleritasConfig.cmake" ]]; then
        echo -e "  ${GREEN}✓${NC} Celeritas already installed"
    else
        echo -e "  ${YELLOW}Note:${NC} Celeritas is best installed via Spack"
        echo ""
        echo -e "  ${BOLD}Recommended installation:${NC}"
        echo -e "    ${CYAN}git clone https://github.com/spack/spack.git${NC}"
        echo -e "    ${CYAN}source spack/share/spack/setup-env.sh${NC}"
        echo -e "    ${CYAN}spack install celeritas +cuda${NC}"
        echo ""

        read -p "  Attempt to build from source? [y/N] " -n 1 -r
        echo
        if [[ $REPLY =~ ^[Yy]$ ]]; then
            echo -e "  Cloning Celeritas..."
            git clone --depth 1 --branch v${CELERITAS_VERSION} \
                https://github.com/celeritas-project/celeritas.git celeritas-src

            mkdir -p celeritas-build && cd celeritas-build

            echo -e "  Configuring..."
            cmake ../celeritas-src \
                -DCMAKE_INSTALL_PREFIX="${CELERITAS_INSTALL}" \
                -DCMAKE_BUILD_TYPE=Release \
                -DCELERITAS_USE_CUDA=ON \
                -DCELERITAS_USE_Geant4=ON \
                -DCELERITAS_USE_VecGeom=OFF \
                > cmake.log 2>&1 || {
                    echo -e "  ${RED}✗${NC} CMake failed. Check dependencies."
                    echo -e "      See: ${PACKAGES_DIR}/celeritas-build/cmake.log"
                    cd "${PACKAGES_DIR}"
                }

            if [[ -f Makefile ]]; then
                echo -e "  Building..."
                make -j${PARALLEL} > build.log 2>&1

                echo -e "  Installing..."
                make install > install.log 2>&1

                cd "${PACKAGES_DIR}"
                echo -e "  ${GREEN}✓${NC} Celeritas installed"
            fi
        else
            echo -e "  ${YELLOW}○${NC} Skipping Celeritas"
        fi
    fi
    echo ""
fi

# ============================================================================
# Install Opticks
# ============================================================================
if [[ "$INSTALL_OPTICKS" == "true" ]]; then
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BOLD}[4/5] Installing Opticks${NC}"
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    OPTICKS_HOME="${PACKAGES_DIR}/opticks"

    if [[ -d "${OPTICKS_HOME}" ]]; then
        echo -e "  ${GREEN}✓${NC} Opticks already cloned"
    else
        echo -e "  ${YELLOW}Note:${NC} Opticks requires OptiX SDK (download from NVIDIA)"
        echo -e "         https://developer.nvidia.com/designworks/optix/download"
        echo ""

        echo -e "  Cloning Opticks..."
        git clone https://bitbucket.org/simoncblyth/opticks.git "${OPTICKS_HOME}"

        echo -e "  ${GREEN}✓${NC} Opticks cloned to ${OPTICKS_HOME}"
        echo ""
        echo -e "  ${BOLD}To complete Opticks installation:${NC}"
        echo -e "    1. Download OptiX 7.x SDK from NVIDIA"
        echo -e "    2. Set environment: ${CYAN}export OPTICKS_HOME=${OPTICKS_HOME}${NC}"
        echo -e "    3. Follow: ${CYAN}https://simoncblyth.bitbucket.io/opticks/${NC}"
    fi
    echo ""
fi

# ============================================================================
# Install Garfield++
# ============================================================================
if [[ "$INSTALL_GARFIELD" == "true" ]]; then
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BOLD}[5/5] Installing Garfield++${NC}"
    echo -e "${MAGENTA}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"

    GARFIELD_INSTALL="${PACKAGES_DIR}/garfield-install"

    if [[ -f "${GARFIELD_INSTALL}/lib/cmake/Garfield/GarfieldConfig.cmake" ]]; then
        echo -e "  ${GREEN}✓${NC} Garfield++ already installed"
    else
        echo -e "  ${YELLOW}Note:${NC} Garfield++ requires ROOT"

        # Check for ROOT
        if which root-config &>/dev/null; then
            ROOT_VER=$(root-config --version)
            echo -e "  ${GREEN}✓${NC} ROOT found: v${ROOT_VER}"

            echo -e "  Cloning Garfield++..."
            git clone https://gitlab.cern.ch/garfield/garfieldpp.git garfield-src

            mkdir -p garfield-build && cd garfield-build

            echo -e "  Configuring..."
            cmake ../garfield-src \
                -DCMAKE_INSTALL_PREFIX="${GARFIELD_INSTALL}" \
                -DCMAKE_BUILD_TYPE=Release \
                -DWITH_EXAMPLES=OFF \
                > cmake.log 2>&1

            echo -e "  Building..."
            make -j${PARALLEL} > build.log 2>&1

            echo -e "  Installing..."
            make install > install.log 2>&1

            cd "${PACKAGES_DIR}"
            rm -rf garfield-build garfield-src

            echo -e "  ${GREEN}✓${NC} Garfield++ installed to ${GARFIELD_INSTALL}"
        else
            echo -e "  ${RED}✗${NC} ROOT not found. Install ROOT first:"
            echo -e "      ${CYAN}sudo apt install root-system${NC}"
            echo -e "      or: ${CYAN}https://root.cern/install/${NC}"
        fi
    fi
    echo ""
fi

# ============================================================================
# Summary
# ============================================================================
echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║${NC}${BOLD}                    Installation Summary                                   ${NC}${CYAN}║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════════════╝${NC}"
echo ""

echo -e "${BOLD}Package locations:${NC}"
echo -e "  GEANT4_Packages: ${CYAN}${PACKAGES_DIR}${NC}"
echo ""

# Check what's installed
[[ -f "${PACKAGES_DIR}/geant4-install/bin/geant4-config" ]] && \
    echo -e "  ${GREEN}✓${NC} Geant4 11.x:  ${PACKAGES_DIR}/geant4-install"

[[ -f "${PACKAGES_DIR}/celeritas-install/lib/cmake/Celeritas/CeleritasConfig.cmake" ]] && \
    echo -e "  ${GREEN}✓${NC} Celeritas:    ${PACKAGES_DIR}/celeritas-install"

[[ -d "${PACKAGES_DIR}/opticks" ]] && \
    echo -e "  ${GREEN}✓${NC} Opticks:      ${PACKAGES_DIR}/opticks"

[[ -f "${PACKAGES_DIR}/garfield-install/lib/cmake/Garfield/GarfieldConfig.cmake" ]] && \
    echo -e "  ${GREEN}✓${NC} Garfield++:   ${PACKAGES_DIR}/garfield-install"

echo ""
echo -e "${BOLD}To use these packages, configure CMake with:${NC}"
echo -e "  ${CYAN}cmake .. -DWITH_CELERITAS=ON -DCeleritas_DIR=${PACKAGES_DIR}/celeritas-install/lib/cmake/Celeritas${NC}"
echo -e "  ${CYAN}         -DWITH_OPTICKS=ON -DOPTICKS_HOME=${PACKAGES_DIR}/opticks${NC}"
echo -e "  ${CYAN}         -DWITH_GARFIELD=ON -DGarfield_DIR=${PACKAGES_DIR}/garfield-install/lib/cmake/Garfield${NC}"
echo ""
