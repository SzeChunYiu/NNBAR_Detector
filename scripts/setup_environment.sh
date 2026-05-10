#!/bin/bash
# ============================================================================
# setup_environment.sh - One-command environment setup for NNBAR simulation
# ============================================================================
#
# This script automatically detects and sets up the correct environment
# for running the NNBAR detector simulation, including GPU packages.
#
# USAGE:
#   source scripts/setup_environment.sh
#
# WHAT IT DOES:
#   1. Detects cluster environment (module system, Spack, etc.)
#   2. Loads required modules or sets paths
#   3. Configures GPU libraries if available
#   4. Exports all necessary environment variables
#
# SUPPORTED ENVIRONMENTS:
#   - Local workstation (Ubuntu/CentOS)
#   - SLURM clusters with Lmod
#   - CVMFS (CERN/HEP clusters)
#   - Spack-based environments
#
# ============================================================================

# Exit on error for script development, but not when sourced
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: This script must be sourced, not executed!"
    echo "Usage: source $0"
    exit 1
fi

# Colors (only if terminal supports them)
if [[ -t 1 ]]; then
    RED='\033[0;31m'
    GREEN='\033[0;32m'
    YELLOW='\033[1;33m'
    BLUE='\033[0;34m'
    CYAN='\033[0;36m'
    NC='\033[0m'
    BOLD='\033[1m'
else
    RED='' GREEN='' YELLOW='' BLUE='' CYAN='' NC='' BOLD=''
fi

# Find script and project directories
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
PACKAGES_DIR="$(dirname "$PROJECT_DIR")/GEANT4_Packages"
OS_NAME="$(uname -s)"

echo ""
echo -e "${CYAN}╔══════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${CYAN}║${NC}${BOLD}              NNBAR Simulation Environment Setup                          ${NC}${CYAN}║${NC}"
echo -e "${CYAN}╚══════════════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# ============================================================================
# Detect Environment Type
# ============================================================================
ENV_TYPE="unknown"
echo -e "${BLUE}[1/4] Detecting environment...${NC}"

# Check for CVMFS (CERN/HEP standard)
if [[ -d "/cvmfs/geant4.cern.ch" ]]; then
    ENV_TYPE="cvmfs"
    echo -e "  ${GREEN}✓${NC} CVMFS detected (CERN/HEP cluster)"

# Check for Lmod module system
elif type module &>/dev/null; then
    ENV_TYPE="lmod"
    echo -e "  ${GREEN}✓${NC} Lmod module system detected"

# Check for Spack
elif type spack &>/dev/null || [[ -d "$HOME/spack" ]]; then
    ENV_TYPE="spack"
    echo -e "  ${GREEN}✓${NC} Spack detected"

# Check for local installation
elif [[ -d "${PACKAGES_DIR}" ]]; then
    ENV_TYPE="local"
    echo -e "  ${GREEN}✓${NC} Local installation detected"

else
    ENV_TYPE="system"
    echo -e "  ${YELLOW}○${NC} Using system packages"
fi

echo ""

# ============================================================================
# Setup Based on Environment Type
# ============================================================================
echo -e "${BLUE}[2/4] Loading packages...${NC}"

case "$ENV_TYPE" in
    # ========================================================================
    # CVMFS (CERN clusters, LHC grid sites)
    # ========================================================================
    cvmfs)
        echo -e "  Loading from CVMFS..."

        # Source the LCG view (includes Geant4, ROOT, etc.)
        LCG_VIEW="/cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc13-opt"
        if [[ -f "${LCG_VIEW}/setup.sh" ]]; then
            source "${LCG_VIEW}/setup.sh"
            echo -e "  ${GREEN}✓${NC} LCG view loaded"
        fi

        # Or source Geant4 directly
        G4_CVMFS="/cvmfs/geant4.cern.ch/geant4/11.2/x86_64-el9-gcc13-optdeb/bin/geant4.sh"
        if [[ -f "${G4_CVMFS}" ]]; then
            source "${G4_CVMFS}"
            echo -e "  ${GREEN}✓${NC} Geant4 loaded from CVMFS"
        fi
        ;;

    # ========================================================================
    # Lmod (SLURM clusters like Compute Canada, NERSC, etc.)
    # ========================================================================
    lmod)
        echo -e "  Loading modules..."

        # Try to load common module names (vary by cluster)
        for mod in gcc cmake cuda geant4 root; do
            if module avail $mod 2>&1 | grep -q $mod; then
                module load $mod 2>/dev/null && echo -e "  ${GREEN}✓${NC} Loaded $mod"
            fi
        done

        # Check for Celeritas module
        if module avail celeritas 2>&1 | grep -q celeritas; then
            module load celeritas && echo -e "  ${GREEN}✓${NC} Loaded celeritas"
        fi
        ;;

    # ========================================================================
    # Spack
    # ========================================================================
    spack)
        echo -e "  Loading from Spack..."

        # Ensure Spack is set up
        if [[ -f "$HOME/spack/share/spack/setup-env.sh" ]]; then
            source "$HOME/spack/share/spack/setup-env.sh"
        fi

        # Load packages if installed
        for pkg in geant4 celeritas root; do
            if spack find $pkg &>/dev/null; then
                spack load $pkg && echo -e "  ${GREEN}✓${NC} Loaded $pkg from Spack"
            fi
        done
        ;;

    # ========================================================================
    # Local installation in GEANT4_Packages
    # ========================================================================
    local)
        echo -e "  Loading from ${PACKAGES_DIR}..."

        # Geant4
        G4_SETUP=""
        if [[ -f "${PACKAGES_DIR}/geant4-install/bin/geant4.sh" ]]; then
            G4_SETUP="${PACKAGES_DIR}/geant4-install/bin/geant4.sh"
        else
            for candidate in "${PACKAGES_DIR}"/install/geant4-*/bin/geant4.sh; do
                if [[ -f "$candidate" ]]; then
                    G4_SETUP="$candidate"
                    break
                fi
            done
        fi

        if [[ -n "$G4_SETUP" ]]; then
            source "$G4_SETUP"
            G4_PREFIX="$(cd "$(dirname "$G4_SETUP")/.." && pwd)"
            export Geant4_DIR="$G4_PREFIX/lib/cmake/Geant4"
            export CMAKE_PREFIX_PATH="$G4_PREFIX:${CMAKE_PREFIX_PATH:-}"
            echo -e "  ${GREEN}✓${NC} Geant4 loaded"
        fi

        # Celeritas
        if [[ -d "${PACKAGES_DIR}/celeritas-install" ]]; then
            export Celeritas_DIR="${PACKAGES_DIR}/celeritas-install/lib/cmake/Celeritas"
            export LD_LIBRARY_PATH="${PACKAGES_DIR}/celeritas-install/lib:$LD_LIBRARY_PATH"
            echo -e "  ${GREEN}✓${NC} Celeritas configured"
        fi

        # Opticks
        if [[ -d "${PACKAGES_DIR}/opticks" ]]; then
            export OPTICKS_HOME="${PACKAGES_DIR}/opticks"
            echo -e "  ${GREEN}✓${NC} Opticks configured"
        fi

        # Garfield++
        if [[ -d "${PACKAGES_DIR}/garfield-install" ]]; then
            export Garfield_DIR="${PACKAGES_DIR}/garfield-install/lib/cmake/Garfield"
            export LD_LIBRARY_PATH="${PACKAGES_DIR}/garfield-install/lib:$LD_LIBRARY_PATH"
            echo -e "  ${GREEN}✓${NC} Garfield++ configured"
        fi
        ;;

    # ========================================================================
    # System packages (Ubuntu apt, etc.)
    # ========================================================================
    system)
        echo -e "  Using system packages..."

        # Check if Geant4 is in system path
        if which geant4-config &>/dev/null; then
            G4_PREFIX=$(geant4-config --prefix)
            if [[ -f "${G4_PREFIX}/bin/geant4.sh" ]]; then
                source "${G4_PREFIX}/bin/geant4.sh"
            fi
            echo -e "  ${GREEN}✓${NC} System Geant4 found"
        else
            echo -e "  ${YELLOW}!${NC} Geant4 not found in PATH"
        fi
        ;;
esac

echo ""

# ============================================================================
# Setup Project-Specific Paths
# ============================================================================
echo -e "${BLUE}[3/4] Setting up project paths...${NC}"

# External dependencies bundled with project
EXTERNAL_DIR="${PROJECT_DIR}/external"

if [[ -d "${EXTERNAL_DIR}/arrow-install" ]]; then
    export Arrow_DIR="${EXTERNAL_DIR}/arrow-install/lib/cmake/Arrow"
    export Parquet_DIR="${EXTERNAL_DIR}/arrow-install/lib/cmake/Parquet"
    export LD_LIBRARY_PATH="${EXTERNAL_DIR}/arrow-install/lib:$LD_LIBRARY_PATH"
    echo -e "  ${GREEN}✓${NC} Arrow/Parquet configured"
elif [[ "$OS_NAME" != "Darwin" && -d "${EXTERNAL_DIR}/arrow-install-linux" ]]; then
    export Arrow_DIR="${EXTERNAL_DIR}/arrow-install-linux/lib/cmake/Arrow"
    export Parquet_DIR="${EXTERNAL_DIR}/arrow-install-linux/lib/cmake/Parquet"
    export LD_LIBRARY_PATH="${EXTERNAL_DIR}/arrow-install-linux/lib:$LD_LIBRARY_PATH"
    echo -e "  ${GREEN}✓${NC} Arrow/Parquet configured"
fi

if [[ -z "${Arrow_DIR:-}" ]]; then
    for prefix in "${CONDA_PREFIX:-}" "$HOME/miniforge3/envs/nnbar_env" "$HOME/miniconda3/envs/nnbar_env"; do
        if [[ -n "$prefix" && -f "$prefix/lib/cmake/Arrow/ArrowConfig.cmake" && -f "$prefix/lib/cmake/Parquet/ParquetConfig.cmake" ]]; then
            export Arrow_DIR="$prefix/lib/cmake/Arrow"
            export Parquet_DIR="$prefix/lib/cmake/Parquet"
            export CMAKE_PREFIX_PATH="$prefix:${CMAKE_PREFIX_PATH:-}"
            echo -e "  ${GREEN}✓${NC} Arrow/Parquet configured from $prefix"
            break
        fi
    done
fi

if [[ -d "${EXTERNAL_DIR}/spdlog-install" ]]; then
    export spdlog_DIR="${EXTERNAL_DIR}/spdlog-install/lib/cmake/spdlog"
    export LD_LIBRARY_PATH="${EXTERNAL_DIR}/spdlog-install/lib:$LD_LIBRARY_PATH"
    echo -e "  ${GREEN}✓${NC} spdlog configured"
elif [[ "$OS_NAME" != "Darwin" && -d "${EXTERNAL_DIR}/spdlog-install-linux" ]]; then
    export spdlog_DIR="${EXTERNAL_DIR}/spdlog-install-linux/lib/cmake/spdlog"
    export LD_LIBRARY_PATH="${EXTERNAL_DIR}/spdlog-install-linux/lib:$LD_LIBRARY_PATH"
    echo -e "  ${GREEN}✓${NC} spdlog configured"
fi

if [[ -z "${spdlog_DIR:-}" ]]; then
    for prefix in "${CONDA_PREFIX:-}" "$HOME/miniforge3/envs/nnbar_env" "$HOME/miniconda3/envs/nnbar_env"; do
        if [[ -n "$prefix" && -f "$prefix/lib/cmake/spdlog/spdlogConfig.cmake" ]]; then
            export spdlog_DIR="$prefix/lib/cmake/spdlog"
            export CMAKE_PREFIX_PATH="$prefix:${CMAKE_PREFIX_PATH:-}"
            echo -e "  ${GREEN}✓${NC} spdlog configured from $prefix"
            break
        fi
    done
fi

if [[ -d "${EXTERNAL_DIR}/json" ]]; then
    export nlohmann_json_DIR="${EXTERNAL_DIR}/json"
    echo -e "  ${GREEN}✓${NC} nlohmann_json configured"
fi

# Export project directory
export NNBAR_PROJECT_DIR="${PROJECT_DIR}"
export NNBAR_DATA_DIR="${PROJECT_DIR}/data"
export NNBAR_CONFIG_DIR="${PROJECT_DIR}/config"

echo ""

# ============================================================================
# Verify Setup
# ============================================================================
echo -e "${BLUE}[4/4] Verifying setup...${NC}"

SETUP_OK=true

# Check Geant4
if which geant4-config &>/dev/null; then
    G4_VER=$(geant4-config --version)
    echo -e "  ${GREEN}✓${NC} Geant4 ${G4_VER}"
else
    echo -e "  ${RED}✗${NC} Geant4 not found"
    SETUP_OK=false
fi

# Check CUDA
if which nvcc &>/dev/null; then
    CUDA_VER=$(nvcc --version | grep release | awk '{print $5}' | tr -d ',')
    echo -e "  ${GREEN}✓${NC} CUDA ${CUDA_VER}"
else
    echo -e "  ${YELLOW}○${NC} CUDA not found (GPU acceleration disabled)"
fi

# Check GPU
if nvidia-smi &>/dev/null; then
    GPU_NAME=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)
    echo -e "  ${GREEN}✓${NC} GPU: ${GPU_NAME}"
else
    echo -e "  ${YELLOW}○${NC} No NVIDIA GPU detected"
fi

echo ""

# ============================================================================
# Summary
# ============================================================================
if [[ "$SETUP_OK" == "true" ]]; then
    echo -e "${GREEN}╔══════════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${GREEN}║${NC}${BOLD}                    Environment Ready!                                    ${NC}${GREEN}║${NC}"
    echo -e "${GREEN}╚══════════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "  ${BOLD}Build:${NC}  ./scripts/build.sh"
    echo -e "  ${BOLD}Run:${NC}    ./build/nnbar-detector-simulation -m macro/signal/run_signal.mac"
    echo ""
else
    echo -e "${RED}╔══════════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${RED}║${NC}${BOLD}                    Setup Incomplete                                      ${NC}${RED}║${NC}"
    echo -e "${RED}╚══════════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
    echo -e "  Some packages are missing. To install:"
    echo -e "    ${CYAN}./scripts/install_gpu_packages.sh${NC}"
    echo ""
fi

# Create a simple cmake configuration helper
export CMAKE_PREFIX_PATH="${EXTERNAL_DIR}/arrow-install:${EXTERNAL_DIR}/spdlog-install:${CMAKE_PREFIX_PATH}"
