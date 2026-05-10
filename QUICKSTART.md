# NNBAR Detector Simulation - Quick Start Guide

## 🚀 One-Command Setup

```bash
# Clone and build
git clone <repository-url> NNBAR_Detector
cd NNBAR_Detector

# Step 1: Setup environment (detects your system automatically)
source scripts/setup_environment.sh

# Step 2: Install dependencies (Arrow, JSON, logging)
./scripts/setup_dependencies.sh

# Step 3: Build
./scripts/build.sh

# Step 4: Run
./build/nnbar-detector-simulation -m macro/signal/run_signal.mac
```

## 📋 Requirements

| Package | Version | Required |
|---------|---------|----------|
| CMake | ≥ 3.16 | ✅ Yes |
| Geant4 | ≥ 10.7 | ✅ Yes |
| GCC/Clang | C++17 | ✅ Yes |
| CUDA | ≥ 11.0 | ❌ Optional (GPU) |

## ⚡ GPU Acceleration (Optional)

For 10-100x faster simulations:

```bash
# Install GPU packages
./scripts/install_gpu_packages.sh --celeritas

# Build with GPU support
cmake .. -DWITH_CELERITAS=ON
make -j$(nproc)
```

## 🖥️ Cluster-Specific Instructions

### SLURM Clusters (Compute Canada, NERSC, etc.)
```bash
# Load modules first
module load gcc cmake cuda geant4

# Then follow standard build steps
source scripts/setup_environment.sh
```

### CVMFS Clusters (CERN, LHC Grid)
```bash
# Environment automatically uses CVMFS
source scripts/setup_environment.sh
# Geant4 loaded from /cvmfs/geant4.cern.ch
```

### Spack Users
```bash
# Install via Spack
spack install geant4 +qt +opengl
spack install celeritas +cuda

# Load and build
spack load geant4 celeritas
source scripts/setup_environment.sh
```

## 🔧 Build Options

| Option | Default | Description |
|--------|---------|-------------|
| `WITH_CELERITAS` | OFF | GPU EM physics acceleration |
| `WITH_OPTICKS` | OFF | GPU optical photon propagation |
| `WITH_GARFIELD` | OFF | Garfield++ TPC simulation |
| `WITH_SCINTILLATION` | OFF | Enable optical photons |
| `WITH_DASHBOARD` | OFF | Qt monitoring interface |
| `MCPL_BUILD` | OFF | Use MCPL particle source |

Example:
```bash
cmake .. -DWITH_CELERITAS=ON -DWITH_SCINTILLATION=ON
```

On macOS, either activate an environment that provides native Arrow/Parquet and
spdlog packages or use `source scripts/setup_environment.sh`; the checked-in
`external/*-linux` bundles are Linux-only.

## 📁 Directory Structure

```
NNBAR_Detector/
├── src/              # Source files
├── include/          # Header files
├── external/         # Bundled dependencies (Arrow, MCPL)
├── config/           # Configuration files
├── macro/            # Geant4 macro files
├── scripts/          # Setup and install scripts
└── build/            # Build output (created by you)
```

## ❓ Troubleshooting

### "Geant4 not found"
```bash
# Check if Geant4 is installed
which geant4-config

# If using modules
module load geant4

# If installed locally
source /path/to/geant4/bin/geant4.sh
```

### "Arrow/Parquet not found"
```bash
# Run the dependency setup
./scripts/setup_dependencies.sh
```

### "CUDA not found" (when using GPU options)
```bash
# Check CUDA installation
which nvcc

# On Ubuntu
sudo apt install nvidia-cuda-toolkit

# On cluster
module load cuda
```

## 📧 Support

Issues: https://github.com/your-repo/issues
