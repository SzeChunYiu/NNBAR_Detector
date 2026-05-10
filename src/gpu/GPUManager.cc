// ============================================================================
// GPUManager.cc
// Central GPU computation manager with automatic CPU fallback
// ============================================================================

#include "gpu/GPUManager.hh"
#include "G4ios.hh"
#include <iomanip>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#define HAS_CUDA 1
#else
#define HAS_CUDA 0
#endif

namespace nnbar {
namespace gpu {

// ============================================================================
// Singleton Instance
// ============================================================================

GPUManager& GPUManager::Instance() {
    static GPUManager instance;
    return instance;
}

GPUManager::~GPUManager() {
    Shutdown();
}

// ============================================================================
// Initialization
// ============================================================================

bool GPUManager::Initialize() {
    std::lock_guard<std::mutex> lock(mutex_);

    if (isInitialized_) {
        return isGPUAvailable_;
    }

    G4cout << "\n╔══════════════════════════════════════════════════════════════════╗" << G4endl;
    G4cout << "║              GPU COMPUTATION MANAGER INITIALIZATION              ║" << G4endl;
    G4cout << "╠══════════════════════════════════════════════════════════════════╣" << G4endl;

#if HAS_CUDA
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);

    if (err != cudaSuccess || deviceCount == 0) {
        G4cout << "║  Status: No CUDA GPU detected                                   ║" << G4endl;
        G4cout << "║  Mode: CPU fallback enabled                                     ║" << G4endl;
        isGPUAvailable_ = false;
    } else {
        // Get device properties
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);

        deviceInfo_.deviceId = 0;
        deviceInfo_.name = prop.name;
        deviceInfo_.totalMemory = prop.totalGlobalMem;
        deviceInfo_.computeCapabilityMajor = prop.major;
        deviceInfo_.computeCapabilityMinor = prop.minor;
        deviceInfo_.multiProcessorCount = prop.multiProcessorCount;
        deviceInfo_.isAvailable = true;

        // Get free memory
        size_t freeMem, totalMem;
        cudaMemGetInfo(&freeMem, &totalMem);
        deviceInfo_.freeMemory = freeMem;

        isGPUAvailable_ = true;

        G4cout << "║  Status: GPU AVAILABLE                                          ║" << G4endl;
        G4cout << "║  Device: " << std::left << std::setw(56) << deviceInfo_.name << "║" << G4endl;
        G4cout << "║  Memory: " << std::setw(6) << (deviceInfo_.totalMemory / (1024*1024*1024))
               << " GB total, " << std::setw(6) << (deviceInfo_.freeMemory / (1024*1024*1024))
               << " GB free" << std::setw(26) << " ║" << G4endl;
        G4cout << "║  Compute: SM " << deviceInfo_.computeCapabilityMajor << "."
               << deviceInfo_.computeCapabilityMinor
               << " (" << deviceInfo_.multiProcessorCount << " multiprocessors)"
               << std::setw(23) << " ║" << G4endl;
    }
#else
    G4cout << "║  Status: CUDA not compiled                                       ║" << G4endl;
    G4cout << "║  Mode: CPU-only mode                                             ║" << G4endl;
    isGPUAvailable_ = false;
#endif

    // Print enabled features
    G4cout << "╠══════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  GPU Features:                                                   ║" << G4endl;
    G4cout << "║    EM Physics (Celeritas): " << std::left << std::setw(38)
           << (enableEMPhysics_ && isGPUAvailable_ ? "ENABLED" : "DISABLED/FALLBACK") << "║" << G4endl;
    G4cout << "║    Optical Photons:        " << std::setw(38)
           << (enableOptical_ && isGPUAvailable_ ? "ENABLED" : "DISABLED/FALLBACK") << "║" << G4endl;
    G4cout << "║    TPC Drift:              " << std::setw(38)
           << (enableTPCDrift_ && isGPUAvailable_ ? "ENABLED" : "DISABLED/FALLBACK") << "║" << G4endl;
    G4cout << std::right;
    G4cout << "╚══════════════════════════════════════════════════════════════════╝\n" << G4endl;

    isInitialized_ = true;
    statistics_.Reset();

    return isGPUAvailable_;
}

void GPUManager::Shutdown() {
    std::lock_guard<std::mutex> lock(mutex_);

    if (!isInitialized_) return;

#if HAS_CUDA
    if (isGPUAvailable_) {
        cudaDeviceReset();
    }
#endif

    isInitialized_ = false;
    isGPUAvailable_ = false;
}

// ============================================================================
// Decision Helpers
// ============================================================================

bool GPUManager::ShouldUseGPU() const {
    switch (computeMode_) {
        case ComputeMode::GPU_ONLY:
            return isGPUAvailable_;
        case ComputeMode::GPU_PREFERRED:
            return isGPUAvailable_;
        case ComputeMode::CPU_ONLY:
            return false;
        case ComputeMode::AUTO:
            return isGPUAvailable_;
        default:
            return isGPUAvailable_;
    }
}

void GPUManager::ReportCPUFallback(const std::string& reason) {
    statistics_.cpuFallbackCount++;
#ifdef DEBUG_VERBOSE
    G4cout << "[GPUManager] CPU fallback: " << reason << G4endl;
#endif
}

// ============================================================================
// Memory Management
// ============================================================================

void* GPUManager::AllocateGPU(size_t bytes) {
#if HAS_CUDA
    if (!isGPUAvailable_) return nullptr;

    void* ptr = nullptr;
    cudaError_t err = cudaMalloc(&ptr, bytes);
    if (err != cudaSuccess) {
        G4cerr << "[GPUManager] cudaMalloc failed: " << cudaGetErrorString(err) << G4endl;
        return nullptr;
    }
    return ptr;
#else
    (void)bytes;
    return nullptr;
#endif
}

void GPUManager::FreeGPU(void* ptr) {
#if HAS_CUDA
    if (ptr && isGPUAvailable_) {
        cudaFree(ptr);
    }
#else
    (void)ptr;
#endif
}

size_t GPUManager::GetAvailableGPUMemory() const {
#if HAS_CUDA
    if (!isGPUAvailable_) return 0;

    size_t freeMem, totalMem;
    cudaMemGetInfo(&freeMem, &totalMem);
    return freeMem;
#else
    return 0;
#endif
}

// ============================================================================
// Synchronization
// ============================================================================

void GPUManager::Synchronize() {
#if HAS_CUDA
    if (isGPUAvailable_) {
        cudaDeviceSynchronize();
    }
#endif
}

// ============================================================================
// Reporting
// ============================================================================

void GPUManager::PrintStatus() const {
    G4cout << "\n[GPUManager] Status:" << G4endl;
    G4cout << "  Initialized: " << (isInitialized_ ? "Yes" : "No") << G4endl;
    G4cout << "  GPU Available: " << (isGPUAvailable_ ? "Yes" : "No") << G4endl;
    if (isGPUAvailable_) {
        G4cout << "  Device: " << deviceInfo_.name << G4endl;
        G4cout << "  Memory: " << (deviceInfo_.freeMemory / (1024*1024)) << " MB free / "
               << (deviceInfo_.totalMemory / (1024*1024)) << " MB total" << G4endl;
    }
    G4cout << "  Compute Mode: ";
    switch (computeMode_) {
        case ComputeMode::GPU_ONLY: G4cout << "GPU_ONLY"; break;
        case ComputeMode::GPU_PREFERRED: G4cout << "GPU_PREFERRED"; break;
        case ComputeMode::CPU_ONLY: G4cout << "CPU_ONLY"; break;
        case ComputeMode::AUTO: G4cout << "AUTO"; break;
    }
    G4cout << G4endl;
}

void GPUManager::PrintStatistics() const {
    G4cout << "\n╔══════════════════════════════════════════════════════════════════╗" << G4endl;
    G4cout << "║                    GPU COMPUTATION STATISTICS                    ║" << G4endl;
    G4cout << "╠══════════════════════════════════════════════════════════════════╣" << G4endl;
    G4cout << "║  EM Tracks (GPU):        " << std::setw(12) << statistics_.gpuEMTracks.load()
           << std::setw(29) << " ║" << G4endl;
    G4cout << "║  Optical Photons (GPU):  " << std::setw(12) << statistics_.gpuOpticalPhotons.load()
           << std::setw(29) << " ║" << G4endl;
    G4cout << "║  TPC Electrons (GPU):    " << std::setw(12) << statistics_.gpuTPCElectrons.load()
           << std::setw(29) << " ║" << G4endl;
    G4cout << "║  Hits Aggregated (GPU):  " << std::setw(12) << statistics_.gpuHits.load()
           << std::setw(29) << " ║" << G4endl;
    G4cout << "╠──────────────────────────────────────────────────────────────────╣" << G4endl;
    G4cout << "║  CPU Fallbacks:          " << std::setw(12) << statistics_.cpuFallbackCount.load()
           << std::setw(29) << " ║" << G4endl;
    G4cout << "║  GPU Time:               " << std::setw(10) << std::fixed << std::setprecision(2)
           << statistics_.gpuTimeMs << " ms" << std::setw(26) << " ║" << G4endl;
    G4cout << "║  CPU Time:               " << std::setw(10) << std::fixed << std::setprecision(2)
           << statistics_.cpuTimeMs << " ms" << std::setw(26) << " ║" << G4endl;
    G4cout << "║  GPU Utilization:        " << std::setw(10) << std::fixed << std::setprecision(1)
           << (statistics_.GetGPUFraction() * 100.0) << " %" << std::setw(27) << " ║" << G4endl;
    G4cout << "╚══════════════════════════════════════════════════════════════════╝\n" << G4endl;
}

}  // namespace gpu
}  // namespace nnbar
