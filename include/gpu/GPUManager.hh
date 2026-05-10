// ============================================================================
// GPUManager.hh
// Central GPU computation manager with automatic CPU fallback
// ============================================================================

#ifndef GPU_MANAGER_HH
#define GPU_MANAGER_HH

#include "config.h"
#include <string>
#include <memory>
#include <atomic>
#include <mutex>

namespace nnbar {
namespace gpu {

// ============================================================================
// GPU Device Information
// ============================================================================
struct GPUDeviceInfo {
    int deviceId = -1;
    std::string name = "No GPU";
    size_t totalMemory = 0;
    size_t freeMemory = 0;
    int computeCapabilityMajor = 0;
    int computeCapabilityMinor = 0;
    int multiProcessorCount = 0;
    bool isAvailable = false;
};

// ============================================================================
// GPU Computation Statistics
// ============================================================================
struct GPUStatistics {
    std::atomic<uint64_t> gpuEMTracks{0};
    std::atomic<uint64_t> gpuOpticalPhotons{0};
    std::atomic<uint64_t> gpuTPCElectrons{0};
    std::atomic<uint64_t> gpuHits{0};
    std::atomic<uint64_t> cpuFallbackCount{0};
    double gpuTimeMs{0};
    double cpuTimeMs{0};
    std::mutex timeMutex;

    void Reset() {
        gpuEMTracks = 0;
        gpuOpticalPhotons = 0;
        gpuTPCElectrons = 0;
        gpuHits = 0;
        cpuFallbackCount = 0;
        std::lock_guard<std::mutex> lock(timeMutex);
        gpuTimeMs = 0;
        cpuTimeMs = 0;
    }

    void AddGPUTime(double ms) {
        std::lock_guard<std::mutex> lock(timeMutex);
        gpuTimeMs += ms;
    }

    void AddCPUTime(double ms) {
        std::lock_guard<std::mutex> lock(timeMutex);
        cpuTimeMs += ms;
    }

    double GetGPUFraction() const {
        double total = gpuTimeMs + cpuTimeMs;
        return total > 0 ? gpuTimeMs / total : 0.0;
    }
};

// ============================================================================
// GPU Computation Mode
// ============================================================================
enum class ComputeMode {
    GPU_ONLY,       // Fail if GPU not available
    GPU_PREFERRED,  // Use GPU if available, fallback to CPU
    CPU_ONLY,       // Force CPU computation
    AUTO            // Automatically choose based on workload
};

// ============================================================================
// GPUManager - Singleton class for managing GPU computation
// ============================================================================
class GPUManager {
public:
    // Singleton access
    static GPUManager& Instance();

    // Delete copy/move
    GPUManager(const GPUManager&) = delete;
    GPUManager& operator=(const GPUManager&) = delete;

    // ========================================================================
    // Initialization
    // ========================================================================

    // Initialize GPU subsystem
    // Returns true if GPU is available
    bool Initialize();

    // Shutdown GPU subsystem
    void Shutdown();

    // ========================================================================
    // Configuration
    // ========================================================================

    // Set computation mode
    void SetComputeMode(ComputeMode mode) { computeMode_ = mode; }
    ComputeMode GetComputeMode() const { return computeMode_; }

    // Enable/disable specific GPU features
    void SetGPUEMPhysics(bool enable) { enableEMPhysics_ = enable; }
    void SetGPUOpticalPhotons(bool enable) { enableOptical_ = enable; }
    void SetGPUTPCDrift(bool enable) { enableTPCDrift_ = enable; }

    bool IsGPUEMPhysicsEnabled() const { return enableEMPhysics_ && isGPUAvailable_; }
    bool IsGPUOpticalEnabled() const { return enableOptical_ && isGPUAvailable_; }
    bool IsGPUTPCDriftEnabled() const { return enableTPCDrift_ && isGPUAvailable_; }

    // ========================================================================
    // Status Queries
    // ========================================================================

    bool IsInitialized() const { return isInitialized_; }
    bool IsGPUAvailable() const { return isGPUAvailable_; }
    const GPUDeviceInfo& GetDeviceInfo() const { return deviceInfo_; }
    const GPUStatistics& GetStatistics() const { return statistics_; }
    GPUStatistics& GetStatistics() { return statistics_; }

    // ========================================================================
    // Decision Helpers
    // ========================================================================

    // Should we use GPU for this computation?
    bool ShouldUseGPU() const;

    // Report fallback to CPU
    void ReportCPUFallback(const std::string& reason);

    // ========================================================================
    // Memory Management
    // ========================================================================

    // Allocate GPU memory (returns nullptr if not available)
    void* AllocateGPU(size_t bytes);

    // Free GPU memory
    void FreeGPU(void* ptr);

    // Get available GPU memory
    size_t GetAvailableGPUMemory() const;

    // ========================================================================
    // Synchronization
    // ========================================================================

    // Synchronize GPU operations
    void Synchronize();

    // ========================================================================
    // Reporting
    // ========================================================================

    // Print GPU status and statistics
    void PrintStatus() const;
    void PrintStatistics() const;

private:
    GPUManager() = default;
    ~GPUManager();

    // State
    bool isInitialized_ = false;
    bool isGPUAvailable_ = false;
    ComputeMode computeMode_ = ComputeMode::GPU_PREFERRED;

    // Feature flags
    bool enableEMPhysics_ = true;
    bool enableOptical_ = true;
    bool enableTPCDrift_ = true;

    // Device info
    GPUDeviceInfo deviceInfo_;

    // Statistics
    GPUStatistics statistics_;

    // Thread safety
    mutable std::mutex mutex_;
};

}  // namespace gpu
}  // namespace nnbar

#endif // GPU_MANAGER_HH
