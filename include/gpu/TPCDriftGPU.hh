// ============================================================================
// TPCDriftGPU.hh
// GPU-accelerated TPC electron drift simulation with CPU fallback
// ============================================================================

#ifndef TPC_DRIFT_GPU_HH
#define TPC_DRIFT_GPU_HH

#include "config.h"
#include <vector>
#include <memory>
#include <mutex>

namespace nnbar {
namespace gpu {

// ============================================================================
// TPC Drift Data Structures
// ============================================================================

struct IonizationCluster {
    float x, y, z;           // Position (mm)
    float time;              // Time (ns)
    int electrons;           // Number of primary electrons
    int eventId;
    int trackId;
};

struct DriftedElectron {
    float x, y, z;           // Final position on readout plane (mm)
    float time;              // Arrival time (ns)
    int clusterId;           // Original cluster ID
    int eventId;
};

struct TPCPadHit {
    int padId;               // Pad identifier
    int moduleId;            // TPC module
    int layerId;             // Layer within module
    float charge;            // Collected charge (arbitrary units)
    float time;              // Hit time (ns)
    int eventId;
};

struct TPCFieldPoint {
    float Ex, Ey, Ez;        // Electric field components (V/cm)
};

// ============================================================================
// GPU TPC Drift Engine
// ============================================================================

class TPCDriftGPU {
public:
    // Singleton access
    static TPCDriftGPU& Instance();

    // Delete copy/move
    TPCDriftGPU(const TPCDriftGPU&) = delete;
    TPCDriftGPU& operator=(const TPCDriftGPU&) = delete;

    // ========================================================================
    // Initialization
    // ========================================================================

    bool Initialize();
    void Shutdown();

    bool IsInitialized() const { return isInitialized_; }
    bool IsGPUAvailable() const { return isGPUEnabled_; }

    // ========================================================================
    // Configuration
    // ========================================================================

    // Set drift properties
    void SetDriftVelocity(float vd) { driftVelocity_ = vd; }  // mm/us
    void SetDiffusionL(float dL) { diffusionL_ = dL; }        // mm/sqrt(cm)
    void SetDiffusionT(float dT) { diffusionT_ = dT; }        // mm/sqrt(cm)
    void SetElectricField(float E) { electricField_ = E; }    // V/cm
    void SetGasGain(float gain) { gasGain_ = gain; }

    // Load electric field map from file
    bool LoadFieldMap(const std::string& filename);

    // Set uniform field (for testing)
    void SetUniformField(float Ex, float Ey, float Ez);

    // ========================================================================
    // Processing
    // ========================================================================

    // Add ionization cluster for processing
    void AddIonizationCluster(const IonizationCluster& cluster);

    // Process all clusters (call at end of event)
    void ProcessEvent();

    // Get pad hits from last processed event
    const std::vector<TPCPadHit>& GetPadHits() const { return padHits_; }

    // Get statistics
    int GetNClusters() const { return nClusters_; }
    int GetTotalInputElectrons() const { return totalInputElectrons_; }
    int GetTotalCollectedElectrons() const { return totalCollectedElectrons_; }
    float GetProcessingTimeMs() const { return processingTimeMs_; }

    // Clear event data
    void ClearEvent();

    // ========================================================================
    // Statistics
    // ========================================================================

    void PrintStatistics() const;

private:
    TPCDriftGPU() = default;
    ~TPCDriftGPU();

    // GPU implementation
    void ProcessEventGPU();

    // CPU fallback implementation
    void ProcessEventCPU();

    // Drift single electron (CPU)
    DriftedElectron DriftElectron(float x, float y, float z, float t, int clusterId, int eventId);

    // Collect electrons on pads
    void CollectOnPads(const std::vector<DriftedElectron>& electrons);

    // State
    bool isInitialized_ = false;
    bool isGPUEnabled_ = false;

    // TPC parameters
    float driftVelocity_ = 5.0f;       // mm/us (typical for Ar/CO2)
    float diffusionL_ = 0.02f;          // mm/sqrt(cm) longitudinal
    float diffusionT_ = 0.02f;          // mm/sqrt(cm) transverse
    float electricField_ = 500.0f;      // V/cm
    float gasGain_ = 1000.0f;           // Amplification factor
    float driftLength_ = 1000.0f;       // mm (max drift distance)

    // Field map (if loaded)
    bool hasFieldMap_ = false;
    std::vector<TPCFieldPoint> fieldMap_;
    int fieldMapNx_ = 0, fieldMapNy_ = 0, fieldMapNz_ = 0;
    float fieldMapDx_ = 0, fieldMapDy_ = 0, fieldMapDz_ = 0;

    // Event data
    std::vector<IonizationCluster> clusters_;
    std::vector<DriftedElectron> driftedElectrons_;
    std::vector<TPCPadHit> padHits_;

    // Statistics
    int nClusters_ = 0;
    int totalInputElectrons_ = 0;
    int totalCollectedElectrons_ = 0;
    float processingTimeMs_ = 0.0f;
    int eventsProcessedGPU_ = 0;
    int eventsProcessedCPU_ = 0;

    // GPU buffers
#if HAS_CUDA
    void* d_clusters_ = nullptr;
    void* d_electrons_ = nullptr;
    void* d_padHits_ = nullptr;
    void* d_fieldMap_ = nullptr;
    void* d_rngStates_ = nullptr;
    size_t maxElectrons_ = 10000000;
#endif

    // Thread safety
    std::mutex mutex_;
};

}  // namespace gpu
}  // namespace nnbar

#endif // TPC_DRIFT_GPU_HH
