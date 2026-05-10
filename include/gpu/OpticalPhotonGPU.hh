// ============================================================================
// OpticalPhotonGPU.hh
// GPU-accelerated optical photon simulation with CPU fallback
// ============================================================================

#ifndef OPTICAL_PHOTON_GPU_HH
#define OPTICAL_PHOTON_GPU_HH

#include "config.h"
#include <vector>
#include <memory>
#include <mutex>

namespace nnbar {
namespace gpu {

// ============================================================================
// Optical Photon Data Structures
// ============================================================================

struct ScintillationPoint {
    float x, y, z;           // Position (mm)
    float energyDeposit;     // Energy deposit (MeV)
    float time;              // Time (ns)
    int materialId;          // Material index
    int eventId;
};

struct OpticalPhoton {
    float x, y, z;           // Position (mm)
    float dx, dy, dz;        // Direction (normalized)
    float wavelength;        // Wavelength (nm)
    float time;              // Time (ns)
    int eventId;
    bool isAlive;
};

struct PMTHit {
    int pmtId;
    int moduleId;
    float time;              // Hit time (ns)
    float wavelength;        // Photon wavelength (nm)
    int eventId;
    int photonCount;
};

// ============================================================================
// GPU Optical Photon Engine
// ============================================================================

class OpticalPhotonGPU {
public:
    // Singleton access
    static OpticalPhotonGPU& Instance();

    // Delete copy/move
    OpticalPhotonGPU(const OpticalPhotonGPU&) = delete;
    OpticalPhotonGPU& operator=(const OpticalPhotonGPU&) = delete;

    // ========================================================================
    // Initialization
    // ========================================================================

    // Initialize the engine
    bool Initialize();

    // Shutdown
    void Shutdown();

    // Check if initialized and GPU available
    bool IsInitialized() const { return isInitialized_; }
    bool IsGPUAvailable() const { return isGPUEnabled_; }

    // ========================================================================
    // Configuration
    // ========================================================================

    // Set scintillator properties
    void SetScintillatorYield(float photonsPerMeV) { scintYield_ = photonsPerMeV; }
    void SetScintillatorDecayTime(float ns) { scintDecayTime_ = ns; }
    void SetScintillatorEmissionSpectrum(const std::vector<float>& wavelengths,
                                          const std::vector<float>& intensities);

    // Set Cherenkov properties
    void SetCherenkovEnabled(bool enable) { cherenkovEnabled_ = enable; }

    // Set material refractive indices
    void SetRefractiveIndex(int materialId, float n) { refractiveIndices_[materialId] = n; }

    // ========================================================================
    // Processing
    // ========================================================================

    // Add scintillation points for processing
    void AddScintillationPoint(const ScintillationPoint& point);

    // Add Cherenkov emission point
    void AddCherenkovPoint(float x, float y, float z,
                           float dx, float dy, float dz,
                           float beta, float trackLength,
                           int materialId, int eventId);

    // Process all pending photons (call at end of event)
    void ProcessEvent();

    // Get PMT hits from last processed event
    const std::vector<PMTHit>& GetPMTHits() const { return pmtHits_; }

    // Get total photon counts
    int GetTotalScintillationPhotons() const { return totalScintPhotons_; }
    int GetTotalCherenkovPhotons() const { return totalCherenkovPhotons_; }
    int GetTotalDetectedPhotons() const { return totalDetectedPhotons_; }

    // Clear event data
    void ClearEvent();

    // ========================================================================
    // Statistics
    // ========================================================================

    void PrintStatistics() const;

private:
    OpticalPhotonGPU() = default;
    ~OpticalPhotonGPU();

    // GPU implementation
    void ProcessEventGPU();

    // CPU fallback implementation
    void ProcessEventCPU();

    // Generate scintillation photons
    void GenerateScintillationPhotons(const ScintillationPoint& point,
                                       std::vector<OpticalPhoton>& photons);

    // Propagate photons to PMTs
    void PropagatePhotons(std::vector<OpticalPhoton>& photons);

    // State
    bool isInitialized_ = false;
    bool isGPUEnabled_ = false;

    // Configuration
    float scintYield_ = 10000.0f;        // photons/MeV
    float scintDecayTime_ = 2.0f;        // ns
    bool cherenkovEnabled_ = true;
    std::vector<float> emissionWavelengths_;
    std::vector<float> emissionIntensities_;
    std::vector<float> refractiveIndices_;

    // Event data
    std::vector<ScintillationPoint> scintPoints_;
    std::vector<OpticalPhoton> cherenkovPhotons_;
    std::vector<PMTHit> pmtHits_;

    // Statistics
    int totalScintPhotons_ = 0;
    int totalCherenkovPhotons_ = 0;
    int totalDetectedPhotons_ = 0;
    int eventsProcessedGPU_ = 0;
    int eventsProcessedCPU_ = 0;

    // GPU buffers (if CUDA available)
#if HAS_CUDA
    void* d_scintPoints_ = nullptr;
    void* d_photons_ = nullptr;
    void* d_pmtHits_ = nullptr;
    size_t maxPhotons_ = 10000000;  // 10M photons max
#endif

    // Thread safety
    std::mutex mutex_;
};

}  // namespace gpu
}  // namespace nnbar

#endif // OPTICAL_PHOTON_GPU_HH
