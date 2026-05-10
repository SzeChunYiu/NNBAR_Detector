// ============================================================================
// OpticksInterface.hh
// Interface for Opticks GPU-accelerated optical photon propagation
// ============================================================================
// Opticks offloads optical photon simulation (scintillation, Cerenkov) to GPU
// using NVIDIA OptiX ray tracing, providing 50-200x speedup.
//
// Requirements:
//   - NVIDIA GPU with RT cores (recommended) or CUDA cores
//   - CUDA Toolkit 11.x+
//   - OptiX 7.x
//   - Opticks library (https://bitbucket.org/simoncblyth/opticks)
//
// Usage:
//   In your PhysicsList or main():
//     #ifdef WITH_OPTICKS
//     OpticksInterface::Instance().Initialize();
//     #endif
// ============================================================================

#ifndef OPTICKS_INTERFACE_HH
#define OPTICKS_INTERFACE_HH

#ifdef WITH_OPTICKS

#include <memory>
#include <vector>
#include <string>

class G4VPhysicalVolume;
class G4Event;
class G4Track;

namespace nnbar {

// ============================================================================
// Photon hit structure returned from GPU
// ============================================================================
struct OpticalHit {
    int detector_id;
    double x, y, z;
    double time;
    double wavelength;
    int flags;
};

// ============================================================================
// OpticksInterface - Singleton for GPU optical photon propagation
// ============================================================================
class OpticksInterface {
public:
    // Singleton access
    static OpticksInterface& Instance();

    // Delete copy/move
    OpticksInterface(const OpticksInterface&) = delete;
    OpticksInterface& operator=(const OpticksInterface&) = delete;

    // Initialize with Geant4 geometry (call after detector construction)
    void Initialize(G4VPhysicalVolume* world);

    // Finalize and cleanup GPU resources
    void Finalize();

    // Check if Opticks is active
    bool IsActive() const { return m_active; }

    // Check if Opticks is enabled (static, checks environment)
    static bool IsEnabled();

    // Enable/disable GPU offloading at runtime
    void SetEnabled(bool enabled);
    bool IsRuntimeEnabled() const { return m_enabled; }

    // Add a scintillation or Cerenkov "genstep" from Geant4
    // These are collected during stepping and processed at end of event
    void AddScintillationGenstep(const G4Track* track,
                                  int numPhotons,
                                  double meanEnergy,
                                  double time);

    void AddCerenkovGenstep(const G4Track* track,
                            int numPhotons,
                            double minWavelength,
                            double maxWavelength);

    // Process all collected gensteps on GPU
    // Call this at EndOfEventAction
    void PropagateOpticalPhotons();

    // Get hits from GPU after propagation
    const std::vector<OpticalHit>& GetHits() const { return m_hits; }

    // Clear hits for next event
    void ClearHits() { m_hits.clear(); }

    // Statistics
    int GetPhotonsGenerated() const { return m_photonsGenerated; }
    int GetPhotonsDetected() const { return m_photonsDetected; }
    double GetGPUTime() const { return m_gpuTime; }

private:
    OpticksInterface();
    ~OpticksInterface();

    // State
    bool m_active = false;
    bool m_enabled = true;
    bool m_initialized = false;

    // Gensteps collected during event
    std::vector<float> m_gensteps;
    int m_numGensteps = 0;

    // Hits returned from GPU
    std::vector<OpticalHit> m_hits;

    // Statistics
    int m_photonsGenerated = 0;
    int m_photonsDetected = 0;
    double m_gpuTime = 0.0;
};

} // namespace nnbar

#endif // WITH_OPTICKS

#endif // OPTICKS_INTERFACE_HH
