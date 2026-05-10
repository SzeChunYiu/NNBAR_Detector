// ============================================================================
// GarfieldGPU.hh
// GPU-accelerated TPC electron drift simulation
// Custom implementation replacing Garfield++ drift with CUDA kernels
// ============================================================================
// Enable with: cmake -DWITH_GARFIELD_GPU=ON ..
// Requires: CUDA Toolkit 11.0+
// ============================================================================

#ifndef GARFIELD_GPU_HH
#define GARFIELD_GPU_HH

#include "config.h"
#include <vector>
#include <memory>
#include <cstdint>

namespace nnbar {

// ============================================================================
// Data Structures for GPU Computation
// ============================================================================

/// Electron state for drift simulation
struct ElectronState {
    float x, y, z;      // Position in cm
    float t;            // Time in ns
    float vx, vy, vz;   // Velocity components
    int32_t status;     // 0=active, 1=collected, 2=absorbed, 3=out-of-bounds
};

/// Ionization cluster from primary track
struct IonizationCluster {
    float x, y, z;      // Position in cm
    float t;            // Time in ns
    int32_t nElectrons; // Number of electrons in cluster
    float energy;       // Energy deposited (eV)
};

/// Gas transport properties (precomputed for efficiency)
struct GasTransportData {
    float driftVelocity;      // cm/ns at reference field
    float diffusionL;         // Longitudinal diffusion (cm/sqrt(cm))
    float diffusionT;         // Transverse diffusion (cm/sqrt(cm))
    float townsendCoeff;      // Townsend coefficient (1/cm)
    float attachmentCoeff;    // Attachment coefficient (1/cm)
    float electricField;      // Field strength (V/cm)
};

/// TPC geometry parameters
struct TPCGeometry {
    float innerRadius;        // Inner radius (beampipe) in cm
    float outerRadius;        // Outer radius in cm
    float halfLength;         // Half-length in z (cm)
    float driftLength;        // Drift distance (cm)
    float cathodeZ;           // Cathode z position
    float anodeZ;             // Anode z position
};

/// Result of drift simulation
struct DriftResult {
    float finalX, finalY, finalZ;  // Final position
    float driftTime;               // Total drift time (ns)
    float pathLength;              // Path length (cm)
    int32_t status;                // Final status
    int32_t nAvalanche;            // Electrons from avalanche (if enabled)
};

// ============================================================================
// GPU Drift Engine Class
// ============================================================================

/**
 * @class GarfieldGPU
 * @brief GPU-accelerated electron drift simulation for TPC
 *
 * Replaces Garfield++ AvalancheMC with CUDA-parallelized drift simulation.
 * Each electron is drifted independently on a GPU thread.
 *
 * Features:
 * - Parallel drift of thousands of electrons simultaneously
 * - Langevin diffusion model with thermal noise
 * - Configurable gas properties from lookup tables
 * - Boundary checking for TPC geometry
 * - Optional avalanche multiplication
 *
 * Performance: 100-1000x speedup vs CPU Garfield++ for >1000 electrons
 */
class GarfieldGPU {
public:
    GarfieldGPU();
    ~GarfieldGPU();

    // ========== Initialization ==========

    /**
     * Initialize GPU resources and check CUDA availability
     * @return true if GPU is available and initialized
     */
    bool Initialize();

    /**
     * Check if GPU acceleration is available
     */
    bool IsGPUAvailable() const { return m_gpuAvailable; }

    /**
     * Get GPU device name
     */
    const char* GetDeviceName() const { return m_deviceName.c_str(); }

    // ========== Configuration ==========

    /**
     * Set TPC geometry
     */
    void SetGeometry(const TPCGeometry& geom);

    /**
     * Set gas transport properties
     * @param arFrac Argon fraction (0-1)
     * @param co2Frac CO2 fraction (0-1)
     * @param pressure Pressure in Torr
     * @param temperature Temperature in Kelvin
     */
    void SetGasProperties(float arFrac, float co2Frac,
                          float pressure = 760.0f, float temperature = 293.15f);

    /**
     * Set electric field strength
     * @param field Field in V/cm
     */
    void SetElectricField(float field);

    /**
     * Enable/disable avalanche multiplication
     */
    void SetAvalancheEnabled(bool enable) { m_avalancheEnabled = enable; }

    /**
     * Set maximum simulation time
     * @param maxTime Maximum drift time in ns
     */
    void SetMaxDriftTime(float maxTime) { m_maxDriftTime = maxTime; }

    // ========== Simulation Interface ==========

    /**
     * Add ionization clusters from primary track
     * @param clusters Vector of ionization clusters
     */
    void AddClusters(const std::vector<IonizationCluster>& clusters);

    /**
     * Add single electron at position
     */
    void AddElectron(float x, float y, float z, float t = 0.0f);

    /**
     * Clear all electrons
     */
    void ClearElectrons();

    /**
     * Perform drift simulation on GPU
     * @return Number of electrons that reached the anode
     */
    int DriftElectrons();

    /**
     * Get drift results after simulation
     */
    const std::vector<DriftResult>& GetResults() const { return m_results; }

    /**
     * Get total number of collected electrons
     */
    int GetCollectedElectrons() const { return m_nCollected; }

    /**
     * Get total drift charge (accounting for avalanche)
     */
    float GetTotalCharge() const { return m_totalCharge; }

    // ========== Statistics ==========

    /**
     * Get GPU kernel execution time
     */
    float GetKernelTimeMs() const { return m_kernelTimeMs; }

    /**
     * Get number of electrons simulated
     */
    int GetNElectrons() const { return static_cast<int>(m_electrons.size()); }

private:
    // ========== GPU Resources ==========
    bool m_gpuAvailable = false;
    bool m_initialized = false;
    std::string m_deviceName;

    // Device pointers (opaque - implemented in .cu file)
    void* m_d_electrons = nullptr;
    void* m_d_results = nullptr;
    void* m_d_gasData = nullptr;
    void* m_d_geometry = nullptr;
    void* m_d_randStates = nullptr;

    // ========== Configuration ==========
    TPCGeometry m_geometry;
    GasTransportData m_gasData;
    bool m_avalancheEnabled = false;
    float m_maxDriftTime = 10000.0f;  // ns

    // ========== Electron Data ==========
    std::vector<ElectronState> m_electrons;
    std::vector<DriftResult> m_results;

    // ========== Statistics ==========
    int m_nCollected = 0;
    float m_totalCharge = 0.0f;
    float m_kernelTimeMs = 0.0f;

    // ========== Internal Methods ==========
    void CalculateGasTransport();
    void AllocateGPUMemory(int nElectrons);
    void FreeGPUMemory();
    void LaunchDriftKernel(int nElectrons);
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * Calculate drift velocity for Ar/CO2 mixture
 * Uses empirical fit to Magboltz data
 */
float CalculateDriftVelocity(float electricField, float pressure,
                              float temperature, float arFrac);

/**
 * Calculate diffusion coefficients
 */
void CalculateDiffusion(float electricField, float pressure,
                        float& diffL, float& diffT, float arFrac);

/**
 * Generate primary ionization clusters using simple model
 * (Alternative to Heed when Garfield++ not available)
 */
std::vector<IonizationCluster> GenerateIonization(
    float x0, float y0, float z0,
    float dx, float dy, float dz,
    float energy,      // Particle energy in MeV
    float pathLength,  // Path length in cm
    float dEdx         // Energy loss in MeV/cm
);

} // namespace nnbar

#endif // GARFIELD_GPU_HH
