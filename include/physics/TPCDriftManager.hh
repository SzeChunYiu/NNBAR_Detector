// ============================================================================
// TPCDriftManager.hh
// Manager class for GPU-accelerated TPC electron drift simulation
// Integrates GarfieldGPU with Geant4 sensitive detectors
// ============================================================================

#ifndef TPC_DRIFT_MANAGER_HH
#define TPC_DRIFT_MANAGER_HH

#include "config.h"
#include <vector>
#include <memory>
#include <mutex>

#if WITH_GARFIELD_GPU
#include "physics/GarfieldGPU.hh"
#endif

namespace nnbar {

/**
 * @struct TPCIonizationData
 * @brief Data for a single ionization deposit in the TPC
 */
struct TPCIonizationData {
    float x, y, z;      // Position in cm
    float t;            // Time in ns
    int nElectrons;     // Number of ionization electrons
    int tpcModule;      // TPC module index
    int tpcLayer;       // TPC layer index
    int trackID;        // Geant4 track ID
};

/**
 * @struct TPCDriftResult
 * @brief Result of drift simulation for collected electrons
 */
struct TPCDriftResult {
    int tpcModule;          // TPC module
    int tpcLayer;           // TPC layer
    float finalX, finalY;   // Final position on readout plane
    float driftTime;        // Drift time (ns)
    int nElectrons;         // Number of electrons (after avalanche)
};

/**
 * @class TPCDriftManager
 * @brief Singleton manager for TPC electron drift simulation
 *
 * Collects ionization data during Geant4 stepping and processes
 * all electrons at end of event using GPU-accelerated drift.
 *
 * Usage:
 *   // In TPCSD::ProcessHits:
 *   TPCDriftManager::Instance()->AddIonization(x, y, z, t, nElectrons, mod, layer, trackID);
 *
 *   // In EventAction::EndOfEventAction:
 *   TPCDriftManager::Instance()->ProcessEvent();
 *   auto results = TPCDriftManager::Instance()->GetResults();
 */
class TPCDriftManager {
public:
    /**
     * Get singleton instance
     */
    static TPCDriftManager* Instance();

    /**
     * Initialize the drift manager
     * @return true if GPU/CPU drift engine initialized successfully
     */
    bool Initialize();

    /**
     * Check if GPU acceleration is available
     */
    bool IsGPUAvailable() const;

    /**
     * Get device name (GPU model or CPU info)
     */
    const char* GetDeviceName() const;

    // ========== Configuration ==========

    /**
     * Set TPC geometry parameters
     */
    void SetTPCGeometry(float innerRadius, float outerRadius,
                        float halfLength, float driftLength);

    /**
     * Set gas properties
     * @param arFrac Argon fraction (0-1)
     * @param co2Frac CO2 fraction (0-1)
     * @param pressure Pressure in Torr
     * @param temperature Temperature in Kelvin
     */
    void SetGasProperties(float arFrac, float co2Frac,
                          float pressure = 760.0f, float temperature = 293.15f);

    /**
     * Enable/disable avalanche multiplication
     */
    void SetAvalancheEnabled(bool enable);

    // ========== Event Processing ==========

    /**
     * Add ionization data from TPC hit
     * Called from TPCSD::ProcessHits
     */
    void AddIonization(float x, float y, float z, float t,
                       int nElectrons, int tpcModule, int tpcLayer, int trackID);

    /**
     * Process all ionization data for current event
     * Drifts all electrons to the readout plane
     * Called at end of event
     */
    void ProcessEvent();

    /**
     * Clear all data for new event
     */
    void ClearEvent();

    /**
     * Get drift results after ProcessEvent()
     */
    const std::vector<TPCDriftResult>& GetResults() const { return m_results; }

    /**
     * Get total number of collected electrons
     */
    int GetTotalCollected() const { return m_totalCollected; }

    /**
     * Get total drift charge (accounting for avalanche)
     */
    float GetTotalCharge() const { return m_totalCharge; }

    // ========== Statistics ==========

    /**
     * Get drift kernel execution time (ms)
     */
    float GetKernelTimeMs() const { return m_kernelTimeMs; }

    /**
     * Get number of ionization clusters in current event
     */
    int GetNClusters() const { return static_cast<int>(m_ionizationData.size()); }

    /**
     * Get total input electrons for current event
     */
    int GetTotalInputElectrons() const { return m_totalInputElectrons; }

    /**
     * Check if drift simulation is enabled
     */
    bool IsEnabled() const { return m_enabled; }

    /**
     * Enable/disable drift simulation
     */
    void SetEnabled(bool enable) { m_enabled = enable; }

private:
    TPCDriftManager();
    ~TPCDriftManager();

    // Singleton instance (thread-safe via std::call_once)
    static TPCDriftManager* s_instance;
    static std::once_flag s_onceFlag;

    // Mutex for thread-safe access to shared data
    mutable std::mutex m_mutex;

    // Configuration
    bool m_initialized = false;
    bool m_enabled = true;

#if WITH_GARFIELD_GPU
    std::unique_ptr<GarfieldGPU> m_driftEngine;
#endif

    // Event data (protected by m_mutex)
    std::vector<TPCIonizationData> m_ionizationData;
    std::vector<TPCDriftResult> m_results;

    // Statistics (protected by m_mutex)
    int m_totalInputElectrons = 0;
    int m_totalCollected = 0;
    float m_totalCharge = 0.0f;
    float m_kernelTimeMs = 0.0f;
};

} // namespace nnbar

#endif // TPC_DRIFT_MANAGER_HH
