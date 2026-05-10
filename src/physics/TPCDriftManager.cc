// ============================================================================
// TPCDriftManager.cc
// Manager class for GPU-accelerated TPC electron drift simulation
// Thread-safe implementation for Geant4 MT mode
// ============================================================================

#include "physics/TPCDriftManager.hh"
#include <iostream>

namespace nnbar {

// Singleton instance with thread-safe initialization
TPCDriftManager* TPCDriftManager::s_instance = nullptr;
std::once_flag TPCDriftManager::s_onceFlag;

TPCDriftManager* TPCDriftManager::Instance() {
    std::call_once(s_onceFlag, []() {
        s_instance = new TPCDriftManager();
    });
    return s_instance;
}

TPCDriftManager::TPCDriftManager() {
#if WITH_GARFIELD_GPU
    m_driftEngine = std::make_unique<GarfieldGPU>();
#endif
}

TPCDriftManager::~TPCDriftManager() {
    s_instance = nullptr;
}

bool TPCDriftManager::Initialize() {
    if (m_initialized) return true;

#if WITH_GARFIELD_GPU
    if (m_driftEngine) {
        m_initialized = m_driftEngine->Initialize();
        if (m_initialized) {
            std::cout << "TPCDriftManager: Initialized with "
                      << GetDeviceName() << std::endl;
        }
        return m_initialized;
    }
    return false;
#else
    std::cout << "TPCDriftManager: GarfieldGPU not enabled in build" << std::endl;
    return false;
#endif
}

bool TPCDriftManager::IsGPUAvailable() const {
#if WITH_GARFIELD_GPU
    return m_driftEngine && m_driftEngine->IsGPUAvailable();
#else
    return false;
#endif
}

const char* TPCDriftManager::GetDeviceName() const {
#if WITH_GARFIELD_GPU
    if (m_driftEngine) {
        return m_driftEngine->GetDeviceName();
    }
#endif
    return "Not available";
}

void TPCDriftManager::SetTPCGeometry(float innerRadius, float outerRadius,
                                      float halfLength, float driftLength) {
#if WITH_GARFIELD_GPU
    if (m_driftEngine) {
        TPCGeometry geom;
        geom.innerRadius = innerRadius;
        geom.outerRadius = outerRadius;
        geom.halfLength = halfLength;
        geom.driftLength = driftLength;
        geom.cathodeZ = 0.0f;
        geom.anodeZ = driftLength;
        m_driftEngine->SetGeometry(geom);
    }
#else
    (void)innerRadius; (void)outerRadius;
    (void)halfLength; (void)driftLength;
#endif
}

void TPCDriftManager::SetGasProperties(float arFrac, float co2Frac,
                                        float pressure, float temperature) {
#if WITH_GARFIELD_GPU
    if (m_driftEngine) {
        m_driftEngine->SetGasProperties(arFrac, co2Frac, pressure, temperature);
    }
#else
    (void)arFrac; (void)co2Frac;
    (void)pressure; (void)temperature;
#endif
}

void TPCDriftManager::SetAvalancheEnabled(bool enable) {
#if WITH_GARFIELD_GPU
    if (m_driftEngine) {
        m_driftEngine->SetAvalancheEnabled(enable);
    }
#else
    (void)enable;
#endif
}

void TPCDriftManager::AddIonization(float x, float y, float z, float t,
                                     int nElectrons, int tpcModule,
                                     int tpcLayer, int trackID) {
    if (!m_enabled) return;

    TPCIonizationData data;
    data.x = x;
    data.y = y;
    data.z = z;
    data.t = t;
    data.nElectrons = nElectrons;
    data.tpcModule = tpcModule;
    data.tpcLayer = tpcLayer;
    data.trackID = trackID;

    // Thread-safe access to shared data
    std::lock_guard<std::mutex> lock(m_mutex);
    m_ionizationData.push_back(data);
    m_totalInputElectrons += nElectrons;
}

void TPCDriftManager::ProcessEvent() {
    // Thread-safe access
    std::lock_guard<std::mutex> lock(m_mutex);

    if (!m_enabled || m_ionizationData.empty()) {
        return;
    }

#if WITH_GARFIELD_GPU
    if (!m_initialized || !m_driftEngine) {
        std::cerr << "TPCDriftManager: Not initialized!" << std::endl;
        return;
    }

    // Clear previous results
    m_driftEngine->ClearElectrons();
    m_results.clear();

    // Add all ionization clusters to drift engine
    for (const auto& ion : m_ionizationData) {
        // Convert position from mm (Geant4) to cm (GarfieldGPU)
        float x_cm = ion.x / 10.0f;
        float y_cm = ion.y / 10.0f;
        float z_cm = ion.z / 10.0f;

        // Add electrons from this cluster
        for (int i = 0; i < ion.nElectrons; i++) {
            m_driftEngine->AddElectron(x_cm, y_cm, z_cm, ion.t);
        }
    }

    // Run drift simulation
    int nCollected = m_driftEngine->DriftElectrons();

    // Get results
    m_totalCollected = nCollected;
    m_totalCharge = m_driftEngine->GetTotalCharge();
    m_kernelTimeMs = m_driftEngine->GetKernelTimeMs();

    // Convert drift results to TPCDriftResult format
    const auto& driftResults = m_driftEngine->GetResults();
    size_t ionIdx = 0;
    size_t electronIdx = 0;

    for (size_t i = 0; i < driftResults.size() && ionIdx < m_ionizationData.size(); i++) {
        const auto& dr = driftResults[i];

        // Map electron back to its ionization cluster
        while (ionIdx < m_ionizationData.size() &&
               electronIdx >= static_cast<size_t>(m_ionizationData[ionIdx].nElectrons)) {
            electronIdx = 0;
            ionIdx++;
        }

        if (dr.status == 1) {  // Collected at anode
            TPCDriftResult result;
            result.tpcModule = m_ionizationData[ionIdx].tpcModule;
            result.tpcLayer = m_ionizationData[ionIdx].tpcLayer;
            result.finalX = dr.finalX * 10.0f;  // cm -> mm
            result.finalY = dr.finalY * 10.0f;
            result.driftTime = dr.driftTime;
            result.nElectrons = dr.nAvalanche;
            m_results.push_back(result);
        }

        electronIdx++;
    }

    std::cout << "TPCDriftManager: Processed " << m_totalInputElectrons
              << " electrons -> " << m_totalCollected << " collected ("
              << m_kernelTimeMs << " ms)" << std::endl;

#else
    std::cout << "TPCDriftManager: GarfieldGPU not available, skipping drift" << std::endl;
#endif
}

void TPCDriftManager::ClearEvent() {
    // Thread-safe access
    std::lock_guard<std::mutex> lock(m_mutex);

    m_ionizationData.clear();
    m_results.clear();
    m_totalInputElectrons = 0;
    m_totalCollected = 0;
    m_totalCharge = 0.0f;
    m_kernelTimeMs = 0.0f;
}

} // namespace nnbar
