// ============================================================================
//
//                 CELERITAS CALORIMETER IMPLEMENTATION
//                 =====================================
//
//  Purpose: Record GPU EM particle energy deposits for physics analysis
//
//  Key Feature: Per-Event Energy Tracking
//  --------------------------------------
//  GeantSimpleCalo accumulates energy over ALL events. To get per-event
//  energy, we store the cumulative total at the START of each event,
//  then subtract it from the current total at the END of the event.
//
// ============================================================================

#include "config.h"

#if WITH_CELERITAS

#include "physics/CeleritasCalorimeter.hh"
#include "physics/CeleritasInterface.hh"

// Celeritas Headers
#include <accel/GeantSimpleCalo.hh>
#include <accel/SharedParams.hh>
#include <accel/TrackingManagerIntegration.hh>
#include <accel/detail/IntegrationSingleton.hh>

// Geant4 Headers
#include <G4LogicalVolume.hh>
#include <G4VSensitiveDetector.hh>
#include <G4SDManager.hh>
#include <G4ios.hh>
#include <G4AutoLock.hh>
#include <G4Threading.hh>
#include <G4RunManager.hh>

namespace {
    G4Mutex calorimetryMutex = G4MUTEX_INITIALIZER;
}

namespace nnbar {

// ============================================================================
// Singleton Instance
// ============================================================================

CeleritasCalorimeter& CeleritasCalorimeter::Instance() {
    static CeleritasCalorimeter instance;
    return instance;
}

// ============================================================================
// Constructor / Destructor
// ============================================================================

CeleritasCalorimeter::CeleritasCalorimeter() {
    G4cout << "[CeleritasCalorimeter] GPU Energy Recording System Created" << G4endl;
}

CeleritasCalorimeter::~CeleritasCalorimeter() {
    G4cout << "[CeleritasCalorimeter] Shutting down" << G4endl;
}

// ============================================================================
// Check If Celeritas Is Active
// ============================================================================

bool CeleritasCalorimeter::IsCeleritasActive() const {
    if (!CeleritasInterface::IsEnabled()) {
        return false;
    }

    try {
        auto& singleton = celeritas::detail::IntegrationSingleton::instance();
        auto mode = singleton.shared_params().mode();
        return mode == celeritas::OffloadMode::enabled;
    } catch (...) {
        return false;
    }
}

// ============================================================================
// Initialize With Detector Volumes
// ============================================================================

void CeleritasCalorimeter::Initialize(G4LogicalVolume* leadGlassVolume,
                                       G4LogicalVolume* scintillatorVolume) {
    G4AutoLock lock(&calorimetryMutex);

    if (isInitialized_) {
        return;
    }

    leadGlassLogicalVolume_ = leadGlassVolume;
    scintillatorLogicalVolume_ = scintillatorVolume;

    G4cout << "[CeleritasCalorimeter] Initialized with:" << G4endl;
    if (leadGlassLogicalVolume_) {
        G4cout << "  - Lead Glass: " << leadGlassLogicalVolume_->GetName() << G4endl;
    }
    if (scintillatorLogicalVolume_) {
        G4cout << "  - Scintillator: " << scintillatorLogicalVolume_->GetName() << G4endl;
    }

    // Reset energy tracking
    {
        std::lock_guard<std::mutex> lock(energyMutex_);
        lastLeadGlass_ = 0.0;
        lastScintillator_ = 0.0;
    }
    eventCount_ = 0;

    isInitialized_ = true;
}

// ============================================================================
// Create Thread-Local Sensitive Detectors
// ============================================================================

void CeleritasCalorimeter::CreateThreadLocalSD() {
    if (!isInitialized_ || !CeleritasInterface::IsEnabled()) {
        return;
    }

    G4AutoLock lock(&calorimetryMutex);

    // If calorimeters already exist, just create thread-local SDs
    if (leadGlassCalorimeter_ || scintillatorCalorimeter_) {
        if (leadGlassCalorimeter_) {
            auto sd = leadGlassCalorimeter_->MakeSensitiveDetector();
            if (sd) threadLocalDetectors_.push_back(std::move(sd));
        }
        if (scintillatorCalorimeter_) {
            auto sd = scintillatorCalorimeter_->MakeSensitiveDetector();
            if (sd) threadLocalDetectors_.push_back(std::move(sd));
        }
        return;
    }

    // First-time setup: Create the calorimeters
    try {
        auto& singleton = celeritas::detail::IntegrationSingleton::instance();
        auto& sharedParams = singleton.shared_params();

        if (sharedParams.mode() != celeritas::OffloadMode::enabled) {
            return;
        }

        // Create non-owning shared_ptr (aliasing constructor)
        auto paramsPtr = std::shared_ptr<celeritas::SharedParams const>(
            std::shared_ptr<void>(),
            &sharedParams
        );

        // Create Lead Glass Calorimeter
        if (leadGlassLogicalVolume_) {
            std::vector<G4LogicalVolume*> volumes = {leadGlassLogicalVolume_};
            leadGlassCalorimeter_ = std::make_shared<celeritas::GeantSimpleCalo>(
                "GPU_LeadGlass", paramsPtr, volumes);
            G4cout << "[CeleritasCalorimeter] Created Lead Glass GPU Calorimeter" << G4endl;

            auto sd = leadGlassCalorimeter_->MakeSensitiveDetector();
            if (sd) threadLocalDetectors_.push_back(std::move(sd));
        }

        // Create Scintillator Calorimeter
        if (scintillatorLogicalVolume_) {
            std::vector<G4LogicalVolume*> volumes = {scintillatorLogicalVolume_};
            scintillatorCalorimeter_ = std::make_shared<celeritas::GeantSimpleCalo>(
                "GPU_Scintillator", paramsPtr, volumes);
            G4cout << "[CeleritasCalorimeter] Created Scintillator GPU Calorimeter" << G4endl;

            auto sd = scintillatorCalorimeter_->MakeSensitiveDetector();
            if (sd) threadLocalDetectors_.push_back(std::move(sd));
        }

    } catch (const std::exception& e) {
        G4cerr << "[CeleritasCalorimeter] ERROR: " << e.what() << G4endl;
    }
}

// ============================================================================
// Begin Event - Increment Event Counter
// ============================================================================
// In multi-threaded mode, we don't snapshot here because events run in parallel.
// Instead, GetPerEventEnergy uses serialized access to track cumulative changes.

void CeleritasCalorimeter::BeginEvent() {
    eventCount_++;
}

// ============================================================================
// Get Per-Event Energy (Serialized for Accuracy)
// ============================================================================
// Returns the energy deposited since the last call to this method.
// Uses mutex to serialize access, ensuring each event gets its correct
// share of GPU energy even in multi-threaded mode.
//
// NOTE: In multi-threaded mode, this may capture energy from partially
// overlapping events. The values are accurate for the run total but may
// have event-to-event variation due to timing of mutex acquisition.

GPUEnergyData CeleritasCalorimeter::GetPerEventEnergy() {
    GPUEnergyData result;

    if (!IsCeleritasActive()) {
        result.isValid = false;
        return result;
    }

    // Serialize access to ensure accurate per-event energy
    std::lock_guard<std::mutex> lock(energyMutex_);

    // Get current cumulative totals
    auto cumulative = GetCumulativeEnergy();
    double currentLG = cumulative["LeadGlass"];
    double currentScint = cumulative["Scintillator"];

    // Calculate energy deposited since last call
    result.leadGlassEnergy = currentLG - lastLeadGlass_;
    result.scintillatorEnergy = currentScint - lastScintillator_;
    result.totalGPUEnergy = result.leadGlassEnergy + result.scintillatorEnergy;
    result.isValid = true;
    result.isApproximate = (G4Threading::GetNumberOfRunningWorkerThreads() > 1);

    // Update last seen values for next call
    lastLeadGlass_ = currentLG;
    lastScintillator_ = currentScint;

    return result;
}

// ============================================================================
// Get Cumulative Energy (Raw Values From GeantSimpleCalo)
// ============================================================================

std::map<std::string, double> CeleritasCalorimeter::GetCumulativeEnergy() const {
    std::map<std::string, double> result;
    result["LeadGlass"] = 0.0;
    result["Scintillator"] = 0.0;

    if (leadGlassCalorimeter_) {
        auto deposits = leadGlassCalorimeter_->CalcTotalEnergyDeposition();
        for (double e : deposits) {
            result["LeadGlass"] += e;
        }
    }

    if (scintillatorCalorimeter_) {
        auto deposits = scintillatorCalorimeter_->CalcTotalEnergyDeposition();
        for (double e : deposits) {
            result["Scintillator"] += e;
        }
    }

    return result;
}

} // namespace nnbar

#endif // WITH_CELERITAS
