// ============================================================================
//
//                      CELERITAS CALORIMETER HEADER FILE
//                      ==================================
//
//  What This File Does:
//  --------------------
//  This file defines a helper class that records energy deposits from
//  particles tracked on the GPU by Celeritas.
//
//  Why We Need This:
//  -----------------
//  When Celeritas runs electromagnetic (EM) particles like electrons,
//  positrons, and gamma rays on the GPU, the normal Geant4 hit recording
//  doesn't work. This class uses Celeritas's "GeantSimpleCalo" feature to
//  capture how much energy these GPU-tracked particles deposit in our
//  Lead Glass and Scintillator detectors.
//
//  How It Works (Simple Explanation):
//  ----------------------------------
//  1. We tell Celeritas which detector volumes we care about
//  2. Celeritas tracks particles on the GPU and remembers energy deposits
//  3. After each event, we calculate PER-EVENT energy by subtracting
//     the previous cumulative total from the current cumulative total
//  4. We return the per-event energy for output and physics analysis
//
// ============================================================================

#ifndef CELERITAS_CALORIMETER_HH
#define CELERITAS_CALORIMETER_HH

#include "config.h"

#if WITH_CELERITAS

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <mutex>
#include <atomic>

class G4LogicalVolume;
class G4VSensitiveDetector;

namespace celeritas {
    class GeantSimpleCalo;
    class SharedParams;
}

namespace nnbar {

// ============================================================================
// GPU Energy Data Structure - Per-Event Energy Information
// ============================================================================
struct GPUEnergyData {
    double leadGlassEnergy = 0.0;      // Energy deposited in Lead Glass (MeV)
    double scintillatorEnergy = 0.0;   // Energy deposited in Scintillator (MeV)
    double totalGPUEnergy = 0.0;       // Total GPU energy (MeV)
    bool isValid = false;              // Whether data is valid for this event
    bool isApproximate = false;        // True if multi-threaded (approximate values)
};

// ============================================================================
// CeleritasCalorimeter - Records EM Energy Deposits From GPU
// ============================================================================
class CeleritasCalorimeter {
public:
    // ------------------------------------------------------------------------
    // Singleton Access
    // ------------------------------------------------------------------------
    static CeleritasCalorimeter& Instance();

    // ------------------------------------------------------------------------
    // Initialization (Call From DetectorConstruction)
    // ------------------------------------------------------------------------
    void Initialize(G4LogicalVolume* leadGlassLV, G4LogicalVolume* scintillatorLV);

    // ------------------------------------------------------------------------
    // Create Thread-Local SDs (Call From RunAction::BeginOfRunAction)
    // ------------------------------------------------------------------------
    void CreateThreadLocalSD();

    // ------------------------------------------------------------------------
    // Get PER-EVENT Energy Deposits (Call From EventAction::EndOfEventAction)
    // ------------------------------------------------------------------------
    // Returns the energy deposited in THIS EVENT ONLY, not cumulative.
    // This is done by tracking previous cumulative values and subtracting.
    GPUEnergyData GetPerEventEnergy();

    // ------------------------------------------------------------------------
    // Get Cumulative Energy (For Debugging)
    // ------------------------------------------------------------------------
    std::map<std::string, double> GetCumulativeEnergy() const;

    // ------------------------------------------------------------------------
    // Begin New Event (Call At Start Of Each Event)
    // ------------------------------------------------------------------------
    // Stores the current cumulative energy so we can calculate per-event later.
    void BeginEvent();

    // ------------------------------------------------------------------------
    // Status Checks
    // ------------------------------------------------------------------------
    bool IsInitialized() const { return isInitialized_; }
    bool IsCeleritasActive() const;

private:
    CeleritasCalorimeter();
    ~CeleritasCalorimeter();

    CeleritasCalorimeter(const CeleritasCalorimeter&) = delete;
    CeleritasCalorimeter& operator=(const CeleritasCalorimeter&) = delete;

    // ------------------------------------------------------------------------
    // State Variables
    // ------------------------------------------------------------------------
    bool isInitialized_ = false;

    // Celeritas Calorimeter Instances
    std::shared_ptr<celeritas::GeantSimpleCalo> leadGlassCalorimeter_;
    std::shared_ptr<celeritas::GeantSimpleCalo> scintillatorCalorimeter_;

    // Geant4 Logical Volumes
    G4LogicalVolume* leadGlassLogicalVolume_ = nullptr;
    G4LogicalVolume* scintillatorLogicalVolume_ = nullptr;

    // Thread-Local Sensitive Detectors
    std::vector<std::unique_ptr<G4VSensitiveDetector>> threadLocalDetectors_;

    // ------------------------------------------------------------------------
    // Per-Event Energy Tracking (Serialized for Accuracy)
    // ------------------------------------------------------------------------
    // Because GeantSimpleCalo accumulates energy globally, we serialize
    // the BeginEvent/GetPerEventEnergy sequence to ensure each event
    // gets its correct share of GPU energy.
    double lastLeadGlass_ = 0.0;
    double lastScintillator_ = 0.0;
    mutable std::mutex energyMutex_;

    // Track number of events for statistics
    std::atomic<int> eventCount_{0};
};

} // namespace nnbar

#endif // WITH_CELERITAS
#endif // CELERITAS_CALORIMETER_HH
