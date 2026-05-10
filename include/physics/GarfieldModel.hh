// ============================================================================
// GarfieldModel.hh
// Garfield++ fast simulation model for realistic TPC electron drift
// ============================================================================
// Enable with: cmake -DWITH_GARFIELD=ON -DGarfield_DIR=/path/to/garfield ..
// ============================================================================

#ifndef GARFIELD_MODEL_HH
#define GARFIELD_MODEL_HH

#include "config.h"

#ifdef WITH_GARFIELD

#include "G4VFastSimulationModel.hh"
#include "G4Region.hh"
#include <memory>

// Forward declarations for Garfield++
namespace Garfield {
    class MediumMagboltz;
    class ComponentAnalyticField;
    class Sensor;
    class TrackHeed;
    class AvalancheMC;
}

namespace nnbar {

/**
 * @class GarfieldModel
 * @brief Fast simulation model using Garfield++ for TPC
 *
 * Replaces Geant4 tracking in TPC gas with Garfield++ simulation:
 * - Primary ionization using Heed (more accurate than Geant4 for gas)
 * - Electron drift with diffusion
 * - Optional avalanche amplification
 * - Signal induction on readout electrodes
 *
 * Based on Geant4GarfieldDegradInterface methodology.
 */
class GarfieldModel : public G4VFastSimulationModel {
public:
    /**
     * Constructor
     * @param modelName Name for this model
     * @param region G4Region where this model applies (TPC gas)
     */
    GarfieldModel(const G4String& modelName, G4Region* region);

    virtual ~GarfieldModel();

    /**
     * Check if this model is applicable to a particle
     * Applicable to charged particles in TPC gas region
     */
    virtual G4bool IsApplicable(const G4ParticleDefinition& particle) override;

    /**
     * Check if model should be invoked for this track
     * Invoked when track is in TPC gas with sufficient energy
     */
    virtual G4bool ModelTrigger(const G4FastTrack& fastTrack) override;

    /**
     * Perform the fast simulation
     * Uses Heed for ionization, Garfield++ for drift
     */
    virtual void DoIt(const G4FastTrack& fastTrack,
                      G4FastStep& fastStep) override;

    /**
     * Set the electric field in the TPC
     * @param field Field strength in V/cm
     */
    void SetElectricField(double field);

    /**
     * Set the TPC gas properties
     * @param arFrac Argon fraction (0-1)
     * @param co2Frac CO2 fraction (0-1)
     * @param pressure Pressure in Torr
     * @param temperature Temperature in Kelvin
     */
    void SetGasProperties(double arFrac, double co2Frac,
                          double pressure = 760., double temperature = 293.15);

    /**
     * Get number of electrons produced in last interaction
     */
    int GetNumberOfElectrons() const { return m_nElectrons; }

    /**
     * Get total number of ions produced
     */
    int GetNumberOfIons() const { return m_nIons; }

private:
    // Initialize Garfield++ components
    void InitializeGarfield();

    // Garfield++ components
    std::unique_ptr<Garfield::MediumMagboltz> m_gas;
    std::unique_ptr<Garfield::ComponentAnalyticField> m_field;
    std::unique_ptr<Garfield::Sensor> m_sensor;
    std::unique_ptr<Garfield::TrackHeed> m_trackHeed;
    std::unique_ptr<Garfield::AvalancheMC> m_avalanche;

    // Configuration
    double m_electricField = 250.0;  // V/cm (typical TPC drift field)
    double m_arFraction = 0.80;
    double m_co2Fraction = 0.20;
    double m_pressure = 760.0;       // Torr
    double m_temperature = 293.15;   // K

    // Statistics
    int m_nElectrons = 0;
    int m_nIons = 0;

    bool m_initialized = false;
};

} // namespace nnbar

#else // !WITH_GARFIELD

// Stub when Garfield++ not available
namespace nnbar {

class GarfieldModel {
public:
    GarfieldModel(const G4String&, G4Region*) {
        G4cout << "GarfieldModel: Garfield++ not compiled - using standard Geant4 TPC simulation" << G4endl;
    }
    void SetElectricField(double) {}
    void SetGasProperties(double, double, double = 760., double = 293.15) {}
    int GetNumberOfElectrons() const { return 0; }
    int GetNumberOfIons() const { return 0; }
};

} // namespace nnbar

#endif // WITH_GARFIELD

#endif // GARFIELD_MODEL_HH
