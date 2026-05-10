// ============================================================================
// GarfieldModel.cc
// Garfield++ fast simulation model implementation for TPC
// ============================================================================

#include "physics/GarfieldModel.hh"
#include "config.h"

#ifdef WITH_GARFIELD

#include "Garfield/MediumMagboltz.hh"
#include "Garfield/ComponentAnalyticField.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4Proton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4SystemOfUnits.hh"
#include "G4FastTrack.hh"
#include "G4FastStep.hh"

namespace nnbar {

GarfieldModel::GarfieldModel(const G4String& modelName, G4Region* region)
    : G4VFastSimulationModel(modelName, region)
{
    G4cout << "=============================================" << G4endl;
    G4cout << "GarfieldModel: Initializing Garfield++ TPC simulation" << G4endl;
    G4cout << "  Region: " << region->GetName() << G4endl;
    G4cout << "=============================================" << G4endl;

    InitializeGarfield();
}

GarfieldModel::~GarfieldModel() {
    G4cout << "GarfieldModel: Cleanup" << G4endl;
}

void GarfieldModel::InitializeGarfield() {
    // Create gas medium (Ar/CO2 mixture)
    m_gas = std::make_unique<Garfield::MediumMagboltz>();
    m_gas->SetComposition("ar", m_arFraction * 100., "co2", m_co2Fraction * 100.);
    m_gas->SetTemperature(m_temperature);
    m_gas->SetPressure(m_pressure);
    m_gas->Initialise(true);  // verbose output

    // Generate gas tables for electron transport
    m_gas->EnableDrift();
    m_gas->EnablePrimaryIonisation();

    // Create uniform electric field component
    m_field = std::make_unique<Garfield::ComponentAnalyticField>();
    m_field->SetMedium(m_gas.get());

    // Set up uniform field (drift direction along z)
    // For a realistic TPC, this would be more complex
    m_field->AddWire(0., 0., 0.01, 1000., "s");  // Simple field setup

    // Create sensor
    m_sensor = std::make_unique<Garfield::Sensor>();
    m_sensor->AddComponent(m_field.get());
    m_sensor->SetArea(-100., -100., -100., 100., 100., 100.);

    // Create Heed track model for primary ionization
    m_trackHeed = std::make_unique<Garfield::TrackHeed>();
    m_trackHeed->SetSensor(m_sensor.get());
    m_trackHeed->EnableDeltaElectronTransport();

    // Create avalanche model (optional, for amplification)
    m_avalanche = std::make_unique<Garfield::AvalancheMC>();
    m_avalanche->SetSensor(m_sensor.get());
    m_avalanche->EnableDiffusion();

    m_initialized = true;
    G4cout << "GarfieldModel: Initialized with Ar/CO2 "
           << (m_arFraction*100) << "/" << (m_co2Fraction*100)
           << " at " << m_pressure << " Torr" << G4endl;
}

G4bool GarfieldModel::IsApplicable(const G4ParticleDefinition& particle) {
    // Apply to charged particles that can ionize gas
    G4String particleName = particle.GetParticleName();
    return (particleName == "e-" || particleName == "e+" ||
            particleName == "mu-" || particleName == "mu+" ||
            particleName == "proton" ||
            particleName == "pi+" || particleName == "pi-");
}

G4bool GarfieldModel::ModelTrigger(const G4FastTrack& fastTrack) {
    // Only trigger if particle has significant energy
    G4double kinEnergy = fastTrack.GetPrimaryTrack()->GetKineticEnergy();

    // Minimum energy for Heed to work properly (~1 keV)
    // For relativistic particles, Heed is most accurate
    if (kinEnergy < 1.0 * keV) return false;

    return m_initialized;
}

void GarfieldModel::DoIt(const G4FastTrack& fastTrack, G4FastStep& fastStep) {
    // Get particle properties
    const G4Track* track = fastTrack.GetPrimaryTrack();
    G4double kinEnergy = track->GetKineticEnergy() / eV;  // Garfield uses eV
    G4double mass = track->GetParticleDefinition()->GetPDGMass() / eV;

    G4ThreeVector position = track->GetPosition();
    G4ThreeVector momentum = track->GetMomentumDirection();

    // Convert to cm (Garfield units)
    double x0 = position.x() / cm;
    double y0 = position.y() / cm;
    double z0 = position.z() / cm;
    double t0 = track->GetGlobalTime() / ns;

    double dx = momentum.x();
    double dy = momentum.y();
    double dz = momentum.z();

    // Set particle type for Heed
    G4String particleName = track->GetParticleDefinition()->GetParticleName();
    if (particleName == "e-" || particleName == "e+") {
        m_trackHeed->SetParticle("electron");
    } else if (particleName == "mu-" || particleName == "mu+") {
        m_trackHeed->SetParticle("muon");
    } else if (particleName == "proton") {
        m_trackHeed->SetParticle("proton");
    } else if (particleName == "pi+" || particleName == "pi-") {
        m_trackHeed->SetParticle("pion");
    }

    // Set kinetic energy
    m_trackHeed->SetKineticEnergy(kinEnergy);

    // Create new track in Garfield++
    m_trackHeed->NewTrack(x0, y0, z0, t0, dx, dy, dz);

    // Collect ionization clusters
    m_nElectrons = 0;
    m_nIons = 0;

    double xc, yc, zc, tc;
    int nc;
    double ec, ekin;

    // Loop over clusters along the track
    while (m_trackHeed->GetCluster(xc, yc, zc, tc, nc, ec, ekin)) {
        m_nElectrons += nc;
        m_nIons += nc;  // One ion per electron-ion pair

        // Optionally: drift electrons to readout
        // m_avalanche->DriftElectron(xc, yc, zc, tc);
    }

    // Kill the Geant4 track (we've handled it)
    fastStep.KillPrimaryTrack();
    fastStep.SetPrimaryTrackPathLength(0.);

    // Set energy deposition
    G4double energyDeposit = kinEnergy * eV;  // Approximate
    fastStep.SetTotalEnergyDeposit(energyDeposit);

    // Store number of electrons for later retrieval
    // (In real implementation, would pass to sensitive detector)
}

void GarfieldModel::SetElectricField(double field) {
    m_electricField = field;
    G4cout << "GarfieldModel: Electric field set to " << field << " V/cm" << G4endl;
}

void GarfieldModel::SetGasProperties(double arFrac, double co2Frac,
                                      double pressure, double temperature) {
    m_arFraction = arFrac;
    m_co2Fraction = co2Frac;
    m_pressure = pressure;
    m_temperature = temperature;

    if (m_gas) {
        m_gas->SetComposition("ar", arFrac * 100., "co2", co2Frac * 100.);
        m_gas->SetTemperature(temperature);
        m_gas->SetPressure(pressure);
        m_gas->Initialise();
    }
}

} // namespace nnbar

#endif // WITH_GARFIELD
