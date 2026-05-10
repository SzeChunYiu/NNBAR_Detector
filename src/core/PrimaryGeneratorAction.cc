// ============================================================================
// PrimaryGeneratorAction.cc
// Primary particle generator for NNBAR detector simulation
// Supports both random signal generation and controlled calibration mode
// ============================================================================

#include "core/PrimaryGeneratorAction.hh"
#if WITH_CRY
#include "core/CRYPrimaryGenerator.hh"
#endif
#include "output/ParquetOutputManager.hh"

#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "Randomize.hh"

#include <cmath>
#include <sstream>
#include <vector>

// Global state
extern G4double event_number_global;
extern G4int run_number;
extern int theta_bin_index;
extern int KE_bin_index;
extern G4int particle_name_input;

namespace {
    G4Mutex PrimaryGeneratorMutex = G4MUTEX_INITIALIZER;
    constexpr double pi = 3.14159265358979323846;

    G4String TrimSignalParticle(const G4String& token) {
        const auto first = token.find_first_not_of(" \t\n\r");
        if (first == G4String::npos) {
            return "";
        }
        const auto last = token.find_last_not_of(" \t\n\r");
        return token.substr(first, last - first + 1);
    }

    std::vector<G4String> SplitSignalParticles(const G4String& configuredParticles) {
        std::vector<G4String> particles;
        std::stringstream stream(configuredParticles);
        std::string token;
        while (std::getline(stream, token, ',')) {
            G4String particle = TrimSignalParticle(token);
            if (!particle.empty()) {
                particles.push_back(particle);
            }
        }
        return particles;
    }
}

// Static member initialization - shared across all threads
G4bool PrimaryGeneratorAction::sCalibrationMode = false;
CalibrationTarget PrimaryGeneratorAction::sCalibTarget = CalibrationTarget::NONE;
G4double PrimaryGeneratorAction::sCalibEnergy = 100.0 * CLHEP::MeV;
G4double PrimaryGeneratorAction::sCalibEnergyMin = 100.0 * CLHEP::MeV;
G4double PrimaryGeneratorAction::sCalibEnergyMax = 100.0 * CLHEP::MeV;
G4String PrimaryGeneratorAction::sCalibParticle = "e-";
G4int PrimaryGeneratorAction::sCalibSurface = 0;
G4double PrimaryGeneratorAction::sCalibSpread = 5.0;

// Signal mode static parameters
G4String PrimaryGeneratorAction::sSignalParticle = "pi+";
G4String PrimaryGeneratorAction::sSignalParticles = "";
G4double PrimaryGeneratorAction::sSignalEnergyMin = 50.0 * CLHEP::MeV;
G4double PrimaryGeneratorAction::sSignalEnergyMax = 500.0 * CLHEP::MeV;

// CRY cosmic mode static parameters
G4bool PrimaryGeneratorAction::sCRYMode = false;
G4String PrimaryGeneratorAction::sCRYParticle = "mu-";
G4double PrimaryGeneratorAction::sCRYEnergyMin = 0.0 * CLHEP::GeV;
G4double PrimaryGeneratorAction::sCRYEnergyMax = 0.5 * CLHEP::GeV;
G4int PrimaryGeneratorAction::sCRYEnergyBinIdx = 0;
G4int PrimaryGeneratorAction::sCRYParticleIdx = 0;
G4String PrimaryGeneratorAction::sCRYDataPath = "";

// Thread-local state
G4ThreadLocal G4int local_event_number;

// Particle types available for signal mode
static const std::vector<G4String> particle_names = {"pi+", "pi-", "e-", "gamma", "pi0"};

// Energy ranges (MeV)
static const std::vector<std::vector<G4double>> Energy_ranges = {
    {0.5, 0.5},
    {2.5, 2.5},
    {100., 100.},
    {400., 400.},
    {50., 100.},
    {100., 200.},
    {200., 500.},
    {500., 700.}
};

// Angle ranges (degrees)
static const std::vector<std::vector<G4double>> angle_ranges = {
    {0., 0.}
};

// Detector geometry constants (from Scintillator_geometry.cc and LeadGlass_geometry.cc)
// These define approximate distances to detector surfaces
namespace DetectorGeometry {
    // Lead glass is positioned beyond the scintillators
    // From lead_glass_position.csv: first block at x=-247.5cm, y=+261.9cm
    // These are the base positions before rotation for each surface
    constexpr G4double lg_x0 = -247.5 * cm;  // x offset in base position
    constexpr G4double lg_y0 = 261.9 * cm;   // y offset (radial distance)
}

// ============================================================================
// Constructor
// ============================================================================

PrimaryGeneratorAction::PrimaryGeneratorAction()
    : fParticleGun(nullptr)
{
    fParticleGun = new G4ParticleGun();
#if WITH_CRY
    fCRYGenerator = new CRYPrimaryGenerator();
#endif
    DefineCommands();
}

// ============================================================================
// Destructor
// ============================================================================

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete fParticleGun;
    delete fMessenger;
    delete fCosmicMessenger;
#if WITH_CRY
    delete fCRYGenerator;
#endif
}

// ============================================================================
// Define Messenger Commands for Calibration Control
// ============================================================================

void PrimaryGeneratorAction::DefineCommands() {
    fMessenger = new G4GenericMessenger(this, "/calibration/",
                                        "Calibration mode controls");

    // Enable/disable calibration mode - use static variable
    fMessenger->DeclareProperty("mode", sCalibrationMode,
        "Enable calibration mode (true/false)");

    // Set calibration target (scintillator or leadglass)
    fMessenger->DeclareMethod("target", &PrimaryGeneratorAction::SetCalibrationTarget,
        "Set calibration target (scintillator/leadglass/tpc_surface)");

    // Set fixed calibration energy - use static variable
    fMessenger->DeclarePropertyWithUnit("energy", "MeV", sCalibEnergy,
        "Set calibration energy (MeV)");

    // Set energy range - use static variables
    fMessenger->DeclarePropertyWithUnit("energy_min", "MeV", sCalibEnergyMin,
        "Set minimum calibration energy (MeV)");
    fMessenger->DeclarePropertyWithUnit("energy_max", "MeV", sCalibEnergyMax,
        "Set maximum calibration energy (MeV)");

    // Set particle type - use static variable
    fMessenger->DeclareProperty("particle", sCalibParticle,
        "Set calibration particle (e-, mu-, pi+, pi-, gamma, proton)");

    // Set target surface - use static variable
    fMessenger->DeclareProperty("surface", sCalibSurface,
        "Set target surface (0=top, 1=right, 2=bottom, 3=left, 4=front, 5=back)");

    // Set angular spread - use static variable
    fMessenger->DeclareProperty("spread", sCalibSpread,
        "Set angular spread in degrees");

    // ========== Signal mode commands ==========
    // Signal mode particle type
    fMessenger->DeclareProperty("signal_particle", sSignalParticle,
        "Set signal mode particle (pi+, pi-, proton, K+, K-)");

    // Signal mode particle list. When set, every listed particle is emitted in
    // the same Geant4 event; otherwise signal_particle preserves legacy mode.
    fMessenger->DeclareProperty("signal_particles", sSignalParticles,
        "Set comma-separated signal particles for one event (pi+,pi-,proton)");

    // Signal mode energy range
    fMessenger->DeclarePropertyWithUnit("signal_energy_min", "MeV", sSignalEnergyMin,
        "Set minimum signal energy (MeV)");
    fMessenger->DeclarePropertyWithUnit("signal_energy_max", "MeV", sSignalEnergyMax,
        "Set maximum signal energy (MeV)");


    // ========== CRY cosmic mode commands ==========
    fCosmicMessenger = new G4GenericMessenger(this, "/cosmic/",
                                             "CRY cosmic ray generator controls");
    fCosmicMessenger->DeclareProperty("mode", sCRYMode,
        "Enable CRY cosmic generator mode (true/false)");
    fCosmicMessenger->DeclareProperty("particle", sCRYParticle,
        "Set CRY particle (mu-, mu+, gamma, e-, neutron, proton)");
    fCosmicMessenger->DeclarePropertyWithUnit("energyMin", "GeV", sCRYEnergyMin,
        "Set minimum CRY kinetic energy (GeV)");
    fCosmicMessenger->DeclarePropertyWithUnit("energyMax", "GeV", sCRYEnergyMax,
        "Set maximum CRY kinetic energy (GeV)");
    fCosmicMessenger->DeclareProperty("energyBin", sCRYEnergyBinIdx,
        "Set CRY energy-bin index (0-5) for weight lookup");
    fCosmicMessenger->DeclareProperty("particleIdx", sCRYParticleIdx,
        "Set CRY particle index (0=mu, 1=gamma, 2=e-, 3=neutron, 4=proton)");
    fCosmicMessenger->DeclareProperty("dataPath", sCRYDataPath,
        "Set CRY data directory path");
}

// ============================================================================
// Set Calibration Target
// ============================================================================

void PrimaryGeneratorAction::SetCalibrationTarget(const G4String& target) {
    if (target == "scintillator" || target == "SCINTILLATOR") {
        sCalibTarget = CalibrationTarget::SCINTILLATOR;
        G4cout << "Calibration target set to: SCINTILLATOR" << G4endl;
    } else if (target == "leadglass" || target == "LEADGLASS" ||
               target == "lead_glass" || target == "LEAD_GLASS") {
        sCalibTarget = CalibrationTarget::LEADGLASS;
        G4cout << "Calibration target set to: LEADGLASS" << G4endl;
    } else if (target == "tpc_surface" || target == "TPC_SURFACE" ||
               target == "tpc" || target == "TPC") {
        sCalibTarget = CalibrationTarget::TPC_SURFACE;
        G4cout << "Calibration target set to: TPC_SURFACE (background mode)" << G4endl;
    } else {
        sCalibTarget = CalibrationTarget::NONE;
        G4cout << "Calibration target set to: NONE (signal mode)" << G4endl;
    }
}



#if WITH_CRY
void PrimaryGeneratorAction::ConfigureCRYGenerator() {
    fCRYGenerator->SetParticleType(sCRYParticle);
    fCRYGenerator->SetEnergyMin(sCRYEnergyMin / CLHEP::GeV);
    fCRYGenerator->SetEnergyMax(sCRYEnergyMax / CLHEP::GeV);
    fCRYGenerator->SetEnergyBinIndex(sCRYEnergyBinIdx);
    fCRYGenerator->SetParticleIndex(sCRYParticleIdx);
    if (!sCRYDataPath.empty()) {
        fCRYGenerator->SetCRYDataPath(sCRYDataPath);
    }
}
#endif

// ============================================================================
// Generate Primaries - Main Entry Point
// ============================================================================

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
#if WITH_CRY
    if (sCRYMode && fCRYGenerator) {
        ConfigureCRYGenerator();
        fCRYGenerator->GenerateCRYPrimaries(anEvent);
        return;
    }
#endif

    // Use static variables for calibration mode (shared across threads)
    if (sCalibrationMode && sCalibTarget != CalibrationTarget::NONE) {
        if (sCalibTarget == CalibrationTarget::TPC_SURFACE) {
            GenerateTPCSurfacePrimaries(anEvent);
        } else {
            GenerateCalibrationPrimaries(anEvent);
        }
    } else {
        GenerateSignalPrimaries(anEvent);
    }
}

// ============================================================================
// Generate Signal Primaries (Original Behavior)
// ============================================================================

void PrimaryGeneratorAction::GenerateSignalPrimaries(G4Event* anEvent) {
    G4AutoLock lock(&PrimaryGeneratorMutex);

    // Set event ID from global counter
    anEvent->SetEventID(static_cast<G4int>(event_number_global));
    local_event_number = anEvent->GetEventID();
    event_number_global++;

    // Generate random position (at origin for signal events)
    G4double x = 0.0;
    G4double y = 0.0;
    G4double z = 0.0;
    G4double t = 0.0;

    std::vector<G4String> signalParticles = SplitSignalParticles(sSignalParticles);
    if (signalParticles.empty()) {
        signalParticles.push_back(sSignalParticle);
    }

    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

    for (const G4String& signalParticle : signalParticles) {
        // Generate random kinetic energy within configured range
        G4double KE;
        if (std::abs(sSignalEnergyMax - sSignalEnergyMin) < 0.001 * MeV) {
            KE = sSignalEnergyMin;
        } else {
            KE = sSignalEnergyMin + G4UniformRand() * (sSignalEnergyMax - sSignalEnergyMin);
        }

        // Generate isotropic momentum direction (uniform on sphere)
        G4double cosTheta = 2.0 * G4UniformRand() - 1.0;  // Uniform in [-1, 1]
        G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
        G4double phi = 2.0 * pi * G4UniformRand();  // Uniform in [0, 2pi]

        G4double px = sinTheta * std::cos(phi);
        G4double py = sinTheta * std::sin(phi);
        G4double pz = cosTheta;

        G4double weight = 1.0;

        // Select particle type from static variable
        G4ParticleDefinition* particle = particleTable->FindParticle(signalParticle);
        if (!particle) {
            G4cerr << "Unknown signal particle: " << signalParticle
                   << ", defaulting to pi+" << G4endl;
            particle = particleTable->FindParticle("pi+");
        }

        // Configure particle gun
        fParticleGun->SetParticleDefinition(particle);
        fParticleGun->SetParticleEnergy(KE);
        fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
        fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz));
        fParticleGun->SetParticleTime(t);

        // Write particle record to Parquet output
        nnbar::ParticleRecord rec;
        rec.event_id = anEvent->GetEventID();
        rec.pid = particle->GetPDGEncoding();
        rec.mass = particle->GetPDGMass() / CLHEP::MeV;
        rec.name = particle->GetParticleName();
        rec.charge = particle->GetPDGCharge();
        rec.ke = KE / CLHEP::MeV;
        rec.angle = 0.0;
        rec.x = x / CLHEP::cm;
        rec.y = y / CLHEP::cm;
        rec.z = z / CLHEP::cm;
        rec.t = t / CLHEP::ms;
        rec.u = px;
        rec.v = py;
        rec.w = pz;
        rec.weight = weight;

        nnbar::ParquetOutputManager::Instance().WriteParticle(rec);

        // Generate one primary vertex for this particle in the shared event.
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}

// ============================================================================
// Generate Calibration Primaries
// ============================================================================

void PrimaryGeneratorAction::GenerateCalibrationPrimaries(G4Event* anEvent) {
    G4AutoLock lock(&PrimaryGeneratorMutex);

    // Set event ID from global counter
    anEvent->SetEventID(static_cast<G4int>(event_number_global));
    local_event_number = anEvent->GetEventID();
    event_number_global++;

    // Start position at origin (inside TPC)
    G4double x = 0.0;
    G4double y = 0.0;
    G4double z = 0.0;
    G4double t = 0.0;

    // Calculate kinetic energy (fixed or random within range) - use static variables
    G4double KE;
    if (std::abs(sCalibEnergyMax - sCalibEnergyMin) < 0.001 * MeV) {
        // Fixed energy mode
        KE = sCalibEnergy;
    } else {
        // Random energy within range (uniform distribution)
        KE = sCalibEnergyMin + G4UniformRand() * (sCalibEnergyMax - sCalibEnergyMin);
    }

    // Get base direction toward target - use static variables
    G4ThreeVector baseDir;
    if (sCalibTarget == CalibrationTarget::SCINTILLATOR) {
        baseDir = GetScintillatorDirection(sCalibSurface);
    } else {
        baseDir = GetLeadGlassDirection(sCalibSurface);
    }

    // Apply angular spread (Gaussian distribution) - use static variable
    G4double spreadRad = sCalibSpread * pi / 180.0;
    G4double theta = G4RandGauss::shoot(0.0, spreadRad);
    G4double phi = G4UniformRand() * 2.0 * pi;

    // Create perpendicular vectors for rotation
    G4ThreeVector perpX, perpY;
    if (std::abs(baseDir.x()) < 0.9) {
        perpX = baseDir.cross(G4ThreeVector(1, 0, 0)).unit();
    } else {
        perpX = baseDir.cross(G4ThreeVector(0, 1, 0)).unit();
    }
    perpY = baseDir.cross(perpX).unit();

    // Apply rotation
    G4ThreeVector dir = baseDir * std::cos(theta)
                      + perpX * std::sin(theta) * std::cos(phi)
                      + perpY * std::sin(theta) * std::sin(phi);
    dir = dir.unit();

    // Get particle definition - use static variable
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(sCalibParticle);
    if (!particle) {
        G4cerr << "Unknown calibration particle: " << sCalibParticle
               << ", defaulting to e-" << G4endl;
        particle = particleTable->FindParticle("e-");
    }

    // Configure particle gun
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(KE);
    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
    fParticleGun->SetParticleMomentumDirection(dir);
    fParticleGun->SetParticleTime(t);

    // Write particle record to Parquet output
    nnbar::ParticleRecord rec;
    rec.event_id = anEvent->GetEventID();
    rec.pid = particle->GetPDGEncoding();
    rec.mass = particle->GetPDGMass() / CLHEP::MeV;
    rec.name = particle->GetParticleName();
    rec.charge = particle->GetPDGCharge();
    rec.ke = KE / CLHEP::MeV;
    rec.angle = 0.0;
    rec.x = x / CLHEP::cm;
    rec.y = y / CLHEP::cm;
    rec.z = z / CLHEP::cm;
    rec.t = t / CLHEP::ms;
    rec.u = dir.x();
    rec.v = dir.y();
    rec.w = dir.z();
    rec.weight = 1.0;

    nnbar::ParquetOutputManager::Instance().WriteParticle(rec);

    // Generate the primary vertex
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // Log first few events for verification with direction details
    if (anEvent->GetEventID() < 5) {
        G4cout << "[Calibration] Event " << anEvent->GetEventID()
               << ": " << sCalibParticle << " @ " << KE/MeV << " MeV"
               << " -> surface " << sCalibSurface
               << " dir=(" << dir.x() << ", " << dir.y() << ", " << dir.z() << ")"
               << G4endl;
    }
}

// ============================================================================
// Get Direction Toward Scintillator Surface
// ============================================================================

G4ThreeVector PrimaryGeneratorAction::GetScintillatorDirection(G4int surface) {
    // Surface indices: 0=top, 1=right, 2=bottom, 3=left, 4=front, 5=back
    switch (surface) {
        case 0: return G4ThreeVector(0, 1, 0);   // Top (+Y)
        case 1: return G4ThreeVector(1, 0, 0);   // Right (+X)
        case 2: return G4ThreeVector(0, -1, 0);  // Bottom (-Y)
        case 3: return G4ThreeVector(-1, 0, 0);  // Left (-X)
        case 4: return G4ThreeVector(0, 0, 1);   // Front (+Z)
        case 5: return G4ThreeVector(0, 0, -1);  // Back (-Z)
        default:
            G4cerr << "Invalid surface index: " << surface << ", using top" << G4endl;
            return G4ThreeVector(0, 1, 0);
    }
}

// ============================================================================
// Get Direction Toward Lead Glass Surface
// ============================================================================

G4ThreeVector PrimaryGeneratorAction::GetLeadGlassDirection(G4int surface) {
    // Lead glass blocks are positioned at specific angles, not along pure axes
    // Base position from CSV: x0=-247.5cm, y0=+261.9cm (center of array)
    // Each surface is rotated by 90 degrees
    using namespace DetectorGeometry;

    G4double angle = surface * 90.0 * CLHEP::deg;
    G4double cos_a = std::cos(angle);
    G4double sin_a = std::sin(angle);

    switch (surface) {
        case 0: case 1: case 2: case 3: {
            // Radial surfaces: calculate actual position after rotation
            // lead_x = lg_x0 * cos(angle) - lg_y0 * sin(angle)
            // lead_y = lg_x0 * sin(angle) + lg_y0 * cos(angle)
            G4double target_x = lg_x0 * cos_a - lg_y0 * sin_a;
            G4double target_y = lg_x0 * sin_a + lg_y0 * cos_a;
            G4ThreeVector dir(target_x, target_y, 0);
            return dir.unit();
        }
        case 4: // Front (+Z)
            return G4ThreeVector(0, 0, 1);
        case 5: // Back (-Z)
            return G4ThreeVector(0, 0, -1);
        default:
            G4cerr << "Invalid lead glass surface: " << surface << ", using surface 0" << G4endl;
            G4double target_x = lg_x0;
            G4double target_y = lg_y0;
            return G4ThreeVector(target_x, target_y, 0).unit();
    }
}

// ============================================================================
// Generate Primaries from TPC Inner Surface (for Compton electron background)
// ============================================================================

void PrimaryGeneratorAction::GenerateTPCSurfacePrimaries(G4Event* anEvent) {
    G4AutoLock lock(&PrimaryGeneratorMutex);

    // Set event ID from global counter
    anEvent->SetEventID(static_cast<G4int>(event_number_global));
    local_event_number = anEvent->GetEventID();
    event_number_global++;

    // TPC geometry: inner radius at beampipe outer surface
    constexpr G4double tpc_inner_radius = 114.0 * cm;  // Just outside beampipe (112 cm)
    constexpr G4double tpc_half_z = 250.0 * cm;        // TPC half-length

    // Generate random position on TPC inner cylindrical surface
    G4double phi = G4UniformRand() * 2.0 * pi;
    G4double z = (G4UniformRand() * 2.0 - 1.0) * tpc_half_z;  // -250 to +250 cm

    G4double x = tpc_inner_radius * std::cos(phi);
    G4double y = tpc_inner_radius * std::sin(phi);
    G4double t = 0.0;

    // Calculate kinetic energy (use calibration energy settings)
    G4double KE;
    if (std::abs(sCalibEnergyMax - sCalibEnergyMin) < 0.001 * MeV) {
        KE = sCalibEnergy;
    } else {
        // Compton electrons typically have energy spectrum, use uniform for now
        KE = sCalibEnergyMin + G4UniformRand() * (sCalibEnergyMax - sCalibEnergyMin);
    }

    // Base direction: OUTWARD from the beampipe surface into TPC gas
    // The TPC gas surrounds the beampipe, so electrons should travel outward
    G4ThreeVector radialDir(std::cos(phi), std::sin(phi), 0.0);

    // Apply angular spread for more realistic Compton scattering directions
    G4double spreadRad = sCalibSpread * pi / 180.0;
    G4double theta_spread = G4RandGauss::shoot(0.0, spreadRad);
    G4double phi_spread = G4UniformRand() * 2.0 * pi;

    // Create perpendicular vectors for spread rotation
    G4ThreeVector zAxis(0, 0, 1);
    G4ThreeVector perpX = radialDir.cross(zAxis).unit();
    G4ThreeVector perpY = radialDir.cross(perpX).unit();

    // Apply spread to create final direction (mostly inward with some spread)
    G4ThreeVector dir = radialDir * std::cos(theta_spread)
                      + perpX * std::sin(theta_spread) * std::cos(phi_spread)
                      + perpY * std::sin(theta_spread) * std::sin(phi_spread);
    dir = dir.unit();

    // Get particle definition (typically e- for Compton electrons)
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle(sCalibParticle);
    if (!particle) {
        G4cerr << "Unknown TPC surface particle: " << sCalibParticle
               << ", defaulting to e-" << G4endl;
        particle = particleTable->FindParticle("e-");
    }

    // Configure particle gun
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticleEnergy(KE);
    fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
    fParticleGun->SetParticleMomentumDirection(dir);
    fParticleGun->SetParticleTime(t);

    // Write particle record to Parquet output
    nnbar::ParticleRecord rec;
    rec.event_id = anEvent->GetEventID();
    rec.pid = particle->GetPDGEncoding();
    rec.mass = particle->GetPDGMass() / CLHEP::MeV;
    rec.name = particle->GetParticleName();
    rec.charge = particle->GetPDGCharge();
    rec.ke = KE / CLHEP::MeV;
    rec.angle = 0.0;
    rec.x = x / CLHEP::cm;
    rec.y = y / CLHEP::cm;
    rec.z = z / CLHEP::cm;
    rec.t = t / CLHEP::ms;
    rec.u = dir.x();
    rec.v = dir.y();
    rec.w = dir.z();
    rec.weight = 1.0;

    nnbar::ParquetOutputManager::Instance().WriteParticle(rec);

    // Generate the primary vertex
    fParticleGun->GeneratePrimaryVertex(anEvent);

    // Log first few events for verification
    if (anEvent->GetEventID() < 5) {
        G4cout << "[TPC Surface] Event " << anEvent->GetEventID()
               << ": " << sCalibParticle << " @ " << KE/MeV << " MeV"
               << " pos=(" << x/cm << ", " << y/cm << ", " << z/cm << ") cm"
               << " dir=(" << dir.x() << ", " << dir.y() << ", " << dir.z() << ")"
               << G4endl;
    }
}
