// ============================================================================
// PrimaryGeneratorAction.hh
// Primary particle generator for NNBAR detector simulation
// Supports both random signal generation and controlled calibration mode
// ============================================================================

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "config.h"

class G4ParticleGun;
class G4Event;
class G4GenericMessenger;
class CRYPrimaryGenerator;

// Calibration target types
enum class CalibrationTarget {
    NONE,           // Normal signal mode
    SCINTILLATOR,   // Target scintillator modules
    LEADGLASS,      // Target lead glass modules
    TPC_SURFACE     // Generate from TPC inner surface (for background)
};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    PrimaryGeneratorAction();
    virtual ~PrimaryGeneratorAction();

    virtual void GeneratePrimaries(G4Event* anEvent) override;

    // Calibration mode setters (called via messenger commands)
    void SetCalibrationMode(G4bool enabled) { sCalibrationMode = enabled; }
    void SetCalibrationTarget(const G4String& target);
    void SetCalibrationEnergy(G4double energy) { sCalibEnergy = energy; }
    void SetCalibrationEnergyMin(G4double energy) { sCalibEnergyMin = energy; }
    void SetCalibrationEnergyMax(G4double energy) { sCalibEnergyMax = energy; }
    void SetCalibrationParticle(const G4String& particle) { sCalibParticle = particle; }
    void SetCalibrationSurface(G4int surface) { sCalibSurface = surface; }
    void SetCalibrationSpread(G4double spread) { sCalibSpread = spread; }

    // CRY cosmic mode setters (called via /cosmic/ messenger commands)
    void SetCRYMode(G4bool enabled) { sCRYMode = enabled; }
    void SetCRYParticle(const G4String& particle) { sCRYParticle = particle; }
    void SetCRYEnergyMin(G4double energy) { sCRYEnergyMin = energy; }
    void SetCRYEnergyMax(G4double energy) { sCRYEnergyMax = energy; }
    void SetCRYEnergyBinIndex(G4int idx) { sCRYEnergyBinIdx = idx; }
    void SetCRYParticleIndex(G4int idx) { sCRYParticleIdx = idx; }
    void SetCRYDataPath(const G4String& path) { sCRYDataPath = path; }

private:
    void DefineCommands();
    void GenerateSignalPrimaries(G4Event* anEvent);
    void GenerateCalibrationPrimaries(G4Event* anEvent);
    void GenerateTPCSurfacePrimaries(G4Event* anEvent);
#if WITH_CRY
    void ConfigureCRYGenerator();
#endif
    G4ThreeVector GetScintillatorDirection(G4int surface);
    G4ThreeVector GetLeadGlassDirection(G4int surface);

    G4ParticleGun* fParticleGun;
    G4GenericMessenger* fMessenger;
    G4GenericMessenger* fCosmicMessenger = nullptr;
#if WITH_CRY
    CRYPrimaryGenerator* fCRYGenerator = nullptr;
#endif

    // Static calibration mode parameters (shared across all threads)
    static G4bool sCalibrationMode;
    static CalibrationTarget sCalibTarget;
    static G4double sCalibEnergy;       // Fixed energy (if min==max)
    static G4double sCalibEnergyMin;    // Minimum energy for range
    static G4double sCalibEnergyMax;    // Maximum energy for range
    static G4String sCalibParticle;     // Particle type for calibration
    static G4int sCalibSurface;         // Target surface (0=top, 1=right, 2=bottom, 3=left, 4=front, 5=back)
    static G4double sCalibSpread;       // Angular spread in degrees

    // Static signal mode parameters (shared across all threads)
    static G4String sSignalParticle;    // Particle type for signal mode
    static G4String sSignalParticles;   // Optional comma-separated signal particles
    static G4double sSignalEnergyMin;   // Minimum energy for signal
    static G4double sSignalEnergyMax;   // Maximum energy for signal

    // Static CRY cosmic mode parameters (shared across all threads)
    static G4bool sCRYMode;
    static G4String sCRYParticle;
    static G4double sCRYEnergyMin;
    static G4double sCRYEnergyMax;
    static G4int sCRYEnergyBinIdx;
    static G4int sCRYParticleIdx;
    static G4String sCRYDataPath;

};

#endif // PrimaryGeneratorAction_h
