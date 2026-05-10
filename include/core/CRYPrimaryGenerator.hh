#pragma once

#ifdef WITH_CRY

#include "G4String.hh"
#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

class CRYSetup;
class CRYGenerator;
class CRYParticle;
class G4Event;
class G4ParticleGun;

class CRYPrimaryGenerator {
public:
    CRYPrimaryGenerator();
    ~CRYPrimaryGenerator();

    CRYPrimaryGenerator(const CRYPrimaryGenerator&) = delete;
    CRYPrimaryGenerator& operator=(const CRYPrimaryGenerator&) = delete;

    void GenerateCRYPrimaries(G4Event* anEvent);

    void SetParticleType(const G4String& name);
    void SetEnergyMin(double emin_GeV);
    void SetEnergyMax(double emax_GeV);
    void SetEnergyBinIndex(int i);
    void SetParticleIndex(int j);
    void SetCRYDataPath(const G4String& path);

private:
    CRYSetup* fSetup = nullptr;
    CRYGenerator* fGenerator = nullptr;
    G4ParticleGun* fGun = nullptr;
    G4String fParticleType = "mu-";
    double fEmin = 0.0;
    double fEmax = 0.5;
    int fEnergyBinIdx = 0;
    int fParticleIdx = 0;
    double fWeight = 1.0;
    G4String fDataPath;
    bool fSetupDirty = true;

    void UpdateCRYSetup();
    double ComputeWeight(int iBin, int jParticle) const;
    int TargetPDG() const;
    bool MatchesRequestedParticle(CRYParticle& particle) const;
    G4ParticleGun* Gun();

    static const double N_ij[6][5];
    static constexpr double S = 1e6;
};

#endif  // WITH_CRY
