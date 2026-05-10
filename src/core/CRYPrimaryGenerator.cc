#include "config.h"

#if WITH_CRY

#include "core/CRYPrimaryGenerator.hh"
#include "output/ParquetOutputManager.hh"

#include "G4AutoLock.hh"
#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYSetup.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <sstream>
#include <vector>

extern G4double event_number_global;

namespace {
G4Mutex CRYPrimaryGeneratorMutex = G4MUTEX_INITIALIZER;

G4String DefaultCRYDataPath() {
    if (const char* cryDir = std::getenv("CRY_DIR")) {
        return G4String(cryDir) + "/data";
    }
    return "./data";
}

double CRYGeant4Random() {
    return G4UniformRand();
}

int ClampIndex(int value, int minValue, int maxValue) {
    return std::max(minValue, std::min(value, maxValue));
}
}  // namespace

const double CRYPrimaryGenerator::N_ij[6][5] = {
    {1.69e11, 2.30e12, 4.02e11, 4.33e11, 2.04e10},
    {1.90e11, 1.09e10, 1.05e10, 1.23e10, 4.34e9},
    {7.69e11, 6.21e9,  5.63e9,  6.03e9,  3.24e9},
    {2.68e11, 7.23e8,  2.24e8,  1.28e8,  1.44e8},
    {2.18e11, 2.30e7,  0.0,     5.92e7,  8.37e7},
    {2.00e11, 0.0,     0.0,     6.25e6,  5.00e6},
};

CRYPrimaryGenerator::CRYPrimaryGenerator()
    : fGun(new G4ParticleGun()), fDataPath(DefaultCRYDataPath()) {
    fWeight = ComputeWeight(fEnergyBinIdx, fParticleIdx);
}

CRYPrimaryGenerator::~CRYPrimaryGenerator() {
    delete fGenerator;
    delete fSetup;
    delete fGun;
}

void CRYPrimaryGenerator::SetParticleType(const G4String& name) {
    if (name == fParticleType) return;
    fParticleType = name;
    fSetupDirty = true;
}

void CRYPrimaryGenerator::SetEnergyMin(double emin_GeV) {
    fEmin = emin_GeV;
}

void CRYPrimaryGenerator::SetEnergyMax(double emax_GeV) {
    fEmax = emax_GeV;
}

void CRYPrimaryGenerator::SetEnergyBinIndex(int i) {
    fEnergyBinIdx = ClampIndex(i, 0, 5);
    fWeight = ComputeWeight(fEnergyBinIdx, fParticleIdx);
}

void CRYPrimaryGenerator::SetParticleIndex(int j) {
    fParticleIdx = ClampIndex(j, 0, 4);
    fWeight = ComputeWeight(fEnergyBinIdx, fParticleIdx);
}

void CRYPrimaryGenerator::SetCRYDataPath(const G4String& path) {
    if (path.empty() || path == fDataPath) return;
    fDataPath = path;
    fSetupDirty = true;
}

void CRYPrimaryGenerator::UpdateCRYSetup() {
    delete fGenerator;
    delete fSetup;
    fGenerator = nullptr;
    fSetup = nullptr;

    const bool wantMuon = (fParticleType == "mu-" || fParticleType == "mu+");
    const bool wantGamma = (fParticleType == "gamma");
    const bool wantElectron = (fParticleType == "e-");
    const bool wantNeutron = (fParticleType == "neutron");
    const bool wantProton = (fParticleType == "proton");

    std::ostringstream config;
    config << "date 1-1-2024 latitude 55.71 altitude 0 subboxLength 2400 "
           << "returnMuons " << (wantMuon ? 1 : 0) << ' '
           << "returnGammas " << (wantGamma ? 1 : 0) << ' '
           << "returnElectrons " << (wantElectron ? 1 : 0) << ' '
           << "returnNeutrons " << (wantNeutron ? 1 : 0) << ' '
           << "returnProtons " << (wantProton ? 1 : 0) << ' '
           << "returnPions 0 returnKaons 0 nParticlesMin 1 nParticlesMax 1000000";

    fSetup = new CRYSetup(config.str(), std::string(fDataPath));
    fSetup->setRandomFunction(CRYGeant4Random);
    fGenerator = new CRYGenerator(fSetup);
    fSetupDirty = false;
}

int CRYPrimaryGenerator::TargetPDG() const {
    if (fParticleType == "mu-") return 13;
    if (fParticleType == "mu+") return -13;
    if (fParticleType == "gamma") return 22;
    if (fParticleType == "e-") return 11;
    if (fParticleType == "neutron") return 2112;
    if (fParticleType == "proton") return 2212;
    return 13;
}

bool CRYPrimaryGenerator::MatchesRequestedParticle(CRYParticle& particle) const {
    return particle.PDGid() == TargetPDG();
}

G4ParticleGun* CRYPrimaryGenerator::Gun() {
    if (!fGun) {
        fGun = new G4ParticleGun();
    }
    return fGun;
}

double CRYPrimaryGenerator::ComputeWeight(int iBin, int jParticle) const {
    iBin = ClampIndex(iBin, 0, 5);
    jParticle = ClampIndex(jParticle, 0, 4);
    const double n_ij = N_ij[iBin][jParticle];
    double sum_i = 0.0;
    for (int i = 0; i < 6; ++i) {
        sum_i += N_ij[i][jParticle];
    }
    if (sum_i <= 0.0 || n_ij <= 0.0) {
        return 0.0;
    }
    return (n_ij / S) * (n_ij / sum_i);
}

void CRYPrimaryGenerator::GenerateCRYPrimaries(G4Event* anEvent) {
    G4AutoLock lock(&CRYPrimaryGeneratorMutex);

    if (!fGenerator || fSetupDirty) {
        UpdateCRYSetup();
    }

    anEvent->SetEventID(static_cast<G4int>(event_number_global));
    event_number_global++;

    CRYParticle* selected = nullptr;
    std::vector<CRYParticle*> particles;
    for (int attempt = 0; attempt < 10000 && selected == nullptr; ++attempt) {
        particles.clear();
        fGenerator->genEvent(&particles);
        for (CRYParticle* particle : particles) {
            if (particle && MatchesRequestedParticle(*particle)) {
                selected = particle;
                break;
            }
        }
        if (!selected) {
            for (CRYParticle* particle : particles) delete particle;
            particles.clear();
        }
    }

    if (!selected) {
        G4cerr << "[CRY] No particle matching " << fParticleType
               << " after 10000 CRY showers; event has no primary" << G4endl;
        return;
    }

    G4ParticleDefinition* definition =
        G4ParticleTable::GetParticleTable()->FindParticle(TargetPDG());
    if (!definition) {
        G4cerr << "[CRY] Geant4 particle not found for " << fParticleType << G4endl;
        for (CRYParticle* particle : particles) delete particle;
        return;
    }

    G4double ke = fEmin * CLHEP::GeV;
    if (fEmax > fEmin) {
        ke = (fEmin + G4UniformRand() * (fEmax - fEmin)) * CLHEP::GeV;
    }
    const G4double x = selected->x() * CLHEP::cm;
    const G4double y = selected->y() * CLHEP::cm;
    const G4double z = 500.0 * CLHEP::cm;
    const G4double t = selected->t() * CLHEP::s;
    G4ThreeVector dir(selected->u(), selected->v(), selected->w());
    if (dir.mag2() == 0.0) {
        dir = G4ThreeVector(0.0, 0.0, -1.0);
    } else {
        dir = dir.unit();
    }

    Gun()->SetParticleDefinition(definition);
    Gun()->SetParticleEnergy(ke);
    Gun()->SetParticlePosition(G4ThreeVector(x, y, z));
    Gun()->SetParticleMomentumDirection(dir);
    Gun()->SetParticleTime(t);

    nnbar::ParticleRecord rec;
    rec.event_id = anEvent->GetEventID();
    rec.pid = definition->GetPDGEncoding();
    rec.mass = definition->GetPDGMass() / CLHEP::MeV;
    rec.name = definition->GetParticleName();
    rec.charge = definition->GetPDGCharge();
    rec.ke = ke / CLHEP::MeV;
    rec.angle = std::acos(std::max(-1.0, std::min(1.0, dir.z())));
    rec.x = x / CLHEP::cm;
    rec.y = y / CLHEP::cm;
    rec.z = z / CLHEP::cm;
    rec.t = t / CLHEP::ms;
    rec.u = dir.x();
    rec.v = dir.y();
    rec.w = dir.z();
    rec.weight = fWeight;
    nnbar::ParquetOutputManager::Instance().WriteParticle(rec);

    Gun()->GeneratePrimaryVertex(anEvent);

    for (CRYParticle* particle : particles) delete particle;
}

#endif  // WITH_CRY
