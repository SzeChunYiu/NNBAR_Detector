// ============================================================================
// G4MCPLGenerator.cc
// MCPL file reader and primary particle generator for NNBAR simulation
// ============================================================================

#include "generator/G4MCPLGenerator.hh"
#include "output/ParquetOutputManager.hh"

#include "G4ParticleGun.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4Exception.hh"
#include "G4ios.hh"

#include <cassert>
#include <filesystem>
#include <sstream>

extern G4double event_number_global;
extern G4int run_number;

mcpl_file_t m_mcplfile;
const mcpl_particle_t* m_p;

namespace {
    G4Mutex MCPLMutex = G4MUTEX_INITIALIZER;
    G4String g_mcplInputFile;

    void CloseMCPLFileIfOpen() {
        if (m_mcplfile.internal) {
            mcpl_close_file(m_mcplfile);
            m_mcplfile.internal = nullptr;
        }
        m_p = nullptr;
    }
}

G4ThreadLocal G4double event_number;
G4ThreadLocal uint32_t local_event_number_MCPL;
G4ThreadLocal int Event_found = 0;

// ============================================================================
// Constructor
// ============================================================================

G4MCPLGenerator::G4MCPLGenerator(const G4String& inputFile)
    : G4VUserPrimaryGeneratorAction(),
      m_currentPDG(0),
      m_currentPartDef(nullptr),
      m_nUnknownPDG(0),
      m_inputFile(inputFile)
{
    G4AutoLock lock(&MCPLMutex);
    m_mcplfile.internal = nullptr;
    if (!m_inputFile.empty()) {
        g_mcplInputFile = m_inputFile;
    }
    m_gun = new G4ParticleGun(1);
}

// ============================================================================
// Destructor
// ============================================================================

G4MCPLGenerator::~G4MCPLGenerator() {
    if (m_nUnknownPDG) {
        std::ostringstream cmt;
        cmt << "Ignored a total of " << m_nUnknownPDG
            << " particles in input due to untranslatable PDG codes";
        G4Exception("G4MCPLGenerator::~G4MCPLGenerator()", "G4MCPLGenerator07",
                    JustWarning, cmt.str().c_str());
    }

    // Only close file from master thread
    if (G4Threading::G4GetThreadId() == 0) {
        G4AutoLock lock(&MCPLMutex);
        CloseMCPLFileIfOpen();
    }
    delete m_gun;
}

void G4MCPLGenerator::SetInputFile(const G4String& inputFile) {
    G4AutoLock lock(&MCPLMutex);
    if (inputFile != g_mcplInputFile) {
        CloseMCPLFileIfOpen();
    }
    g_mcplInputFile = inputFile;
    G4cout << "MCPL input file set to: " << g_mcplInputFile << G4endl;
}

G4String G4MCPLGenerator::GetInputFile() {
    G4AutoLock lock(&MCPLMutex);
    return g_mcplInputFile;
}

// ============================================================================
// Use Particle Filter
// ============================================================================

bool G4MCPLGenerator::UseParticle(const mcpl_particle_t*) const {
    return true;
}

// ============================================================================
// Modify Particle (hook for subclasses)
// ============================================================================

void G4MCPLGenerator::ModifyParticle(G4ThreeVector&, G4ThreeVector&,
                                      G4ThreeVector&, G4double&, G4double&) const {
    // Default: no modification
}

// ============================================================================
// Generate Primaries
// ============================================================================

void G4MCPLGenerator::GeneratePrimaries(G4Event* evt) {
    G4AutoLock lock(&MCPLMutex);

    if (!EnsureOpen()) {
        G4RunManager::GetRunManager()->AbortRun(true);
        return;
    }

    auto particleTable = G4ParticleTable::GetParticleTable();

    local_event_number_MCPL = static_cast<uint32_t>(event_number_global);
    evt->SetEventID(static_cast<G4int>(event_number_global));
    event_number_global++;

    Event_found = 0;

    // Process particles with matching event ID
    while (m_p && m_p->userflags == local_event_number_MCPL) {
        Event_found = 1;

        // Look up particle definition
        auto PartDef = particleTable->FindParticle(m_p->pdgcode);
        if (!PartDef && (m_p->pdgcode / 100000000 == 10)) {
            PartDef = G4IonTable::GetIonTable()->GetIon(m_p->pdgcode);
        }

        if (!PartDef) {
            m_p = mcpl_read(m_mcplfile);
            continue;
        }

        // Set up particle properties
        m_gun->SetParticleDefinition(PartDef);

        G4ThreeVector pos(m_p->position[0], m_p->position[1], m_p->position[2]);
        pos *= CLHEP::cm;

        G4ThreeVector dir(m_p->direction[0], m_p->direction[1], m_p->direction[2]);
        G4ThreeVector pol(m_p->polarisation[0], m_p->polarisation[1], m_p->polarisation[2]);

        G4double KE = m_p->ekin;  // Already in MeV
        G4double time = 0.0;
        G4double weight = m_p->weight;

        ModifyParticle(pos, dir, pol, time, weight);

        m_gun->SetParticleMomentumDirection(dir);
        m_gun->SetParticlePosition(pos);
        m_gun->SetParticleEnergy(KE);
        m_gun->SetParticleTime(time);
        m_gun->SetParticlePolarization(pol);

        // Write particle record to Parquet
        nnbar::ParticleRecord rec;
        rec.event_id = evt->GetEventID();
        rec.pid = m_p->pdgcode;
        rec.mass = PartDef->GetPDGMass() / CLHEP::MeV;
        rec.name = PartDef->GetParticleName();
        rec.charge = PartDef->GetPDGCharge();
        rec.ke = KE;
        rec.angle = 0.0;
        rec.x = pos.x() / CLHEP::cm;
        rec.y = pos.y() / CLHEP::cm;
        rec.z = pos.z() / CLHEP::cm;
        rec.t = time / CLHEP::ms;
        rec.u = dir.x();
        rec.v = dir.y();
        rec.w = dir.z();
        rec.weight = weight;

        nnbar::ParquetOutputManager::Instance().WriteParticle(rec);

        m_gun->GeneratePrimaryVertex(evt);

        // Read next particle
        m_p = mcpl_read(m_mcplfile);
    }

    // If current particle has different event ID, we're done with this event
    // Skip to next matching event for future calls
    if (m_p && m_p->userflags != local_event_number_MCPL && Event_found == 1) {
        // Done - particle belongs to next event
        return;
    }

    // Search for particles matching our event ID
    while (m_p && m_p->userflags != local_event_number_MCPL) {
        m_p = mcpl_read(m_mcplfile);
    }

    // Process found particles
    while (m_p && m_p->userflags == local_event_number_MCPL) {
        Event_found = 1;

        auto PartDef = particleTable->FindParticle(m_p->pdgcode);
        if (!PartDef && (m_p->pdgcode / 100000000 == 10)) {
            PartDef = G4IonTable::GetIonTable()->GetIon(m_p->pdgcode);
        }

        if (!PartDef) {
            m_p = mcpl_read(m_mcplfile);
            continue;
        }

        m_gun->SetParticleDefinition(PartDef);

        G4ThreeVector pos(m_p->position[0], m_p->position[1], m_p->position[2]);
        pos *= CLHEP::cm;

        G4ThreeVector dir(m_p->direction[0], m_p->direction[1], m_p->direction[2]);
        G4ThreeVector pol(m_p->polarisation[0], m_p->polarisation[1], m_p->polarisation[2]);

        G4double KE = m_p->ekin;
        G4double time = 0.0;
        G4double weight = m_p->weight;

        ModifyParticle(pos, dir, pol, time, weight);

        m_gun->SetParticleMomentumDirection(dir);
        m_gun->SetParticlePosition(pos);
        m_gun->SetParticleEnergy(KE);
        m_gun->SetParticleTime(time);
        m_gun->SetParticlePolarization(pol);

        // Write particle record
        nnbar::ParticleRecord rec;
        rec.event_id = evt->GetEventID();
        rec.pid = m_p->pdgcode;
        rec.mass = PartDef->GetPDGMass() / CLHEP::MeV;
        rec.name = PartDef->GetParticleName();
        rec.charge = PartDef->GetPDGCharge();
        rec.ke = KE;
        rec.angle = 0.0;
        rec.x = pos.x() / CLHEP::cm;
        rec.y = pos.y() / CLHEP::cm;
        rec.z = pos.z() / CLHEP::cm;
        rec.t = time / CLHEP::ms;
        rec.u = dir.x();
        rec.v = dir.y();
        rec.w = dir.z();
        rec.weight = weight;

        nnbar::ParquetOutputManager::Instance().WriteParticle(rec);

        m_gun->GeneratePrimaryVertex(evt);

        m_p = mcpl_read(m_mcplfile);
    }

    if (Event_found == 0) {
        G4cout << "MCPL file exhausted - aborting run" << G4endl;
        G4RunManager::GetRunManager()->AbortRun(true);
    }
}

// ============================================================================
// Find Next Valid Particle
// ============================================================================

void G4MCPLGenerator::FindNext() {
    while ((m_p = mcpl_read(m_mcplfile))) {
        if (!UseParticle(m_p)) continue;
        if (!(m_p->weight > 0.0)) continue;
        if (m_p->pdgcode == 0) continue;

        if (!LookupPDG(m_p->pdgcode)) {
            ++m_nUnknownPDG;
            if (m_nUnknownPDG <= 100) {
                std::ostringstream cmt;
                cmt << "Ignoring particle with untranslatable PDG code ("
                    << m_p->pdgcode << ")";
                G4Exception("G4MCPLGenerator::GeneratePrimaries()",
                           "G4MCPLGenerator05", JustWarning, cmt.str().c_str());
                if (m_nUnknownPDG == 100) {
                    G4Exception("G4MCPLGenerator::GeneratePrimaries()",
                               "G4MCPLGenerator06", JustWarning,
                               "Limit reached. Suppressing further warnings.");
                }
            }
            continue;
        }
        break;
    }
}

bool G4MCPLGenerator::EnsureOpen() {
    if (m_mcplfile.internal && m_p) {
        return true;
    }

    if (g_mcplInputFile.empty()) {
        G4Exception("G4MCPLGenerator::EnsureOpen()", "G4MCPLGeneratorNoInput",
                    FatalException,
                    "MCPL mode is enabled, but no MCPL input file was configured. "
                    "Use /particle_generator/set_mcpl_file before /run/beamOn, "
                    "or run the executable with -g/--gun for particle-gun mode.");
        return false;
    }

    const std::filesystem::path inputPath(g_mcplInputFile.c_str());
    if (!std::filesystem::exists(inputPath)) {
        std::ostringstream cmt;
        cmt << "Configured MCPL input file does not exist: "
            << g_mcplInputFile
            << "\nUse /particle_generator/set_mcpl_file with a valid path, "
               "or run particle-gun macros with -g/--gun.";
        G4Exception("G4MCPLGenerator::EnsureOpen()", "G4MCPLGeneratorMissingInput",
                    FatalException, cmt.str().c_str());
        return false;
    }

    m_mcplfile = mcpl_open_file(g_mcplInputFile.c_str());
    m_p = nullptr;
    FindNext();

    if (!m_p) {
        std::ostringstream cmt;
        cmt << "MCPL input file contains no usable particles: "
            << g_mcplInputFile;
        G4Exception("G4MCPLGenerator::EnsureOpen()", "G4MCPLGeneratorEmptyInput",
                    FatalException, cmt.str().c_str());
        return false;
    }

    return true;
}

// ============================================================================
// PDG Code Lookup with Cache
// ============================================================================

G4ParticleDefinition* G4MCPLGenerator::LookupPDG(G4int pdgcode) {
    if (m_currentPDG == pdgcode) {
        return m_currentPartDef;
    }

    m_currentPDG = pdgcode;

    auto it = m_pdg2pdef.find(pdgcode);
    if (it != m_pdg2pdef.end()) {
        m_currentPartDef = it->second;
    } else {
        m_currentPartDef = G4ParticleTable::GetParticleTable()->FindParticle(pdgcode);
        if (!m_currentPartDef && (pdgcode / 100000000 == 10)) {
            m_currentPartDef = G4IonTable::GetIonTable()->GetIon(pdgcode);
        }
        m_pdg2pdef[pdgcode] = m_currentPartDef;
    }

    return m_currentPartDef;
}
