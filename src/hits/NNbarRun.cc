// ============================================================================
// NNbarRun.cc
// Custom G4Run class for recording detector hits to Parquet output
// ============================================================================

#include "hits/NNbarRun.hh"
#include "hits/NNbarHit.hh"
#include "output/ParquetOutputManager.hh"

#include "G4VHitsCollection.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Run.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include <algorithm>
#include <vector>

#include "config.h"

namespace {
    G4Mutex RunActionMutex = G4MUTEX_INITIALIZER;
}

// ============================================================================
// Constructor
// ============================================================================

NNbarRun::NNbarRun() : G4Run(),
    fCarbonHCID(-1),
    fSiliconHCID(-1),
    fBeampipeHCID(-1),
    fTPCHCID(-1),
    fScintHCID(-1),
    fLeadGlassHCID(-1),
    fPMTHCID(-1)
{
    // Get hit collection IDs
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();

    fCarbonHCID = sdManager->GetCollectionID("CarbonHitCollection");
    fSiliconHCID = sdManager->GetCollectionID("SiliconHitCollection");
    fBeampipeHCID = sdManager->GetCollectionID("TubeHitCollection");
    fTPCHCID = sdManager->GetCollectionID("TPCHitCollection");
    fScintHCID = sdManager->GetCollectionID("ScintillatorHitCollection");
    fLeadGlassHCID = sdManager->GetCollectionID("LeadGlassHitCollection");
    fPMTHCID = sdManager->GetCollectionID("PMTHitCollection");
}

// ============================================================================
// Destructor
// ============================================================================

NNbarRun::~NNbarRun() {}

// ============================================================================
// Record Event - Main hit recording logic
// ============================================================================

void NNbarRun::RecordEvent(const G4Event* event) {
    if (fCarbonHCID < 0) {
        return;
    }

    G4HCofThisEvent* hce = event->GetHCofThisEvent();
    if (!hce) {
        return;
    }

    // Lock for thread safety
    G4AutoLock lock(&RunActionMutex);

    auto& output = nnbar::ParquetOutputManager::Instance();
    const int eventId = event->GetEventID();

    // Process Carbon hits
    ProcessCarbonHits(hce, eventId, output);

    // Process Beampipe hits
    ProcessBeampipeHits(hce, eventId, output);

    // Process Silicon hits
    ProcessSiliconHits(hce, eventId, output);

    // Process TPC hits
    ProcessTPCHits(hce, eventId, output);

    // Process Scintillator hits
    ProcessScintillatorHits(hce, eventId, output);

    // Process LeadGlass hits
    ProcessLeadGlassHits(hce, eventId, output);

    // Process PMT hits
    ProcessPMTHits(hce, eventId, output);
}

// ============================================================================
// Process Carbon Target Hits
// ============================================================================

void NNbarRun::ProcessCarbonHits(G4HCofThisEvent* hce, int eventId,
                                  nnbar::ParquetOutputManager& output) {
    auto* hits = static_cast<NNbarHitsCollection*>(hce->GetHC(fCarbonHCID));
    if (!hits) return;

    for (size_t h = 0; h < hits->entries(); ++h) {
        auto* hit = (*hits)[h];

        nnbar::CarbonRecord rec;
        rec.event_id = eventId;
        rec.track_id = hit->GetTrackID();
        rec.parent_id = hit->GetParentID();
        rec.name = hit->GetName();
        rec.process = hit->GetProcess();
        rec.step_info = hit->GetStepInfo();
        rec.origin = hit->GetOriginVolName();
        rec.x = hit->GetPosX() / CLHEP::cm;
        rec.y = hit->GetPosY() / CLHEP::cm;
        rec.z = hit->GetPosZ() / CLHEP::cm;
        rec.px = hit->GetPX();
        rec.py = hit->GetPY();
        rec.pz = hit->GetPZ();
        rec.t = hit->GetTime() / CLHEP::ns;
        rec.ke = hit->GetKinEn() / CLHEP::MeV;
        rec.edep = hit->GetEdep() / CLHEP::MeV;

        output.WriteCarbon(rec);
    }
}

// ============================================================================
// Process Beampipe Hits
// ============================================================================

void NNbarRun::ProcessBeampipeHits(G4HCofThisEvent* hce, int eventId,
                                    nnbar::ParquetOutputManager& output) {
    auto* hits = static_cast<NNbarHitsCollection*>(hce->GetHC(fBeampipeHCID));
    if (!hits) return;

    for (size_t h = 0; h < hits->entries(); ++h) {
        auto* hit = (*hits)[h];

        nnbar::BeampipeRecord rec;
        rec.event_id = eventId;
        rec.track_id = hit->GetTrackID();
        rec.parent_id = hit->GetParentID();
        rec.name = hit->GetName();
        rec.process = hit->GetProcess();
        rec.step_info = hit->GetStepInfo();
        rec.current_vol = hit->GetVolName();
        rec.origin = hit->GetOriginVolName();
        rec.x = hit->GetPosX() / CLHEP::cm;
        rec.y = hit->GetPosY() / CLHEP::cm;
        rec.z = hit->GetPosZ() / CLHEP::cm;
        rec.px = hit->GetPX();
        rec.py = hit->GetPY();
        rec.pz = hit->GetPZ();
        rec.t = hit->GetTime() / CLHEP::ns;
        rec.ke = hit->GetKinEn() / CLHEP::MeV;
        rec.edep = hit->GetEdep() / CLHEP::MeV;

        output.WriteBeampipe(rec);
    }
}

// ============================================================================
// Process Silicon Tracker Hits
// ============================================================================

void NNbarRun::ProcessSiliconHits(G4HCofThisEvent* hce, int eventId,
                                   nnbar::ParquetOutputManager& output) {
    auto* hits = static_cast<NNbarHitsCollection*>(hce->GetHC(fSiliconHCID));
    if (!hits) return;

    for (size_t h = 0; h < hits->entries(); ++h) {
        auto* hit = (*hits)[h];

        // Only record hits with energy deposit
        G4double edep = hit->GetEdep();
        if (edep <= 0.0) continue;

        nnbar::SiliconRecord rec;
        rec.event_id = eventId;
        rec.track_id = hit->GetTrackID();
        rec.parent_id = hit->GetParentID();
        rec.name = hit->GetName();
        rec.process = hit->GetProcess();
        rec.step_info = hit->GetStepInfo();
        rec.origin = hit->GetOriginVolName();
        rec.layer_id = hit->GetXID();
        rec.x = hit->GetPosX() / CLHEP::cm;
        rec.y = hit->GetPosY() / CLHEP::cm;
        rec.z = hit->GetPosZ() / CLHEP::cm;
        rec.px = hit->GetPX();
        rec.py = hit->GetPY();
        rec.pz = hit->GetPZ();
        rec.t = hit->GetTime() / CLHEP::ns;
        rec.ke = hit->GetKinEn() / CLHEP::MeV;
        rec.edep = edep / CLHEP::MeV;

        output.WriteSilicon(rec);
    }
}

// ============================================================================
// Process TPC Hits
// ============================================================================

void NNbarRun::ProcessTPCHits(G4HCofThisEvent* hce, int eventId,
                               nnbar::ParquetOutputManager& output) {
    auto* hits = static_cast<NNbarHitsCollection*>(hce->GetHC(fTPCHCID));
    if (!hits) return;

    for (size_t h = 0; h < hits->entries(); ++h) {
        auto* hit = (*hits)[h];

        // Skip ionization-only processes
        G4String proc = hit->GetProcess();
        if (proc == "eIoni" || proc == "hIoni" || proc == "muIoni") {
            continue;
        }

        // Only record hits with electrons
        G4double electrons = hit->GetPhotons();  // Reused for electron count
        if (electrons <= 0) continue;

        nnbar::TPCRecord rec;
        rec.event_id = eventId;
        rec.track_id = hit->GetTrackID();
        rec.parent_id = hit->GetParentID();
        rec.name = hit->GetName();
        rec.process = proc;
        rec.step_info = hit->GetStepInfo();
        rec.origin = hit->GetOriginVolName();
        rec.current_vol = hit->GetVolName();
        rec.module_id = hit->GetMod_ID();
        rec.layer_id = hit->GetXID();
        rec.x = hit->GetPosX() / CLHEP::cm;
        rec.y = hit->GetPosY() / CLHEP::cm;
        rec.z = hit->GetPosZ() / CLHEP::cm;
        rec.px = hit->GetPX();
        rec.py = hit->GetPY();
        rec.pz = hit->GetPZ();
        rec.t = hit->GetTime() / CLHEP::ns;
        rec.ke = hit->GetKinEn() / CLHEP::MeV;
        rec.edep = hit->GetEdep() / CLHEP::MeV;
        rec.track_length = hit->GetTrackLength() / CLHEP::cm;
        rec.electrons = electrons;

        output.WriteTPC(rec);
    }
}

// ============================================================================
// Process Scintillator Hits
// ============================================================================

void NNbarRun::ProcessScintillatorHits(G4HCofThisEvent* hce, int eventId,
                                        nnbar::ParquetOutputManager& output) {
    auto* hits = static_cast<NNbarHitsCollection*>(hce->GetHC(fScintHCID));
    if (!hits) return;

    for (size_t h = 0; h < hits->entries(); ++h) {
        auto* hit = (*hits)[h];

        // Skip optical photons
        if (hit->GetName() == "opticalphoton") continue;

        // Only record hits with scintillation photons
        G4int photons = hit->GetPhotons();
        if (photons <= 0) continue;

        nnbar::ScintillatorRecord rec;
        rec.event_id = eventId;
        rec.track_id = hit->GetTrackID();
        rec.parent_id = hit->GetParentID();
        rec.name = hit->GetName();
        rec.process = hit->GetProcess();
        rec.step_info = hit->GetStepInfo();
        rec.origin = hit->GetOriginVolName();
        rec.volume = hit->GetVolName();
        rec.module_id = hit->GetMod_ID();
        rec.layer_id = hit->GetXID();
        rec.stave_id = hit->GetStave_ID();
        rec.x = hit->GetPosX() / CLHEP::cm;
        rec.y = hit->GetPosY() / CLHEP::cm;
        rec.z = hit->GetPosZ() / CLHEP::cm;
        rec.particle_x = hit->GetPosX_particle() / CLHEP::cm;
        rec.particle_y = hit->GetPosY_particle() / CLHEP::cm;
        rec.particle_z = hit->GetPosZ_particle() / CLHEP::cm;
        rec.x_local = hit->GetLocalPosX() / CLHEP::cm;
        rec.y_local = hit->GetLocalPosY() / CLHEP::cm;
        rec.z_local = hit->GetLocalPosZ() / CLHEP::cm;
        rec.t = hit->GetTime() / CLHEP::ns;
        rec.ke = hit->GetKinEn() / CLHEP::MeV;
        rec.edep = hit->GetEdep() / CLHEP::MeV;
        rec.photons = photons;

        output.WriteScintillator(rec);
    }
}

// ============================================================================
// Process Lead Glass Hits
// ============================================================================

void NNbarRun::ProcessLeadGlassHits(G4HCofThisEvent* hce, int eventId,
                                     nnbar::ParquetOutputManager& output) {
    auto* hits = static_cast<NNbarHitsCollection*>(hce->GetHC(fLeadGlassHCID));
    if (!hits) return;

    for (size_t h = 0; h < hits->entries(); ++h) {
        auto* hit = (*hits)[h];

        // Skip optical photons
        if (hit->GetName() == "opticalphoton") continue;

        // Record all hits (previously filtered by Cherenkov photons)
        G4int photons = hit->GetPhotons();

        nnbar::LeadGlassRecord rec;
        rec.event_id = eventId;
        rec.track_id = hit->GetTrackID();
        rec.parent_id = hit->GetParentID();
        rec.name = hit->GetName();
        rec.process = hit->GetProcess();
        rec.step_info = hit->GetStepInfo();
        rec.origin = hit->GetOriginVolName();
        rec.module_id = hit->GetXID();
        rec.x = hit->GetPosX() / CLHEP::cm;
        rec.y = hit->GetPosY() / CLHEP::cm;
        rec.z = hit->GetPosZ() / CLHEP::cm;
        rec.t = hit->GetTime() / CLHEP::ns;
        rec.ke = hit->GetKinEn() / CLHEP::MeV;
        rec.edep = hit->GetEdep() / CLHEP::MeV;
        rec.photons = photons;

        output.WriteLeadGlass(rec);
    }
}

// ============================================================================
// Process PMT Hits
// ============================================================================

void NNbarRun::ProcessPMTHits(G4HCofThisEvent* hce, int eventId,
                               nnbar::ParquetOutputManager& output) {
    auto* hits = static_cast<NNbarHitsCollection*>(hce->GetHC(fPMTHCID));
    if (!hits) return;

    // Collect photon data per PMT module
    struct PMTData {
        G4double x = 0, y = 0, z = 0;
        std::vector<G4double> times;
        std::vector<G4double> energies;
    };
    std::map<G4int, PMTData> pmtMap;

    for (size_t h = 0; h < hits->entries(); ++h) {
        auto* hit = (*hits)[h];

        // Only store optical photon hits
        if (hit->GetName() != "opticalphoton") continue;

        G4int moduleId = hit->GetXID();
        auto& pmt = pmtMap[moduleId];

        pmt.x = hit->GetPosX() / CLHEP::cm;
        pmt.y = hit->GetPosY() / CLHEP::cm;
        pmt.z = hit->GetPosZ() / CLHEP::cm;
        pmt.times.push_back(hit->GetTime() / CLHEP::ns);
        pmt.energies.push_back(hit->GetKinEn() / CLHEP::MeV);
    }

    // Write aggregated PMT data
    for (const auto& [moduleId, pmt] : pmtMap) {
        nnbar::PMTRecord rec;
        rec.event_id = eventId;
        rec.module_id = moduleId;
        rec.x = pmt.x;
        rec.y = pmt.y;
        rec.z = pmt.z;
        rec.photons = static_cast<int32_t>(pmt.times.size());
        rec.ke = pmt.energies;
        rec.t = pmt.times;

        output.WritePMT(rec);
    }
}

// ============================================================================
// Merge Runs (for multi-threaded execution)
// ============================================================================

void NNbarRun::Merge(const G4Run* aRun) {
    G4Run::Merge(aRun);
}
