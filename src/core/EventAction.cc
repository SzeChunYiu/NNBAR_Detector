// ============================================================================
// EventAction.cc
// Event-level action for NNBAR detector simulation
// Collects track/hit data and sends to EventDisplay
// ============================================================================

#include "core/EventAction.hh"
#include "core/SteppingAction.hh"
#include "config.h"

#include <iomanip>

#include "G4Event.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"

#if WITH_DASHBOARD
#include "gui/DashboardWindow.hh"
#include "gui/EventDisplay.hh"
#endif

#if WITH_GARFIELD_GPU
#include "physics/TPCDriftManager.hh"
#endif

#if WITH_CELERITAS
#include "physics/CeleritasCalorimeter.hh"
#endif

#include "output/ParquetOutputManager.hh"

// ============================================================================
// Constructor
// ============================================================================

EventAction::EventAction(SteppingAction* steppingAction)
    : G4UserEventAction(), fSteppingAction(steppingAction) {}

// ============================================================================
// Destructor
// ============================================================================

EventAction::~EventAction() {}

// ============================================================================
// Begin of Event Action
// ============================================================================

void EventAction::BeginOfEventAction(const G4Event* /*event*/) {
#if WITH_DASHBOARD
    // Clear previous event data from display
    auto& dashboard = nnbar::DashboardWindow::Instance();
    if (dashboard.GetEventDisplay()) {
        dashboard.GetEventDisplay()->ClearEvent();
    }
#endif

#if WITH_GARFIELD_GPU
    // Clear TPC drift data for new event
    nnbar::TPCDriftManager::Instance()->ClearEvent();
#endif

#if WITH_CELERITAS
    // Record current cumulative GPU energy so we can calculate per-event later
    auto& celCalo = nnbar::CeleritasCalorimeter::Instance();
    if (celCalo.IsInitialized()) {
        celCalo.BeginEvent();
    }
#endif
}

// ============================================================================
// End of Event Action
// ============================================================================

void EventAction::EndOfEventAction(const G4Event* event) {
    // Print progress periodically
    int eventID = event->GetEventID();
    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();

    if (printModulo > 0 && eventID % printModulo == 0) {
#if DEBUG_VERBOSE
        G4cout << "---> End of event: " << eventID << G4endl;
#endif
    }

    // ======== GPU/Accelerator Status Printout ========
    G4cout << "\n╔══════════════════════════════════════════════════════════════════╗" << G4endl;
    G4cout << "║  EVENT " << std::setw(6) << eventID << " SUMMARY"
           << std::setw(40) << " ║" << G4endl;
    G4cout << "╠══════════════════════════════════════════════════════════════════╣" << G4endl;

#if WITH_GARFIELD_GPU
    // Process TPC electron drift using GPU/OpenMP
    auto* driftMgr = nnbar::TPCDriftManager::Instance();
    int nClusters = driftMgr->GetNClusters();
    int nInputElectrons = driftMgr->GetTotalInputElectrons();

    G4cout << "║  [TPC Drift - GarfieldGPU]" << std::setw(42) << " ║" << G4endl;
    G4cout << "║    Device: " << std::left << std::setw(54)
           << driftMgr->GetDeviceName() << " ║" << G4endl;
    G4cout << std::right;

    if (driftMgr->IsEnabled() && nClusters > 0) {
        G4cout << "║    Ionization clusters: " << std::setw(8) << nClusters
               << std::setw(33) << " ║" << G4endl;
        G4cout << "║    Input electrons:     " << std::setw(8) << nInputElectrons
               << std::setw(33) << " ║" << G4endl;

        driftMgr->ProcessEvent();

        int nCollected = driftMgr->GetTotalCollected();
        float efficiency = nInputElectrons > 0 ?
            100.0f * nCollected / nInputElectrons : 0.0f;
        float kernelTime = driftMgr->GetKernelTimeMs();

        G4cout << "║    Collected electrons: " << std::setw(8) << nCollected
               << " (" << std::fixed << std::setprecision(1) << efficiency << "%)"
               << std::setw(20) << " ║" << G4endl;
        G4cout << "║    Drift time:          " << std::setw(8) << std::setprecision(2)
               << kernelTime << " ms" << std::setw(28) << " ║" << G4endl;
    } else {
        G4cout << "║    Status: " << std::left << std::setw(54)
               << (driftMgr->IsEnabled() ? "No TPC ionization this event" : "DISABLED")
               << " ║" << G4endl;
        G4cout << std::right;
    }
#else
    G4cout << "║  [TPC Drift] NOT COMPILED (WITH_GARFIELD_GPU=OFF)"
           << std::setw(15) << " ║" << G4endl;
#endif

#if WITH_CELERITAS
    G4cout << "╠──────────────────────────────────────────────────────────────────╣" << G4endl;
    G4cout << "║  [EM Physics - Celeritas GPU]" << std::setw(38) << " ║" << G4endl;
    G4cout << "║    Status: " << std::left << std::setw(54)
           << "ENABLED - e-/e+/gamma offloaded to GPU" << " ║" << G4endl;
    G4cout << std::right;

    // Get PER-EVENT energy deposits from GPU via GeantSimpleCalo
    auto& celCalo = nnbar::CeleritasCalorimeter::Instance();
    fCurrentGPUEnergy = nnbar::GPUEnergyData();  // Reset for this event
    if (celCalo.IsInitialized() && celCalo.IsCeleritasActive()) {
        fCurrentGPUEnergy = celCalo.GetPerEventEnergy();
        if (fCurrentGPUEnergy.isValid) {
            G4cout << "║    " << std::left << std::setw(18) << "LeadGlass Edep:"
                   << std::right << std::fixed << std::setprecision(2) << std::setw(10)
                   << fCurrentGPUEnergy.leadGlassEnergy / MeV << " MeV" << std::setw(25) << " ║" << G4endl;
            G4cout << "║    " << std::left << std::setw(18) << "Scintillator Edep:"
                   << std::right << std::fixed << std::setprecision(2) << std::setw(10)
                   << fCurrentGPUEnergy.scintillatorEnergy / MeV << " MeV" << std::setw(25) << " ║" << G4endl;
            G4cout << "║    " << std::left << std::setw(18) << "Total GPU Edep:"
                   << std::right << std::fixed << std::setprecision(2) << std::setw(10)
                   << fCurrentGPUEnergy.totalGPUEnergy / MeV << " MeV" << std::setw(25) << " ║" << G4endl;

            // Write GPU energy to Parquet output (skip zero-energy records)
            if (fCurrentGPUEnergy.totalGPUEnergy > 0.0) {
                auto& parquetMgr = nnbar::ParquetOutputManager::Instance();
                if (parquetMgr.IsInitialized()) {
                    nnbar::GPUEnergyRecord record;
                    record.event_id = eventID;
                    record.lead_glass_energy = fCurrentGPUEnergy.leadGlassEnergy / MeV;
                    record.scintillator_energy = fCurrentGPUEnergy.scintillatorEnergy / MeV;
                    record.total_gpu_energy = fCurrentGPUEnergy.totalGPUEnergy / MeV;
                    parquetMgr.WriteGPUEnergy(record);
                }
            }
        }
    }
#else
    G4cout << "╠──────────────────────────────────────────────────────────────────╣" << G4endl;
    G4cout << "║  [EM Physics] Geant4 CPU (Celeritas not compiled)"
           << std::setw(16) << " ║" << G4endl;
#endif

    G4cout << "╚══════════════════════════════════════════════════════════════════╝\n" << G4endl;

    // Send event data to the display
    SendEventToDisplay(event);
}

// ============================================================================
// Send Event Data to Dashboard Display
// ============================================================================

void EventAction::SendEventToDisplay(const G4Event* event) {
#if WITH_DASHBOARD
    auto& dashboard = nnbar::DashboardWindow::Instance();
    nnbar::EventDisplay* display = dashboard.GetEventDisplay();
    if (!display) return;

    // Get trajectory container
    G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
    if (!trajectoryContainer) return;

    int runId = G4RunManager::GetRunManager()->GetCurrentRun()->GetRunID();
    int eventId = event->GetEventID();

    // Create event display data
    nnbar::EventDisplayData eventData;
    eventData.runId = runId;
    eventData.eventId = eventId;
    eventData.totalEnergy = 0.0;

    // Process trajectories
    size_t nTraj = trajectoryContainer->entries();
    for (size_t i = 0; i < nTraj; i++) {
        G4VTrajectory* vTraj = (*trajectoryContainer)[i];
        G4Trajectory* traj = dynamic_cast<G4Trajectory*>(vTraj);
        if (!traj) continue;

        nnbar::Track track;
        track.trackId = traj->GetTrackID();
        track.parentId = traj->GetParentID();
        track.pdg = traj->GetPDGEncoding();
        track.initialEnergy = traj->GetInitialKineticEnergy() / MeV;
        track.particleName = traj->GetParticleName();
        track.creatorProcess = "";  // Not easily available from trajectory

        // Skip optical photons (PDG code 0 or -22 typically used)
        if (track.pdg == 0 || track.pdg == -22) continue;

        // Skip very low energy tracks
        if (track.initialEnergy < 1.0) continue;  // Skip tracks < 1 MeV

        eventData.totalEnergy += track.initialEnergy;

        // Extract trajectory points
        int nPoints = traj->GetPointEntries();
        for (int j = 0; j < nPoints; j++) {
            G4VTrajectoryPoint* point = traj->GetPoint(j);
            if (!point) continue;

            G4ThreeVector pos = point->GetPosition();
            nnbar::TrackPoint pt;
            pt.x = pos.x() / cm;  // Convert to cm
            pt.y = pos.y() / cm;
            pt.z = pos.z() / cm;
            pt.t = 0;  // Time not available from trajectory point
            pt.ke = track.initialEnergy;  // Approximate

            track.points.push_back(pt);
        }

        // Only add tracks with at least 2 points
        if (track.points.size() >= 2) {
            eventData.tracks.push_back(track);

            // Count particle statistics
            if (track.parentId == 0) {
                eventData.nPrimaries++;
            } else {
                eventData.nSecondaries++;
            }

            int absPdg = std::abs(track.pdg);
            switch (absPdg) {
                case 22: eventData.nGammas++; break;
                case 11: eventData.nElectrons++; break;
                case 2112: eventData.nNeutrons++; break;
                case 2212: eventData.nProtons++; break;
                case 211: case 111: eventData.nPions++; break;
                case 13: eventData.nMuons++; break;
                default: eventData.nOther++; break;
            }
        }
    }

    // Get energy deposition data from SteppingAction
    if (fSteppingAction) {
        eventData.edepTPC = fSteppingAction->GetEdepTPC();
        eventData.edepScintillator = fSteppingAction->GetEdepScintillator();
        eventData.edepLeadGlass = fSteppingAction->GetEdepLeadGlass();
        eventData.edepOther = fSteppingAction->GetEdepOther();
        eventData.photonsTPC = fSteppingAction->GetPhotonsTPC();
        eventData.photonsScintillator = fSteppingAction->GetPhotonsScintillator();
        eventData.photonsLeadGlass = fSteppingAction->GetPhotonsLeadGlass();
        eventData.scintPhotonCount = fSteppingAction->GetScintPhotonCount();
        eventData.cerenkovPhotonCount = fSteppingAction->GetCerenkovPhotonCount();

#if WITH_CELERITAS
        // Add GPU energy deposits from Celeritas GeantSimpleCalo
        // Use the per-event energy stored earlier in EndOfEventAction
        if (fCurrentGPUEnergy.isValid) {
            eventData.edepLeadGlass += fCurrentGPUEnergy.leadGlassEnergy / MeV;
            eventData.edepScintillator += fCurrentGPUEnergy.scintillatorEnergy / MeV;
        }
#endif

        eventData.totalEnergy = eventData.edepTPC + eventData.edepScintillator +
                                eventData.edepLeadGlass + eventData.edepOther;
    }

    // Update event display (thread-safe)
    display->SetEventData(eventData);

    // Update dashboard statistics
    nnbar::EventStatistics stats;
    stats.event_id = eventId;
    stats.n_primaries = eventData.nPrimaries;
    stats.n_secondaries = eventData.nSecondaries;
    stats.total_edep = eventData.totalEnergy;
    stats.n_tpc_hits = eventData.photonsTPC;
    stats.n_scintillator_hits = eventData.photonsScintillator;
    stats.n_leadglass_hits = eventData.photonsLeadGlass;
    stats.n_optical_photons = eventData.scintPhotonCount + eventData.cerenkovPhotonCount;
    dashboard.UpdateEventStats(stats);
    dashboard.UpdateEventCount(eventId + 1);

#if DEBUG_VERBOSE
    G4cout << "EventDisplay: Event " << eventId
           << " - " << eventData.tracks.size() << " tracks"
           << ", E = " << eventData.totalEnergy << " MeV" << G4endl;
#endif

#else
    (void)event;
#endif // WITH_DASHBOARD
}
