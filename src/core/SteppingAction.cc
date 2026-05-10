// ============================================================================
// SteppingAction - Records particle interactions at each simulation step
// Outputs interaction data to Parquet format for analysis
// ============================================================================

#include "core/SteppingAction.hh"
#include "core/EventAction.hh"
#include "core/Analysis.hh"
#include "output/ParquetOutputManager.hh"
#include "util/GeometryManager.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4DynamicParticle.hh"
#include "G4Threading.hh"

namespace {G4Mutex SteppingMutex = G4MUTEX_INITIALIZER;}

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;
  fEdepTPC = 0.0;
  fEdepScintillator = 0.0;
  fEdepLeadGlass = 0.0;
  fEdepOther = 0.0;
  fPhotonsTPC = 0;
  fPhotonsScintillator = 0;
  fPhotonsLeadGlass = 0;
}

void SteppingAction::ResetEdep() {
  fEdepTPC = 0.0;
  fEdepScintillator = 0.0;
  fEdepLeadGlass = 0.0;
  fEdepOther = 0.0;
  fPhotonsTPC = 0;
  fPhotonsScintillator = 0;
  fPhotonsLeadGlass = 0;
  fScintillationCounter = 0;
  fCerenkovCounter = 0;
  // Clear per-volume energy data in GeometryManager
  nnbar::GeometryManager::Instance().ClearEventData();
}

//....

SteppingAction::~SteppingAction()
{ 
}

//....

void SteppingAction::UserSteppingAction(const G4Step* step)
{

  G4AutoLock lock(&SteppingMutex);
  G4int eventNumber = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();

  // Reset counters at start of new event
  if (eventNumber != fEventNumber) {
    fEventNumber = eventNumber;
    ResetEdep();
  }

  G4Track* track = step->GetTrack();
  G4String name = track->GetParticleDefinition()->GetParticleName();
  G4String current_vol = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();

  // Accumulate energy deposition by detector type AND per-volume
  G4double edep = step->GetTotalEnergyDeposit();
  if (edep > 0) {
    G4double edepMeV = edep / CLHEP::MeV;
    auto& geoMgr = nnbar::GeometryManager::Instance();

    // Determine detector type from volume name
    if (current_vol.find("TPC") != std::string::npos ||
        current_vol.find("tpc") != std::string::npos) {
      fEdepTPC += edepMeV;
    } else if (current_vol.find("Scint") != std::string::npos ||
               current_vol.find("scint") != std::string::npos) {
      fEdepScintillator += edepMeV;
      // Track per-bar energy in scintillator
      // Get module copy number (parent volume) and bar copy number
      auto* touchable = step->GetPreStepPoint()->GetTouchable();
      G4int barCopyNo = touchable->GetReplicaNumber(0);  // Bar level
      G4int layerCopyNo = touchable->GetReplicaNumber(1); // Layer level
      G4int moduleCopyNo = touchable->GetReplicaNumber(2); // Module level
      (void)layerCopyNo; // Layer not used for now
      geoMgr.AddScintillatorEnergyDeposit(moduleCopyNo, barCopyNo, edepMeV);
    } else if (current_vol.find("LeadGlass") != std::string::npos ||
               current_vol.find("leadglass") != std::string::npos ||
               current_vol.find("LG") != std::string::npos) {
      fEdepLeadGlass += edepMeV;
      // Track per-block energy in lead glass
      G4int copyNo = step->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(0);
      geoMgr.AddLeadGlassEnergyDeposit(copyNo, edepMeV);
    } else {
      fEdepOther += edepMeV;
    }
  }

  // Count optical photons by detector
  if (name == "opticalphoton") {
    // Count by production process
    if (track->GetParentID() > 0 && track->GetCreatorProcess()) {
      G4String procName = track->GetCreatorProcess()->GetProcessName();
      if (procName == "Scintillation") {
        fScintillationCounter++;
      } else if (procName == "Cerenkov") {
        fCerenkovCounter++;
      }
    }

    // Count by detector location
    if (current_vol.find("TPC") != std::string::npos) {
      fPhotonsTPC++;
    } else if (current_vol.find("Scint") != std::string::npos) {
      fPhotonsScintillator++;
    } else if (current_vol.find("LeadGlass") != std::string::npos ||
               current_vol.find("LG") != std::string::npos) {
      fPhotonsLeadGlass++;
    }
    return;  // Don't process further for optical photons
  }

  G4int parentID = track->GetParentID();
  G4String origin = track->GetOriginTouchable()->GetVolume()->GetName();
  G4int current_replica = step->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(0);
  G4int origin_replica = track->GetOriginTouchable()->GetReplicaNumber(0);

  G4double KE = track->GetKineticEnergy();

  G4String proc = "primary";
  if (parentID>0){proc = track ->GetCreatorProcess()->GetProcessName();}

  std::vector<G4String> proc_name_ban_list = {"hIoni","eIoni","ionIoni","muIoni"};
  std::vector<G4String> vol_name_ban_list = {"TPC_PV"}; //"TPC_PV","B4C_PV"


  //std::cout << " interaction output " << std::endl;

  // //else if
  // if it is not a primary particle， we are interested if the following criteria is met
  // 1. Parent ID > 0
  // 2. it is in its origin volume
  // 3. it is not optical photon
  // 4. it is in the origin replica
  // ###########################################
  // Need to get the following quantities:
  // Track ID , parent ID, name << proc << KE << time << xyz << volume << momentum dir

  if(parentID>0 && step->IsFirstStepInVolume() && current_vol == origin && KE > 1.0*CLHEP::MeV && current_replica == origin_replica){ //

      G4int track_ID = track->GetTrackID();
      G4ThreeVector momentum = track->GetMomentumDirection();
      G4double mass = track->GetParticleDefinition()->GetPDGMass();
      G4String proc = track ->GetCreatorProcess()->GetProcessName();
      G4double time = track->GetGlobalTime();
      G4ThreeVector vertex = track->GetVertexPosition();

      //std::cout << " interaction output " << std::endl;

      // Skip banned processes/volumes combinations
      bool inBannedVol = std::find(vol_name_ban_list.begin(), vol_name_ban_list.end(), current_vol) != vol_name_ban_list.end();
      bool isBannedProc = std::find(proc_name_ban_list.begin(), proc_name_ban_list.end(), proc) != proc_name_ban_list.end();

      if (!(inBannedVol && isBannedProc)) {
          // Write interaction record to Parquet
          nnbar::InteractionRecord rec;
          rec.event_id = eventNumber;
          rec.track_id = track_ID;
          rec.parent_id = parentID;
          rec.name = name;
          rec.process = proc;
          rec.current_vol = current_vol;
          rec.origin = origin;
          rec.mass = mass;
          rec.ke = KE;
          rec.t = time;
          rec.x = vertex[0] / CLHEP::cm;
          rec.y = vertex[1] / CLHEP::cm;
          rec.z = vertex[2] / CLHEP::cm;
          rec.px = momentum[0];
          rec.py = momentum[1];
          rec.pz = momentum[2];

          nnbar::ParquetOutputManager::Instance().WriteInteraction(rec);
      }
  }
}
