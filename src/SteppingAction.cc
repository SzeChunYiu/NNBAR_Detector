#include "SteppingAction.hh"
#include "EventAction.hh"
#include "Analysis.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4DynamicParticle.hh"
#include "G4Threading.hh"
#include <arrow/io/file.h>
#include <parquet/stream_writer.h>
#include "parquet_writer.h"


// //....
// extern G4ThreadLocal G4int local_event_number;
// extern G4ThreadLocal G4int local_event_number_MCPL;
// extern G4int event_number;

extern parquet::StreamWriter Interaction_os;
extern std::ofstream calorimeter_photon_outFile;

namespace {G4Mutex SteppingMutex = G4MUTEX_INITIALIZER;}

SteppingAction::SteppingAction()
: G4UserSteppingAction()
{ 
  fScintillationCounter = 0;
  fCerenkovCounter      = 0;
  fEventNumber = -1;
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

  G4Track* track = step->GetTrack();
  G4String name = track->GetParticleDefinition()->GetParticleName();

  if (name!="opticalphoton"){

    G4int parentID = track->GetParentID();
    G4int ID = track->GetTrackID();
    G4String current_vol = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();
    G4String origin = track->GetOriginTouchable()->GetVolume()->GetName();
    G4String current_replica = step->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(0);
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

    if(parentID>0 & step->IsFirstStepInVolume() & current_vol == origin & KE > 1.0*CLHEP::MeV & current_replica == origin_replica){ //
      
        G4int track_ID = track->GetTrackID(); 
        G4ThreeVector momentum = track->GetMomentumDirection();
        G4double mass = track->GetParticleDefinition()->GetPDGMass();
        G4String proc = track ->GetCreatorProcess()->GetProcessName();
        G4double time = track->GetGlobalTime();
        G4ThreeVector vertex = track->GetVertexPosition();
        G4double weight = track->GetWeight ();

        //std::cout << " interaction output " << std::endl; 

        //std::find(vol_name_ban_list.begin(), vol_name_ban_list.end(), current_vol) != vol_name_ban_list.end() & 
        if (std::find(vol_name_ban_list.begin(), vol_name_ban_list.end(), current_vol) != vol_name_ban_list.end()& 
            std::find(proc_name_ban_list.begin(), proc_name_ban_list.end(), proc) != proc_name_ban_list.end()){}
        else{
              Interaction_os << eventNumber << track_ID << parentID << name << proc 
                  << current_vol << origin << mass <<  KE 
                  << time << vertex[0]/CLHEP::cm << vertex[1]/CLHEP::cm << vertex[2]/CLHEP::cm 
                  << momentum[0] <<  momentum[1] <<  momentum[2]
                  << parquet::EndRow;
                  //Interaction_os.EndRowGroup();
        }
    }
  }
}
