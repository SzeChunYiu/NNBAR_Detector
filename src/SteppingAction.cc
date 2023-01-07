//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

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
//....
extern G4ThreadLocal G4int local_event_number;
extern G4ThreadLocal G4int local_event_number_MCPL;
extern G4int event_number;
extern std::ofstream pi0_outFile;
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
  G4int parentID = track->GetParentID();
  G4int ID = track->GetTrackID();
  G4String current_vol = step->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();
  G4String origin = track->GetOriginTouchable()->GetVolume()->GetName();
  G4String current_replica = step->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(0);
  G4int origin_replica = track->GetOriginTouchable()->GetReplicaNumber(0);
  G4String name = track->GetParticleDefinition()->GetParticleName();//PDGEncoding(); //(*secondary)[j]
  
  G4String proc = "primary";
  if (parentID>0){proc = track ->GetCreatorProcess()->GetProcessName();}


  std::vector<G4String> proc_name_ban_list = {"hIoni","eIoni","ionIoni","muIoni"}; 
  std::vector<G4String> vol_name_ban_list = {"TPC_PV","TPCPV_blocks"}; //"TPC_PV","B4C_PV"
  
  // if (parentID==0){
  //   G4StepPoint* PreStep = step->GetPreStepPoint();
  //   //Position
  //   G4ThreeVector pos1 = PreStep->GetPosition();

  //   G4int track_ID = track->GetTrackID(); 
  //   G4double KE = track->GetKineticEnergy();
  //   G4ThreeVector momentum = track->GetMomentumDirection();
  //   G4double mass = track->GetParticleDefinition()->GetPDGMass();
  //   G4String proc = "primary";
  //   G4double time = track->GetGlobalTime();
  //   G4ThreeVector vertex = track->GetVertexPosition();
  //   pi0_outFile << eventNumber << ","  << track_ID <<"," << parentID << "," << name  <<"," << proc << "," << current_vol << ","<< origin << ","  << mass << "," <<  KE << "," << std::setprecision(13)<< time <<  std::setprecision(5) << "," << pos1[0]/CLHEP::m << "," << pos1[1]/CLHEP::m << "," << pos1[2]/CLHEP::m << "," << momentum[0] << "," << momentum[1] << "," << momentum[2] <<G4endl;
  // }

  // else if (parentID==1 & proc == "neutronInelastic"){
  //     G4StepPoint* PreStep = step->GetPreStepPoint();
  //   //Position
  //   G4ThreeVector pos1 = PreStep->GetPosition();
  //   G4int track_ID = track->GetTrackID(); 
  //   G4double KE = track->GetKineticEnergy();
  //   G4ThreeVector momentum = track->GetMomentumDirection();
  //   G4double mass = track->GetParticleDefinition()->GetPDGMass();
  //   G4double time = track->GetGlobalTime();
  //   G4ThreeVector vertex = track->GetVertexPosition();
  //   pi0_outFile << eventNumber << ","  << track_ID <<"," << parentID << "," << name  <<"," << proc << "," << current_vol << ","<< origin << ","  << mass << "," <<  KE << "," << std::setprecision(13)<< time <<  std::setprecision(5) << "," << pos1[0]/CLHEP::m << "," << pos1[1]/CLHEP::m << "," << pos1[2]/CLHEP::m << "," << momentum[0] << "," << momentum[1] << "," << momentum[2] <<G4endl;
  // }

  //else if
  if(parentID>0 & step->IsFirstStepInVolume() & current_vol == origin & name!="opticalphoton"& current_replica == origin_replica){ //
    
    G4int track_ID = track->GetTrackID(); 
    G4double KE = track->GetKineticEnergy();
    G4ThreeVector momentum = track->GetMomentumDirection();
    G4double mass = track->GetParticleDefinition()->GetPDGMass();
    G4String proc = track ->GetCreatorProcess()->GetProcessName();
    G4double time = track->GetGlobalTime();
    G4ThreeVector vertex = track->GetVertexPosition();
    G4double weight = track->GetWeight ();

    //std::find(vol_name_ban_list.begin(), vol_name_ban_list.end(), current_vol) != vol_name_ban_list.end() & 
    if (std::find(vol_name_ban_list.begin(), vol_name_ban_list.end(), current_vol) != vol_name_ban_list.end()& 
        std::find(proc_name_ban_list.begin(), proc_name_ban_list.end(), proc) != proc_name_ban_list.end()){}
    else{pi0_outFile << eventNumber << ","  << track_ID <<"," << parentID << "," << name  <<"," << proc << "," << current_vol << ","<< origin << ","  << mass << "," <<  KE << "," << std::setprecision(13)<< time <<  std::setprecision(5) << "," << vertex[0]/CLHEP::m << "," << vertex[1]/CLHEP::m << "," << vertex[2]/CLHEP::m << "," << momentum[0] << "," << momentum[1] << "," << momentum[2] << "," << std::setprecision(10) << weight <<G4endl;}
  }

  // if(parentID>0 & step->IsFirstStepInVolume() & current_vol == origin & name=="opticalphoton"& current_replica == origin_replica){

  //   G4int module_ID = track->GetOriginTouchable()->GetReplicaNumber(2);
  //   G4int layer_ID = track->GetOriginTouchable()->GetReplicaNumber(1);
  //   G4int stave_ID = track->GetOriginTouchable()->GetReplicaNumber(0);
  //   G4int track_ID = track->GetTrackID(); 
  //   G4double KE = track->GetKineticEnergy();
  //   G4ThreeVector momentum = track->GetMomentumDirection();
  //   G4double mass = track->GetParticleDefinition()->GetPDGMass();
  //   G4String proc = track ->GetCreatorProcess()->GetProcessName();
  //   G4double time = track->GetGlobalTime();
  //   G4ThreeVector vertex = track->GetVertexPosition();
  //   G4double weight = track->GetWeight ();

  //   calorimeter_photon_outFile << eventNumber <<","<<track_ID<<","<<parentID<<","<< name << "," << proc<<","
  //                   << current_vol << "," << origin << ","
  //                   << module_ID<<","<<layer_ID<< "," << stave_ID << ","
  //                   << std::setprecision(13) << time << std::setprecision(5) <<"," << KE/CLHEP::eV << "," 
  //                   << vertex[0]/CLHEP::m << "," << vertex[1]/CLHEP::m << "," << vertex[2]/CLHEP::m << "," <<  std::setprecision(10) << weight << G4endl;
  // }







}
//}
//} 