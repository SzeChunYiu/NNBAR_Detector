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


#include "EventAction.hh"
#include "Analysis.hh"
#include "ScintillatorSD.hh"
#include "AbsorberSD.hh"
#include "TubeSD.hh"
#include "PMTSD.hh"
#include "TPCSD.hh"
#include "Scint_DetSD.hh"
#include "ShieldSD.hh"
#include "NNbarHit.hh"

#include "G4VHitsCollection.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include "config.h"
//....

extern G4double event_number; 
//extern G4double silicon_time_trigger;

extern std::vector<std::vector<G4double>> PMT_record;
extern std::vector<std::vector<G4double>> scint_record;
//extern std::ofstream PMT_outFile; 
extern std::ofstream Silicon_outFile;
extern std::ofstream Tube_outFile;
extern std::ofstream TPC_outFile;
extern std::ofstream Scint_layer_outFile; 
extern std::ofstream Abs_outFile;

#if VERSION_SHIELD==1
    extern std::ofstream Shield_outFile;
#endif

EventAction::EventAction(): 
    G4UserEventAction(),
    fAbsoEdepHCID(-1),
    fGapEdepHCID(-1),
    fAbsoTrackLengthHCID(-1),
    fCerenkovHCID(-1),
    fScintTrackLengthHCID(-1),
    scintHitsCollectionID(-1),
    absHitsCollectionID(-1),
    tubeHitsCollectionID(-1),
    TPCHitsCollectionID(-1),
    PMTHitsCollectionID(-1),
    ShieldHitsCollectionID(-1),
    SiliconHitsCollectionID(-1)
{}

//....

EventAction::~EventAction()
{}

//....

G4THitsMap<G4double>* 
EventAction::GetHitsCollection(G4int hcID,
                                  const G4Event* event) const
{
  auto hitsCollection 
    = static_cast<G4THitsMap<G4double>*>(
        event->GetHCofThisEvent()->GetHC(hcID));
  
  if ( ! hitsCollection ) {
    G4ExceptionDescription msg;
    msg << "Cannot access hitsCollection ID " << hcID; 
    G4Exception("EventAction::GetHitsCollection()",
      "MyCode0003", FatalException, msg);
  }         

  return hitsCollection;
}    

//....

G4double EventAction::GetSum(G4THitsMap<G4double>* hitsMap) const
{
  G4double sumValue = 0.;
  for ( auto it : *hitsMap->GetMap() ) {
    // hitsMap->GetMap() returns the map of std::map<G4int, G4double*>
    sumValue += *(it.second);
  }
  return sumValue;  
}  

//.....

void EventAction::BeginOfEventAction(const G4Event* /*event*/)
{
    G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
    if(scintHitsCollectionID == -1) {
       scintHitsCollectionID = pSDManager->GetCollectionID("ScintillatorHitCollection");
       absHitsCollectionID = pSDManager->GetCollectionID("AbsorberHitCollection");
       tubeHitsCollectionID = pSDManager->GetCollectionID("TubeHitCollection");
       TPCHitsCollectionID = pSDManager->GetCollectionID("TPCHitCollection");
       PMTHitsCollectionID = pSDManager->GetCollectionID("PMTHitCollection");
       SiliconHitsCollectionID = pSDManager->GetCollectionID("SiliconHitCollection");
       #if VERSION_SHIELD==1
           ShieldHitsCollectionID = pSDManager->GetCollectionID("ShieldHitCollection");
       #endif
  }
}

//....

void EventAction::EndOfEventAction(const G4Event* event)
{  
 
    if(scintHitsCollectionID  < 0) {return;}

	G4HCofThisEvent* HCE = event->GetHCofThisEvent();

    int CHCID = -1;
    if (CHCID<0) {CHCID = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitCollection");}
    NNbarHitsCollection* ScintHits = 0;

    int CHCID2 = -1;
    if (CHCID2<0) {CHCID2 = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitCollection");}
    NNbarHitsCollection* AbsHits   = 0;
    
    int CHCID3 = -1;
    if (CHCID3<0) {CHCID3 = G4SDManager::GetSDMpointer()->GetCollectionID("TubeHitCollection");}
    NNbarHitsCollection* TubeHits  = 0;
    
    int CHCID4 = -1;
    if (CHCID4<0) {CHCID4 = G4SDManager::GetSDMpointer()->GetCollectionID("TPCHitCollection");}
    NNbarHitsCollection* TPCHits  = 0;
    
    int CHCID5 = -1;
    if (CHCID5<0) {CHCID5 = G4SDManager::GetSDMpointer()->GetCollectionID("PMTHitCollection");}
    NNbarHitsCollection* PMTHits  = 0;
    
    int CHCID6 = -1;
    if (CHCID6<0) {CHCID6 = G4SDManager::GetSDMpointer()->GetCollectionID("SiliconHitCollection");}
    NNbarHitsCollection* SiliconHits  = 0;
    
    #if VERSION_SHIELD==1    
        int CHCID7 = -1;
        if (CHCID7<0) {CHCID7 = G4SDManager::GetSDMpointer()->GetCollectionID("ShieldHitCollection");}
        NNbarHitsCollection* ShieldHits  = 0;
    #endif

    if (HCE) {

        G4AnalysisManager* analysis = G4AnalysisManager::Instance();

        G4int b = 1;
        G4int ltime     = 0.;
        G4int parentID  = 0;
        G4String proc   = "";
        G4String name   = "";
        G4double time   = 0.;
        G4int trID      = 0;
        G4int i         = 0;
        G4double kinEn  = 0.;
        G4double eDep   = 0.;
        G4double trackl = 0.;	
        G4int hitCount  = 0;
        G4int org_replica = 99;
        G4int group_ID = 999;
        G4int module_ID = 999;
        
        // Book vector to keep track of Edep in each Scintillator Sheet
        G4double EdepPerSheet[10] = {0., 0., 0., 0., 0.,0., 0., 0., 0., 0.};
        G4int ScintPerSheet[10] = { 0,0,0,0,0,0,0,0,0,0 };
        G4int ScintPerSheet_new[10] = { 0,0,0,0,0,0,0,0,0,0 };
        G4double totEdep   = 0.;     
        G4double eDepScint = 0.;
        G4double eDepAbs   = 0.;
        G4double eDepTube  = 0.; 
        G4double extraEdep = 0.;
        G4double eDepCompt = 0.;
        G4double eDepInelastic= 0.;
        G4double eDephIoni = 0.;
        G4double eDepHadElas = 0.;
        G4double eDepPrimary = 0.;
        G4double eDepOther = 0.;
        G4int cerenkovCounter = 0;
        G4int scint_photons = 0;
        G4int scint_photons_check = 0;
        G4int PMT_photons = 0;



        ScintHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID));
        
        std::vector<G4int> Scint_Group_index;
        std::vector<G4int> Scint_Group_index_reduced;
        std::vector<G4int> Scint_Module_index;
        std::vector<G4int> Scint_Module_index_reduced;
        std::vector<G4int> Scint_Layer_index;
        std::vector<G4int> Scint_photon_array;
        std::vector<G4int> Scint_time_array;
        std::vector<G4int> event_ID_array;

        if (ScintHits) {
           hitCount = ScintHits->entries();
           for (G4int h=0; h<hitCount; h++) { 
               name     = ((*ScintHits)[h]) -> GetName();
               if (name != "opticalphoton"){
                    //ltime    = ((*ScintHits)[h]) -> GetLocalTime();
                    parentID = ((*ScintHits)[h]) -> GetParentID();
                    trID     = ((*ScintHits)[h]) -> GetTrackID();
                    proc     = ((*ScintHits)[h]) -> GetProcess();
                    x = ((*ScintHits)[h]) -> GetPosX();
                    y = ((*ScintHits)[h]) -> GetPosY();
                    z = ((*ScintHits)[h]) -> GetPosZ();
                    
                    G4double particle_x = ((*ScintHits)[h]) -> GetPosX_particle();
                    G4double particle_y = ((*ScintHits)[h]) -> GetPosY_particle();
                    G4double particle_z = ((*ScintHits)[h]) -> GetPosZ_particle();
                    
                    time     = ((*ScintHits)[h]) -> GetTime();

                    auto stave_ID = ((*ScintHits)[h]) -> GetStave_ID(); 
                    i        = ((*ScintHits)[h]) -> GetXID();
                    //group_ID = ((*ScintHits)[h]) -> GetGroup_ID();
                    module_ID = ((*ScintHits)[h]) -> GetMod_ID();

                    kinEn    = ((*ScintHits)[h]) -> GetKinEn();
                    eDep     = ((*ScintHits)[h]) -> GetEdep();
                    trackl   = ((*ScintHits)[h]) -> GetTrackLength();
                    G4int scint_photons_per_hit = ((*ScintHits)[h])->GetPhotons();
                    // writing output

                    if (scint_photons_per_hit>0){
                    Scint_layer_outFile << event_number<<","<<trID<<","<<parentID<<","<<name<<","<<proc<<","<< module_ID<<","<<i<< "," << stave_ID << ","
                    << time <<","<<kinEn<<","<<eDep<<","<< scint_photons_per_hit << "," << x <<","<< y << ","<< z << "," << particle_x << "," << particle_y << "," << particle_z  <<G4endl;
                    }
                    //std::cout<< group_ID<<","<<module_ID<<","<<i<< "," << stave_ID << std::endl;

               }

            }
        }
 
        AbsHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID2));
        if(AbsHits) {
            for (G4int h=0; h<AbsHits->entries(); h++) {
                name     = ((*AbsHits)[h]) -> GetName();
                if (name != "opticalphoton"){
                    ltime           = ((*AbsHits)[h]) -> GetLocalTime();
                    parentID	= ((*AbsHits)[h]) -> GetParentID();
                    proc            = ((*AbsHits)[h]) -> GetProcess();
                    G4String name   = ((*AbsHits)[h]) -> GetName();
                    G4double time   = ((*AbsHits)[h]) -> GetTime();
                    G4int trID      = ((*AbsHits)[h]) -> GetTrackID();
                    G4int i         = ((*AbsHits)[h]) -> GetXID();
                    G4double kinEn  = ((*AbsHits)[h]) -> GetKinEn();
                    G4double eDep   = ((*AbsHits)[h]) -> GetEdep();
                    G4double trackl = ((*AbsHits)[h]) -> GetTrackLength();
                    G4double photons_cerenkov = ((*AbsHits)[h])->GetPhotons();
                    x = ((*AbsHits)[h]) -> GetPosX();
                    y = ((*AbsHits)[h]) -> GetPosY();
                    z = ((*AbsHits)[h]) -> GetPosZ();

                    if (photons_cerenkov>0){
                    Abs_outFile << event_number<<","<<trID << ","<<parentID<<","<<name<<","<<proc<<","<<i<<","
                    << time<<","<<kinEn<<","<<eDep<<","<< trackl <<","<< photons_cerenkov <<  "," << x <<","<< y << ","<<z<<G4endl;
                    }
                }
                //cerenkovCounter = cerenkovCounter + photons_cerenkov;
            }
        }
        
        TubeHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID3));
        if (TubeHits) {
	    hitCount = TubeHits->entries();
            for (G4int h=0; h<hitCount; h++) {
                ltime    = ((*TubeHits)[h]) -> GetLocalTime();
                parentID = ((*TubeHits)[h]) -> GetParentID();
                proc     = ((*TubeHits)[h]) -> GetProcess();
                name     = ((*TubeHits)[h]) -> GetName();
                time     = ((*TubeHits)[h]) -> GetTime(); 
                trID     = ((*TubeHits)[h]) -> GetTrackID();
                kinEn    = ((*TubeHits)[h]) -> GetKinEn();
                eDep     = ((*TubeHits)[h]) -> GetEdep();
                trackl   = ((*TubeHits)[h]) -> GetTrackLength();	
                G4double x = ((*TubeHits)[h]) -> GetPosX();
                G4double y = ((*TubeHits)[h]) -> GetPosY();
                G4double z = ((*TubeHits)[h]) -> GetPosZ();
            
                Tube_outFile << event_number << "," << trID << "," << parentID << "," << name << "," << x <<","<<y<<","<<z<<","<<
                time<< "," <<kinEn<<","<<eDep << "," << trackl <<G4endl;
            
            }         
        }

        TPCHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID4));
        if (TPCHits) {
            hitCount = TPCHits->entries();
            for (G4int h=0; h<hitCount; h++) {
                
                ltime    = ((*TPCHits)[h]) -> GetLocalTime();
                parentID = ((*TPCHits)[h]) -> GetParentID();
                proc     = ((*TPCHits)[h]) -> GetProcess();
                name     = ((*TPCHits)[h]) -> GetName();
                time     = ((*TPCHits)[h]) -> GetTime();
                trID     = ((*TPCHits)[h]) -> GetTrackID();
                kinEn    = ((*TPCHits)[h]) -> GetKinEn();
                eDep     = ((*TPCHits)[h]) -> GetEdep();
                module_ID = ((*TPCHits)[h]) -> GetMod_ID();
                i = ((*TPCHits)[h]) -> GetXID();
                x = ((*TPCHits)[h]) -> GetPosX();
                y = ((*TPCHits)[h]) -> GetPosY();
                z = ((*TPCHits)[h]) -> GetPosZ();
                G4double electrons = ((*TPCHits)[h]) -> GetPhotons(); 
                G4double trackl = ((*TPCHits)[h]) -> GetTrackLength();

                //SD_outFile << event_number<<","<<trID<<","<<parentID<<","<<name<<","<<proc<<","<< "TPC"<<","<< 0 <<","<< 0 <<","<<i<<","<<x<<","<<y<<","<<z<<","<<
                //time<<","<<kinEn<<","<<eDep<<","<< photons <<G4endl;
                if (electrons>0){
                    if (proc !="eIoni" && proc!= "hIoni") // dont want to see the electrons ... 
                    TPC_outFile << event_number<< "," << module_ID << "," << i << "," << trID << "," << parentID << "," << name << ","<< proc <<","  << x <<","<<y<<","<<z<<","<<
                    time<< "," <<kinEn<<","<<eDep << "," << electrons << "," << trackl <<G4endl;
                    
                }
                //std::cout << " TPC hit " << module_ID << ", layer: " << i << ", " << eDep/CLHEP::MeV  << std::endl;
            }
        }

        // PMT Sensitive volume
        std::vector<G4int> PMT_index;
        std::vector<G4int> PMT_index_reduced;
        std::vector<G4double> PMT_time;
        std::vector<G4double> PMT_KE;

        PMTHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID5));
        
        if (PMTHits) {
            hitCount = PMTHits->entries();
            for (G4int h=0; h<hitCount; h++) {
                if (((*PMTHits)[h]) -> GetName() == "opticalphoton"){ // only store photon hit
                    G4double time = ((*PMTHits)[h]) -> GetTime(); 
                    G4double trID = ((*PMTHits)[h]) -> GetTrackID(); 
                    G4int i = ((*PMTHits)[h]) -> GetXID();
                    PMT_index.push_back(i);
                    PMT_time.push_back(time);
                    PMT_KE.push_back(((*PMTHits)[h]) -> GetKinEn());
                }
            }
        }

        #if VERSION_SHIELD==1    
            ShieldHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID7));
            if (ShieldHits) {
            hitCount = ShieldHits->entries();
                for (G4int h=0; h<hitCount; h++) {
                    name     = ((*ShieldHits)[h]) -> GetName();
                    if (name != "opticalphoton"){
                        ltime    = ((*ShieldHits)[h]) -> GetLocalTime();
                        parentID = ((*ShieldHits)[h]) -> GetParentID();
                        proc     = ((*ShieldHits)[h]) -> GetProcess();
                        name     = ((*ShieldHits)[h]) -> GetName();
                        time     = ((*ShieldHits)[h]) -> GetTime(); 
                        trID     = ((*ShieldHits)[h]) -> GetTrackID();
                        kinEn    = ((*ShieldHits)[h]) -> GetKinEn();
                        eDep     = ((*ShieldHits)[h]) -> GetEdep();
                        trackl   = ((*ShieldHits)[h]) -> GetTrackLength();	
                        G4double x = ((*ShieldHits)[h]) -> GetPosX();
                        G4double y = ((*ShieldHits)[h]) -> GetPosY();
                        G4double z = ((*ShieldHits)[h]) -> GetPosZ();
                        G4String vol_name = ((*ShieldHits)[h]) -> GetVolName();
                        G4double pX = ((*ShieldHits)[h]) -> GetPX();
                        G4double pY = ((*ShieldHits)[h]) -> GetPY();
                        G4double pZ = ((*ShieldHits)[h]) -> GetPZ();
                        
                        if (vol_name != "LeadShield"){
                        Shield_outFile << event_number << "," << trID << "," << parentID << "," << name << "," << proc << "," << vol_name << "," << x <<","<<y<<","<<z<<","
                        <<time<< "," <<kinEn<<","<< pX <<","<< pY <<","<< pZ <<","<< eDep << "," << trackl <<G4endl;
                        }
                    }
                }         
            }
        #endif

        /**
        PMT_index_reduced = PMT_index;
        sort(PMT_index_reduced.begin(), PMT_index_reduced.end());
        PMT_index_reduced.erase(unique(PMT_index_reduced.begin(), PMT_index_reduced.end()), PMT_index_reduced.end() );

        // writing the PMT outputs
        //std::cout  << PMT_index_reduced.size() << std::endl;
        PMT_outFile << "{" << "\"Event_ID\":" << event_number << ",";
        PMT_outFile << "\"PMT_signal\":[";

        for (int i = 0; i < PMT_index_reduced.size();i++){
            
            if (i > 0){PMT_outFile << ',';}
            PMT_outFile << "{\"PMT_ID\":" << PMT_index_reduced[i] << ",";
            PMT_outFile << "\"Time\":[" ;
            int count = 0; 
            for (int j=0; j < PMT_index.size();j++){  //
                if (PMT_index[j] == PMT_index_reduced[i]){
                    if (count > 0) {PMT_outFile <<",";}
                    PMT_outFile << PMT_time[j];
                    count++;
                }
            }
            PMT_outFile << "],";

            PMT_outFile << "\"KE\":[" ;
            count = 0;
            for (int j=0; j < PMT_index.size();j++){  //
                if (PMT_index[j] == PMT_index_reduced[i]){
                    if (count > 0) {PMT_outFile <<",";}
                    PMT_outFile << PMT_KE[j];
                    count++;
                }
            }

            PMT_outFile << "]}";
        }
        PMT_outFile << "]";
        PMT_outFile << "}" <<G4endl;
        ***/

        SiliconHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID6));        
        if (SiliconHits) {
            hitCount = SiliconHits->entries();
            for (G4int h=0; h<hitCount; h++) {
                G4double time = ((*SiliconHits)[h]) -> GetTime(); 
                G4double trID = ((*SiliconHits)[h]) -> GetTrackID(); 
                G4int i = ((*SiliconHits)[h]) -> GetXID();
                G4String name     = ((*SiliconHits)[h]) -> GetName();
                G4int parentID = ((*SiliconHits)[h]) -> GetParentID();
                G4String proc = ((*SiliconHits)[h]) -> GetProcess();
                G4double kinEn    = ((*SiliconHits)[h]) -> GetKinEn();
                eDep     = ((*SiliconHits)[h]) -> GetEdep();
                module_ID = ((*SiliconHits)[h]) -> GetMod_ID();
                x = ((*SiliconHits)[h]) -> GetPosX();
                y = ((*SiliconHits)[h]) -> GetPosY();
                z = ((*SiliconHits)[h]) -> GetPosZ();
                G4double trackl = ((*SiliconHits)[h]) -> GetTrackLength();

                Silicon_outFile << event_number<<"," << i << "," << trID << "," << parentID << "," << name << "," << x <<","<<y<<","<<z<<","<<
                time<< "," <<kinEn<<","<<eDep << "," << trackl <<G4endl;
            }
        }
    }
    
    else {G4cout << "No HCE" << G4endl;}
    auto eventID = event->GetEventID();
    auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {G4cout << "---> End of event: " << eventID << G4endl;}
    
    event_number ++;
    
}  
