#include "EventAction.hh"
#include "Analysis.hh"

#include "Sensitive_Detector/ScintillatorSD.hh"
#include "Sensitive_Detector/LeadGlassSD.hh"
#include "Sensitive_Detector/TubeSD.hh"
#include "Sensitive_Detector/PMTSD.hh"
#include "Sensitive_Detector/TPCSD.hh"
#include "Sensitive_Detector/Scint_DetSD.hh"
#include "Sensitive_Detector/ShieldSD.hh"

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


EventAction::EventAction(): 
    G4UserEventAction(),
    fAbsoEdepHCID(-1),
    fGapEdepHCID(-1),
    fAbsoTrackLengthHCID(-1),
    fCerenkovHCID(-1),
    fScintTrackLengthHCID(-1),
    //DetArea_HitsCollectionID(-1),
    CarbonHitsCollectionID(-1),
    SiliconHitsCollectionID(-1),
    tubeHitsCollectionID(-1),
    TPCHitsCollectionID(-1),
    scintHitsCollectionID(-1),
    absHitsCollectionID(-1),
    PMTHitsCollectionID(-1),
    ShieldHitsCollectionID(-1)
{}

//....

EventAction::~EventAction()
{}

//....

G4THitsMap<G4double>* 
EventAction::GetHitsCollection(G4int hcID,const G4Event* event) const
{
  auto hitsCollection = static_cast<G4THitsMap<G4double>*>(event->GetHCofThisEvent()->GetHC(hcID));
  
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
    // std::cout<<"BEGIN of EVENT ACTION" << std::endl;
    // G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
    // pSDManager->ListTree();
    // std::cout<<"************" << std::endl;

    // if(CarbonHitsCollectionID == -1) {
    //    CarbonHitsCollectionID = pSDManager->GetCollectionID("CarbonHitCollection");
    //    SiliconHitsCollectionID = pSDManager->GetCollectionID("SiliconHitCollection");
    //    tubeHitsCollectionID = pSDManager->GetCollectionID("TubeHitCollection");
    //    TPCHitsCollectionID = pSDManager->GetCollectionID("TPCHitCollection");
       
    //    scintHitsCollectionID = pSDManager->GetCollectionID("ScintillatorHitCollection");
    //    absHitsCollectionID = pSDManager->GetCollectionID("AbsorberHitCollection");
    //    PMTHitsCollectionID = pSDManager->GetCollectionID("PMTHitCollection");
       
    //    #if VERSION_SHIELD==1
    //       ShieldHitsCollectionID = pSDManager->GetCollectionID("ShieldHitCollection");
    //    #endif
       
       
    // }
}

void EventAction::EndOfEventAction(const G4Event* event)
{  
 
    // if(CarbonHitsCollectionID  < 0) {return;}

	// G4HCofThisEvent* HCE = event->GetHCofThisEvent();

    // int CHCID_tube = -1;
    // if (CHCID_tube<0) {CHCID_tube = G4SDManager::GetSDMpointer()->GetCollectionID("TubeHitCollection");}
    // NNbarHitsCollection* TubeHits  = 0;

    // int CHCID_carbon = -1;
    // if (CHCID_carbon<0) {CHCID_carbon = G4SDManager::GetSDMpointer()->GetCollectionID("CarbonHitCollection");}
    // NNbarHitsCollection* CarbonHits  = 0;
    
    // int CHCID_scint = -1;
    // if (CHCID_scint<0) {CHCID_scint = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitCollection");}
    // NNbarHitsCollection* ScintHits = 0;

    // int CHCID_abs = -1;
    // if (CHCID_abs<0) {CHCID_abs = G4SDManager::GetSDMpointer()->GetCollectionID("AbsorberHitCollection");}
    // NNbarHitsCollection* AbsHits   = 0;
    
    // int CHCID_TPC = -1;
    // if (CHCID_TPC<0) {CHCID_TPC = G4SDManager::GetSDMpointer()->GetCollectionID("TPCHitCollection");}
    // NNbarHitsCollection* TPCHits  = 0;
    
    // int CHCID_PMT = -1;
    // if (CHCID_PMT<0) {CHCID_PMT = G4SDManager::GetSDMpointer()->GetCollectionID("PMTHitCollection");}
    // NNbarHitsCollection* PMTHits  = 0;
    
    // int CHCID_Si = -1;
    // if (CHCID_Si<0) {CHCID_Si = G4SDManager::GetSDMpointer()->GetCollectionID("SiliconHitCollection");}
    // NNbarHitsCollection* SiliconHits  = 0;
    
    // // #if VERSION_SHIELD==1    
    // //     int CHCID_shield = -1;
    // //     if (CHCID_shield<0) {CHCID_shield = G4SDManager::GetSDMpointer()->GetCollectionID("ShieldHitCollection");}
    // //     NNbarHitsCollection* ShieldHits  = 0;
    // // #endif


    // //std::cout << "Event Action Starts " << std::endl;
    // if (HCE) {

    //     G4AnalysisManager* analysis = G4AnalysisManager::Instance();
    //     //std::cout << "Event Action Starts Carbon" << std::endl;
    //     CarbonHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_carbon));
    //     if (CarbonHits) {
    // 	    hitCount = CarbonHits->entries();
    //         for (G4int h=0; h<hitCount; h++) {
    //             ltime    = ((*CarbonHits)[h]) -> GetLocalTime();
    //             parentID = ((*CarbonHits)[h]) -> GetParentID();
    //             proc     = ((*CarbonHits)[h]) -> GetProcess();
    //             name     = ((*CarbonHits)[h]) -> GetName();
    //             time     = ((*CarbonHits)[h]) -> GetTime(); 
    //             trID     = ((*CarbonHits)[h]) -> GetTrackID();
    //             kinEn    = ((*CarbonHits)[h]) -> GetKinEn();
    //             eDep     = ((*CarbonHits)[h]) -> GetEdep();
    //             trackl   = ((*CarbonHits)[h]) -> GetTrackLength();	
    //             G4double x = ((*CarbonHits)[h]) -> GetPosX();
    //             G4double y = ((*CarbonHits)[h]) -> GetPosY();
    //             G4double z = ((*CarbonHits)[h]) -> GetPosZ();
    //             G4double px = ((*CarbonHits)[h]) -> GetPX();
    //             G4double py = ((*CarbonHits)[h]) -> GetPY();
    //             G4double pz = ((*CarbonHits)[h]) -> GetPZ();
    //             int step_info = ((*CarbonHits)[h]) -> GetStepInfo();
    //             G4String origin_volume = ((*CarbonHits)[h]) -> GetOriginVolName();
                
    //             Carbon_outFile << int(event_number) << "," << trID << "," << parentID << "," << name << "," << proc << "," << step_info << "," << origin_volume
    //             <<","<< x <<","<<y<<","<<z << "," << px <<","<<py<<","<<pz<<"," << std::setprecision(13) << time << std::setprecision(4) << "," <<kinEn<<","<<eDep <<G4endl;
            
    //         }         
    //     }

    //     SiliconHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_Si));        
    //     if (SiliconHits) {
    //         hitCount = SiliconHits->entries();
    //         for (G4int h=0; h<hitCount; h++) {
    //             G4double time = ((*SiliconHits)[h]) -> GetTime(); 
    //             G4double trID = ((*SiliconHits)[h]) -> GetTrackID(); 
    //             G4int i = ((*SiliconHits)[h]) -> GetXID();
    //             G4String name     = ((*SiliconHits)[h]) -> GetName();
    //             G4int parentID = ((*SiliconHits)[h]) -> GetParentID();
    //             G4String proc = ((*SiliconHits)[h]) -> GetProcess();
    //             G4double kinEn    = ((*SiliconHits)[h]) -> GetKinEn();
    //             eDep     = ((*SiliconHits)[h]) -> GetEdep();
    //             module_ID = ((*SiliconHits)[h]) -> GetMod_ID();
    //             x = ((*SiliconHits)[h]) -> GetPosX();
    //             y = ((*SiliconHits)[h]) -> GetPosY();
    //             z = ((*SiliconHits)[h]) -> GetPosZ();
    //             G4double trackl = ((*SiliconHits)[h]) -> GetTrackLength();

    //             Silicon_outFile << int(event_number)<<"," << i << "," << trID << "," << parentID << "," << name << "," << x <<","<<y<<","<<z<<","<<
    //             std::setprecision(13) << time<< "," <<kinEn<<","<<eDep << "," << trackl <<G4endl;
    //         }
    //     }

    //     TubeHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_tube));
    //     //std::cout << "Event Action Starts Tube" << std::endl;
    //     if (TubeHits) {
	//     hitCount = TubeHits->entries();
    //         for (G4int h=0; h<hitCount; h++) {
    //             ltime    = ((*TubeHits)[h]) -> GetLocalTime();
    //             parentID = ((*TubeHits)[h]) -> GetParentID();
    //             proc     = ((*TubeHits)[h]) -> GetProcess();
    //             name     = ((*TubeHits)[h]) -> GetName();
    //             time     = ((*TubeHits)[h]) -> GetTime(); 
    //             trID     = ((*TubeHits)[h]) -> GetTrackID();
    //             kinEn    = ((*TubeHits)[h]) -> GetKinEn();
    //             eDep     = ((*TubeHits)[h]) -> GetEdep();
    //             trackl   = ((*TubeHits)[h]) -> GetTrackLength();	
    //             G4double x = ((*TubeHits)[h]) -> GetPosX();
    //             G4double y = ((*TubeHits)[h]) -> GetPosY();
    //             G4double z = ((*TubeHits)[h]) -> GetPosZ();
    //             G4String current_vol = ((*TubeHits)[h]) -> GetVolName();
    //             int step_info = ((*TubeHits)[h]) -> GetStepInfo();
    //             G4String origin_volume = ((*TubeHits)[h]) -> GetOriginVolName();

    //             Tube_outFile << int(event_number) << "," << trID << "," << parentID << "," << name << "," << proc << "," << current_vol << "," << step_info << "," << origin_volume
    //             <<"," << x <<","<<y<<","<<z<<","
    //             << std::setprecision(13) << time << std::setprecision(4) << "," <<kinEn<<","<<eDep << "," << trackl <<G4endl;
            
    //         }         
    //     }

    //     TPCHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_TPC));
    //     if (TPCHits) {
    //         hitCount = TPCHits->entries();
    //         for (G4int h=0; h<hitCount; h++) {
                
    //             ltime    = ((*TPCHits)[h]) -> GetLocalTime();
    //             parentID = ((*TPCHits)[h]) -> GetParentID();
    //             proc     = ((*TPCHits)[h]) -> GetProcess();
    //             name     = ((*TPCHits)[h]) -> GetName();
    //             time     = ((*TPCHits)[h]) -> GetTime();
    //             trID     = ((*TPCHits)[h]) -> GetTrackID();
    //             kinEn    = ((*TPCHits)[h]) -> GetKinEn();
    //             eDep     = ((*TPCHits)[h]) -> GetEdep();
    //             module_ID = ((*TPCHits)[h]) -> GetMod_ID();
    //             i = ((*TPCHits)[h]) -> GetXID();
    //             x = ((*TPCHits)[h]) -> GetPosX();
    //             y = ((*TPCHits)[h]) -> GetPosY();
    //             z = ((*TPCHits)[h]) -> GetPosZ();
    //             G4double electrons = ((*TPCHits)[h]) -> GetPhotons(); 
    //             G4double trackl = ((*TPCHits)[h]) -> GetTrackLength();

    //             //SD_outFile << event_number<<","<<trID<<","<<parentID<<","<<name<<","<<proc<<","<< "TPC"<<","<< 0 <<","<< 0 <<","<<i<<","<<x<<","<<y<<","<<z<<","<<
    //             //time<<","<<kinEn<<","<<eDep<<","<< photons <<G4endl;
    //             if (electrons>0){
    //                 if (proc !="eIoni" && proc!= "hIoni") {// dont want to see the electrons ... 
    //                 TPC_outFile << int(event_number)<< "," << module_ID << "," << i << "," << trID << "," << parentID << "," << name << ","<< proc <<","  << x <<","<<y<<","<<z<<","
    //                 << std::setprecision(13) << time<< "," <<kinEn<<","<<eDep << "," << electrons << "," << trackl <<G4endl;}
                    
    //             }
    //             //std::cout << " TPC hit " << module_ID << ", layer: " << i << ", " << eDep/CLHEP::MeV  << std::endl;
    //         }
    //     }
    //     ScintHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_scint));
    //     if (ScintHits) {
    //        hitCount = ScintHits->entries();
    //        for (G4int h=0; h<hitCount; h++) { 
    //            name     = ((*ScintHits)[h]) -> GetName();
    //            if (name != "opticalphoton"){
    //                 //ltime    = ((*ScintHits)[h]) -> GetLocalTime();
    //                 parentID = ((*ScintHits)[h]) -> GetParentID();
    //                 trID     = ((*ScintHits)[h]) -> GetTrackID();
    //                 proc     = ((*ScintHits)[h]) -> GetProcess();
    //                 x = ((*ScintHits)[h]) -> GetPosX();
    //                 y = ((*ScintHits)[h]) -> GetPosY();
    //                 z = ((*ScintHits)[h]) -> GetPosZ();
                    
    //                 G4double particle_x = ((*ScintHits)[h]) -> GetPosX_particle();
    //                 G4double particle_y = ((*ScintHits)[h]) -> GetPosY_particle();
    //                 G4double particle_z = ((*ScintHits)[h]) -> GetPosZ_particle();
                    
    //                 time     = ((*ScintHits)[h]) -> GetTime();

    //                 auto stave_ID = ((*ScintHits)[h]) -> GetStave_ID(); 
    //                 i        = ((*ScintHits)[h]) -> GetXID();
    //                 //group_ID = ((*ScintHits)[h]) -> GetGroup_ID();
    //                 module_ID = ((*ScintHits)[h]) -> GetMod_ID();

    //                 kinEn    = ((*ScintHits)[h]) -> GetKinEn();
    //                 eDep     = ((*ScintHits)[h]) -> GetEdep();
    //                 trackl   = ((*ScintHits)[h]) -> GetTrackLength();
    //                 G4int scint_photons_per_hit = ((*ScintHits)[h])->GetPhotons();
                    
    //                 int step_info = ((*ScintHits)[h]) -> GetStepInfo();
    //                 G4String origin_volume = ((*ScintHits)[h]) -> GetOriginVolName();
                    
    //                 // writing output

    //                 if (scint_photons_per_hit>0){
    //                 Scint_layer_outFile << int(event_number)<<","<<trID<<","<<parentID<<","<<name<<","<<proc<<","
    //                 << origin_volume << "," << step_info << ","
    //                 << module_ID<<","<<i<< "," << stave_ID << ","
    //                 << std::setprecision(13) << time <<","<<kinEn<<","<<eDep<<","<< scint_photons_per_hit << "," << x <<","<< y << ","<< z << "," << particle_x << "," << particle_y << "," << particle_z  <<G4endl;
    //                 }
    //                 //std::cout<< group_ID<<","<<module_ID<<","<<i<< "," << stave_ID << std::endl;

    //            }

    //         }
    //     }

    //     AbsHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_abs));
    //     if(AbsHits) {
    //         for (G4int h=0; h<AbsHits->entries(); h++) {
    //             name     = ((*AbsHits)[h]) -> GetName();
    //             if (name != "opticalphoton"){
    //                 ltime           = ((*AbsHits)[h]) -> GetLocalTime();
    //                 parentID	= ((*AbsHits)[h]) -> GetParentID();
    //                 proc            = ((*AbsHits)[h]) -> GetProcess();
    //                 G4String name   = ((*AbsHits)[h]) -> GetName();
    //                 G4double time   = ((*AbsHits)[h]) -> GetTime();
    //                 G4int trID      = ((*AbsHits)[h]) -> GetTrackID();
    //                 G4int i         = ((*AbsHits)[h]) -> GetXID();
    //                 G4double kinEn  = ((*AbsHits)[h]) -> GetKinEn();
    //                 G4double eDep   = ((*AbsHits)[h]) -> GetEdep();
    //                 G4double trackl = ((*AbsHits)[h]) -> GetTrackLength();
    //                 G4double photons_cerenkov = ((*AbsHits)[h])->GetPhotons();
    //                 x = ((*AbsHits)[h]) -> GetPosX();
    //                 y = ((*AbsHits)[h]) -> GetPosY();
    //                 z = ((*AbsHits)[h]) -> GetPosZ();

    //                 int step_info = ((*AbsHits)[h]) -> GetStepInfo();
    //                 G4String origin_volume = ((*AbsHits)[h]) -> GetOriginVolName();

    //                 if (photons_cerenkov>0){
    //                 Abs_outFile << int(event_number)<<","<<trID << ","<<parentID<<","<<name<<","<<proc<<","
    //                 << origin_volume << "," << step_info << ","
    //                 <<i<<","
    //                 << std::setprecision(13)<< time<<","<<kinEn<<","<<eDep<<","<< trackl <<","<< photons_cerenkov <<  "," << x <<","<< y << ","<<z<<G4endl;
    //                 }
    //             }
    //             //cerenkovCounter = cerenkovCounter + photons_cerenkov;
    //         }
    //     }

        
    //     // PMT Sensitive volume
    //     std::vector<G4int> PMT_index;
    //     std::vector<G4int> PMT_index_reduced;
    //     std::vector<G4double> PMT_time;
    //     std::vector<G4double> PMT_KE;

    //     PMTHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_PMT));
        
    //     if (PMTHits) {
    //         hitCount = PMTHits->entries();
    //         for (G4int h=0; h<hitCount; h++) {
    //             if (((*PMTHits)[h]) -> GetName() == "opticalphoton"){ // only store photon hit
    //                 G4double time = ((*PMTHits)[h]) -> GetTime(); 
    //                 G4double trID = ((*PMTHits)[h]) -> GetTrackID(); 
    //                 G4int i = ((*PMTHits)[h]) -> GetXID();
    //                 PMT_index.push_back(i);
    //                 PMT_time.push_back(time);
    //                 PMT_KE.push_back(((*PMTHits)[h]) -> GetKinEn());
    //             }
    //         }
    //     }

    //     // #if VERSION_SHIELD==1    
    //     //     ShieldHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_shield));
    //     //     if (ShieldHits) {
    //     //     hitCount = ShieldHits->entries();
    //     //         for (G4int h=0; h<hitCount; h++) {
    //     //             name     = ((*ShieldHits)[h]) -> GetName();
    //     //             if (name != "opticalphoton"){
    //     //                 ltime    = ((*ShieldHits)[h]) -> GetLocalTime();
    //     //                 parentID = ((*ShieldHits)[h]) -> GetParentID();
    //     //                 proc     = ((*ShieldHits)[h]) -> GetProcess();
    //     //                 name     = ((*ShieldHits)[h]) -> GetName();
    //     //                 time     = ((*ShieldHits)[h]) -> GetTime(); 
    //     //                 trID     = ((*ShieldHits)[h]) -> GetTrackID();
    //     //                 kinEn    = ((*ShieldHits)[h]) -> GetKinEn();
    //     //                 eDep     = ((*ShieldHits)[h]) -> GetEdep();
    //     //                 trackl   = ((*ShieldHits)[h]) -> GetTrackLength();	
    //     //                 G4double x = ((*ShieldHits)[h]) -> GetPosX();
    //     //                 G4double y = ((*ShieldHits)[h]) -> GetPosY();
    //     //                 G4double z = ((*ShieldHits)[h]) -> GetPosZ();
    //     //                 G4String vol_name = ((*ShieldHits)[h]) -> GetVolName();
    //     //                 G4double pX = ((*ShieldHits)[h]) -> GetPX();
    //     //                 G4double pY = ((*ShieldHits)[h]) -> GetPY();
    //     //                 G4double pZ = ((*ShieldHits)[h]) -> GetPZ();
                        
    //     //                 if (vol_name != "LeadShield"){
    //     //                 Shield_outFile << event_number << "," << trID << "," << parentID << "," << name << "," << proc << "," << vol_name << "," << x <<","<<y<<","<<z<<","
    //     //                 <<time<< "," <<kinEn<<","<< pX <<","<< pY <<","<< pZ <<","<< eDep << "," << trackl <<G4endl;
    //     //                 }
    //     //             }
    //     //         }         
    //     //     }
    //     // #endif
       


    // }
    // //std::cout << "Event Action Ends===" << std::endl;
    // // else {G4cout << "No HCE" << G4endl;}
    // int eventID = event->GetEventID();
    // auto printModulo = G4RunManager::GetRunManager()->GetPrintProgress();
    // if ( ( printModulo > 0 ) && ( eventID % printModulo == 0 ) ) {G4cout << "---> End of event: " << event_number << " local event ID: "<< eventID << G4endl;}
    
    // event_number ++;
}  
