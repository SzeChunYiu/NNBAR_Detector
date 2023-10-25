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
#include "NNbarRun.hh"

#include <arrow/io/file.h>
#include <parquet/stream_writer.h>
#include "parquet_writer.h"

#include "G4VHitsCollection.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"
#include "G4Run.hh"

#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "Randomize.hh"
#include <iomanip>
#include "config.h"
//....

extern G4ThreadLocal G4double event_number; 
extern G4ThreadLocal G4int local_event_number;

extern G4ThreadLocal G4int local_event_number_MCPL;
extern G4ThreadLocal G4double event_number_event_action;

extern parquet::StreamWriter Carbon_os;
extern parquet::StreamWriter Silicon_os;
extern parquet::StreamWriter Beampipe_os;
extern parquet::StreamWriter TPC_os;
extern parquet::StreamWriter Scintillator_os;
extern parquet::StreamWriter LeadGlass_os;
extern parquet::StreamWriter Shield_os;
extern parquet::StreamWriter PMT_os;

extern parquetwriter::Writer PMT_output;

namespace {G4Mutex RunActionMutex = G4MUTEX_INITIALIZER;}

NNbarRun::NNbarRun():
	G4Run(),
    fAbsoEdepHCID(-1),
    fGapEdepHCID(-1),
    fAbsoTrackLengthHCID(-1),
    fCerenkovHCID(-1),
    fScintTrackLengthHCID(-1),
    CarbonHitsCollectionID(-1),
    SiliconHitsCollectionID(-1),
    BeampipeHitsCollectionID(-1),
    TPCHitsCollectionID(-1),
    scintHitsCollectionID(-1),
    LeadGlassHitsCollectionID(-1),
    PMTHitsCollectionID(-1),
    ShieldHitsCollectionID(-1)

{
	G4SDManager* pSDManager = G4SDManager::GetSDMpointer();
    if(CarbonHitsCollectionID == -1) {
       CarbonHitsCollectionID = pSDManager->GetCollectionID("CarbonHitCollection");
       SiliconHitsCollectionID = pSDManager->GetCollectionID("SiliconHitCollection");
       BeampipeHitsCollectionID = pSDManager->GetCollectionID("TubeHitCollection");
       TPCHitsCollectionID = pSDManager->GetCollectionID("TPCHitCollection");
       scintHitsCollectionID = pSDManager->GetCollectionID("ScintillatorHitCollection");
       LeadGlassHitsCollectionID = pSDManager->GetCollectionID("LeadGlassHitCollection");
       PMTHitsCollectionID = pSDManager->GetCollectionID("PMTHitCollection");
    }       
}

NNbarRun::~NNbarRun() {}

void NNbarRun::RecordEvent(const G4Event* event)
{
    if(CarbonHitsCollectionID  < 0) {return;}

	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
    int CHCID_tube = -1;
    if (CHCID_tube<0) {CHCID_tube = G4SDManager::GetSDMpointer()->GetCollectionID("TubeHitCollection");}
    NNbarHitsCollection* BeampipeHits  = 0;

    int CHCID_carbon = -1;
    if (CHCID_carbon<0) {CHCID_carbon = G4SDManager::GetSDMpointer()->GetCollectionID("CarbonHitCollection");}
    NNbarHitsCollection* CarbonHits  = 0;
    
    int CHCID_scint = -1;
    if (CHCID_scint<0) {CHCID_scint = G4SDManager::GetSDMpointer()->GetCollectionID("ScintillatorHitCollection");}
    NNbarHitsCollection* ScintHits = 0;

    int CHCID_LeadGlass = -1;
    if (CHCID_LeadGlass<0) {CHCID_LeadGlass = G4SDManager::GetSDMpointer()->GetCollectionID("LeadGlassHitCollection");}
    NNbarHitsCollection* LeadGlassHits   = 0;
    
    int CHCID_TPC = -1;
    if (CHCID_TPC<0) {CHCID_TPC = G4SDManager::GetSDMpointer()->GetCollectionID("TPCHitCollection");}
    NNbarHitsCollection* TPCHits  = 0;
    
    int CHCID_PMT = -1;
    if (CHCID_PMT<0) {CHCID_PMT = G4SDManager::GetSDMpointer()->GetCollectionID("PMTHitCollection");}
    NNbarHitsCollection* PMTHits  = 0;
    
    int CHCID_Si = -1;
    if (CHCID_Si<0) {CHCID_Si = G4SDManager::GetSDMpointer()->GetCollectionID("SiliconHitCollection");}
    NNbarHitsCollection* SiliconHits  = 0;
    
    // #if VERSION_SHIELD==1    
    //     int CHCID_shield = -1;
    //     if (CHCID_shield<0) {CHCID_shield = G4SDManager::GetSDMpointer()->GetCollectionID("ShieldHitCollection");}
    //     NNbarHitsCollection* ShieldHits  = 0;
    // #endif

    if (HCE){

        G4AutoLock lock(&RunActionMutex);
        
        CarbonHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_carbon));
        if (CarbonHits) {
    	    hitCount = CarbonHits->entries();
            for (G4int h=0; h<hitCount; h++) {
                 
                G4int trID     = ((*CarbonHits)[h]) -> GetTrackID();
                G4int parentID = ((*CarbonHits)[h]) -> GetParentID();
                G4String name     = ((*CarbonHits)[h]) -> GetName();
                G4String proc     = ((*CarbonHits)[h]) -> GetProcess();
                int step_info = ((*CarbonHits)[h]) -> GetStepInfo();
                G4String origin_volume = ((*CarbonHits)[h]) -> GetOriginVolName();
                
                G4double x = ((*CarbonHits)[h]) -> GetPosX();
                G4double y = ((*CarbonHits)[h]) -> GetPosY();
                G4double z = ((*CarbonHits)[h]) -> GetPosZ();
                G4double px = ((*CarbonHits)[h]) -> GetPX();
                G4double py = ((*CarbonHits)[h]) -> GetPY();
                G4double pz = ((*CarbonHits)[h]) -> GetPZ();
                G4double time     = ((*CarbonHits)[h]) -> GetTime(); 
                G4double KE    = ((*CarbonHits)[h]) -> GetKinEn();
                G4double eDep     = ((*CarbonHits)[h]) -> GetEdep();
            
                Carbon_os << event->GetEventID() << trID << parentID << name 
                          << proc << step_info << origin_volume 
                          << x/CLHEP::cm << y/CLHEP::cm << z/CLHEP::cm
                          << px << py << pz
                          << time/CLHEP::ns << KE/CLHEP::MeV << eDep/CLHEP::MeV
                          << parquet::EndRow;
                
                // Carbon_outfile << event->GetEventID() << "," << trID << "," << parentID << "," << name << "," 
                //                << proc << "," << step_info << "," << origin_volume <<","
                //                << x/CLHEP::cm <<","<<y/CLHEP::cm<<","<<z/CLHEP::cm << "," 
                //                << px <<","<<py<<","<<pz<<","
                //                << std::setprecision(13) << time/CLHEP::ns << std::setprecision(4) << "," 
                //                << KE/CLHEP::MeV << "," << eDep/CLHEP::MeV << G4endl;

            }        
        }
        
        
        BeampipeHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_tube));
        if (BeampipeHits) {
	    hitCount = BeampipeHits->entries();
            for (G4int h=0; h<hitCount; h++) {
                
                G4int trID     = ((*BeampipeHits)[h]) -> GetTrackID();
                G4int parentID = ((*BeampipeHits)[h]) -> GetParentID();
                G4String name     = ((*BeampipeHits)[h]) -> GetName();
                G4String proc     = ((*BeampipeHits)[h]) -> GetProcess();
                int step_info = ((*BeampipeHits)[h]) -> GetStepInfo();
                G4String origin_volume = ((*BeampipeHits)[h]) -> GetOriginVolName();
                G4String current_vol = ((*BeampipeHits)[h]) -> GetVolName();
                
                G4double x = ((*BeampipeHits)[h]) -> GetPosX();
                G4double y = ((*BeampipeHits)[h]) -> GetPosY();
                G4double z = ((*BeampipeHits)[h]) -> GetPosZ();
                G4double px = ((*BeampipeHits)[h]) -> GetPX();
                G4double py = ((*BeampipeHits)[h]) -> GetPY();
                G4double pz = ((*BeampipeHits)[h]) -> GetPZ();
                G4double time     = ((*BeampipeHits)[h]) -> GetTime(); 
                G4double KE    = ((*BeampipeHits)[h]) -> GetKinEn();
                G4double eDep     = ((*BeampipeHits)[h]) -> GetEdep();
                
                Beampipe_os << event->GetEventID() << trID << parentID << name 
                          << proc << step_info << current_vol << origin_volume
                          << x/CLHEP::cm << y/CLHEP::cm << z/CLHEP::cm
                          << px << py << pz
                          << time/CLHEP::ns << KE/CLHEP::MeV << eDep/CLHEP::MeV
                          << parquet::EndRow;

                // Beampipe_outfile << event->GetEventID() << "," << trID << "," << parentID << "," << name << "," << proc << "," 
                // << step_info << "," << current_vol << "," << origin_volume <<"," 
                // << x/CLHEP::cm <<","<<y/CLHEP::cm<<","<<z/CLHEP::cm<<","<< px <<"," << py << "," << pz <<","
                // << std::setprecision(13) << time/CLHEP::ns << std::setprecision(4) << "," 
                // << KE/CLHEP::MeV<<","<< eDep/CLHEP::MeV <<G4endl;
            }  
        }

        SiliconHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_Si));        
        if (SiliconHits) {
            hitCount = SiliconHits->entries();
            for (G4int h=0; h<hitCount; h++) {

                G4int trID     = ((*SiliconHits)[h]) -> GetTrackID();
                G4int parentID = ((*SiliconHits)[h]) -> GetParentID();
                G4String name     = ((*SiliconHits)[h]) -> GetName();
                G4String proc     = ((*SiliconHits)[h]) -> GetProcess();
                int step_info = ((*SiliconHits)[h]) -> GetStepInfo();
                G4String origin_volume = ((*SiliconHits)[h]) -> GetOriginVolName();
                G4int layer_ID = ((*SiliconHits)[h]) -> GetXID();

                G4double x = ((*SiliconHits)[h]) -> GetPosX();
                G4double y = ((*SiliconHits)[h]) -> GetPosY();
                G4double z = ((*SiliconHits)[h]) -> GetPosZ();
                G4double px = ((*SiliconHits)[h]) -> GetPX();
                G4double py = ((*SiliconHits)[h]) -> GetPY();
                G4double pz = ((*SiliconHits)[h]) -> GetPZ();
                G4double time     = ((*SiliconHits)[h]) -> GetTime(); 
                G4double KE    = ((*SiliconHits)[h]) -> GetKinEn();
                G4double eDep     = ((*SiliconHits)[h]) -> GetEdep();
                
                //std::cout << event->GetEventID() << "," << name << "," << proc << ","<< "KE:" << KE << " edep::"<<eDep << std::endl;

                if (eDep>0.){

                    Silicon_os << event->GetEventID() << trID << parentID << name 
                            << proc << step_info << origin_volume << layer_ID
                            << x/CLHEP::cm << y/CLHEP::cm << z/CLHEP::cm
                            << px << py << pz
                            << time/CLHEP::ns << KE/CLHEP::MeV << eDep/CLHEP::MeV
                            << parquet::EndRow;
                    

                    //  Silicon_outfile << event->GetEventID() <<"," << i << "," << trID << "," << parentID << "," << name << "," 
                    //                  << proc << "," << step_info << "," << origin_volume << "," << layer_ID << ","
                    //                  << x/CLHEP::cm << "," <<y/CLHEP::cm << "," << z/CLHEP::cm <<"," 
                    //                  << px <<"," << py << "," << pz <<","
                    //                  << std::setprecision(13) << time/CLHEP::ns << std::setprecision(4) << "," 
                    //                  << KE/CLHEP::MeV << "," << eDep/CLHEP::MeV <<G4endl;

                }
            }  
        }
        
        TPCHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_TPC));
        if (TPCHits) {
            hitCount = TPCHits->entries();
            for (G4int h=0; h<hitCount; h++) {
            
                G4int trID     = ((*TPCHits)[h]) -> GetTrackID();
                G4int parentID = ((*TPCHits)[h]) -> GetParentID();
                G4String name     = ((*TPCHits)[h]) -> GetName();
                G4String proc     = ((*TPCHits)[h]) -> GetProcess();
                int step_info = ((*TPCHits)[h]) -> GetStepInfo();
                G4String origin_volume = ((*TPCHits)[h]) -> GetOriginVolName();
                G4String current_vol = ((*TPCHits)[h]) -> GetVolName();

                G4double x = ((*TPCHits)[h]) -> GetPosX();
                G4double y = ((*TPCHits)[h]) -> GetPosY();
                G4double z = ((*TPCHits)[h]) -> GetPosZ();
                G4double px = ((*TPCHits)[h]) -> GetPX();
                G4double py = ((*TPCHits)[h]) -> GetPY();
                G4double pz = ((*TPCHits)[h]) -> GetPZ();
                G4double time     = ((*TPCHits)[h]) -> GetTime(); 
                G4double KE    = ((*TPCHits)[h]) -> GetKinEn();
                G4double eDep     = ((*TPCHits)[h]) -> GetEdep();
                G4double trackl = ((*TPCHits)[h]) -> GetTrackLength();
                
                //TPC unique information
                G4int module_ID = ((*TPCHits)[h]) -> GetMod_ID();
                G4int layer = ((*TPCHits)[h]) -> GetXID(); // TPC layer 
                G4double electrons = ((*TPCHits)[h]) -> GetPhotons(); 

                if (electrons>0){
                if (proc !="eIoni" && proc!= "hIoni"&& proc!= "muIoni") {// dont want to see the electrons ... 
                    
                    TPC_os<< event->GetEventID() << trID << parentID << name 
                          << proc << step_info << origin_volume << current_vol
                          << module_ID << layer
                          << x/CLHEP::cm << y/CLHEP::cm << z/CLHEP::cm
                          << px << py << pz
                          << time/CLHEP::ns << KE/CLHEP::MeV << eDep/CLHEP::MeV << trackl/CLHEP::cm << electrons
                          << parquet::EndRow;                    
                    
                    // TPC_outfile << event->GetEventID() << "," <<  trID << "," << parentID << "," << name << ","
                    //             << proc <<"," << step_info <<","<< origin_volume << "," << current_vol << ","
                    //             << module_ID << "," << layer << ","
                    //             << x/CLHEP::cm << "," << y/CLHEP::cm << "," << z/CLHEP::cm <<","
                    //             << px << ","<< py << ","<< pz<< ","
                    //             << std::setprecision(13) << time/CLHEP::ns << "," 
                    //             << KE/CLHEP::MeV << "," << eDep/CLHEP::MeV << "," << trackl/CLHEP::cm << "," << electrons <<G4endl;
                }   
                }   
            }
        }

        ScintHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_scint));
        if (ScintHits) {
           hitCount = ScintHits->entries();
           for (G4int h=0; h<hitCount; h++) { 
               G4String name = ((*ScintHits)[h]) -> GetName();
               if (name != "opticalphoton"){
                    
                    G4int trID     = ((*ScintHits)[h]) -> GetTrackID();
                    G4int parentID = ((*ScintHits)[h]) -> GetParentID();
                    G4String name     = ((*ScintHits)[h]) -> GetName();
                    G4String proc     = ((*ScintHits)[h]) -> GetProcess();
                    int step_info = ((*ScintHits)[h]) -> GetStepInfo();
                    G4String current_volume = ((*ScintHits)[h]) -> GetVolName();
                    G4String origin_volume = ((*ScintHits)[h]) -> GetOriginVolName();
                    G4double particle_x = ((*ScintHits)[h]) -> GetPosX_particle();
                    G4double particle_y = ((*ScintHits)[h]) -> GetPosY_particle();
                    G4double particle_z = ((*ScintHits)[h]) -> GetPosZ_particle();
                    G4double px = ((*ScintHits)[h]) -> GetPX();
                    G4double py = ((*ScintHits)[h]) -> GetPY();
                    G4double pz = ((*ScintHits)[h]) -> GetPZ();
                    G4double time     = ((*ScintHits)[h]) -> GetTime(); 
                    G4double KE    = ((*ScintHits)[h]) -> GetKinEn();
                    G4double eDep     = ((*ScintHits)[h]) -> GetEdep();

                    //Scintillator information
                    auto stave_ID  = ((*ScintHits)[h]) -> GetStave_ID();
                    auto layer_ID  = ((*ScintHits)[h]) -> GetXID();
                    auto module_ID = ((*ScintHits)[h]) -> GetMod_ID();
                    G4double x = ((*ScintHits)[h]) -> GetPosX();
                    G4double y = ((*ScintHits)[h]) -> GetPosY();
                    G4double z = ((*ScintHits)[h]) -> GetPosZ();

                    G4double x_local = ((*ScintHits)[h]) -> GetLocalPosX();
                    G4double y_local = ((*ScintHits)[h]) -> GetLocalPosY();
                    G4double z_local = ((*ScintHits)[h]) -> GetLocalPosZ();

                    G4int scint_photons_per_hit = ((*ScintHits)[h])->GetPhotons();
                    
                    // WARNING::Origin volume has bug, it only returns the last volume the particle was in...
                    if (scint_photons_per_hit>0){

                        Scintillator_os<< event->GetEventID() << trID << parentID << name 
                        << proc << step_info << origin_volume << current_volume
                        << module_ID <<  layer_ID << stave_ID << x/CLHEP::cm << y/CLHEP::cm << z/CLHEP::cm
                        << particle_x/CLHEP::cm << particle_y/CLHEP::cm << particle_z/CLHEP::cm
                        << x_local/CLHEP::cm << y_local/CLHEP::cm << z_local/CLHEP::cm
                        << time/CLHEP::ns << KE/CLHEP::MeV << eDep/CLHEP::MeV << scint_photons_per_hit
                        << parquet::EndRow;                        
                        

                        // Scintillator_outfile<< event->GetEventID()<< "," << trID<< "," << parentID<< "," << name << ","
                        // << proc << ","<< step_info<< "," << origin_volume<< ","
                        // << module_ID<< "," <<  layer_ID<< "," << stave_ID<< "," 
                        // << x/CLHEP::cm << ","<< y/CLHEP::cm << ","<< z/CLHEP::cm<< ","
                        // << particle_x/CLHEP::cm << ","<< particle_y/CLHEP::cm << ","<< particle_z/CLHEP::cm<< ","
                        // << time/CLHEP::ns << ","<< KE/CLHEP::MeV << ","<< eDep/CLHEP::MeV << ","<< scint_photons_per_hit
                        // << G4endl;
                    }
                }
            }
        }

        LeadGlassHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_LeadGlass));
        if(LeadGlassHits) {
            for (G4int h=0; h<LeadGlassHits->entries(); h++) {
                G4String name = ((*LeadGlassHits)[h]) -> GetName();
                if (name != "opticalphoton"){
                    
                    // WARNING::Origin volume has bug, it only returns the last volume the particle was in...
                    G4int trID     = ((*LeadGlassHits)[h]) -> GetTrackID();
                    G4int parentID = ((*LeadGlassHits)[h]) -> GetParentID();
                    G4String name     = ((*LeadGlassHits)[h]) -> GetName();
                    G4String proc     = ((*LeadGlassHits)[h]) -> GetProcess();
                    int step_info = ((*LeadGlassHits)[h]) -> GetStepInfo();
                    G4String origin_volume = ((*LeadGlassHits)[h]) -> GetOriginVolName();

                    G4double x = ((*LeadGlassHits)[h]) -> GetPosX();
                    G4double y = ((*LeadGlassHits)[h]) -> GetPosY();
                    G4double z = ((*LeadGlassHits)[h]) -> GetPosZ();
                    G4double time     = ((*LeadGlassHits)[h]) -> GetTime(); 
                    G4double KE    = ((*LeadGlassHits)[h]) -> GetKinEn();
                    G4double eDep     = ((*LeadGlassHits)[h]) -> GetEdep();
                    G4double trackl = ((*LeadGlassHits)[h]) -> GetTrackLength();

                    //LeadGlass information
                    G4int module_ID = ((*LeadGlassHits)[h]) -> GetXID();
                    G4int LeadGlass_photons_per_hit = ((*LeadGlassHits)[h])->GetPhotons();

                    if (LeadGlass_photons_per_hit>0){
                        
                        LeadGlass_os<< event->GetEventID() << trID << parentID << name 
                        << proc 
                        << step_info 
                        << origin_volume
                        << module_ID << x/CLHEP::cm << y/CLHEP::cm << z/CLHEP::cm
                        << time/CLHEP::ns << KE/CLHEP::MeV << eDep/CLHEP::MeV << LeadGlass_photons_per_hit
                        << parquet::EndRow;                    
                        
                        // LeadGlass_outfile<< event->GetEventID() << ","<< trID << ","<< parentID << ","<< name << ","
                        // << proc << ","<< step_info<< "," << origin_volume<< ","
                        // << module_ID<< "," << x/CLHEP::cm<< "," << y/CLHEP::cm<< "," << z/CLHEP::cm<< ","
                        // << time/CLHEP::ns << "," << KE/CLHEP::MeV << ","<< eDep/CLHEP::MeV<< "," << LeadGlass_photons_per_hit
                        // << G4endl;                      
                    }
                }
            }
        }

        // PMT Sensitive volume
        std::vector<G4int> PMT_index;
        std::vector<G4int> PMT_index_reduced;
        std::vector<G4double> PMT_time;
        std::vector<G4double> PMT_KE;
        std::vector<G4double> PMT_x;
        std::vector<G4double> PMT_y;
        std::vector<G4double> PMT_z;

        PMTHits = (NNbarHitsCollection*)(HCE->GetHC(CHCID_PMT));
        
        if (PMTHits) {
            hitCount = PMTHits->entries();
            for (G4int h=0; h<hitCount; h++) {

                if (((*PMTHits)[h]) -> GetName() == "opticalphoton"){ // only store photon hit
                    
                    G4double time = ((*PMTHits)[h]) -> GetTime(); 
                    G4int module_ID = ((*PMTHits)[h]) -> GetXID();
                    PMT_x.push_back(((*PMTHits)[h]) -> GetPosX());
                    PMT_y.push_back(((*PMTHits)[h]) -> GetPosY());
                    PMT_z.push_back(((*PMTHits)[h]) -> GetPosZ());
                    PMT_index.push_back(module_ID);
                    PMT_time.push_back(time);
                    PMT_KE.push_back(((*PMTHits)[h]) -> GetKinEn());
                }
            }
        }
        
        PMT_index_reduced = PMT_index;
        sort(PMT_index_reduced.begin(), PMT_index_reduced.end());
        PMT_index_reduced.erase(unique(PMT_index_reduced.begin(), PMT_index_reduced.end()), PMT_index_reduced.end() );
        
        // writing the PMT outputs
        for (int i = 0; i < PMT_index_reduced.size();i++){
            
            int count = 0;
            G4double x_pmt = 0;
            G4double y_pmt = 0;
            G4double z_pmt = 0;
            std::vector<G4double> PMT_photon_time;
            std::vector<G4double> PMT_photon_KE;

            for (int j=0; j < PMT_index.size();j++){ // go through all pmt entries
                if (PMT_index[j] == PMT_index_reduced[i]){ 
                    x_pmt = PMT_x[j]/CLHEP::cm;
                    y_pmt = PMT_y[j]/CLHEP::cm;
                    z_pmt = PMT_z[j]/CLHEP::cm;
                    G4double t_pmt = PMT_time[j];
                    G4double KE_pmt = PMT_KE[j];
                    PMT_photon_time.push_back(t_pmt);
                    PMT_photon_KE.push_back(KE_pmt);
                    count++;
                }
            }

            PMT_output.fill("Event_ID", event->GetEventID());
            PMT_output.fill("Module_ID", PMT_index_reduced[i]);
            PMT_output.fill("x", x_pmt);
            PMT_output.fill("y", y_pmt);
            PMT_output.fill("z", z_pmt);
            PMT_output.fill("photons", count);
            PMT_output.fill("KE", PMT_photon_KE);
            PMT_output.fill("t", PMT_photon_time);
            PMT_output.end_row();
        }
        
    }
}

void NNbarRun::Merge(const G4Run* aRun)
{
  const NNbarRun* localRun = static_cast<const NNbarRun*>(aRun);
  G4Run::Merge(aRun);
} 
