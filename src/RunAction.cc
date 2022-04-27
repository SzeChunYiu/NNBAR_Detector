#include "RunAction.hh"
#include "Analysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include <iostream>
#include <fstream>
#include <string>
#include "config.h"
//....

extern G4int run_number; 
extern int particle_name_file_index;
extern G4double event_number;
extern std::vector<std::vector<G4double>> particle_gun_record;
extern std::vector<std::vector<G4double>> PMT_record;
extern std::vector<std::vector<G4double>> scint_record;
//extern std::ofstream PMT_outFile;
extern std::ofstream Silicon_outFile;
extern std::ofstream Shield_outFile;
extern std::ofstream TPC_outFile;
extern std::ofstream Scint_layer_outFile;
extern std::ofstream Abs_outFile; 
extern std::ofstream pi0_outFile;
extern std::ofstream Tube_outFile;
extern std::ofstream Particle_outFile; 
extern std::ofstream Lead_glass_outFile;
extern std::ofstream Scint_outFile;


string particle_name_list[7] = { "", "_neutron", "_proton", "_gamma","_electron","_muon","_pion" };

RunAction::RunAction()
 : G4UserRunAction()
{ 
    fMessenger = new G4GenericMessenger(this, "/particle_generator/", "Name the particle for the file name");
    //G4GenericMessenger::Command& filenameCMD = fMessenger->DeclareProperty("Particle_index", particle_name_file_index, "Index of the particle in the order");
    //filenameCMD.SetParameterName("Particle Index", true);
    //filenameCMD.SetDefaultValue("99");

  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);     

  // Create analysis manager, ROOT specified in Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories 
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);

  // Book histograms

  G4cout << "Booking histograms " << G4endl;

  // add histogram here if needed // 


  analysisManager->FinishNtuple();
}

//....

RunAction::~RunAction()
{
  delete G4AnalysisManager::Instance();  
}

//.....

void RunAction::BeginOfRunAction(const G4Run*)
{ 
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  
  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  
  G4String fileName = "./output/calo-sim.root";
  analysisManager->OpenFile(fileName);

  TPC_outFile.open("./output/TPC_output_"+std::to_string(run_number)+".txt");
  Scint_layer_outFile.open("./output/Scintillator_output_"+std::to_string(run_number)+".txt");
  Abs_outFile.open("./output/Lead_Glass_output_"+std::to_string(run_number)+".txt");
  Silicon_outFile.open("./output/Silicon_output_"+std::to_string(run_number)+".txt");
  pi0_outFile.open("./output/pi0_gamma_"+std::to_string(run_number)+".txt");
  Tube_outFile.open("./output/Tube_output"+std::to_string(run_number)+".txt");
  Particle_outFile.open("./output/particle_output_"+std::to_string(run_number)+".txt");
  
  #if VERSION_SHIELD==1      
    Shield_outFile.open("./output/Shield_output_"+std::to_string(run_number)+".txt");
  #endif
  
  Particle_outFile << "Event_ID,Run_ID,PID,Mass,Charge,KE,x,y,z,t,u,v,w"<< G4endl;
  Scint_layer_outFile << "Event_ID,Track_ID,Parent_ID,Name,Proc,module_ID,layer,index,t,KE,eDep,photons,x,y,z,particle_x,particle_y,particle_z"<<G4endl;
  Abs_outFile<< "Event_ID,Track_ID,Parent_ID,Name,Proc,index,t,KE,eDep,trackl,photons,x,y,z"<<G4endl;
  TPC_outFile << "Event_ID,module_ID,Layer,Track_ID,Parent_ID,Name,proc,x,y,z,t,KE,eDep,electrons,trackl"<< G4endl; //x,y,z,
  Silicon_outFile << "Event_ID,layer,Track_ID,Parent_ID,Name,x,y,z,t,KE,eDep,trackl" <<G4endl;
  pi0_outFile << "Event_ID,Track_ID,KE,px,py,pz"<<G4endl; 
  Tube_outFile << "Event_ID,Track_ID,Parent_ID,Name,proc,x,y,z,t,KE,eDep,trackl" <<G4endl;
  
  #if VERSION_SHIELD==1    
    Shield_outFile<<"Event_ID,Track_ID,Parent_ID,Name,proc,Volume,x,y,z,t,KE,px,py,pz,eDep,trackl"<<G4endl;
  #endif 
}

//....

void RunAction::EndOfRunAction(const G4Run*)
{
  
  auto analysisManager = G4AnalysisManager::Instance();
  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();

  // write the PMT records 

  //PMT_outFile.close();
  Silicon_outFile.close();
  TPC_outFile.close();
  Scint_layer_outFile.close();
  Abs_outFile.close();
  pi0_outFile.close();
  Tube_outFile.close();
  Scint_outFile.close();
  Particle_outFile.close();
  Lead_glass_outFile.close();

  #if VERSION_SHIELD==1    
  Shield_outFile.close();
  #endif

  // clear all stored data 
  G4double event_number = 0.0;
  particle_gun_record.clear(); PMT_record.clear(); scint_record.clear();

  run_number++;

}

//....
