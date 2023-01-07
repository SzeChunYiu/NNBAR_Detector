#include "RunAction.hh"
#include "Analysis.hh"
#include "NNbarRun.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include <iostream>
#include <fstream>
#include <string>
#include "config.h"


#include "mcpl.h"

//....

extern G4int run_number; 
extern int particle_name_file_index;
extern G4double event_number;
extern std::vector<std::vector<G4double>> particle_gun_record;
extern std::vector<std::vector<G4double>> PMT_record;
extern std::vector<std::vector<G4double>> scint_record;
//extern std::ofstream PMT_outFile;
extern std::ofstream pi0_outFile;
extern std::ofstream calorimeter_photon_outFile ;
extern std::ofstream Particle_outFile; 

extern std::ofstream  Silicon_outFile;
extern std::ofstream  Carbon_outFile;
extern std::ofstream  Tube_outFile;
extern std::ofstream  TPC_outFile;
extern std::ofstream  Scint_layer_outFile;
extern std::ofstream  Abs_outFile;

#if VERSION_SHIELD==1   
extern std::ofstream Shield_outFile;
#endif

extern int theta_bin_index;
extern int KE_bin_index;
extern int particle_name_input;

extern mcpl_outfile_t f;
extern mcpl_particle_t *p;

string particle_name_list[7] = { "", "_neutron", "_proton", "_gamma","_electron","_muon","_pion" };

RunAction::RunAction(): G4UserRunAction()
{ 
  fMessenger = new G4GenericMessenger(this, "/particle_generator/", "Name the particle for the file name");
  G4GenericMessenger::Command& angleCMD = fMessenger->DeclareProperty("angle_index",theta_bin_index,"Angle index of Event");
	angleCMD.SetParameterName("angle_index",true);
	angleCMD.SetRange("angle_index>=0");
	angleCMD.SetDefaultValue("0");

	G4GenericMessenger::Command& KE_CMD = fMessenger->DeclareProperty("KE_index",KE_bin_index,"KE index of Event");
	KE_CMD.SetParameterName("KE_index",true);
	KE_CMD.SetRange("KE_index>=0");
	KE_CMD.SetDefaultValue("0");

	//Just need a control to the run nummber  
	G4GenericMessenger::Command& Run_number_CMD = fMessenger->DeclareProperty("set_run_number",run_number,"Run number of the simulation");
	Run_number_CMD.SetParameterName("run_number",true);
	Run_number_CMD.SetRange("run_number>=0");
	Run_number_CMD.SetDefaultValue("0");    

  G4GenericMessenger::Command& Particle_name_CMD = fMessenger->DeclareProperty("set_particle",particle_name_input,"Please give particle name [Default: pi+]");
	Particle_name_CMD.SetParameterName("particle_name_input",true);
	//Particle_name_CMD.SetDefaultValue("211");    
	
  G4RunManager::GetRunManager()->SetPrintProgress(1);
  // empty? //
}

RunAction::~RunAction()
{}

G4Run* RunAction::GenerateRun()
{return new NNbarRun;}

void RunAction::BeginOfRunAction(const G4Run*aRun)
{ 

  if(IsMaster()){  // if this is master, then open the new files!
  Carbon_outFile.open("./output/Carbon_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  Silicon_outFile.open("./output/Silicon_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  Tube_outFile.open("./output/Tube_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  TPC_outFile.open("./output/TPC_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  Scint_layer_outFile.open("./output/Scintillator_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  Abs_outFile.open("./output/LeadGlass_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  
  Particle_outFile.open("./output/particle_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  pi0_outFile.open("./output/interaction_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  calorimeter_photon_outFile.open("./output/calorimeter_photon_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  
  #if VERSION_SHIELD==1      
    Shield_outFile.open("./output/Shield_output_"+std::to_string(theta_bin_index)+"_"+std::to_string(KE_bin_index)+"_"+std::to_string(run_number)+".txt");
  #endif
  
  Particle_outFile << "Event_ID,Run_ID,PID,Mass,Charge,KE,angle,x,y,z,t,u,v,w,weight"<< G4endl;
  Silicon_outFile << "Event_ID,layer,Track_ID,Parent_ID,Name,proc,step,Origin,x,y,z,u,v,w,t,KE,eDep,trackl" <<G4endl;
  Carbon_outFile << "Event_ID,Track_ID,Parent_ID,Name,proc,step,Origin,x,y,z,px,py,pz,t,KE,eDep" <<G4endl;
  Tube_outFile << "Event_ID,Track_ID,Parent_ID,Name,proc,volume,step,Origin,x,y,z,u,v,w,t,KE,eDep,trackl" <<G4endl;
  TPC_outFile << "Event_ID,module_ID,Layer,Track_ID,Parent_ID,Name,proc,step,Origin,x,y,z,t,KE,eDep,electrons,trackl"<< G4endl; //x,y,z,
  Scint_layer_outFile << "Event_ID,Track_ID,Parent_ID,Name,Proc,Volume,origin,step,module_ID,layer,index,t,KE,eDep,photons,x,y,z,particle_x,particle_y,particle_z"<<G4endl;
  Abs_outFile<< "Event_ID,Track_ID,Parent_ID,Name,Proc,Volume,origin,step,index,t,KE,eDep,trackl,photons,x,y,z"<<G4endl;
  pi0_outFile << "Event_ID,Track_ID,Parent_ID,Name,Proc,Volume,Origin,mass,KE,t,x,y,z,px,py,pz,weight" <<G4endl;
  calorimeter_photon_outFile << "Event_ID,Track_ID,Parent_ID,Name,Proc,Volume,Origin,Module,Layer,Index,t,KE,x,y,z,weight" <<G4endl;

  #if VERSION_SHIELD==1    
    Shield_outFile<<"Event_ID,Track_ID,Parent_ID,Name,proc,Volume,x,y,z,t,KE,px,py,pz,eDep,trackl"<<G4endl;
  #endif 

  }
}

void RunAction::EndOfRunAction(const G4Run*run)
{
  
  const NNbarRun* myrun = dynamic_cast<const NNbarRun*>(run);  

  if(IsMaster()){
  pi0_outFile.close();  
  calorimeter_photon_outFile.close();  
  Particle_outFile.close();
  Silicon_outFile.close();
  Carbon_outFile.close();
  Tube_outFile.close();
  TPC_outFile.close();
  Scint_layer_outFile.close();
  Abs_outFile.close();
  //PMT_outFile.close();
  #if VERSION_SHIELD==1    
  Shield_outFile.close();
  #endif

  // clear all stored data 
  G4double event_number = 0.0;
  particle_gun_record.clear(); PMT_record.clear(); scint_record.clear();
  }

  if(IsMaster()){run_number++;}
}

//....
