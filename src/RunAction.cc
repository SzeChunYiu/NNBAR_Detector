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
#include "G4Allocator.hh"
#include "NNbarHit.hh"
#include "G4MCPLGenerator.hh"
#include "mcpl.h"
#include <filesystem> // or #include <filesystem> for C++17 and up
#define BOOST_NO_SCOPED_ENUMS
#include <boost/filesystem.hpp>
#undef BOOST_NO_SCOPED_ENUMS
#include <arrow/io/file.h>
#include <parquet/stream_writer.h>
#include "parquet_writer.h"

//Simulation settings
extern G4int run_number; 
extern int particle_name_file_index;
extern G4double event_number;

extern int theta_bin_index;
extern int KE_bin_index;
extern int particle_name_input;
extern G4double event_number_global;

G4String folder_name = "";
string particle_name_list[7] = { "", "_neutron", "_proton", "_gamma","_electron","_muon","_pion" };

//MCPL file setting among runs
extern mcpl_file_t m_mcplfile;
extern const mcpl_particle_t * m_p;

//output file for each run
std::shared_ptr<arrow::io::FileOutputStream> Particle_outFile;
std::shared_ptr<arrow::io::FileOutputStream> Interaction_outFile;
std::shared_ptr<arrow::io::FileOutputStream> Carbon_outFile;
std::shared_ptr<arrow::io::FileOutputStream> Silicon_outFile;
std::shared_ptr<arrow::io::FileOutputStream> Beampipe_outFile;
std::shared_ptr<arrow::io::FileOutputStream> TPC_outFile;
std::shared_ptr<arrow::io::FileOutputStream> Scintillator_outFile;
std::shared_ptr<arrow::io::FileOutputStream> LeadGlass_outFile;

parquet::StreamWriter Particle_os;
parquet::StreamWriter Interaction_os;
parquet::StreamWriter Carbon_os;
parquet::StreamWriter Silicon_os;
parquet::StreamWriter Beampipe_os;
parquet::StreamWriter TPC_os;
parquet::StreamWriter Scintillator_os;
parquet::StreamWriter LeadGlass_os;

parquetwriter::Writer PMT_output;

extern G4ThreadLocal G4Allocator<NNbarHit>* NNbarHitAllocator;

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

  //Just need a control to the global event nummber  
	G4GenericMessenger::Command& Event_number_CMD = fMessenger->DeclareProperty("set_event_number",event_number_global,"Event number of the simulation");
	Event_number_CMD.SetParameterName("Event_number",true);
	Event_number_CMD.SetRange("Event_number>=0");
	Event_number_CMD.SetDefaultValue("0.");    

  //Just need a control to the output folder name  
	G4GenericMessenger::Command& Folder_CMD = fMessenger->DeclareMethod("set_folder_name",&RunAction::Set_output_folder,"Output folder name");
	Folder_CMD.SetParameterName("Folder_name",true);

  //Just need a control to the mcpl file name  
	G4GenericMessenger::Command& MCPL_File_CMD = fMessenger->DeclareMethod("set_mcpl_file", &RunAction::Set_MCPL_File, "Set and read the mcpl file");

  G4GenericMessenger::Command& Particle_name_CMD = fMessenger->DeclareProperty("set_particle",particle_name_input,"Please give particle name [Default: pi+]");
	Particle_name_CMD.SetParameterName("particle_name_input",true);   
	
  G4RunManager::GetRunManager()->SetPrintProgress(1);
}

RunAction::~RunAction()
{}

G4Run* RunAction::GenerateRun()
{return new NNbarRun;}

void RunAction::BeginOfRunAction(const G4Run*aRun)
{ 
  if(!IsMaster()){
  //std::cout << " = = = = = = = Memory = = = = = = = = = = " << std::endl;
  if (NNbarHitAllocator!=0){
    std::cout<<NNbarHitAllocator->GetAllocatedSize ()<<std::endl;
    NNbarHitAllocator->ResetStorage();
    //std::cout<< "After reset " << NNbarHitAllocator->GetAllocatedSize ()<<std::endl;
  }
  //std::cout << " = = = = = = = Memory check end = = = = = = = = = = " << std::endl;
  }
  
  if(IsMaster()){  // if this is master, then open the new files!
  
  G4String output_directory = ("./output/"+(folder_name));
  std::cout << "creating directory :: " << output_directory <<" :: Status : " << boost::filesystem::create_directories(output_directory) << std::endl;

    //Particle output file
    PARQUET_ASSIGN_OR_THROW(Particle_outFile, arrow::io::FileOutputStream::Open("./output/"+(folder_name)+"/Particle_output_"+std::to_string(run_number)+".parquet"));
    parquet::WriterProperties::Builder Particle_builder; Particle_builder.compression(parquet::Compression::SNAPPY); parquet::schema::NodeVector Particle_fields;
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("Event_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32,-1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("PID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("Mass", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("Name", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("Charge", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE,-1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("KE", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("angle", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE,-1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("x", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("t", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("u", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("v", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("w", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Particle_fields.push_back(parquet::schema::PrimitiveNode::Make("weight", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
    std::shared_ptr<parquet::schema::GroupNode> Particle_schema = std::static_pointer_cast<parquet::schema::GroupNode>(parquet::schema::GroupNode::Make("schema", parquet::Repetition::REQUIRED, Particle_fields));
    Particle_os = parquet::StreamWriter {parquet::ParquetFileWriter::Open(Particle_outFile, Particle_schema, Particle_builder.build())};
    
    //Interaction output file =>> Caused a seg fault at the end of the simulation (maybe because it is empty?)
    PARQUET_ASSIGN_OR_THROW(Interaction_outFile, arrow::io::FileOutputStream::Open("./output/"+(folder_name)+"/Interaction_output_"+std::to_string(run_number)+".parquet"));
    parquet::WriterProperties::Builder Interaction_builder; Interaction_builder.compression(parquet::Compression::SNAPPY); parquet::schema::NodeVector Interaction_fields;
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("Event_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32,-1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("Track_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("Parent_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("Name", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("Proc", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("Current_Vol", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("Origin", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("m", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("KE", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("t", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("x", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("px", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("py", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Interaction_fields.push_back(parquet::schema::PrimitiveNode::Make("pz", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
    std::shared_ptr<parquet::schema::GroupNode> Interaction_schema = std::static_pointer_cast<parquet::schema::GroupNode>(parquet::schema::GroupNode::Make("schema", parquet::Repetition::OPTIONAL, Interaction_fields));
    Interaction_os = parquet::StreamWriter {parquet::ParquetFileWriter::Open(Interaction_outFile, Interaction_schema, Interaction_builder.build())};
 
  // Carbon output file
    PARQUET_ASSIGN_OR_THROW(Carbon_outFile, arrow::io::FileOutputStream::Open("./output/"+(folder_name)+"/Carbon_output_"+std::to_string(run_number)+".parquet"));
    parquet::WriterProperties::Builder Carbon_SD_builder; Carbon_SD_builder.compression(parquet::Compression::SNAPPY);parquet::schema::NodeVector Carbon_fields;
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("Event_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32,-1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("Track_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("Parent_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("Name", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("Proc", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("Step_info", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("Origin", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("x", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("px", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("py", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("pz", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("t", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("KE", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Carbon_fields.push_back(parquet::schema::PrimitiveNode::Make("eDep", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
    std::shared_ptr<parquet::schema::GroupNode> Carbon_SD_schema = std::static_pointer_cast<parquet::schema::GroupNode>(parquet::schema::GroupNode::Make("schema", parquet::Repetition::REQUIRED, Carbon_fields));
    Carbon_os = parquet::StreamWriter {parquet::ParquetFileWriter::Open(Carbon_outFile, Carbon_SD_schema, Carbon_SD_builder.build())};

  // Silicon output file
    PARQUET_ASSIGN_OR_THROW(Silicon_outFile, arrow::io::FileOutputStream::Open("./output/"+(folder_name)+"/Silicon_output_"+std::to_string(run_number)+".parquet"));
    parquet::WriterProperties::Builder Silicon_SD_builder; Silicon_SD_builder.compression(parquet::Compression::SNAPPY); parquet::schema::NodeVector Silicon_fields;
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("Event_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32,-1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("Track_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("Parent_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("Name", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("Proc", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("Step_info", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("Origin", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("Layer_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("x", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("px", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("py", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("pz", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("t", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("KE", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Silicon_fields.push_back(parquet::schema::PrimitiveNode::Make("eDep", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
    std::shared_ptr<parquet::schema::GroupNode> Silicon_SD_schema = std::static_pointer_cast<parquet::schema::GroupNode>(parquet::schema::GroupNode::Make("schema", parquet::Repetition::REQUIRED, Silicon_fields));
    Silicon_os = parquet::StreamWriter {parquet::ParquetFileWriter::Open(Silicon_outFile, Silicon_SD_schema, Silicon_SD_builder.build())};
  
  // Beampipe output file
    PARQUET_ASSIGN_OR_THROW(Beampipe_outFile, arrow::io::FileOutputStream::Open("./output/"+(folder_name)+"/Beampipe_output_"+std::to_string(run_number)+".parquet"));
    parquet::WriterProperties::Builder Beampipe_SD_builder; Beampipe_SD_builder.compression(parquet::Compression::SNAPPY); parquet::schema::NodeVector Beampipe_fields;
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("Event_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32,-1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("Track_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("Parent_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("Name", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("Proc", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("Step_info", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("Current_Vol", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("Origin", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("x", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("px", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("py", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("pz", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("t", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("KE", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Beampipe_fields.push_back(parquet::schema::PrimitiveNode::Make("eDep", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
    std::shared_ptr<parquet::schema::GroupNode> Beampipe_SD_schema = std::static_pointer_cast<parquet::schema::GroupNode>(parquet::schema::GroupNode::Make("schema", parquet::Repetition::REQUIRED, Beampipe_fields));
    Beampipe_os = parquet::StreamWriter {parquet::ParquetFileWriter::Open(Beampipe_outFile, Beampipe_SD_schema, Beampipe_SD_builder.build())};
  
  // TPC output file
    PARQUET_ASSIGN_OR_THROW(TPC_outFile, arrow::io::FileOutputStream::Open("./output/"+(folder_name)+"/TPC_output_"+std::to_string(run_number)+".parquet"));
    parquet::WriterProperties::Builder TPC_SD_builder; TPC_SD_builder.compression(parquet::Compression::SNAPPY); parquet::schema::NodeVector TPC_fields;
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Event_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32,-1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Track_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Parent_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Name", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Proc", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Step_info", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Origin", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Current_Vol", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Module_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("Layer_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("x", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("px", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("py", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("pz", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("t", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("KE", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("eDep", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("trackl", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        TPC_fields.push_back(parquet::schema::PrimitiveNode::Make("electrons", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
    std::shared_ptr<parquet::schema::GroupNode> TPC_SD_schema = std::static_pointer_cast<parquet::schema::GroupNode>(parquet::schema::GroupNode::Make("schema", parquet::Repetition::REQUIRED, TPC_fields));
    TPC_os = parquet::StreamWriter {parquet::ParquetFileWriter::Open(TPC_outFile, TPC_SD_schema, TPC_SD_builder.build())};
  
  // Scintillator output file
    PARQUET_ASSIGN_OR_THROW(Scintillator_outFile, arrow::io::FileOutputStream::Open("./output/"+(folder_name)+"/Scintillator_output_"+std::to_string(run_number)+".parquet"));
    parquet::WriterProperties::Builder Scintillator_SD_builder; Scintillator_SD_builder.compression(parquet::Compression::SNAPPY); parquet::schema::NodeVector Scintillator_fields;
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Event_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Track_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Parent_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Name", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Proc", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Step_info", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Origin", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Volume", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Module_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Layer_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("Stave_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("x", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));    
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("particle_x", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));    
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("particle_y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("particle_z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("x_local", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));    
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("y_local", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("z_local", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("t", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("KE", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("eDep", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        Scintillator_fields.push_back(parquet::schema::PrimitiveNode::Make("photons", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
    std::shared_ptr<parquet::schema::GroupNode> Scintillator_SD_schema = std::static_pointer_cast<parquet::schema::GroupNode>(parquet::schema::GroupNode::Make("schema", parquet::Repetition::REQUIRED, Scintillator_fields));
    Scintillator_os = parquet::StreamWriter {parquet::ParquetFileWriter::Open(Scintillator_outFile, Scintillator_SD_schema, Scintillator_SD_builder.build())};
  
  // LeadGlass output file
    PARQUET_ASSIGN_OR_THROW(LeadGlass_outFile, arrow::io::FileOutputStream::Open("./output/"+(folder_name)+"/LeadGlass_output_"+std::to_string(run_number)+".parquet"));
    parquet::WriterProperties::Builder LeadGlass_SD_builder; LeadGlass_SD_builder.compression(parquet::Compression::SNAPPY); parquet::schema::NodeVector LeadGlass_fields;
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("Event_ID",parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));//
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("Track_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("Parent_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("Name", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("Proc", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("Step_info", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("Origin", parquet::Repetition::REQUIRED, parquet::Type::BYTE_ARRAY,parquet::ConvertedType::UTF8));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("Module_ID", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("x", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));    
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("y", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("z", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("t", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("KE", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("eDep", parquet::Repetition::REQUIRED, parquet::Type::DOUBLE,parquet::ConvertedType::NONE, -1));
        LeadGlass_fields.push_back(parquet::schema::PrimitiveNode::Make("photons", parquet::Repetition::REQUIRED, parquet::Type::INT32,parquet::ConvertedType::INT_32, -1));
    
    std::shared_ptr<parquet::schema::GroupNode> LeadGlass_SD_schema = std::static_pointer_cast<parquet::schema::GroupNode>(parquet::schema::GroupNode::Make("schema", parquet::Repetition::REQUIRED, LeadGlass_fields));
    LeadGlass_os = parquet::StreamWriter {parquet::ParquetFileWriter::Open(LeadGlass_outFile, LeadGlass_SD_schema, LeadGlass_SD_builder.build())};
  
    // PMT output file
        std::ifstream layout_file("./config/PMT_output_layout.json");
        PMT_output.set_layout(layout_file);
        PMT_output.set_dataset_name("./output/"+(folder_name)+"/PMT_output_"+std::to_string(run_number)); // must give a name to the output
        PMT_output.initialize();
  }
}

void RunAction::EndOfRunAction(const G4Run*run)
{
  
  const NNbarRun* myrun = dynamic_cast<const NNbarRun*>(run);  

  if(IsMaster()){
  // clear all stored data 
  G4double event_number = 0.0;
  }

  if(IsMaster()){
    run_number++;
    Particle_os.EndRowGroup();
    Interaction_os.EndRowGroup();
    Carbon_os.EndRowGroup(); 
    Silicon_os.EndRowGroup();
    Beampipe_os.EndRowGroup();
    TPC_os.EndRowGroup();
    Scintillator_os.EndRowGroup();
    LeadGlass_os.EndRowGroup();
    PMT_output.finish();
  }


}


void RunAction::Set_MCPL_File(G4String& filename){
  std::cout<< " Changing to mcpl file " << filename<<std::endl;
    m_mcplfile = mcpl_open_file(filename.c_str());
    m_p = NULL;
    m_p = mcpl_read(m_mcplfile);
}
//....

void RunAction::Set_output_folder(G4String& foldername_){
  if(IsMaster()){  // if this is master, then open the new files!
    folder_name = foldername_;
    G4String output_directory = ("./output/"+(folder_name));
    std::cout << "creating directory :: " << output_directory <<" :: Status : " << boost::filesystem::create_directories(output_directory) << std::endl;
  }
}