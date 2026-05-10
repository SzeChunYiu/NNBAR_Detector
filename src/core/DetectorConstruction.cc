#include "core/DetectorConstruction.hh"
#include "config.h"
#include "detector/beampipe_geometry.hh"
#include "detector/Silicon_geometry.hh"
#include "detector/TPC_geometry.hh"
#include "detector/Scintillator_geometry.hh"
#include "detector/LeadGlass_geometry.hh"
#include "detector/Cosmic_Shielding_geometry.hh"
#include "detector/beampipe_shielding_geometry.hh"

#include "sensitive/SiliconSD.hh"
#include "sensitive/TPCSD.hh"
#include "sensitive/TubeSD.hh"

#include "sensitive/CarbonSD.hh"
#include "sensitive/ScintillatorSD.hh"
#include "sensitive/LeadGlassSD.hh"
#include "sensitive/PMTSD.hh"

#if WITH_CELERITAS
#include "physics/CeleritasCalorimeter.hh"
#endif

#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AutoDelete.hh"
#include "G4RotationMatrix.hh"

#include "G4UniformElectricField.hh"
#include "G4EqMagElectricField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4CashKarpRKF45.hh"
#include "G4DormandPrince745.hh"
#include "G4ChordFinder.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4PropagatorInField.hh"
#include "util/ElectricField.hh"

#include "G4SDManager.hh"
#include "G4SDChargedFilter.hh"
#include "G4PSPopulation.hh"
#include "G4SDParticleFilter.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSTrackLength.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PhysicalConstants.hh"

#include <G4ProductionCuts.hh>
#include "G4RegionStore.hh"
#include "G4Region.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include<string>
//#include "G4GDMLParser.hh"
#include <boost/lexical_cast.hpp>
#include "config.h"
#include "util/GeometryParameters.hh"
#include "util/GeometryManager.hh"

using boost::lexical_cast;

std::vector<G4LogicalVolume*> Beampipe_output;
std::vector<G4LogicalVolume*> Silicon_output;
std::vector<G4LogicalVolume*> TPC_output;
std::vector<G4LogicalVolume*> Scintillator_output;
std::vector<G4LogicalVolume*> LeadGlass_output;
std::vector<G4LogicalVolume*> CosmicShielding_output;
std::vector<G4LogicalVolume*> Beampipe_Shielding_output;

DetectorConstruction::DetectorConstruction():G4VUserDetectorConstruction(),fCheckOverlaps(false) {}
DetectorConstruction::~DetectorConstruction(){}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  DefineMaterials();
  return DefineVolumes();
}

void DetectorConstruction::DefineMaterials()
{ 
  
  G4double a; G4double z; G4double density;
  auto nistManager = G4NistManager::Instance();
  // Vacuum
  G4Material* vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole, density= universe_mean_density,kStateGas, 2.73*kelvin, 3.e-18*pascal);
  // -- optical properties of vacuum
  G4double vacuum_Energy[] = { 2.0 * eV,7.0 * eV, 7.14 * eV }; const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double); G4double vacuum_RIND[] = { 1.,1.,1. };
  G4MaterialPropertiesTable * vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum); vacuum->SetMaterialPropertiesTable(vacuum_mt);

  //Carbon Target
  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Material* Carbon_target = new G4Material("Carbon_target" , density=3.52*g/cm3, 1);
  Carbon_target->AddElement(elC, 1);

}

//....
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto carbonMaterial = G4Material::GetMaterial("Carbon_target");

  // World
  auto worldSizeXY = 20.0 * m; auto worldSizeZ = 450.0 * m; //1 * calorThickness; makes it longer because I don't want to shift the coordinates
  auto worldS = new G4Box("WorldS",worldSizeXY/2., worldSizeXY/2., worldSizeZ/2.);
  auto worldLV = new G4LogicalVolume(worldS,defaultMaterial,"World");
  G4VPhysicalVolume* worldPV = new G4PVPlacement(0,G4ThreeVector(),worldLV,"WorldPV",0,false,0,fCheckOverlaps);
  
  G4OpticalSurface* op_World= new G4OpticalSurface("LeadGlass");
  op_World->SetType(dielectric_metal);
  op_World->SetFinish(polished);
  op_World->SetModel(unified);
  // basically it acts effectively as a non-reflective coating [should work on the bottom surface of the Lead Glass]
  G4double pp_World[] = { 2.0 * eV, 3.5 * eV }; const G4int num_World = sizeof(pp_World) / sizeof(G4double);
  G4double reflectivity_World[] = { 0.0, 0.0}; G4double efficiency_World[] = { 0.0, 0.0 };
  G4MaterialPropertiesTable* LeadGlassProperty = new G4MaterialPropertiesTable();
  LeadGlassProperty->AddProperty("REFLECTIVITY", pp_World, reflectivity_World, num_World);
  LeadGlassProperty->AddProperty("EFFICIENCY", pp_World, efficiency_World, num_World);
  op_World ->SetMaterialPropertiesTable(LeadGlassProperty);
  new G4LogicalSkinSurface("name",worldLV,op_World);
  
  // Carbon
  G4double carbon_radius = 30*cm; G4double carbon_len = 0.01*cm; G4double carbon_angle = 360. * deg; //default, carbon 100 um
  auto carbonS = new G4Cons("CarbonS", 0.,carbon_radius,0.,carbon_radius,carbon_len,0.,carbon_angle);
  
    #if TARGET_BUILD==1
      auto carbonLV = new G4LogicalVolume(carbonS,carbonMaterial,"CarbonLV");
    #else
      auto carbonLV = new G4LogicalVolume(carbonS,defaultMaterial,"CarbonLV");
    #endif 
    
  // Place carbon target (unused PV pointer - placement still creates the volume)
  new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),carbonLV,"CarbonPV",worldLV,false,0,fCheckOverlaps);

  //Silicon construction
  auto siliconBuilder = new ::Silicon();
  Silicon_output = siliconBuilder->Construct_Volumes(worldLV);

  // Beampipe Construction
  auto beampipeBuilder = new ::Beampipe();
  Beampipe_output = beampipeBuilder->Construct_Volumes(worldLV);

  // TPC Construction
  auto tpcBuilder = new ::TPC();
  TPC_output = tpcBuilder->Construct_Volumes(worldLV);

  //Scintillator Construction
  auto scintillatorBuilder = new ::Scintillator();
  Scintillator_output = scintillatorBuilder->Construct_Volumes(worldLV);

  // Lead glass Construction
  auto leadGlassBuilder = new ::LeadGlass();
  LeadGlass_output = leadGlassBuilder->Construct_Volumes(worldLV);

  /////////////////////////////////////////////////////////
  ////////////// Cosmic shielding construction ////////////
  /////////////////////////////////////////////////////////
  auto cosmicShieldingBuilder = new ::CosmicShielding();
  CosmicShielding_output = cosmicShieldingBuilder->Construct_Volumes(worldLV);


  /////////////////////////////////////////////////////////
  ////////////// Beampipe shielding construction ////////////
  /////////////////////////////////////////////////////////
  auto beampipeShieldingBuilder = new ::Beampipe_Shielding();
  Beampipe_Shielding_output = beampipeShieldingBuilder->Construct_Volumes(worldLV);

  // ========================================================================
  // Initialize GeometryManager for volume lookup and visualization
  // ========================================================================
  nnbar::GeometryManager::Instance().Initialize();
  RegisterTPCGeometry();

  // ========================================================================
  // Register geometry parameters for EventDisplay (dynamic geometry sync)
  // ========================================================================
  RegisterGeometryParameters();

  return worldPV;
}

void DetectorConstruction::RegisterGeometryParameters()
{
  // Access global geometry values from beampipe_geometry.cc and TPC_geometry.cc
  extern G4double Beampipe_5_radius_1;
  extern G4double Beampipe_5_radius_2;
  extern G4double Beampipe_thickness;
  extern G4double Beampipe_5_len;
  extern G4double Beampipe_6_len;
  extern G4double TPC_drift_len;
  extern G4double TPC_wall_thickness;

  auto& geom = nnbar::GeometryParameters::Instance();

  // Beampipe (convert from Geant4 units to cm)
  geom.Set("beampipe.inner_radius", Beampipe_5_radius_1 / cm);
  geom.Set("beampipe.outer_radius", Beampipe_5_radius_2 / cm);
  geom.Set("beampipe.thickness", Beampipe_thickness / cm);
  geom.Set("beampipe.half_z", Beampipe_5_len / 2.0 / cm);

  // TPC common
  G4double TPC_z = (Beampipe_5_len + 2.0 * Beampipe_6_len) / 2.0;
  geom.Set("tpc.drift_length", TPC_drift_len / cm);
  geom.Set("tpc.half_z", TPC_z / 2.0 / cm);

  // TPC Type II (horizontal)
  G4double TPC_x_2 = 2.0 * Beampipe_5_radius_2 + 2.0 * TPC_wall_thickness;
  G4double TPC_y_2 = TPC_drift_len + 2.0 * TPC_wall_thickness;
  geom.Set("tpc.type2.width", TPC_x_2 / cm);
  geom.Set("tpc.type2.height", TPC_y_2 / cm);
  geom.Set("tpc.type2.y", (Beampipe_5_radius_2 + TPC_y_2 / 2.0) / cm);

  // TPC Type I (vertical)
  G4double TPC_x_1 = TPC_drift_len + 2.0 * TPC_wall_thickness;
  G4double TPC_y_1 = (2.0 * Beampipe_5_radius_2 + 2.0 * TPC_y_2) / 2.0;
  geom.Set("tpc.type1.width", TPC_x_1 / cm);
  geom.Set("tpc.type1.height", TPC_y_1 / cm);
  geom.Set("tpc.type1.x", (TPC_x_2 / 2.0 + TPC_x_1 / 2.0) / cm);
  geom.Set("tpc.type1.y", (TPC_y_1 / 2.0) / cm);

  G4cout << "========================================" << G4endl;
  G4cout << "Geometry registered for EventDisplay:" << G4endl;
  G4cout << "  Beampipe inner: " << geom.BeampipeInnerRadius() << " cm" << G4endl;
  G4cout << "  Beampipe outer: " << geom.BeampipeOuterRadius() << " cm" << G4endl;
  G4cout << "  TPC Type II: " << geom.TPCTypeIIWidth() << " x " << geom.TPCTypeIIHeight() << " cm" << G4endl;
  G4cout << "  TPC Type I:  " << geom.TPCTypeIWidth() << " x " << geom.TPCTypeIHeight() << " cm" << G4endl;
  G4cout << "========================================" << G4endl;
}

void DetectorConstruction::RegisterTPCGeometry()
{
  // Access global geometry values from beampipe_geometry.cc and TPC_geometry.cc
  extern G4double Beampipe_5_radius_2;
  extern G4double Beampipe_5_len;
  extern G4double Beampipe_6_len;
  extern G4double TPC_drift_len;
  extern G4double TPC_wall_thickness;

  auto& geoMgr = nnbar::GeometryManager::Instance();

  // Calculate TPC dimensions
  G4double TPC_z = (Beampipe_5_len + 2.0 * Beampipe_6_len) / 2.0;
  G4double TPC_y_2 = TPC_drift_len + 2.0 * TPC_wall_thickness;
  G4double TPC_x_2 = 2.0 * Beampipe_5_radius_2 + 2.0 * TPC_wall_thickness;
  G4double TPC_x_1 = TPC_drift_len + 2.0 * TPC_wall_thickness;
  G4double TPC_y_1 = (2.0 * Beampipe_5_radius_2 + 2.0 * TPC_y_2) / 2.0;

  // Module positions (from TPC_geometry.cc)
  G4double pos_y_type2 = Beampipe_5_radius_2 + TPC_y_2 / 2.0;
  G4double pos_x_type1 = TPC_x_2 / 2.0 + TPC_x_1 / 2.0;
  G4double pos_y_type1 = TPC_y_1 / 2.0;
  G4double pos_z_front = -TPC_z / 2.0;
  G4double pos_z_back = TPC_z / 2.0;

  // Register all 12 TPC modules
  // Front modules (z < 0)
  geoMgr.RegisterTPCModule(0, 2, 0, pos_y_type2, pos_z_front,
                           TPC_x_2/mm, TPC_y_2/mm, TPC_z/mm, 1, -1);  // Type II, drift -Y
  geoMgr.RegisterTPCModule(1, 1, -pos_x_type1, pos_y_type1, pos_z_front,
                           TPC_x_1/mm, TPC_y_1/mm, TPC_z/mm, 0, 1);   // Type I, drift +X
  geoMgr.RegisterTPCModule(2, 1, -pos_x_type1, -pos_y_type1, pos_z_front,
                           TPC_x_1/mm, TPC_y_1/mm, TPC_z/mm, 0, 1);   // Type I, drift +X
  geoMgr.RegisterTPCModule(3, 2, 0, -pos_y_type2, pos_z_front,
                           TPC_x_2/mm, TPC_y_2/mm, TPC_z/mm, 1, 1);   // Type II, drift +Y
  geoMgr.RegisterTPCModule(4, 1, pos_x_type1, -pos_y_type1, pos_z_front,
                           TPC_x_1/mm, TPC_y_1/mm, TPC_z/mm, 0, -1);  // Type I, drift -X
  geoMgr.RegisterTPCModule(5, 1, pos_x_type1, pos_y_type1, pos_z_front,
                           TPC_x_1/mm, TPC_y_1/mm, TPC_z/mm, 0, -1);  // Type I, drift -X

  // Back modules (z > 0)
  geoMgr.RegisterTPCModule(6, 2, 0, pos_y_type2, pos_z_back,
                           TPC_x_2/mm, TPC_y_2/mm, TPC_z/mm, 1, -1);  // Type II, drift -Y
  geoMgr.RegisterTPCModule(7, 1, pos_x_type1, pos_y_type1, pos_z_back,
                           TPC_x_1/mm, TPC_y_1/mm, TPC_z/mm, 0, -1);  // Type I, drift -X
  geoMgr.RegisterTPCModule(8, 1, pos_x_type1, -pos_y_type1, pos_z_back,
                           TPC_x_1/mm, TPC_y_1/mm, TPC_z/mm, 0, -1);  // Type I, drift -X
  geoMgr.RegisterTPCModule(9, 2, 0, -pos_y_type2, pos_z_back,
                           TPC_x_2/mm, TPC_y_2/mm, TPC_z/mm, 1, 1);   // Type II, drift +Y
  geoMgr.RegisterTPCModule(10, 1, -pos_x_type1, -pos_y_type1, pos_z_back,
                           TPC_x_1/mm, TPC_y_1/mm, TPC_z/mm, 0, 1);   // Type I, drift +X
  geoMgr.RegisterTPCModule(11, 1, -pos_x_type1, pos_y_type1, pos_z_back,
                           TPC_x_1/mm, TPC_y_1/mm, TPC_z/mm, 0, 1);   // Type I, drift +X

  G4cout << "GeometryManager: Registered " << geoMgr.GetNumTPCModules() << " TPC modules" << G4endl;
}

void DetectorConstruction::ConstructSDandField()
{ 
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  G4String Carbon_DetectorName = "CarbonLV";
  CarbonSD* Carbon_Detector = new CarbonSD(Carbon_DetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(Carbon_Detector);
  SetSensitiveDetector("CarbonLV", Carbon_Detector);

  G4String siliconDetectorName = "SiliconLV" ;
  SiliconSD* siliconDetector = new SiliconSD(siliconDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(siliconDetector);
  for (size_t i=0;i<Silicon_output.size();i++){Silicon_output[i] -> SetSensitiveDetector(siliconDetector);}

  G4String tubeDetectorName = "BeampipeLV" ;
  TubeSD* tubeDetector = new TubeSD(tubeDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(tubeDetector);
  for (size_t i=0;i<Beampipe_output.size();i++){Beampipe_output[i] -> SetSensitiveDetector(tubeDetector);}

  G4String TPCDetectorName = "TPCLV";
  TPCSD* TPCDetector = new TPCSD(TPCDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(TPCDetector);
  for (size_t i=0;i<TPC_output.size();i++){TPC_output[i] -> SetSensitiveDetector(TPCDetector);}

  G4String scintDetectorName = "ScintLV" ;
  ScintillatorSD* scintDetector = new ScintillatorSD(scintDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(scintDetector);
  for (size_t i=0;i<Scintillator_output.size();i++){Scintillator_output[i] -> SetSensitiveDetector(scintDetector);}

  G4String LeadGlassDetectorName = "LeadGlassLV" ;
  LeadGlassSD* LeadGlassDetector = new LeadGlassSD(LeadGlassDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(LeadGlassDetector);
  LeadGlass_output[0] -> SetSensitiveDetector(LeadGlassDetector);
  
  G4String PMTDetectorName = "PMTLV" ;
  PMTSD* PMTDetector = new PMTSD(PMTDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(PMTDetector);
  LeadGlass_output[1] -> SetSensitiveDetector(PMTDetector);

#if WITH_CELERITAS
  // Initialize Celeritas calorimeters for GPU EM energy recording
  // This records energy deposits from GPU-tracked e-, e+, gamma in LeadGlass and Scintillator
  G4LogicalVolume* lgLV = LeadGlass_output.size() > 0 ? LeadGlass_output[0] : nullptr;
  G4LogicalVolume* scintLV = Scintillator_output.size() > 0 ? Scintillator_output[0] : nullptr;
  nnbar::CeleritasCalorimeter::Instance().Initialize(lgLV, scintLV);
  G4cout << "[DetectorConstruction] Celeritas calorimeters initialized for GPU EM recording" << G4endl;
#endif

  // Electric field
  G4EqMagElectricField* fEquation;
  G4MagIntegratorStepper* fStepper;
  G4FieldManager* fFieldMgr;
  G4double fMinStep;
  G4ChordFinder* fChordFinder;

  auto fElectricField = new ElectricField();
	fEquation = new G4EqMagElectricField(fElectricField);
	G4int nvar = 8;
	fStepper = new G4DormandPrince745(fEquation,nvar);

	// Field Manager
	fFieldMgr = new G4FieldManager();
	fFieldMgr->SetDetectorField(fElectricField);

	fMinStep = 1* mm; // minimal step

	G4TransportationManager* transportMgr = G4TransportationManager::GetTransportationManager();

  std::cout << "Electric field ...  " << std::endl;

	fFieldMgr->SetDeltaOneStep(1*mm);

	G4MagInt_Driver* fIntgrDriver = new G4MagInt_Driver(fMinStep,fStepper,fStepper->GetNumberOfVariables());
	fChordFinder = new G4ChordFinder(fIntgrDriver);
	fFieldMgr->SetChordFinder(fChordFinder);
  TPC_output[0]->SetFieldManager(fFieldMgr, true);
  TPC_output[1]->SetFieldManager(fFieldMgr, true);	
  transportMgr->GetPropagatorInField()->SetLargestAcceptableStep(1.*cm);

}