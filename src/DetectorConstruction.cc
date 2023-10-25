#include "DetectorConstruction.hh"
#include "Detector_Module/beampipe_geometry.hh"
#include "Detector_Module/Silicon_geometry.hh"
#include "Detector_Module/TPC_geometry.hh"
#include "Detector_Module/Scintillator_geometry.hh"
#include "Detector_Module/LeadGlass_geometry.hh"
#include "Detector_Module/Cosmic_Shielding_geometry.hh"
#include "Detector_Module/beampipe_shielding_geometry.hh"

#include "Sensitive_Detector/SiliconSD.hh"
#include "Sensitive_Detector/TPCSD.hh"
#include "Sensitive_Detector/TubeSD.hh"

#include "Sensitive_Detector/CarbonSD.hh"
#include "Sensitive_Detector/ScintillatorSD.hh"
#include "Sensitive_Detector/LeadGlassSD.hh"
#include "Sensitive_Detector/PMTSD.hh"

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
#include "ElectricField.hh"

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

using boost::lexical_cast;

std::vector<G4LogicalVolume*> Beampipe_output;
std::vector<G4LogicalVolume*> Silicon_output;
std::vector<G4LogicalVolume*> TPC_output;
std::vector<G4LogicalVolume*> Scintillator_output;
std::vector<G4LogicalVolume*> LeadGlass_output;
std::vector<G4LogicalVolume*> CosmicShielding_output;
std::vector<G4LogicalVolume*> Beampipe_Shielding_output;

DetectorConstruction::DetectorConstruction():G4VUserDetectorConstruction(),fCheckOverlaps(true) {}
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
    
  auto carbonPV = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),carbonLV,"CarbonPV",worldLV,false,0,fCheckOverlaps);

  //Silicon constrcution !!! NO INNER TRACKER FOR HIBEAM
  Silicon* Silicon; 
  auto Silicon_LV_list = Silicon -> Construct_Volumes(worldLV);
  Silicon_output = Silicon_LV_list;

  // Beampipe Construction
  Beampipe* beampipe;
  auto Beampipe_LV_list = beampipe->Construct_Volumes(worldLV);
  Beampipe_output = Beampipe_LV_list;
  
  // TPC Construction
  TPC* TPC; 
  TPC_output = TPC -> Construct_Volumes(worldLV);

  //Scintillator Construction
  Scintillator*Scintillator; 
  Scintillator_output = Scintillator-> Construct_Volumes(worldLV);
  
  // Lead glass Construction 
  LeadGlass*LeadGlass;
  LeadGlass_output = LeadGlass->Construct_Volumes(worldLV);

  /////////////////////////////////////////////////////////
  ////////////// Cosmic shielding construction ////////////
  /////////////////////////////////////////////////////////
  CosmicShielding*CosmicShielding;
  CosmicShielding_output = CosmicShielding->Construct_Volumes(worldLV);


  /////////////////////////////////////////////////////////
  ////////////// Beampipe shielding construction ////////////
  /////////////////////////////////////////////////////////
  Beampipe_Shielding*Beampipe_Shielding;
  Beampipe_Shielding_output = Beampipe_Shielding->Construct_Volumes(worldLV);


  return worldPV;
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
  for (int i=0;i<Silicon_output.size();i++){Silicon_output[i] -> SetSensitiveDetector(siliconDetector);}

  G4String tubeDetectorName = "BeampipeLV" ;
  TubeSD* tubeDetector = new TubeSD(tubeDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(tubeDetector);
  for (int i=0;i<Beampipe_output.size();i++){Beampipe_output[i] -> SetSensitiveDetector(tubeDetector);}

  G4String TPCDetectorName = "TPCLV";
  TPCSD* TPCDetector = new TPCSD(TPCDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(TPCDetector);
  for (int i=0;i<TPC_output.size();i++){TPC_output[i] -> SetSensitiveDetector(TPCDetector);}

  G4String scintDetectorName = "ScintLV" ;
  ScintillatorSD* scintDetector = new ScintillatorSD(scintDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(scintDetector);
  for (int i=0;i<Scintillator_output.size();i++){Scintillator_output[i] -> SetSensitiveDetector(scintDetector);}

  G4String LeadGlassDetectorName = "LeadGlassLV" ;
  LeadGlassSD* LeadGlassDetector = new LeadGlassSD(LeadGlassDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(LeadGlassDetector);
  LeadGlass_output[0] -> SetSensitiveDetector(LeadGlassDetector);
  
  G4String PMTDetectorName = "PMTLV" ;
  PMTSD* PMTDetector = new PMTSD(PMTDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(PMTDetector);
  LeadGlass_output[1] -> SetSensitiveDetector(PMTDetector);

  // Electric field 

  G4ElectricField* fEMfield;
  G4EqMagElectricField* fEquation;
  G4MagIntegratorStepper* fStepper;
  G4FieldManager* fFieldMgr;
  G4double fMinStep;
  G4ChordFinder* fChordFinder;

  auto fElectricField = new ElectricField();
	fEquation = new G4EqMagElectricField(fElectricField);
	G4int nvar = 8;
	fStepper = new G4DormandPrince745(fEquation,nvar);

	//global Feild Manager
	fFieldMgr = new G4FieldManager();
	fFieldMgr->SetDetectorField(fElectricField); //fEMfield<->fElectricField

	fMinStep = 1* mm; // minimal step of 10 microns //1*mm  0.0005?!

	G4FieldManager* globalFieldManager;
	G4TransportationManager* transportMgr = G4TransportationManager::GetTransportationManager();

  std::cout << "Electric field ...  " << std::endl;
	double MaxTrackingStep = 0.01;
	globalFieldManager = transportMgr->GetFieldManager();

	// Relative accuracy values:
	G4double minEps = 1 * mm;  //   Minimum & value for smallest steps  (1e-5)
	G4double maxEps = 100.0 *mm;  //   Maximum & value for largest steps  (1e-4)

	fFieldMgr->SetDeltaOneStep(1*mm);  // 0.5 micrometer

	G4MagInt_Driver* fIntgrDriver = new G4MagInt_Driver(fMinStep,fStepper,fStepper->GetNumberOfVariables());
	fChordFinder = new G4ChordFinder(fIntgrDriver);
	fFieldMgr->SetChordFinder(fChordFinder);
  TPC_output[0]->SetFieldManager(fFieldMgr, true);
  TPC_output[1]->SetFieldManager(fFieldMgr, true);	
  transportMgr->GetPropagatorInField()->SetLargestAcceptableStep(1.*cm);

}