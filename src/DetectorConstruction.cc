#include "DetectorConstruction.hh"
#include "ScintillatorSD.hh"
#include "AbsorberSD.hh"
#include "TubeSD.hh"
#include "TPCSD.hh"
#include "PMTSD.hh"
#include "SiliconSD.hh"
#include "ShieldSD.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
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

G4ElectricField* fEMfield;
G4EqMagElectricField* fEquation;
G4MagIntegratorStepper* fStepper;
G4FieldManager* fFieldMgr;
G4double fMinStep;
G4ChordFinder* fChordFinder;

//for the file output
extern std::ofstream Lead_glass_outFile;
extern std::ofstream Scint_outFile;

// read the position of the lead glass blocks stored in csv file
std::string filename_data_lead_glass_pos = "./lead_glass_position/lead_glass_position.csv";
std::string filename_data_lead_glass_pos_fb = "./lead_glass_position/lead_glass_position_fb.csv";
std::string filename_data_pmt_pos = "./lead_glass_position/pmt_pos.csv";
std::vector<std::vector<double>> data_lead_glass_pos;
std::vector<std::vector<double>> data_lead_glass_pos_fb;
std::vector<std::vector<double>> data_pmt;

void import_lead_glass_pos(std::string file_name, std::vector<std::vector<double> >& data) {
	
	std::string row;
	std::ifstream init_file(file_name.c_str());

	// open file
	if (init_file.is_open()) {
		std::cerr << "Opening Position file : "<< file_name << " ... " << std::endl;
		// loop in each line in the file
		int count_line = 0;
		while (getline(init_file, row)) {
			count_line++;
			std::istringstream iss(row);
			// initialize a vector to store the row elements
			std::vector<double> row;
			std::string token;
			// get each element by splitting the row string by commas
			while (std::getline(iss, token, ',')) {
				// convert the string to int (or use stof for floats)
				row.push_back(boost::lexical_cast<double>(token.c_str()));
			}
			data.push_back(row);
			//std::cerr << "Reading the " << count_line << " th line in file " << file_name <<" ... "<< std::endl;
		}
		init_file.close();
		std::cerr << "Lead_glass position loaded " << std::endl;
	}
	else
		std::cerr << "ERROR: Unable to open file" << std::endl;
	return;
}

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction(), fCheckOverlaps(true)
{
  import_lead_glass_pos(filename_data_lead_glass_pos, data_lead_glass_pos);
  import_lead_glass_pos(filename_data_lead_glass_pos_fb, data_lead_glass_pos_fb);
  import_lead_glass_pos(filename_data_pmt_pos, data_pmt);
}
//....
DetectorConstruction::~DetectorConstruction()
{}
//....
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  // Define volumes
  return DefineVolumes();
}
//....
void DetectorConstruction::DefineMaterials()
{ 

  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  G4double a; G4double z; G4double density;

  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elTi = nistManager->FindOrBuildElement("Ti");
  G4Element* elAs = nistManager->FindOrBuildElement("As");
  G4Element* elPb = nistManager->FindOrBuildElement("Pb");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elSi = nistManager->FindOrBuildElement("Si");
  G4Element* elN = nistManager->FindOrBuildElement("N");
  G4Element* elAr = nistManager->FindOrBuildElement("Ar");
  G4Element* elAl = nistManager->FindOrBuildElement("Al");
  G4Element* elMg = nistManager->FindOrBuildElement("Mg");
  G4Element* elMn = nistManager->FindOrBuildElement("Mn");
  G4Element* elNi = nistManager->FindOrBuildElement("Ni");
  G4Element* elCr = nistManager->FindOrBuildElement("Cr");
  G4Element* elF = nistManager->FindOrBuildElement("F");
  G4Element* elFe = nistManager->FindOrBuildElement("Fe");

  // Vacuum
  G4Material* vacuum = new G4Material("Galactic", z=1., a=1.01*g/mole, density= universe_mean_density,kStateGas, 2.73*kelvin, 3.e-18*pascal);
  // -- optical properties of vacuum
  G4double vacuum_Energy[] = { 2.0 * eV,7.0 * eV, 7.14 * eV }; const G4int vacnum = sizeof(vacuum_Energy) / sizeof(G4double); G4double vacuum_RIND[] = { 1.,1.,1. };
  G4MaterialPropertiesTable * vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND, vacnum); vacuum->SetMaterialPropertiesTable(vacuum_mt);
  
  // Tube
  G4Material* Al = new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);

  G4Material* StainlessSteel = new G4Material("StainlessSteel",   density = 8.02*g/cm3,5);
  StainlessSteel->AddElement(elMn, 0.02);
  StainlessSteel->AddElement(elSi, 0.01);
  StainlessSteel->AddElement(elCr, 0.19);
  StainlessSteel->AddElement(elNi, 0.10);
  StainlessSteel->AddElement(elFe, 0.68);

  // Silicon
  G4Material* Silicon = new G4Material("Silicon", z=14., a= 28.0855*g/mole, density = 2.33*g/cm3);

  // --------Air
  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, 2);
  Air->AddElement(elN, 0.7); Air->AddElement(elO, 0.3);

  // ---------FR4
  //----- Epoxy
  G4Material* Epoxy = new G4Material("Epoxy" , density=1.2*g/cm3, 2);
  Epoxy->AddElement(elH, 2); Epoxy->AddElement(elC, 2); 
  //----- SiO2 (Quarz)
  G4Material* SiO2 = new G4Material("SiO2",density= 2.200*g/cm3, 2);
  SiO2->AddElement(elSi, 1); SiO2->AddElement(elO , 2); 
  //FR4 (Glass + Epoxy)
  G4Material* FR4 = new G4Material("FR4" , density=1.86*g/cm3, 2);
  FR4->AddMaterial(Epoxy, 0.472); FR4->AddMaterial(SiO2, 0.528);

  //Carbon Target
  G4Material* Carbon_target = new G4Material("Carbon_target" , density=3.52*g/cm3, 1);
  Carbon_target->AddElement(elC, 1);

  // ----------TPC
  // CO2
  G4Material* CO2 = new G4Material("CO2", density= 1.98*g/cm3, 2);
  CO2->AddElement(elO, 2); CO2->AddElement(elC, 1);
  // Ar/CO2 80/20
  G4Material* Gas = new G4Material("Gas", density=1.823*mg/cm3, 2);
  Gas->AddElement(elAr, .8); Gas->AddMaterial(CO2, .2);

  // BC-408 taken from datasheet
  G4Material* Scint = new G4Material("Scint", 1.023*g/cm3, 2);
  Scint->AddElement(elH, 0.524573); Scint->AddElement(elC, 1 - 0.524573);

  // Lead-glass (taken from PDG)
  G4Material* Abs = new G4Material("Abs", 6.22*g/cm3, 5);
  Abs->AddElement(elO, 0.156453); Abs->AddElement(elSi, 0.080866); Abs->AddElement(elTi, 0.008092); Abs->AddElement(elAs, .002651); Abs->AddElement(elPb, 0.751938);

  // MgF2 coating of the lead glass
  G4Material *AlMgF2 = new G4Material("AlMgF2", density = 2.9007*g/cm3, 3);
  AlMgF2->AddElement(elAl, 0.331); AlMgF2->AddElement(elF, 0.408); AlMgF2->AddElement(elMg, 0.261);
  G4Material* Lead = new G4Material("Lead", z=82., a= 207.2*g/mole, density = 11.29*g/cm3);


  // PMT window quartz actually
  G4Material* PMT_window = new G4Material("PMT_window_mat",density= 2.200*g/cm3, 2);
  PMT_window->AddElement(elSi, 1); PMT_window->AddElement(elO , 2); 
  // -- optical properties of PMT
  G4double pmt_window_Energy[] = { 2.0 * eV, 7.0 * eV}; const G4int pmt_window_num = sizeof(pmt_window_Energy) / sizeof(G4double); G4double pmt_window_RIND[] = { 1.53,1.53 };
  G4MaterialPropertiesTable * pmt_window_mt = new G4MaterialPropertiesTable();
  pmt_window_mt->AddProperty("RINDEX", pmt_window_Energy, pmt_window_RIND, pmt_window_num); PMT_window->SetMaterialPropertiesTable(pmt_window_mt);
 
  
 // ----------------- Generate and Add Material Properties Table ----------------
 G4double PhotonWavelength[] =
      { 2325.4*nm, 1970.1*nm, 1529.6*nm, 1060.0*nm,
        1014.0*nm, 852.10*nm, 706.50*nm, 656.30*nm,
        643.80*nm, 632.80*nm, 589.30*nm, 587.60*nm,
        546.10*nm, 486.10*nm, 480.00*nm, 435.80*nm,
        404.70*nm, 365.00*nm};

 const G4int nEntries = sizeof(PhotonWavelength)/sizeof(G4double);
 G4double PhotonEnergy[nEntries];
 for (int i=0; i < nEntries; ++i) {PhotonEnergy[i] = (1240.*nm/PhotonWavelength[i])*eV;};

 //// Lead Glass Schott SF5
 G4double refractiveIndex[] =
        { 1.63289, 1.63785, 1.64359, 1.65104,
          1.65206, 1.65664, 1.66327, 1.66661,
          1.66756, 1.66846, 1.67252, 1.67270,
          1.67764, 1.68750, 1.68876, 1.69986,
          1.71069, 1.73056};

  G4MaterialPropertiesTable* absMPT = new G4MaterialPropertiesTable();
  absMPT->AddProperty("RINDEX", PhotonEnergy, refractiveIndex, nEntries)->SetSpline(true);
  Abs->SetMaterialPropertiesTable(absMPT);
  //G4cout << "Absorber Properties -------" << G4endl;
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  absMPT->DumpTable();

//Scintillator Optical Properties

  const G4int nEntries2 = 12;

  G4double ScintPhotonEnergy[nEntries2] =

  { 2.08*eV, 2.38*eV, 2.58*eV, 2.7*eV, 2.76*eV,
    2.82*eV, 2.92*eV, 2.95*eV, 3.02*eV, 3.1*eV,
    3.26*eV, 3.44*eV};

  G4double rindex_scint[nEntries2] =
    {1.58, 1.58, 1.58, 1.58, 1.58,
     1.58, 1.58, 1.58, 1.58, 1.58,
     1.58, 1.58};

  G4double atten_scint[nEntries2] =
    {210*cm, 210*cm, 210*cm, 210*cm, 210*cm,
     210*cm, 210*cm, 210*cm, 210*cm, 210*cm,
     210*cm, 210*cm};

  G4double scintilFast[nEntries2] = {.04, .07, .20, .49, .84, 1.0, .83, .55, .40, .17,.03, 0.};
  G4double scintilSlow[nEntries2] = {.04, .07, .20, .49, .84, 1.0, .83, .55, .40, .17,.03, 0.};

  G4MaterialPropertiesTable *scintMPT = new G4MaterialPropertiesTable();
  scintMPT->AddProperty("RINDEX", ScintPhotonEnergy, rindex_scint, nEntries2)->SetSpline(true);
  scintMPT->AddProperty("ABSLENGTH", ScintPhotonEnergy, atten_scint, nEntries2)->SetSpline(true);
  scintMPT->AddProperty("FASTCOMPONENT", ScintPhotonEnergy, scintilFast, nEntries2)->SetSpline(true);
  scintMPT->AddProperty("SLOWCOMPONENT", ScintPhotonEnergy, scintilSlow, nEntries2)->SetSpline(true);
  // 64% of Antracene: 17400
  scintMPT->AddConstProperty("SCINTILLATIONYIELD", 100./ MeV); //original 11136000.
  scintMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  scintMPT->AddConstProperty("FASTTIMECONSTANT", 1.0*ns); // org: 0.9
  scintMPT->AddConstProperty("SLOWTIMECONSTANT", 1.0*ns); // org: 2.1
  scintMPT->AddConstProperty("YIELDRATIO", 1.);
  Scint->SetMaterialPropertiesTable(scintMPT);
  scintMPT->DumpTable();
}

//....
G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  auto absorberMaterial = G4Material::GetMaterial("Abs");
  auto scintMaterial = G4Material::GetMaterial("Scint");
  auto tubeMaterial = G4Material::GetMaterial("Aluminum"); // or stainless steel // StainlessSteel,Aluminum
  auto FR4Material = G4Material::GetMaterial("FR4");
  auto TPCMaterial = G4Material::GetMaterial("Gas");
  auto SiliconMaterial = G4Material::GetMaterial("Silicon");
  auto pmtMaterial = G4Material::GetMaterial("PMT_window_mat");
  auto leadglasscoatMaterial =  G4Material::GetMaterial("AlMgF2");
  auto shieldMaterial = G4Material::GetMaterial("Lead");
  auto carbonMaterial = G4Material::GetMaterial("Carbon_target");

  if ( ! defaultMaterial || ! absorberMaterial || ! scintMaterial) {
    G4ExceptionDescription msg; msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()","MyCode0001", FatalException, msg);}

  // World
  auto worldSizeXY = 10.0 * m; auto worldSizeZ = 10.0 * m; //1 * calorThickness;
  auto worldS = new G4Box("WorldS",worldSizeXY/2., worldSizeXY/2., worldSizeZ/2.);
  auto worldLV = new G4LogicalVolume(worldS,defaultMaterial,"World");
  auto worldPV = new G4PVPlacement(0,G4ThreeVector(),worldLV,"WorldPV",0,false,0,fCheckOverlaps);

  G4double carbon_radius = 60.0*cm; G4double carbon_len = 0.01*cm; G4double carbon_angle = 360. * deg;
  auto carbonS = new G4Cons("CarbonS", 0.,carbon_radius,0.,carbon_radius,carbon_len,0.,carbon_angle);
  auto carbonLV = new G4LogicalVolume(carbonS,carbonMaterial,"CarbonLV");
  auto carbonPV = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),carbonLV,"CarbonPV",worldLV,false,0,fCheckOverlaps);


  // Silicon
  G4double silicon_radius_1 = 80.0*cm;//87.97*cm; 
  G4double silicon_radius_2 = 90.0*cm;//97.97*cm;
  G4double silicon_len = 3.0*m; G4double silicon_thickness = 0.03*cm; G4double silicon_angle = 360. * deg;
  auto siliconS_1 = new G4Cons("siliconS_1", silicon_radius_1,silicon_radius_1+silicon_thickness,silicon_radius_1,silicon_radius_1+silicon_thickness,silicon_len,0.,silicon_angle);
  auto siliconLV_1 = new G4LogicalVolume(siliconS_1,SiliconMaterial,"siliconLV_1");
  auto siliconPV_1 = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),siliconLV_1,"siliconPV_1",worldLV,false,0,fCheckOverlaps);
  
  auto siliconS_2 = new G4Cons("siliconS_2", silicon_radius_2,silicon_radius_2+silicon_thickness,silicon_radius_2,silicon_radius_2+silicon_thickness,silicon_len,0.,silicon_angle);
  auto siliconLV_2 = new G4LogicalVolume(siliconS_2,SiliconMaterial,"siliconLV_2");
  auto siliconPV_2 = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),siliconLV_2,"siliconPV_2",worldLV,false,1,fCheckOverlaps);

  // Aluminum Tube
  G4double tube_radius = 1.0*m; G4double tube_len = 3.0*m; G4double tube_thickness = 2.0*cm; G4double tube_angle = 360. * deg;
  auto tubeS = new G4Cons("TubeS", tube_radius-tube_thickness,tube_radius,tube_radius-tube_thickness,tube_radius,tube_len,0.,tube_angle);
  auto tubeLV = new G4LogicalVolume(tubeS,tubeMaterial,"TubeLV");
  auto tubePV = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),tubeLV,"TubePV",worldLV,false,0,fCheckOverlaps);
  
  G4Region* tube_region = new G4Region("Tube_region");
  tube_region->AddRootLogicalVolume(tubeLV); 
  G4Region* aRegion = G4RegionStore::GetInstance()->GetRegion("Tube_region");
  G4ProductionCuts* tubecut = new G4ProductionCuts();
  tubecut->SetProductionCut(100*cm,"gamma");
  tubecut->SetProductionCut(1000*cm,"e-");
  tubecut->SetProductionCut(1000*cm,"e+");
  tubecut->SetProductionCut(10000*cm,"proton");
  aRegion->SetProductionCuts(tubecut);
  
  // TPC
  G4double FR4Thickness_1 = 0.16*cm; G4double FR4Thickness_2 = 0.02*cm;
  G4double gasThickness_1 = 2.0 *cm; G4double gasThickness_2 = 50.*cm;
  G4double airThickness = 1.34 *cm;
  G4double TPCThickness = FR4Thickness_1*2.+gasThickness_2+airThickness;

  // --- TPC container
  G4double TPC_x_1 = 0.85*m; G4double TPC_y_1 = 1.87*m; 
  G4double TPC_x_2 = 2.04*m; G4double TPC_y_2 = 0.85*m; G4double TPC_z = 2.0*m;
  G4double TPC_total_length = 2.0 * TPC_x_1 + TPC_x_2; // beware of this value, may affect the lead glass and scintillators

  auto TPCS_1 = new G4Box("TPCS_1",TPC_x_1/2.,TPC_y_1/2.,TPC_z/2.);
  auto TPCS_2 = new G4Box("TPCS_2",TPC_x_2/2.,TPC_y_2/2.,TPC_z/2.);
  
  G4double TPC_layer_x = 1.0*cm; G4double TPC_layer_y = 1.0*cm; G4double TPC_layer_z = 1.0*cm; 

  auto TPCS_layer = new G4Box("TPCS_layer",TPC_layer_x/2.,TPC_layer_y/2.,TPC_z/2.); //TPC long bar
  auto TPCS_block = new G4Box("TPCS_block",TPC_layer_x/2.,TPC_layer_y/2.,TPC_layer_z/2.); //TPC blocks along z dir in the long bar

  auto TPCLV_1 = new G4LogicalVolume(TPCS_1,defaultMaterial,"TPCLV_1");
  auto TPCLV_2 = new G4LogicalVolume(TPCS_2,defaultMaterial,"TPCLV_2");
  TPCLV = new G4LogicalVolume(TPCS_layer,defaultMaterial,"TPCLV_bar");
  auto TPCLV_block = new G4LogicalVolume(TPCS_block,TPCMaterial,"TPCLV");

  G4Region* TPC_region = new G4Region("TPC_region");
  TPCLV->SetRegion(TPC_region);
  TPC_region->AddRootLogicalVolume(TPCLV_block);
  G4Region* TPC_Region_ = G4RegionStore::GetInstance()->GetRegion("TPC_region");
  G4ProductionCuts* TPCcut = new G4ProductionCuts();
  TPCcut->SetProductionCut(0.0001*cm,"gamma");
  TPCcut->SetProductionCut(0.00001*nm,"e-"); // 130cm -> 1.3 MeV   220cm -> 3.4 MeV 400*cm -> 50 MeV
  TPC_Region_->SetProductionCuts(TPCcut);


  // installing the small TPC blocks in the TPC long bars
  int TPC_z_index = 0;

  int TPC_number_z = TPC_z / TPC_layer_z;
  for (int j = 0; j<TPC_number_z;j++){
    G4ThreeVector TPC_block_pos = G4ThreeVector(0.,0.,-TPC_z/2.+(2.*j+1.)/2.*TPC_layer_z); 
    new G4PVPlacement(0,TPC_block_pos,TPCLV_block,"TPCPV_blocks",TPCLV,false,TPC_z_index,fCheckOverlaps);
    std::cout << j << std::endl;
    TPC_z_index++;
  }
  
  // installing small TPC virtual volumes for type I TPC
  int TPC_index = 0;
  int TPC_number_x1 = TPC_x_1 / TPC_layer_x;
  int TPC_number_y1 = TPC_y_1 / TPC_layer_y;
  for (int j = 0; j<TPC_number_y1;j++){ 
    for (int i = 0; i<TPC_number_x1; i++){ 
      G4ThreeVector TPC_layer_pos = G4ThreeVector(-TPC_x_1/2.+(2.0*i+1.0)/2.*TPC_layer_x,-TPC_y_1/2.+(2.*j+1.)/2.*TPC_layer_y,0.); 
      new G4PVPlacement(0,TPC_layer_pos,TPCLV,"TPCPV_layer",TPCLV_1,false,TPC_index,false);
      TPC_index++;
    }
  }

  // installing small TPC virtual volumes for type II TPC
  TPC_index = 0;
  int TPC_number_x2 = TPC_x_2 / TPC_layer_x;
  int TPC_number_y2 = TPC_y_2 / TPC_layer_y;
  for (int j = 0; j<TPC_number_y2;j++){
    for (int i = 0; i<TPC_number_x2; i++){
      G4ThreeVector TPC_layer_pos = G4ThreeVector(-TPC_x_2/2.+(2.0*i+1.0)/2.*TPC_layer_x,-TPC_y_2/2.+(2.*j+1.)/2.*TPC_layer_y,0.); 
      new G4PVPlacement(0,TPC_layer_pos,TPCLV,"TPCPV_layer",TPCLV_2,false,TPC_index,false);
      TPC_index++;
    }
  }

  //------ front TPCs
  auto TPC_pos1 = G4ThreeVector(0.,tube_radius+tube_thickness+TPC_y_2/2.,-TPC_z/2.);
  auto TPC_pos2 = G4ThreeVector(-(tube_radius+tube_thickness+TPC_y_2/2.),TPC_y_1/2.,-TPC_z/2.);
  auto TPC_pos3 = G4ThreeVector(-(tube_radius+tube_thickness+TPC_y_2/2.),-TPC_y_1/2.,-TPC_z/2.);  
  auto TPC_pos4 = G4ThreeVector(0.,-(tube_radius+tube_thickness+TPC_y_2/2.),-TPC_z/2.);
  auto TPC_pos5 = G4ThreeVector((tube_radius+tube_thickness+TPC_y_2/2.),-TPC_y_1/2.,-TPC_z/2.);
  auto TPC_pos6 = G4ThreeVector((tube_radius+tube_thickness+TPC_y_2/2.),TPC_y_1/2.,-TPC_z/2.);

  //------ back TPCs
  auto TPC_pos7 = G4ThreeVector(0.,tube_radius+tube_thickness+TPC_y_2/2.,TPC_z/2.);
  auto TPC_pos8 = G4ThreeVector((tube_radius+tube_thickness+TPC_y_2/2.),TPC_y_1/2.,TPC_z/2.);
  auto TPC_pos9 = G4ThreeVector((tube_radius+tube_thickness+TPC_y_2/2.),-TPC_y_1/2.,TPC_z/2.);
  auto TPC_pos10 = G4ThreeVector(0.,-(tube_radius+tube_thickness+TPC_y_2/2.),TPC_z/2.);
  auto TPC_pos11 = G4ThreeVector(-(tube_radius+tube_thickness+TPC_y_2/2.),-TPC_y_1/2.,TPC_z/2.);
  auto TPC_pos12 = G4ThreeVector(-(tube_radius+tube_thickness+TPC_y_2/2.),TPC_y_1/2.,TPC_z/2.);
  

  new G4PVPlacement(0,TPC_pos1,TPCLV_2,"TPCPV1",worldLV,false,0,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos2,TPCLV_1,"TPCPV2",worldLV,false,1,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos3,TPCLV_1,"TPCPV3",worldLV,false,2,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos4,TPCLV_2,"TPCPV4",worldLV,false,3,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos5,TPCLV_1,"TPCPV5",worldLV,false,4,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos6,TPCLV_1,"TPCPV6",worldLV,false,5,fCheckOverlaps);

  new G4PVPlacement(0,TPC_pos7,TPCLV_2,"TPCPV7",worldLV,false,6,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos8,TPCLV_1,"TPCPV8",worldLV,false,7,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos9,TPCLV_1,"TPCPV9",worldLV,false,8,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos10,TPCLV_2,"TPCPV10",worldLV,false,9,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos11,TPCLV_1,"TPCPV11",worldLV,false,10,fCheckOverlaps);
  new G4PVPlacement(0,TPC_pos12,TPCLV_1,"TPCPV12",worldLV,false,11,fCheckOverlaps);


  // Scintillator
  
  G4double scint_layer_t = 3.*cm; G4double scint_layers = 10.; int nofLayers=10;
  G4double scint_t = scint_layers * scint_layer_t; G4double scint_length = 1.6*m;
  double n_scint_module = 3.0; G4double dx = 5.0*cm; G4double dy=5.0*cm; G4double scint_w = (TPC_total_length + dy + scint_t - 2.*dx)/n_scint_module;
  G4double dz = 0.05*m;

  G4double scint_light_detector_w = 1.0 *cm; G4double scint_light_detector_length = 1.0 *cm;   
  
  auto scintS = new G4Box("ScintS",scint_w/2., scint_t/2., scint_length/2.);
  auto scintLV = new G4LogicalVolume(scintS,defaultMaterial,"ScintLV");
  
  std::vector<G4VPhysicalVolume *> scint_PV_array;
  std::vector<G4VPhysicalVolume *> scint_Module_PV_array;
  std::vector<G4VPhysicalVolume *> scint_detector_PV_array; 

  // scintillator layers 
  // layer 1 and 2: (1 = horizontal, 2 = vertical staves)

  auto scint_layerH_S = new G4Box("ScintS",scint_w/2., scint_layer_t/2., scint_length/2.);
  auto scint_layerH_LV = new G4LogicalVolume(scint_layerH_S, defaultMaterial,"Scint_layerH_LV");
  
  auto scint_layerV_S = new G4Box("Scint2S",scint_w/2., scint_layer_t/2., scint_length/2.);
  auto scint_layerV_LV = new G4LogicalVolume(scint_layerV_S, defaultMaterial,"Scint_layerV_LV");
  
  // Two types of staves here, horizontal and vertical
  // --- horizontal staves (along x)
  auto scint_StaveH_x = scint_w/8.;
  auto scint_StaveH_z = scint_length;
  auto scint_StaveH_S = new G4Box("Scint_StaveS",scint_StaveH_x/2., scint_layer_t/2., scint_StaveH_z/2.);
  auto scint_StaveH_LV = new G4LogicalVolume(scint_StaveH_S,scintMaterial,"Scint_StaveH_LV");
  
  //--- vertical staves (divided along z)
  auto scint_StaveV_x = scint_w;
  auto scint_StaveV_z = scint_length/8.;
  auto scint_StaveV_S = new G4Box("Scint_StaveS",scint_StaveV_x/2., scint_layer_t/2.,scint_StaveV_z/2.);
  auto scint_StaveV_LV = new G4LogicalVolume(scint_StaveV_S,scintMaterial,"Scint_StaveV_LV");
 
  for (int i; i<10; i++){
    if (i%2==0){new G4PVPlacement(0, G4ThreeVector(0.,-scint_t/2.+(2*i+1)/2.*scint_layer_t,0.),scint_layerH_LV,"Scint_layerPV",scintLV,false,i,fCheckOverlaps);}
    else{new G4PVPlacement(0, G4ThreeVector(0.,-scint_t/2.+(2*i+1)/2.*scint_layer_t,0.),scint_layerV_LV,"Scint_layerPV",scintLV,false,i,fCheckOverlaps);}
  }

  // placing the staves 
  for (int i; i<8; i++){
    new G4PVPlacement(0,G4ThreeVector(-scint_w/2.+(2*i+1)/2.*scint_StaveH_x,0.,0.),scint_StaveH_LV,"Scint_layerPV",scint_layerH_LV,false,i,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,0.,-scint_length/2.+(2*i+1)/2.*scint_StaveV_z),scint_StaveV_LV,"Scint_layerPV",scint_layerV_LV,false,i,fCheckOverlaps);
  }

  // scint_PV_array.push_back(new G4PVPlacement(0, G4ThreeVector(0.,0.,-scint_Stave_2_w/2.+(2*i+1)/2.*scint_Stave_2_w),scint_layer2S,"Scint_layerPV",scint_layer2LV,false,i,fCheckOverlaps));

  auto scint_pos11 = G4ThreeVector(-(3*scint_w+2*dx)/2. + scint_w/2., 0. , -(3*scint_length+2*dz)/2.+scint_length/2.);
  auto scint_pos21 = scint_pos11 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos31 = scint_pos21 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos41 = scint_pos11 + G4ThreeVector(0.,0., dz + scint_length);
  auto scint_pos51 = scint_pos41 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos61 = scint_pos51 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos71 = scint_pos41 + G4ThreeVector(0.,0., dz + scint_length);
  auto scint_pos81 = scint_pos71 + G4ThreeVector(dx + scint_w,0.,0.);
  auto scint_pos91 = scint_pos81 + G4ThreeVector(dx + scint_w,0.,0.);


  auto scint_pos12 = G4ThreeVector(0.,-(3*scint_w+2*dx)/2. + scint_w/2., -(3*scint_length+2*dz)/2.+scint_length/2.);
  auto scint_pos22 = scint_pos12 + G4ThreeVector(0.,dx + scint_w,0.);
  auto scint_pos32 = scint_pos22 + G4ThreeVector(0.,dx + scint_w,0.);
  auto scint_pos42 = scint_pos12 + G4ThreeVector(0.,0., dz + scint_length);
  auto scint_pos52 = scint_pos42 + G4ThreeVector(0.,dx + scint_w,0.);
  auto scint_pos62 = scint_pos52 + G4ThreeVector(0.,dx + scint_w,0.);
  auto scint_pos72 = scint_pos42 + G4ThreeVector(0.,0., dz + scint_length);
  auto scint_pos82 = scint_pos72 + G4ThreeVector(0.,dx + scint_w,0.);
  auto scint_pos92 = scint_pos82 + G4ThreeVector(0.,dx + scint_w,0.);

  std::vector<G4ThreeVector> scint_pos_array;
  auto scint_group_pos1 = G4ThreeVector((3.*scint_w+2.*dx-1.0*TPC_total_length)/2.,TPC_total_length/2.+scint_t/2.+dy,0.);
  auto scint_group_pos2 = G4ThreeVector(-(scint_t/2.+dy+TPC_total_length/2.),(-TPC_total_length/2.+(3.*scint_w+2.*dx)/2.),0.);
  auto scint_group_pos3 = G4ThreeVector(-(3.*scint_w+2*dx-1.0*TPC_total_length)/2.,-(TPC_total_length/2.+scint_t/2.+dy),0.);
  auto scint_group_pos4 = G4ThreeVector((scint_t/2.+dy+TPC_total_length/2.),-(-TPC_total_length/2.+(3.*scint_w+2.*dx)/2.),0.);
  
  
  scint_pos_array.push_back(scint_group_pos1);
  scint_pos_array.push_back(scint_group_pos2);
  scint_pos_array.push_back(scint_group_pos3);
  scint_pos_array.push_back(scint_group_pos4);

  std::vector<G4RotationMatrix *> scint_rot_array;
  G4RotationMatrix* zRot1 = new G4RotationMatrix; zRot1 -> rotateZ(0.*deg);
  G4RotationMatrix* zRot2 = new G4RotationMatrix; zRot2 -> rotateZ(270.*deg);
  G4RotationMatrix* zRot3 = new G4RotationMatrix; zRot3 -> rotateZ(180.*deg);
  G4RotationMatrix* zRot4 = new G4RotationMatrix; zRot4 -> rotateZ(90.*deg);

  scint_rot_array.push_back(zRot1);
  scint_rot_array.push_back(zRot2);
  scint_rot_array.push_back(zRot3);
  scint_rot_array.push_back(zRot4);

  int scint_mod_index_count = 0;

  for (int i = 0; i<4;i++){
    std::cout << "in the loop!" << std::endl;
    if (i%2 == 0){ // if it is the 
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos11+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+0,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos21+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+1,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos31+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+2,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos41+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+3,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos51+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+4,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos61+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+5,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos71+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+6,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos81+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+7,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos91+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+8,fCheckOverlaps));
    }
    else{
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos12+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+0,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos22+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+1,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos32+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+2,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos42+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+3,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos52+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+4,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos62+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+5,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos72+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+6,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos82+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+7,fCheckOverlaps));
      scint_Module_PV_array.push_back(new G4PVPlacement(scint_rot_array[i],scint_pos92+scint_pos_array[i],scintLV,"ScintPV",worldLV,false,scint_mod_index_count+8,fCheckOverlaps));
    }
    scint_mod_index_count = scint_mod_index_count+9;
  }

  
  // - - - front and back scintillators
  // Type I dimensions
  G4double scint_x_fb1 = TPC_total_length/2. + dy + scint_t;
  G4double scint_y_fb1 = scint_x_fb1 - tube_thickness - tube_radius;
  // Type II dimensions
  G4double scint_x_fb2 = scint_y_fb1; // not sure why they are the same but it works 
  G4double scint_y_fb2 = 2.*scint_x_fb1 - 2.*scint_y_fb1;

  // virtual volume for the layers
  auto scint_fb1_S = new G4Box("ScintS",scint_x_fb1/2.,scint_y_fb1/2.,scint_t/2.);
  auto scint_fb1_LV = new G4LogicalVolume(scint_fb1_S,defaultMaterial,"Scint_fb_LV1");
  
  auto scint_fb2_S = new G4Box("ScintS",scint_x_fb2/2.,scint_y_fb2/2.,scint_t/2.);
  auto scint_fb2_LV = new G4LogicalVolume(scint_fb2_S,defaultMaterial,"Scint_fb_LV2");
  
  // defining the layers
  // type 1 layer 1 and 2: (1 = horizontal, 2 = vertical staves)
  auto scint_layer_fb11_S = new G4Box("ScintS", scint_x_fb1/2.,scint_y_fb1/2., scint_layer_t/2.);
  auto scint_layer_fb11_LV = new G4LogicalVolume(scint_layer_fb11_S,defaultMaterial,"Scint_Stave_LV11");
  auto scint_layer_fb12_S= new G4Box("ScintS", scint_x_fb1/2.,scint_y_fb1/2., scint_layer_t/2.);
  auto scint_layer_fb12_LV = new G4LogicalVolume(scint_layer_fb12_S,defaultMaterial,"Scint_Stave_LV12");
  
  auto scint_layer_fb21_S = new G4Box("ScintS", scint_x_fb2/2.,scint_y_fb2/2., scint_layer_t/2.);
  auto scint_layer_fb21_LV = new G4LogicalVolume(scint_layer_fb21_S,defaultMaterial,"Scint_Stave_LV21");
  auto scint_layer_fb22_S = new G4Box("ScintS", scint_x_fb2/2.,scint_y_fb2/2., scint_layer_t/2.);
  auto scint_layer_fb22_LV = new G4LogicalVolume(scint_layer_fb22_S,defaultMaterial,"Scint_Stave_LV22");

  // // staves for the front and back layers
  // Staves for type 1 fb scint
  // --- vertical staves (divide x)
  G4double scint_stave_xV_fb1 = scint_x_fb1/8.; 
  G4double scint_stave_yV_fb1 = scint_y_fb1;
  auto scint_stave_fb1V_S = new G4Box("ScintS", scint_stave_xV_fb1/2.,scint_stave_yV_fb1/2., scint_layer_t/2.);
  auto scint_stave_fb1V_LV = new G4LogicalVolume(scint_stave_fb1V_S,scintMaterial,"Scint_fb_Stave1V_LV");

 // --- horizontal staves (divide y)
  G4double scint_stave_xH_fb1 = scint_x_fb1;
  G4double scint_stave_yH_fb1 = scint_y_fb1/8.;
  auto scint_stave_fb1H_S = new G4Box("ScintS", scint_stave_xH_fb1/2.,scint_stave_yH_fb1/2., scint_layer_t/2.);
  auto scint_stave_fb1H_LV = new G4LogicalVolume(scint_stave_fb1H_S,scintMaterial,"Scint_fb_Stave1H_LV");
  
  // Staves for type 2 fb scint
  // --- vertical staves (divide x)
  G4double scint_stave_xV_fb2 = scint_x_fb2/8.; 
  G4double scint_stave_yV_fb2 = scint_y_fb2;
  auto scint_stave_fb2V_S = new G4Box("ScintS", scint_stave_xV_fb2/2.,scint_stave_yV_fb2/2., scint_layer_t/2.);
  auto scint_stave_fb2V_LV = new G4LogicalVolume(scint_stave_fb2V_S,scintMaterial,"Scint_fb_Stave2V_LV");

  // --- horizontal staves (divide y)
  G4double scint_stave_xH_fb2 = scint_x_fb2;
  G4double scint_stave_yH_fb2 = scint_y_fb2/8.;
  auto scint_stave_fb2H_S = new G4Box("ScintS", scint_stave_xH_fb2/2.,scint_stave_yH_fb2/2., scint_layer_t/2.);
  auto scint_stave_fb2H_LV = new G4LogicalVolume(scint_stave_fb2H_S,scintMaterial,"Scint_fb_Stave2H_LV");

  for (int i; i<10; i++){
    if (i%2==0){new G4PVPlacement(0, G4ThreeVector(0.,0.,-scint_t/2.+(2*i+1)/2.*scint_layer_t),scint_layer_fb11_LV,"Scint_layerPV",scint_fb1_LV,false,i,fCheckOverlaps);}
    else{new G4PVPlacement(0, G4ThreeVector(0.,0.,-scint_t/2.+(2*i+1)/2.*scint_layer_t),scint_layer_fb12_LV,"Scint_layerPV",scint_fb1_LV,false,i,fCheckOverlaps);}
    
    if (i%2==0){new G4PVPlacement(0, G4ThreeVector(0.,0.,-scint_t/2.+(2*i+1)/2.*scint_layer_t),scint_layer_fb21_LV,"Scint_layerPV",scint_fb2_LV,false,i,fCheckOverlaps);}
    else{new G4PVPlacement(0, G4ThreeVector(0.,0.,-scint_t/2.+(2*i+1)/2.*scint_layer_t),scint_layer_fb22_LV,"Scint_layerPV",scint_fb2_LV,false,i,fCheckOverlaps);}
  }

  // placing the staves 
  for (int i; i<8; i++){
    // Type I scintillator 
    // --- horizontal ones for odd number layers
    new G4PVPlacement(0,G4ThreeVector(0.,-scint_y_fb1/2.+(2*i+1)/2.*scint_stave_yH_fb1,0.),scint_stave_fb1H_LV,"Scint_layerPV",scint_layer_fb11_LV,false,i,fCheckOverlaps);
    // --- vertical ones for even number layers
    new G4PVPlacement(0,G4ThreeVector(-scint_x_fb1/2.+(2*i+1)/2.*scint_stave_xV_fb1,0.,0.),scint_stave_fb1V_LV,"Scint_layerPV",scint_layer_fb12_LV,false,i,fCheckOverlaps);
    
    // Type II scintillator 
    // --- horizontal ones for odd number layers
    new G4PVPlacement(0,G4ThreeVector(0.,-scint_y_fb2/2.+(2*i+1)/2.*scint_stave_yH_fb2,0.),scint_stave_fb2H_LV,"Scint_layerPV",scint_layer_fb21_LV,false,i,fCheckOverlaps);
    // --- vertical ones for even number layers
    new G4PVPlacement(0,G4ThreeVector(-scint_x_fb2/2.+(2*i+1)/2.*scint_stave_xV_fb2,0.,0.),scint_stave_fb2V_LV,"Scint_layerPV",scint_layer_fb22_LV,false,i,fCheckOverlaps); 
  }
  
  // rotation is introduced to make a better indexing of the layers
  G4RotationMatrix * scint_fb_rot =  new G4RotationMatrix; scint_fb_rot -> rotateX(180.0*deg);
  
  auto scint_fb_pos1 = G4ThreeVector(scint_x_fb1/2.,tube_radius+tube_thickness+scint_y_fb1/2.,-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos2 = G4ThreeVector(-scint_x_fb1/2.,tube_radius+tube_thickness+scint_y_fb1/2.,-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos3 = G4ThreeVector(scint_x_fb1/2.,-(tube_radius+tube_thickness+scint_y_fb1/2.),-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos4 = G4ThreeVector(-scint_x_fb1/2.,-(tube_radius+tube_thickness+scint_y_fb1/2.),-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos5 = G4ThreeVector(-(tube_radius+tube_thickness+scint_y_fb1/2.),0.,-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos6 = G4ThreeVector((tube_radius+tube_thickness+scint_y_fb1/2.),0.,-(3*scint_length+2*dz)/2.-scint_t/2.);
  auto scint_fb_pos7 = G4ThreeVector(scint_x_fb1/2.,tube_radius+tube_thickness+scint_y_fb1/2.,(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos8 = G4ThreeVector(-scint_x_fb1/2.,tube_radius+tube_thickness+scint_y_fb1/2.,(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos9 = G4ThreeVector(scint_x_fb1/2.,-(tube_radius+tube_thickness+scint_y_fb1/2.),(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos10 = G4ThreeVector(-scint_x_fb1/2.,-(tube_radius+tube_thickness+scint_y_fb1/2.),(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos11 = G4ThreeVector(-(tube_radius+tube_thickness+scint_y_fb1/2.),0.,(3*scint_length+2*dz)/2.+scint_t/2.);
  auto scint_fb_pos12 = G4ThreeVector((tube_radius+tube_thickness+scint_y_fb1/2.),0.,(3*scint_length+2*dz)/2.+scint_t/2.);
  
  auto scint_fb_PV1 = new G4PVPlacement(scint_fb_rot,scint_fb_pos2,scint_fb1_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+0,fCheckOverlaps); 
  auto scint_fb_PV2 = new G4PVPlacement(scint_fb_rot,scint_fb_pos5,scint_fb2_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+1,fCheckOverlaps); 
  auto scint_fb_PV3 = new G4PVPlacement(scint_fb_rot,scint_fb_pos4,scint_fb1_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+2,fCheckOverlaps); 
  auto scint_fb_PV4 = new G4PVPlacement(scint_fb_rot,scint_fb_pos3,scint_fb1_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+3,fCheckOverlaps); 
  auto scint_fb_PV5 = new G4PVPlacement(scint_fb_rot,scint_fb_pos6,scint_fb2_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+4,fCheckOverlaps); 
  auto scint_fb_PV6 = new G4PVPlacement(scint_fb_rot,scint_fb_pos1,scint_fb1_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+5,fCheckOverlaps); 

  auto scint_fb_PV7 = new G4PVPlacement(0,scint_fb_pos7,scint_fb1_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+6,fCheckOverlaps);
  auto scint_fb_PV8 = new G4PVPlacement(0,scint_fb_pos12,scint_fb2_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+7,fCheckOverlaps);
  auto scint_fb_PV9 = new G4PVPlacement(0,scint_fb_pos10,scint_fb1_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+8,fCheckOverlaps);
  auto scint_fb_PV10 = new G4PVPlacement(0,scint_fb_pos9,scint_fb1_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+9,fCheckOverlaps);
  auto scint_fb_PV11 = new G4PVPlacement(0,scint_fb_pos11,scint_fb2_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+10,fCheckOverlaps);
  auto scint_fb_PV12 = new G4PVPlacement(0,scint_fb_pos8,scint_fb1_LV,"ScintPVfb",worldLV,false,scint_mod_index_count+11,fCheckOverlaps);
  
  // cuts for scintillator
  G4Region* scint_region = new G4Region("Scint_region");
  scint_region->AddRootLogicalVolume(scint_StaveH_LV);
  scint_region->AddRootLogicalVolume(scint_StaveV_LV);
  scint_region->AddRootLogicalVolume(scint_stave_fb1V_LV);
  scint_region->AddRootLogicalVolume(scint_stave_fb1H_LV);
  scint_region->AddRootLogicalVolume(scint_stave_fb2V_LV);
  scint_region->AddRootLogicalVolume(scint_stave_fb2H_LV); 
  G4Region* scint_Region = G4RegionStore::GetInstance()->GetRegion("Scint_region");
  G4ProductionCuts* scintcut = new G4ProductionCuts();
  scintcut->SetProductionCut(500.0*cm,"gamma"); // 130cm -> 1.3 MeV   220cm -> 3.4 MeV 400*cm -> 50 MeV
  scintcut->SetProductionCut(15.0*mm,"e-");
  scintcut->SetProductionCut(15.0*mm,"e+");
  scintcut->SetProductionCut(15.0*mm,"proton");
  scint_Region->SetProductionCuts(scintcut);
  
  std::cout<<"LEAD GLASS"<<std::endl;

  // Lead Glass
  G4double lead_glass_xy = 8.*cm; G4double lead_glass_z = 25.*cm;
  G4double coating_thickness = 0.01*mm;
  auto lead_glass_y_level = TPC_total_length/2. +dy + scint_t;
  
  G4double offset_lead_glass = 0.0*mm;
  auto absorberS = new G4Box("AbsoS",(lead_glass_xy+coating_thickness)/2., lead_glass_z/2., (lead_glass_xy+coating_thickness)/2.); // is the whole lead glass module including the coating
  auto absorberLV = new G4LogicalVolume(absorberS,defaultMaterial,"AbsoLV");
  int lead_index = 0;
  auto leadglassS = new G4Box("LeadGlassS",(lead_glass_xy)/2., (lead_glass_z-coating_thickness)/2., (lead_glass_xy)/2.); // just the lead glass itself
  auto leadglassLV = new G4LogicalVolume(leadglassS,absorberMaterial,"LeadGlassLV");
  new G4PVPlacement(0,G4ThreeVector(0., -lead_glass_z/2.+(lead_glass_z-coating_thickness)/2., 0.),leadglassLV,"LeadGlassPV",absorberLV,false,0,fCheckOverlaps);

  auto hallowS = new G4Box("HallowS",(lead_glass_xy)/2., (lead_glass_z)/2., (lead_glass_xy)/2.); // just the lead glass itself
  
  // Coating of the lead glass 
  auto coatingS = new G4SubtractionSolid("Coating", absorberS, hallowS,0, G4ThreeVector(0.,0.,0.));
  auto coatingLV = new G4LogicalVolume(coatingS,leadglasscoatMaterial,"CoatingLV");                              
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),coatingLV,"Coating",absorberLV,false,0,fCheckOverlaps);  // checking overlaps 
  
  // PMT window
  auto PMTS = new G4Box("PMTS",(lead_glass_xy)/2., coating_thickness/2., (lead_glass_xy)/2.);
  auto PMTLV = new G4LogicalVolume(PMTS,defaultMaterial,"PMTLV");            
  new G4PVPlacement(0,G4ThreeVector(0., -lead_glass_z/2.+(lead_glass_z-coating_thickness)+coating_thickness/2., 0.),PMTLV,"PMTPV",absorberLV,false,0,fCheckOverlaps);

  G4Region* abs_region = new G4Region("Abs_region");
  abs_region->AddRootLogicalVolume(leadglassLV); 
  G4Region* absRegion = G4RegionStore::GetInstance()->GetRegion("Abs_region");
  G4ProductionCuts* abscut = new G4ProductionCuts();
  abscut->SetProductionCut(15.0 *cm,"gamma"); 
  abscut->SetProductionCut(1 *mm,"e-");
  abscut->SetProductionCut(1 *mm,"e+");
  abscut->SetProductionCut(1 *mm,"proton");
  absRegion->SetProductionCuts(abscut);

  ///***calculate the position directly from python
  
  std::cout << " Constructing the lead glass blocks from data" << std::endl;
  std::vector<std::vector<G4RotationMatrix *>> rot_array_dir;
  std::vector<std::vector<G4VPhysicalVolume *>> lg_PV_array;
  G4double rot_angle = 90.0*deg;
 
  //(disabled)

  for (int i = 0; i < 4; i++){ 
    std::vector<G4RotationMatrix *>temp_rot_array_dir;
    std::vector<G4VPhysicalVolume*>temp_PV_array;
    for (int j = 0 ; j < data_lead_glass_pos.size();j++){temp_rot_array_dir.push_back(new G4RotationMatrix);}
    rot_array_dir.push_back(temp_rot_array_dir);
    lg_PV_array.push_back(temp_PV_array);
  }
  
  // Placing the lead glass blocks on the 4 surfaces 
  for (int i = 0 ; i <4;i++){
    for (int j = 0 ; j <data_lead_glass_pos.size();j++){ // data_lead_glass_pos.size()

      G4double lead_x = 0.; G4double lead_y=0.;
      G4double lead_x0 = 0. ; G4double lead_y0 = 0.*cm;
      lead_x0 = data_lead_glass_pos[j][0]*cm; lead_y0 = data_lead_glass_pos[j][1]*cm; 

      double i_angle = (double) i;
      lead_x = lead_x0*cos(i_angle*rot_angle) - (lead_y0)*sin(i_angle*rot_angle);
      lead_y = lead_x0*sin(i_angle*rot_angle) + (lead_y0)*cos(i_angle*rot_angle);
      rot_array_dir[i][j]->rotate(-data_lead_glass_pos[j][3]*deg, G4ThreeVector(cos(i_angle*rot_angle),sin(i_angle*rot_angle),0.));
      rot_array_dir[i][j]-> rotateZ(data_lead_glass_pos[j][5]*deg -i_angle*rot_angle); //data_lead_glass_pos[i][5]*deg 
      lg_PV_array[i].push_back(new G4PVPlacement(rot_array_dir[i][j],G4ThreeVector(lead_x, lead_y,data_lead_glass_pos[j][2]*cm) ,absorberLV,"AbsoPV",worldLV,false,lead_index,1));
      
      Lead_glass_outFile << lead_index << "," << lead_x/cm << "," << lead_y/cm << "," << data_lead_glass_pos[j][2] << G4endl;
      //std::cout << lead_index << "," << lead_x/cm << "," << lead_y/cm << "," << data_lead_glass_pos[j][2] << std::endl; 
      lead_index ++ ;
    }
  }
  
  // for the Front and back lead glass
  G4double lead_glass_y_level_fb = (3*scint_length+2*dz)/2.+scint_t;
  std::cout << "Full Detector half length " << lead_glass_y_level_fb << std::endl;
  std::vector<G4RotationMatrix *> rot_array_dir111;
  std::vector<G4RotationMatrix *> rot_array_dir211;

  for (int i = 0 ; i < data_lead_glass_pos_fb.size();i++){rot_array_dir111.push_back(new G4RotationMatrix);}
  for (int i = 0 ; i < data_lead_glass_pos_fb.size();i++){rot_array_dir211.push_back(new G4RotationMatrix);}
  
  for (int i = 0 ; i < data_lead_glass_pos_fb.size();i++){
        rot_array_dir111[i]-> rotateX(-data_lead_glass_pos_fb[i][3]*deg+90.*deg); rot_array_dir111[i]-> rotateZ(data_lead_glass_pos_fb[i][5]*deg); //rot_array_dir111[i]-> rotateY(-data_lead_glass_pos_fb[i][5]*deg);
        new G4PVPlacement(rot_array_dir111[i],
        G4ThreeVector(data_lead_glass_pos_fb[i][0]*cm,data_lead_glass_pos_fb[i][2]*cm,-(data_lead_glass_pos_fb[i][1]*cm))
        ,absorberLV,"AbsoPV",worldLV,false,lead_index,0);
        lead_index ++ ;        
  }
  
  for (int i = 0 ; i < data_lead_glass_pos_fb.size();i++){
        rot_array_dir211[i]-> rotateX(data_lead_glass_pos_fb[i][3]*deg-90.*deg); rot_array_dir211[i]-> rotateZ(data_lead_glass_pos_fb[i][5]*deg); //rot_array_dir111[i]-> rotateY(-data_lead_glass_pos_fb[i][5]*deg);
        new G4PVPlacement(rot_array_dir211[i],
        G4ThreeVector(data_lead_glass_pos_fb[i][0]*cm,data_lead_glass_pos_fb[i][2]*cm,(data_lead_glass_pos_fb[i][1]*cm))
        ,absorberLV,"AbsoPV",worldLV,false,lead_index,0);
        lead_index ++ ;        
  }


  /////Shield building 

  #if VERSION_SHIELD == 1
    // Lead Shield  
    G4double shield_offset = 50.0*cm;
    
    // ===== TOP and BOTTOM SURFACE ===== //
    G4double lead_top_x = 2.0*TPC_x_1+TPC_x_2+2.0*dz+2.0*scint_t+2.0*lead_glass_z+2.0*shield_offset; 
    G4double lead_top_y = 100.*cm; 
    G4double lead_top_z = 2.0*TPC_z+2.0*scint_t+2.0*dz+2.0*lead_glass_z+4.0*shield_offset;

    // Lead Shield top & bottom
    auto LeadShield_top_S = new G4Box("LeadShield_top_S",lead_top_x/2., lead_top_y/2., lead_top_z/2.);
    auto LeadShield_top_LV = new G4LogicalVolume(LeadShield_top_S,shieldMaterial,"LeadShield_top_LV");
    new G4PVPlacement(0,G4ThreeVector(0.,lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y/2.,0.),LeadShield_top_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,-(lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y/2.),0.),LeadShield_top_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    

    //measuring plane 0.1 cm thick 
    auto Measuring_plane_top_S = new G4Box("Measuring_plane_top_S",lead_top_x/2., 0.1*cm/2., lead_top_z/2.);
    auto Measuring_plane_top_LV = new G4LogicalVolume(Measuring_plane_top_S,defaultMaterial,"Measuring_plane_top_LV");

    new G4PVPlacement(0,G4ThreeVector(0.,(lead_glass_y_level+lead_glass_z+15.0*cm),0.),Measuring_plane_top_LV,"MeasurePlane",worldLV,false,0,fCheckOverlaps);
    
    // ===== 2 Sides ===== //
    // Lead shield Sides  
    G4double lead_side_y = 2.0*(lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y); 
    G4double lead_side_x = 100.*cm; 
    G4double lead_side_z = 2.0*TPC_z+2.0*scint_t+2.0*dz+2.0*lead_glass_z+4.0*shield_offset;

    auto LeadShield_side_S = new G4Box("LeadShield_side_S",lead_side_x/2., lead_side_y/2., lead_side_z/2.);
    auto LeadShield_side_LV = new G4LogicalVolume(LeadShield_side_S,shieldMaterial,"LeadShield_side_LV");
    new G4PVPlacement(0,G4ThreeVector(lead_top_x/2.+lead_side_x/2.,0.,0.),LeadShield_side_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(-(lead_top_x/2.+lead_side_x/2.),0.,0.),LeadShield_side_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    // lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y-lead_side_y/2.
    
    // Lead Shield FB 
    G4double lead_fb_x = 2.0*TPC_x_1+TPC_x_2+2.0*dz+2.0*scint_t+2.0*lead_glass_z+2.0*shield_offset+2.0*lead_side_x;
    G4double lead_fb_y = lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y-tube_thickness-tube_radius;
    G4double lead_fb_z = 100.0*cm;

    auto LeadShield_fb_S = new G4Box("LeadShield_fb_S",lead_fb_x/2., lead_fb_y/2., lead_fb_z/2.);
    auto LeadShield_fb_LV = new G4LogicalVolume(LeadShield_fb_S,shieldMaterial,"LeadShield_fb_LV");
    new G4PVPlacement(0,G4ThreeVector(0.,lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y-lead_fb_y/2.,lead_top_z/2.+lead_fb_z/2.),LeadShield_fb_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y-lead_fb_y/2.,-(lead_top_z/2.+lead_fb_z/2.)),LeadShield_fb_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,-(lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y-lead_fb_y/2.),lead_top_z/2.+lead_fb_z/2.),LeadShield_fb_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,-(lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y-lead_fb_y/2.),-(lead_top_z/2.+lead_fb_z/2.)),LeadShield_fb_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);

    // G4double lead_fb_bot_x = 2.0*TPC_x_1+TPC_x_2+2.0*dz+2.0*scint_t+2.0*lead_glass_z+2.0*shield_offset+2.0*lead_side_x;
    // G4double lead_fb_bot_y = lead_glass_y_level+lead_glass_z-tube_thickness-tube_radius;
    // G4double lead_fb_bot_z = 100.0*cm;
    
    // auto LeadShield_fb_bot_S = new G4Box("LeadShield_fb_bot_S",lead_fb_bot_x/2., lead_fb_bot_y/2., lead_fb_bot_z/2.);
    // auto LeadShield_fb_bot_LV = new G4LogicalVolume(LeadShield_fb_bot_S,shieldMaterial,"LeadShield_fb_bot_LV");
    // new G4PVPlacement(0,G4ThreeVector(0.,-(lead_glass_y_level+lead_glass_z)+lead_fb_bot_y/2.,lead_top_z/2.+lead_fb_z/2.),LeadShield_fb_bot_LV,"LeadShield_fb_bot_PV",worldLV,false,0,fCheckOverlaps);
    // new G4PVPlacement(0,G4ThreeVector(0.,-(lead_glass_y_level+lead_glass_z)+lead_fb_bot_y/2.,-(lead_top_z/2.+lead_fb_z/2.)),LeadShield_fb_bot_LV,"LeadShield_fb_bot_PV",worldLV,false,0,fCheckOverlaps);

    G4double lead_fb_mid_x = (lead_fb_x-2.0*tube_radius-2.0*tube_thickness)/2.0;
    G4double lead_fb_mid_y = 2.0*tube_radius+2.0*tube_thickness;
    G4double lead_fb_mid_z = 100.0*cm;
    
    auto LeadShield_fb_mid_S = new G4Box("LeadShield_fb_mid_S",lead_fb_mid_x/2., lead_fb_mid_y/2., lead_fb_mid_z/2.);
    auto LeadShield_fb_mid_LV = new G4LogicalVolume(LeadShield_fb_mid_S,shieldMaterial,"LeadShield_fb_mid_LV");
    new G4PVPlacement(0,G4ThreeVector(lead_fb_x/2.-lead_fb_mid_x/2.,0.,lead_top_z/2.+lead_fb_z/2.),LeadShield_fb_mid_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(-(lead_fb_x/2.-lead_fb_mid_x/2.),0.,lead_top_z/2.+lead_fb_z/2.),LeadShield_fb_mid_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);

    new G4PVPlacement(0,G4ThreeVector(lead_fb_x/2.-lead_fb_mid_x/2.,0.,-(lead_top_z/2.+lead_fb_z/2.)),LeadShield_fb_mid_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(-(lead_fb_x/2.-lead_fb_mid_x/2.),0.,-(lead_top_z/2.+lead_fb_z/2.)),LeadShield_fb_mid_LV,"LeadShield",worldLV,false,0,fCheckOverlaps);
    //new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),LeadShield_fb_mid_LV,"LeadShield_fb_mid_PV",worldLV,false,0,fCheckOverlaps);
    
    // ***********************************************
    // ======= Scintillator Shield building  =========
    // ***********************************************
    // Scint Shield top 
    G4double scint_top_x = 2.0*TPC_x_1+TPC_x_2+2.0*dz+2.0*scint_t+2.0*lead_glass_z+2.0*shield_offset+2.0*lead_side_x; 
    G4double scint_top_y = 5.*cm; 
    G4double scint_top_z = 2.0*TPC_z+2.0*scint_t+2.0*dz+2.0*lead_glass_z+4.0*shield_offset+2.0*lead_fb_z;

    auto ScintShield_top_1_S = new G4Box("ScintShield_top_S",scint_top_x/2., scint_top_y/2., scint_top_z/2.);
    auto ScintShield_top_1_LV = new G4LogicalVolume(ScintShield_top_1_S,defaultMaterial,"ScintShield_top_1_LV");

    auto ScintShield_top_2_S = new G4Box("ScintShield_top_S",scint_top_x/2., scint_top_y/2., scint_top_z/2.);
    auto ScintShield_top_2_LV = new G4LogicalVolume(ScintShield_top_2_S,defaultMaterial,"ScintShield_top_2_LV");
    
    new G4PVPlacement(0,G4ThreeVector(0.,lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y+scint_top_y/2.,0.),ScintShield_top_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);  
    new G4PVPlacement(0,G4ThreeVector(0.,lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y+3.*scint_top_y/2.,0.),ScintShield_top_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,-(lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y+scint_top_y/2.),0.),ScintShield_top_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);  
    new G4PVPlacement(0,G4ThreeVector(0.,-(lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y+3.*scint_top_y/2.),0.),ScintShield_top_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);

    // ===> Filling Scint Shield top with scint staves
    int num_x_top_staves_1 = 20;
    int num_z_top_staves_1 = 40;
    
    G4double x_top_staves_1 = scint_top_x/num_x_top_staves_1*mm;
    G4double z_top_staves_1 = scint_top_z/num_z_top_staves_1*mm;


    auto ScintShield_stave_top_1_S = new G4Box("ScintShield_top_S",x_top_staves_1/2., scint_top_y/2., z_top_staves_1/2.);
    auto ScintShield_stave_top_1_LV = new G4LogicalVolume(ScintShield_stave_top_1_S,scintMaterial,"ScintShield_stave_top_1_LV");

    for (int i = 0;i<num_x_top_staves_1;i++){
      for (int j=0;j<num_z_top_staves_1;j++){
        new G4PVPlacement(0,G4ThreeVector(-scint_top_x/2.+((2.0*i+1)*x_top_staves_1/2.),0.,-scint_top_z/2.+((2.0*j+1)*z_top_staves_1/2.)),ScintShield_stave_top_1_LV,"ScintShield",ScintShield_top_1_LV,false,0,fCheckOverlaps);
        new G4PVPlacement(0,G4ThreeVector(-scint_top_x/2.+((2.0*i+1)*x_top_staves_1/2.),0.,-scint_top_z/2.+((2.0*j+1)*z_top_staves_1/2.)),ScintShield_stave_top_1_LV,"ScintShield",ScintShield_top_2_LV,false,0,fCheckOverlaps);
      }
    }

    // Scint Shield sides
    G4double scint_sides_y = 2.*(lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y+2.*scint_top_y); 
    G4double scint_sides_x = 5.*cm; 
    G4double scint_sides_z = 2.0*TPC_z+2.0*scint_t+2.0*dz+2.0*lead_glass_z+4.0*shield_offset+2.0*lead_fb_z;

    auto ScintShield_sides_1_S = new G4Box("ScintShield_sides_S",scint_sides_x/2., scint_sides_y/2., scint_sides_z/2.);
    auto ScintShield_sides_1_LV = new G4LogicalVolume(ScintShield_sides_1_S,defaultMaterial,"ScintShield_sides_1_LV");

    auto ScintShield_sides_2_S = new G4Box("ScintShield_sides_S",scint_sides_x/2., scint_sides_y/2., scint_sides_z/2.);
    auto ScintShield_sides_2_LV = new G4LogicalVolume(ScintShield_sides_2_S,defaultMaterial,"ScintShield_sides_2_LV");

    new G4PVPlacement(0,G4ThreeVector(lead_top_x/2.+lead_side_x+scint_sides_x/2.,0.,0.),ScintShield_sides_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(lead_top_x/2.+lead_side_x+3.*scint_sides_x/2.,0.,0.),ScintShield_sides_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(-(lead_top_x/2.+lead_side_x+scint_sides_x/2.),0.,0.),ScintShield_sides_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(-(lead_top_x/2.+lead_side_x+3.*scint_sides_x/2.),0.,0.),ScintShield_sides_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    
    // ===> Filling Scint Shield top with scint staves
    int num_y_sides_staves_1 = 20;
    int num_z_sides_staves_1 = 40;
    
    G4double y_sides_staves_1 = scint_sides_y/num_y_sides_staves_1*mm;
    G4double z_sides_staves_1 = scint_sides_z/num_z_sides_staves_1*mm;

    auto ScintShield_stave_sides_1_S = new G4Box("ScintShield_top_S",scint_sides_x/2., y_sides_staves_1/2., z_sides_staves_1/2.);
    auto ScintShield_stave_sides_1_LV = new G4LogicalVolume(ScintShield_stave_sides_1_S,scintMaterial,"ScintShield_stave_sides_1_LV");

    for (int i = 0;i<num_y_sides_staves_1;i++){
      for (int j=0;j<num_z_sides_staves_1;j++){
        new G4PVPlacement(0,G4ThreeVector(0.,-scint_sides_y/2.+((2.0*i+1)*y_sides_staves_1/2.),-scint_sides_z/2.+((2.0*j+1)*z_sides_staves_1/2.)),ScintShield_stave_sides_1_LV,"ScintShield",ScintShield_sides_1_LV,false,0,fCheckOverlaps);
        new G4PVPlacement(0,G4ThreeVector(0.,-scint_sides_y/2.+((2.0*i+1)*y_sides_staves_1/2.),-scint_sides_z/2.+((2.0*j+1)*z_sides_staves_1/2.)),ScintShield_stave_sides_1_LV,"ScintShield",ScintShield_sides_2_LV,false,0,fCheckOverlaps);
      }
    }

    //Scint Shield Front-Back 
    G4double scint_fb_x = 2.0*TPC_x_1+TPC_x_2+2.0*dz+2.0*scint_t+2.0*lead_glass_z+2.0*shield_offset+2.0*lead_side_x + 4.*scint_sides_x;
    G4double scint_fb_y = lead_glass_y_level+lead_glass_z+shield_offset+lead_top_y-tube_thickness-tube_radius+2.*scint_sides_x;
    G4double scint_fb_z = 5.0*cm;

    auto ScintShield_fb_1_S = new G4Box("ScintShield_fb_S",scint_fb_x/2., scint_fb_y/2., scint_fb_z/2.);
    auto ScintShield_fb_1_LV = new G4LogicalVolume(ScintShield_fb_1_S,defaultMaterial,"ScintShield_fb_1_LV");

    auto ScintShield_fb_2_S = new G4Box("ScintShield_fb_S",scint_fb_x/2., scint_fb_y/2., scint_fb_z/2.);
    auto ScintShield_fb_2_LV = new G4LogicalVolume(ScintShield_fb_2_S,defaultMaterial,"ScintShield_fb_2_LV");

    // front top
    new G4PVPlacement(0,G4ThreeVector(0.,lead_top_x/2.+lead_side_x+2.*scint_sides_x-scint_fb_y/2.,lead_top_z/2.+lead_fb_z+scint_fb_z/2.),ScintShield_fb_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,lead_top_x/2.+lead_side_x+2.*scint_sides_x-scint_fb_y/2.,lead_top_z/2.+lead_fb_z+3.*scint_fb_z/2.),ScintShield_fb_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    // front bottom
    new G4PVPlacement(0,G4ThreeVector(0.,-(lead_top_x/2.+lead_side_x+2.*scint_sides_x-scint_fb_y/2.),lead_top_z/2.+lead_fb_z+scint_fb_z/2.),ScintShield_fb_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,-(lead_top_x/2.+lead_side_x+2.*scint_sides_x-scint_fb_y/2.),lead_top_z/2.+lead_fb_z+3.*scint_fb_z/2.),ScintShield_fb_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    
    // back top
    new G4PVPlacement(0,G4ThreeVector(0.,lead_top_x/2.+lead_side_x+2.*scint_sides_x-scint_fb_y/2.,-(lead_top_z/2.+lead_fb_z+scint_fb_z/2.)),ScintShield_fb_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,lead_top_x/2.+lead_side_x+2.*scint_sides_x-scint_fb_y/2.,-(lead_top_z/2.+lead_fb_z+3.*scint_fb_z/2.)),ScintShield_fb_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    
    //back bottom
    new G4PVPlacement(0,G4ThreeVector(0.,-(lead_top_x/2.+lead_side_x+2.*scint_sides_x-scint_fb_y/2.),-(lead_top_z/2.+lead_fb_z+scint_fb_z/2.)),ScintShield_fb_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(0.,-(lead_top_x/2.+lead_side_x+2.*scint_sides_x-scint_fb_y/2.),-(lead_top_z/2.+lead_fb_z+3.*scint_fb_z/2.)),ScintShield_fb_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    
    // ===> Filling Scint Shield FB with scint staves
    int num_x_fb_staves_1 = 20;
    int num_y_fb_staves_1 = 10;

    G4double x_fb_staves_1 = scint_fb_x/num_x_fb_staves_1*mm;
    G4double y_fb_staves_1 = scint_fb_y/num_y_fb_staves_1*mm;

    auto ScintShield_stave_fb_1_S = new G4Box("ScintShield_fb_S",x_fb_staves_1/2.,y_fb_staves_1/2.,scint_fb_z/2.);
    auto ScintShield_stave_fb_1_LV = new G4LogicalVolume(ScintShield_stave_fb_1_S,scintMaterial,"ScintShield_stave_fb_1_LV");

    for (int i = 0;i<num_x_fb_staves_1;i++){
      for (int j=0;j<num_y_fb_staves_1;j++){
        new G4PVPlacement(0,G4ThreeVector(-scint_fb_x/2.+((2.0*i+1)*x_fb_staves_1/2.),-scint_fb_y/2.+((2.0*j+1)*y_fb_staves_1/2.),0.),ScintShield_stave_fb_1_LV,"ScintShield",ScintShield_fb_1_LV,false,0,fCheckOverlaps);
        new G4PVPlacement(0,G4ThreeVector(-scint_fb_x/2.+((2.0*i+1)*x_fb_staves_1/2.),-scint_fb_y/2.+((2.0*j+1)*y_fb_staves_1/2.),0.),ScintShield_stave_fb_1_LV,"ScintShield",ScintShield_fb_2_LV,false,0,fCheckOverlaps);
      }
    }

    // Scint Shield Front Back mid
    G4double scint_fb_mid_x = (lead_fb_x+2.*scint_sides_x-2.0*tube_radius-2.0*tube_thickness)/2.0;
    G4double scint_fb_mid_y = 2.0*tube_radius+2.0*tube_thickness;
    G4double scint_fb_mid_z = 5.0*cm;
    
    auto ScintShield_fb_mid_1_S = new G4Box("ScintShield_fb_mid_S",scint_fb_mid_x/2., scint_fb_mid_y/2., scint_fb_mid_z/2.);
    auto ScintShield_fb_mid_1_LV = new G4LogicalVolume(ScintShield_fb_mid_1_S,defaultMaterial,"ScintShield_fb_mid_1_LV");
    
    auto ScintShield_fb_mid_2_S = new G4Box("ScintShield_fb_mid_S",scint_fb_mid_x/2., scint_fb_mid_y/2., scint_fb_mid_z/2.);
    auto ScintShield_fb_mid_2_LV = new G4LogicalVolume(ScintShield_fb_mid_2_S,defaultMaterial,"ScintShield_fb_mid_2_LV");

    //front right
    new G4PVPlacement(0,G4ThreeVector((lead_fb_x+4.*scint_sides_x)/2.-scint_fb_mid_x/2.,0.,lead_top_z/2.+lead_fb_mid_z+scint_fb_mid_z/2.),ScintShield_fb_mid_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector((lead_fb_x+4.*scint_sides_x)/2.-scint_fb_mid_x/2.,0.,lead_top_z/2.+lead_fb_mid_z+3.*scint_fb_mid_z/2.),ScintShield_fb_mid_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    //front left
    new G4PVPlacement(0,G4ThreeVector(-((lead_fb_x+4.*scint_sides_x)/2.-scint_fb_mid_x/2.),0.,lead_top_z/2.+lead_fb_mid_z+scint_fb_mid_z/2.),ScintShield_fb_mid_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(-((lead_fb_x+4.*scint_sides_x)/2.-scint_fb_mid_x/2.),0.,lead_top_z/2.+lead_fb_mid_z+3.*scint_fb_mid_z/2.),ScintShield_fb_mid_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);

    //back right
    new G4PVPlacement(0,G4ThreeVector((lead_fb_x+4.*scint_sides_x)/2.-scint_fb_mid_x/2.,0.,-(lead_top_z/2.+lead_fb_mid_z+scint_fb_mid_z/2.)),ScintShield_fb_mid_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector((lead_fb_x+4.*scint_sides_x)/2.-scint_fb_mid_x/2.,0.,-(lead_top_z/2.+lead_fb_mid_z+3.*scint_fb_mid_z/2.)),ScintShield_fb_mid_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);

    //back left
    new G4PVPlacement(0,G4ThreeVector(-((lead_fb_x+4.*scint_sides_x)/2.-scint_fb_mid_x/2.),0.,-(lead_top_z/2.+lead_fb_mid_z+scint_fb_mid_z/2.)),ScintShield_fb_mid_1_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);
    new G4PVPlacement(0,G4ThreeVector(-((lead_fb_x+4.*scint_sides_x)/2.-scint_fb_mid_x/2.),0.,-(lead_top_z/2.+lead_fb_mid_z+3.*scint_fb_mid_z/2.)),ScintShield_fb_mid_2_LV,"ScintShield",worldLV,false,0,fCheckOverlaps);

    // ===> Filling Scint Shield FB with scint staves
    int num_x_fb_mid_staves_1 = 4;
    int num_y_fb_mid_staves_1 = 4;

    G4double x_fb_mid_staves_1 = scint_fb_mid_x/num_x_fb_mid_staves_1*mm;
    G4double y_fb_mid_staves_1 = scint_fb_mid_y/num_y_fb_mid_staves_1*mm;

    auto ScintShield_stave_fb_mid_1_S = new G4Box("ScintShield_fb_mid_S",x_fb_mid_staves_1/2.,y_fb_mid_staves_1/2.,scint_fb_z/2.);
    auto ScintShield_stave_fb_mid_1_LV = new G4LogicalVolume(ScintShield_stave_fb_mid_1_S,scintMaterial,"ScintShield_stave_fb_mid_1_LV");
    
    for (int i = 0;i<num_x_fb_mid_staves_1;i++){
      for (int j=0;j<num_y_fb_mid_staves_1;j++){
        new G4PVPlacement(0,G4ThreeVector(-scint_fb_mid_x/2.+((2.0*i+1)*x_fb_mid_staves_1/2.),-scint_fb_mid_y/2.+((2.0*j+1)*y_fb_mid_staves_1/2.),0.),ScintShield_stave_fb_mid_1_LV,"ScintShield",ScintShield_fb_mid_1_LV,false,0,fCheckOverlaps);
        new G4PVPlacement(0,G4ThreeVector(-scint_fb_mid_x/2.+((2.0*i+1)*x_fb_mid_staves_1/2.),-scint_fb_mid_y/2.+((2.0*j+1)*y_fb_mid_staves_1/2.),0.),ScintShield_stave_fb_mid_1_LV,"ScintShield",ScintShield_fb_mid_2_LV,false,0,fCheckOverlaps);
      }
    }

    // cuts for shield
    G4Region* shield_region = new G4Region("Shield_region");
    shield_region->AddRootLogicalVolume(LeadShield_top_LV);
    shield_region->AddRootLogicalVolume(LeadShield_side_LV);
    shield_region->AddRootLogicalVolume(LeadShield_fb_LV);
    shield_region->AddRootLogicalVolume(LeadShield_fb_mid_LV);

    shield_region->AddRootLogicalVolume(ScintShield_stave_top_1_LV);
    shield_region->AddRootLogicalVolume(ScintShield_stave_sides_1_LV);
    shield_region->AddRootLogicalVolume(ScintShield_stave_fb_1_LV);
    shield_region->AddRootLogicalVolume(ScintShield_stave_fb_mid_1_LV );
    G4Region* shield_Region = G4RegionStore::GetInstance()->GetRegion("Shield_region");
    G4ProductionCuts* shieldcut = new G4ProductionCuts();
    shieldcut->SetProductionCut(8.0*cm,"gamma"); // 130cm -> 1.3 MeV   220cm -> 3.4 MeV 400*cm -> 50 MeV
    shieldcut->SetProductionCut(5.0*mm,"e-");
    shieldcut->SetProductionCut(5.0*mm,"e+");
    shieldcut->SetProductionCut(15.0*mm,"proton");
    shield_Region->SetProductionCuts(shieldcut);
  #endif
  
  // **************************************************
  // ******** optics of the lead glass blocks  ********
  // **************************************************

  // --- the surface between the lead glass and the world 
  // --- Surface 1: should be highly reflective from glass to world    
  
  
  std::cout << "Setting up the optics ..." << std::endl;

  G4OpticalSurface* op_glass_world= new G4OpticalSurface("glass_World");
  op_glass_world->SetType(dielectric_metal);
  op_glass_world->SetFinish(polished);
  op_glass_world->SetModel(unified);
  
  G4double pp_glass_world[] = { 2.0 * eV, 3.5 * eV }; const G4int num_glass_world = sizeof(pp_glass_world) / sizeof(G4double);
  G4double reflectivity_glass_world[] = { 0.0, 0.0}; G4double efficiency_glass_world[] = { 0.0, 0.0 };
  G4MaterialPropertiesTable* GlassWorldProperty = new G4MaterialPropertiesTable();
  GlassWorldProperty->AddProperty("REFLECTIVITY", pp_glass_world, reflectivity_glass_world, num_glass_world);
  GlassWorldProperty->AddProperty("EFFICIENCY", pp_glass_world, efficiency_glass_world, num_glass_world);
  op_glass_world ->SetMaterialPropertiesTable(GlassWorldProperty);
  
  new G4LogicalSkinSurface("name",worldLV,op_glass_world);

  G4OpticalSurface* op_coating= new G4OpticalSurface("coating");
  op_coating->SetType(dielectric_metal);
  op_coating->SetFinish(polished);
  op_coating->SetModel(unified);
  
  G4double pp_coating[] = { 2.0 * eV, 3.5 * eV }; const G4int num_coating = sizeof(pp_coating) / sizeof(G4double);
  G4double reflectivity_coating[] = { 0.1, 0.1}; G4double efficiency_coating[] = { 0.0, 0.0 };
  G4MaterialPropertiesTable* CoatingProperty = new G4MaterialPropertiesTable();
  CoatingProperty->AddProperty("REFLECTIVITY", pp_coating, reflectivity_coating, num_coating);
  CoatingProperty->AddProperty("EFFICIENCY", pp_coating, efficiency_coating, num_coating);
  op_coating ->SetMaterialPropertiesTable(CoatingProperty);
  
  new G4LogicalSkinSurface("name",coatingLV,op_coating);
  //new G4LogicalSkinSurface("name",PMTLV,op_glass_world);
  
  //scintillators
  new G4LogicalSkinSurface("name",scintLV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_layerH_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_layerV_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_StaveH_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_StaveV_LV,op_glass_world);

  // scintillators front and back
  new G4LogicalSkinSurface("name",scint_fb1_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_fb2_LV,op_glass_world);

  new G4LogicalSkinSurface("name",scint_layer_fb11_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_layer_fb12_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_layer_fb21_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_layer_fb22_LV,op_glass_world);
  
  new G4LogicalSkinSurface("name",scint_stave_fb1H_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_stave_fb1V_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_stave_fb2H_LV,op_glass_world);
  new G4LogicalSkinSurface("name",scint_stave_fb2V_LV,op_glass_world);

  new G4LogicalSkinSurface("name",absorberLV,op_glass_world);
  
  #if VERSION_SHIELD == 1
    new G4LogicalSkinSurface("name",ScintShield_stave_top_1_LV,op_glass_world);
    new G4LogicalSkinSurface("name",worldLV,op_glass_world);
  #endif


  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  TPCLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  TPCLV_block->SetVisAttributes (G4VisAttributes::GetInvisible());
  leadglassLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  coatingLV->SetVisAttributes (G4VisAttributes::GetInvisible());
  auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true);
  auto green_color= new G4VisAttributes(G4Colour(0.788235,0.890196,0.741176)); green_color->SetVisibility(true);
  // original green: 0.517647,0.772549,0.556863
  auto orange_color= new G4VisAttributes(G4Colour(0.988235,0.686275,0.243137)); orange_color->SetVisibility(true);
  auto red_color= new G4VisAttributes(G4Colour(0.956863,0.0901961,0.494118)); red_color->SetVisibility(true); 
  auto blue_color= new G4VisAttributes(G4Colour(0.447059,0.623529,0.811765)); blue_color->SetVisibility(true);
  auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549)); black_color->SetVisibility(true);

  absorberLV->SetVisAttributes(black_color); // green
  scint_StaveH_LV->SetVisAttributes(blue_color); //blue
  scint_StaveV_LV->SetVisAttributes(blue_color); //blue
  scint_stave_fb1V_LV->SetVisAttributes(blue_color); //blue
  scint_stave_fb1H_LV->SetVisAttributes(blue_color); //blue
  scint_stave_fb2V_LV->SetVisAttributes(blue_color); //blue
  scint_stave_fb2H_LV->SetVisAttributes(blue_color); //blue
  TPCLV_1->SetVisAttributes(red_color); // red
  TPCLV_2->SetVisAttributes(red_color); // red
  siliconLV_1->SetVisAttributes(orange_color); //orange
  siliconLV_2->SetVisAttributes(orange_color);
  tubeLV -> SetVisAttributes(grey_color);
  carbonLV -> SetVisAttributes(grey_color);

  #if VERSION_SHIELD ==1
    // visual attributes for the shields 
    Measuring_plane_top_LV -> SetVisAttributes(green_color);
    LeadShield_top_LV->SetVisAttributes(black_color);
    LeadShield_side_LV->SetVisAttributes(black_color);
    LeadShield_fb_LV->SetVisAttributes(black_color);
    LeadShield_fb_mid_LV->SetVisAttributes(black_color);
    //LeadShield_fb_bot_LV->SetVisAttributes(black_color);

    ScintShield_top_1_LV->SetVisAttributes(blue_color); 
    ScintShield_top_2_LV->SetVisAttributes(blue_color);
    ScintShield_sides_1_LV->SetVisAttributes(blue_color);
    ScintShield_sides_2_LV->SetVisAttributes(blue_color);
    ScintShield_fb_1_LV->SetVisAttributes(blue_color);
    ScintShield_fb_2_LV->SetVisAttributes(blue_color);
    ScintShield_fb_mid_1_LV->SetVisAttributes(blue_color);
    ScintShield_fb_mid_2_LV->SetVisAttributes(blue_color);
    
    ScintShield_stave_top_1_LV->SetVisAttributes(red_color);
    ScintShield_stave_sides_1_LV->SetVisAttributes(red_color);
    ScintShield_stave_fb_1_LV->SetVisAttributes(red_color);
    ScintShield_stave_fb_mid_1_LV->SetVisAttributes(red_color);
  #endif

  return worldPV;
}

//....

void DetectorConstruction::ConstructSDandField()
{ 

  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare Scintillator as SinctillatorSD
  G4String scintDetectorName = "ScintLV" ;
  ScintillatorSD* scintDetector = new ScintillatorSD(scintDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(scintDetector);
  SetSensitiveDetector("Scint_StaveV_LV", scintDetector);
  SetSensitiveDetector("Scint_StaveH_LV", scintDetector);
  SetSensitiveDetector("Scint_fb_Stave1H_LV", scintDetector);
  SetSensitiveDetector("Scint_fb_Stave1V_LV", scintDetector);
  SetSensitiveDetector("Scint_fb_Stave2H_LV", scintDetector);
  SetSensitiveDetector("Scint_fb_Stave2V_LV", scintDetector);
  
  // declare absorber as AbsorberSD
  G4String absorberDetectorName = "AbsoLV" ;
  AbsorberSD* absorberDetector = new AbsorberSD(absorberDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(absorberDetector);
  SetSensitiveDetector("LeadGlassLV", absorberDetector);
  
  // declare vacuum as TubeSD
  G4String tubeDetectorName = "TubeLV" ;
  TubeSD* tubeDetector = new TubeSD(tubeDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(tubeDetector);
  SetSensitiveDetector("TubeLV", tubeDetector);
  
  // declare vacuum as TPCSD
  G4String TPCDetectorName = "TPCLV" ;
  TPCSD* TPCDetector = new TPCSD(TPCDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(TPCDetector);
  SetSensitiveDetector("TPCLV", TPCDetector);

  G4String Silicon_DetectorName = "siliconLV";
  SiliconSD* Silicon_Detector = new SiliconSD(Silicon_DetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(Silicon_Detector);
  SetSensitiveDetector("siliconLV_1", Silicon_Detector);
  SetSensitiveDetector("siliconLV_2", Silicon_Detector);

  G4String PMTDetectorName = "PMTLV" ;
  PMTSD* PMTDetector = new PMTSD(PMTDetectorName);
  G4SDManager::GetSDMpointer()->AddNewDetector(PMTDetector);
  //SetSensitiveDetector("PMTLV", PMTDetector);

  #if VERSION_SHIELD ==1
    G4String ShieldDetectorName = "ShieldLV" ;
    ShieldSD* shieldDetector = new ShieldSD(ShieldDetectorName);
    G4SDManager::GetSDMpointer()->AddNewDetector(shieldDetector);
    SetSensitiveDetector("LeadShield_top_LV", shieldDetector);
    SetSensitiveDetector("LeadShield_side_LV", shieldDetector);
    SetSensitiveDetector("LeadShield_fb_LV", shieldDetector);
    SetSensitiveDetector("LeadShield_fb_mid_LV", shieldDetector);
    SetSensitiveDetector("ScintShield_stave_top_1_LV", shieldDetector);
    //SetSensitiveDetector("ScintShield_top_2_LV", shieldDetector);
    SetSensitiveDetector("ScintShield_stave_sides_1_LV", shieldDetector);
    //SetSensitiveDetector("ScintShield_sides_2_LV", shieldDetector);  
    SetSensitiveDetector("ScintShield_stave_fb_1_LV", shieldDetector);
    //SetSensitiveDetector("ScintShield_fb_2_LV", shieldDetector);
    SetSensitiveDetector("ScintShield_stave_fb_mid_1_LV", shieldDetector);
    SetSensitiveDetector("Measuring_plane_top_LV", shieldDetector);
  #endif

  //SetSensitiveDetector("ScintShield_fb_mid_2_LV",shieldDetector);

  auto fElectricField = new ElectricField(); //new G4UniformElectricField(G4ThreeVector(0.0,0.0,400*volt/cm));
	G4EqMagElectricField*fEquation = new G4EqMagElectricField(fElectricField);  //fEMfield<->fElectricField
  
	G4int nvar = 8;
	fStepper = new G4DormandPrince745(fEquation,nvar);

	//global Feild Manager
	fFieldMgr = new G4FieldManager();
	fFieldMgr->SetDetectorField(fElectricField); //fEMfield<->fElectricField

	fMinStep = 0.01* nm; // minimal step of 10 microns //1*mm  0.0005?!

	G4FieldManager* globalFieldManager;
	G4TransportationManager* transportMgr = G4TransportationManager::GetTransportationManager();

  std::cout << "Electric field ...  " << std::endl;
	double MaxTrackingStep = 0.01;
	globalFieldManager = transportMgr->GetFieldManager();

	// Relative accuracy values:
	G4double minEps = 0.1 * nm;  //   Minimum & value for smallest steps  (1e-5)
	G4double maxEps = 10.0 *mm;  //   Maximum & value for largest steps  (1e-4)

	//fFieldMgr->SetDeltaOneStep(0.05*mm);  // 0.5 micrometer

	G4MagInt_Driver* fIntgrDriver = new G4MagInt_Driver(fMinStep,fStepper,fStepper->GetNumberOfVariables());
	fChordFinder = new G4ChordFinder(fIntgrDriver);
	fFieldMgr->SetChordFinder(fChordFinder);
  //TPCLV->SetFieldManager(fFieldMgr, true);
  //TPCLV_2->SetFieldManager(fFieldMgr, true);	
  transportMgr->GetPropagatorInField()->SetLargestAcceptableStep(1.*cm);


}

//....

