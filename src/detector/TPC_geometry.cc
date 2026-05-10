// TPC
#include "core/DetectorConstruction.hh"
#include "detector/TPC_geometry.hh"
#include "detector/beampipe_geometry.hh"

#include "G4SDManager.hh"
#include "sensitive/TPCSD.hh"

#include "G4UnionSolid.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4NistManager.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <G4ProductionCuts.hh>
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4UserLimits.hh"


TPC::TPC():G4VUserDetectorConstruction(){} //:G4VUserDetectorConstruction()
TPC::~TPC(){}

std::vector<G4LogicalVolume*> TPC_Construction_list;

// TPC dimensions
extern G4double Beampipe_5_len;
extern G4double Beampipe_5_radius_2;
extern G4double Beampipe_6_len;

G4double TPC_drift_len = 0.85*m;
G4double TPC_wall_thickness = 2.0*mm;
G4double TPC_z;// length of TPC in Z direction


void TPC::DefineMaterials()
{
  auto nistManager = G4NistManager::Instance();
  G4double density;

  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elAr = nistManager->FindOrBuildElement("Ar");

  // ============================================================================
  // Ar/CO2 80/20 TPC Gas Mixture
  // ============================================================================
  // Properties at STP (20°C, 1 atm):
  //   - Ar density: 1.66 mg/cm³
  //   - CO2 density: 1.84 mg/cm³
  //   - 80/20 mixture density: ~1.70 mg/cm³ (by volume fraction)
  //   - W-value (energy per ion pair): ~26 eV for Ar/CO2
  //   - Primary ionization: ~25-30 clusters/cm for MIP
  //   - Drift velocity: ~4 cm/µs at 400 V/cm
  // ============================================================================

  // CO2 at gas phase density (1.84 mg/cm³ at STP)
  G4Material* CO2 = new G4Material("CO2", density=1.84*mg/cm3, 2, kStateGas, 293.15*kelvin, 1.0*atmosphere);
  CO2->AddElement(elC, 1);
  CO2->AddElement(elO, 2);

  // Ar/CO2 80/20 by volume (common TPC gas mixture)
  // Volume fractions convert to mass fractions:
  //   - Ar: 80% volume -> ~78% by mass (Ar: 40 g/mol)
  //   - CO2: 20% volume -> ~22% by mass (CO2: 44 g/mol)
  // Effective density at STP: 1.70 mg/cm³
  G4Material* Gas = new G4Material("Gas", density=1.70*mg/cm3, 2, kStateGas, 293.15*kelvin, 1.0*atmosphere);
  Gas->AddElement(elAr, 0.78);     // ~78% by mass
  Gas->AddMaterial(CO2, 0.22);    // ~22% by mass

  // Set mean excitation energy for accurate dE/dx calculation
  // Ar: I = 188 eV, CO2: I = 85 eV
  // Weighted average for 80/20: ~167 eV
  Gas->GetIonisation()->SetMeanExcitationEnergy(167.0*eV);

  G4cout << "========================================" << G4endl;
  G4cout << "TPC Gas: Ar/CO2 (80/20) configured" << G4endl;
  G4cout << "  Density: " << Gas->GetDensity()/(mg/cm3) << " mg/cm³" << G4endl;
  G4cout << "  Mean excitation energy: 167 eV" << G4endl;
  G4cout << "  W-value: ~26 eV/ion-pair (physics)" << G4endl;
  G4cout << "========================================" << G4endl;
}

std::vector<G4LogicalVolume*> TPC::Construct_Volumes(G4LogicalVolume* mother)
{
  // OUR TPC GEOMETRY OUTPUT SENT BACK TO DETECTOR_CONSTRUCTION.CC //
  //std::vector<G4LogicalVolume*> TPC_Construction_list;
  //####################################################################//

  ////////////////////////////////////////////////////
  ///////////// Defining materials ///////////////////
  ////////////////////////////////////////////////////
    DefineMaterials();
    auto TPCMaterial = G4Material::GetMaterial("Gas");
    auto TPC_wall_Material = G4Material::GetMaterial("Aluminum");
  
  /////////////////////////////////////////////////////
  //////////// Numbers for TPC construction ///////////
  /////////////////////////////////////////////////////

    G4double TPC_drift_len = 0.85*m;
    G4double TPC_wall_thickness = 2.0*mm;
    G4double TPC_layer_thickness = 1.0*cm; // For TPC dEdx calculation
    G4double TPC_z = (Beampipe_5_len+2.*Beampipe_6_len)/2.; // 2.45*m; // length of TPC in Z direction

    //TPC type II -- the horizontal ones
    G4double TPC_x_2 = 2.0*Beampipe_5_radius_2+ 2.*TPC_wall_thickness; 
    G4double TPC_y_2 = TPC_drift_len + 2.*TPC_wall_thickness; 
    // TPC type II layers
    G4double TPC_layer_x_2 = TPC_x_2 - 2.*TPC_wall_thickness; //TPC_x_2 - 2.*TPC_wall_thickness
    G4double TPC_layer_y_2 = TPC_layer_thickness;
    G4double TPC_layer_z_2 = TPC_z - 2.0*TPC_wall_thickness;
    G4int TPC_layer_2_number = TPC_drift_len/TPC_layer_thickness;


    //TPC type I -- the vertical ones
    G4double TPC_x_1 = TPC_drift_len + 2.0*TPC_wall_thickness;
    G4double TPC_y_1 = (2.0*Beampipe_5_radius_2+2.0*TPC_y_2)/2.;
    // TPC type I layers
    G4double TPC_layer_x_1 = TPC_layer_thickness;
    G4double TPC_layer_y_1 = TPC_y_1 - 2.0*TPC_wall_thickness;
    G4double TPC_layer_z_1 = TPC_z - 2.0*TPC_wall_thickness;
    G4int TPC_layer_1_number = TPC_drift_len/TPC_layer_thickness;

  ///////////////////////////////////////////////////////////
  ///////////// TPC geometry construction ///////////////////
  ///////////////////////////////////////////////////////////

    //vectrical TPCs
    auto TPCS_1 = new G4Box("TPCS_1",TPC_x_1/2.,TPC_y_1/2.,TPC_z/2.); 
    auto TPCLV_1 = new G4LogicalVolume(TPCS_1,TPC_wall_Material,"TPCLV_1");
        auto TPC_1_layer_S = new G4Box("TPCS_1_layer",TPC_layer_x_1/2.,TPC_layer_y_1/2.,TPC_layer_z_1/2.);
        auto TPC_1_layer_LV = new G4LogicalVolume(TPC_1_layer_S,TPCMaterial,"TPC_1_layer_LV");
        for (int i=0;i<TPC_layer_1_number;i++){
          G4ThreeVector TPC_layer_pos = G4ThreeVector(-TPC_drift_len/2.+(double(i)+0.5)*TPC_layer_thickness,0.,0.); 
          new G4PVPlacement(0,TPC_layer_pos,TPC_1_layer_LV,"TPC_1_layer_PV",TPCLV_1,false,i,false);
        }

    std::cout << Beampipe_5_radius_2 <<"," <<TPC_x_2/2.<<" "<< TPC_y_2/2.<<" "<<TPC_z/2.<< "TPC_layer_x_2::"<<TPC_layer_x_2<<std::endl;
    //horizontal TPCs
    auto TPCS_2 = new G4Box("TPCS_2",TPC_x_2/2.,TPC_y_2/2.,TPC_z/2.); 
    auto TPCLV_2 = new G4LogicalVolume(TPCS_2,TPC_wall_Material,"TPCLV_2");
        auto TPC_2_layer_S = new G4Box("TPCS_2_layer",TPC_layer_x_2/2.,TPC_layer_y_2/2.,TPC_layer_z_2/2.);
        auto TPC_2_layer_LV = new G4LogicalVolume(TPC_2_layer_S,TPCMaterial,"TPC_2_layer_LV");
        for (int i=0;i<TPC_layer_2_number;i++){
          G4ThreeVector TPC_layer_pos = G4ThreeVector(0.,-TPC_drift_len/2.+(double(i)+0.5)*TPC_layer_thickness,0.); 
          new G4PVPlacement(0,TPC_layer_pos,TPC_2_layer_LV,"TPC_2_layer_PV",TPCLV_2,false,i,false);
        }
    

    //------ front TPCs
    auto TPC_pos1 = G4ThreeVector(0.,Beampipe_5_radius_2+TPC_y_2/2.,-TPC_z/2.);
    auto TPC_pos2 = G4ThreeVector(-(TPC_x_2/2.+TPC_x_1/2.),TPC_y_1/2.,-TPC_z/2.);
    auto TPC_pos3 = G4ThreeVector(-(TPC_x_2/2.+TPC_x_1/2.),-TPC_y_1/2.,-TPC_z/2.);  
    auto TPC_pos4 = G4ThreeVector(0.,-(Beampipe_5_radius_2+TPC_y_2/2.),-TPC_z/2.);
    auto TPC_pos5 = G4ThreeVector((TPC_x_2/2.+TPC_x_1/2.),-TPC_y_1/2.,-TPC_z/2.);
    auto TPC_pos6 = G4ThreeVector((TPC_x_2/2.+TPC_x_1/2.),TPC_y_1/2.,-TPC_z/2.);

    //------ back TPCs
    auto TPC_pos7 = G4ThreeVector(0.,Beampipe_5_radius_2+TPC_y_2/2.,TPC_z/2.);
    auto TPC_pos8 = G4ThreeVector((TPC_x_2/2.+TPC_x_1/2.),TPC_y_1/2.,TPC_z/2.);
    auto TPC_pos9 = G4ThreeVector((TPC_x_2/2.+TPC_x_1/2.),-TPC_y_1/2.,TPC_z/2.);
    auto TPC_pos10 = G4ThreeVector(0.,-(Beampipe_5_radius_2+TPC_y_2/2.),TPC_z/2.);
    auto TPC_pos11 = G4ThreeVector(-(TPC_x_2/2.+TPC_x_1/2.),-TPC_y_1/2.,TPC_z/2.);
    auto TPC_pos12 = G4ThreeVector(-(TPC_x_2/2.+TPC_x_1/2.),TPC_y_1/2.,TPC_z/2.);
    
    new G4PVPlacement(0,TPC_pos1,TPCLV_2,"TPCPV",mother,false,0,false);
    new G4PVPlacement(0,TPC_pos2,TPCLV_1,"TPCPV",mother,false,1,false);
    new G4PVPlacement(0,TPC_pos3,TPCLV_1,"TPCPV",mother,false,2,false);
    new G4PVPlacement(0,TPC_pos4,TPCLV_2,"TPCPV",mother,false,3,false);
    new G4PVPlacement(0,TPC_pos5,TPCLV_1,"TPCPV",mother,false,4,false);
    new G4PVPlacement(0,TPC_pos6,TPCLV_1,"TPCPV",mother,false,5,false);

    new G4PVPlacement(0,TPC_pos7,TPCLV_2,"TPCPV",mother,false,6,false);
    new G4PVPlacement(0,TPC_pos8,TPCLV_1,"TPCPV",mother,false,7,false);
    new G4PVPlacement(0,TPC_pos9,TPCLV_1,"TPCPV",mother,false,8,false);
    new G4PVPlacement(0,TPC_pos10,TPCLV_2,"TPCPV",mother,false,9,false);
    new G4PVPlacement(0,TPC_pos11,TPCLV_1,"TPCPV",mother,false,10,false);
    new G4PVPlacement(0,TPC_pos12,TPCLV_1,"TPCPV",mother,false,11,false);

  //////////////////////////////////////////////////////////////////////
  //////////////////// Writing the outputs /////////////////////////////
  //////////////////////////////////////////////////////////////////////

    TPC_Construction_list.push_back(TPCLV_1);
    TPC_Construction_list.push_back(TPCLV_2);
    TPC_Construction_list.push_back(TPC_1_layer_LV);
    TPC_Construction_list.push_back(TPC_2_layer_LV);
  
  //////////////////////////////////////////////////////////////////////
  //////////////////// Region Settings /////////////////////////////////
  //////////////////////////////////////////////////////////////////////

    //Region for the TPC wall
    G4Region* TPC_wall_region = new G4Region("TPC_wall_region");
    TPC_wall_region->AddRootLogicalVolume(TPCLV_1);
    TPC_wall_region->AddRootLogicalVolume(TPCLV_2);
    G4Region* TPC_wall_region_ = G4RegionStore::GetInstance()->GetRegion("TPC_wall_region");
    G4ProductionCuts* TPC_wall_cut = new G4ProductionCuts();
    TPC_wall_cut->SetProductionCut(0.1*cm,"gamma");
    TPC_wall_cut->SetProductionCut(0.1*mm,"e-");
    TPC_wall_region_->SetProductionCuts(TPC_wall_cut);

    // User limits for TPC wall to kill very low energy particles
    G4UserLimits* tpcWallLimits = new G4UserLimits();
    tpcWallLimits->SetUserMinEkine(1.0*keV);  // Minimum kinetic energy
    tpcWallLimits->SetUserMaxTime(10.0*ms);    // Max tracking time
    TPCLV_1->SetUserLimits(tpcWallLimits);
    TPCLV_2->SetUserLimits(tpcWallLimits);

    // Region for the TPC gas volume (Ar/CO2)
    // Production cuts define the minimum range for secondary particle production
    G4Region* TPC_region = new G4Region("TPC_region");
    TPC_region->AddRootLogicalVolume(TPC_1_layer_LV);
    TPC_region->AddRootLogicalVolume(TPC_2_layer_LV);

    G4Region* TPC_Region_ = G4RegionStore::GetInstance()->GetRegion("TPC_region");
    G4ProductionCuts* TPCcut = new G4ProductionCuts();
    TPCcut->SetProductionCut(0.1*mm, "gamma");  // Reasonable for gas detector
    TPCcut->SetProductionCut(0.1*mm, "e-");     // Fixed: was 0.00001nm (unphysical)
    TPCcut->SetProductionCut(0.1*mm, "e+");
    TPCcut->SetProductionCut(0.1*mm, "proton");
    TPC_Region_->SetProductionCuts(TPCcut);

    G4cout << "TPC: Production cuts set to 0.1 mm for gas region" << G4endl;

    // User limits to prevent tracking of very low energy particles
    // This kills particles below 1 keV in TPC gas volume to avoid infinite tracking
    // at sub-keV energies where physics is not meaningful
    G4UserLimits* tpcLimits = new G4UserLimits();
    tpcLimits->SetUserMinEkine(1.0*keV);  // Minimum kinetic energy
    tpcLimits->SetUserMaxTime(10.0*ms);   // Max tracking time per track
    TPC_1_layer_LV->SetUserLimits(tpcLimits);
    TPC_2_layer_LV->SetUserLimits(tpcLimits);
    G4cout << "TPC: User limits set - min KE 1 keV, max time 10 ms" << G4endl;

  ////////////////////////////////////////////////////
  ///////////////// Color settings ///////////////////
  ////////////////////////////////////////////////////

    auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549));black_color->SetVisibility(true);
    auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true);
    auto green_color= new G4VisAttributes(G4Colour(0.788235,0.890196,0.741176)); green_color->SetVisibility(true);
    auto orange_color= new G4VisAttributes(G4Colour(0.988235,0.686275,0.243137)); orange_color->SetVisibility(true);
    auto purple_color= new G4VisAttributes(G4Colour(0.48,0.27,0.833)); purple_color->SetVisibility(true);
    
    TPCLV_1 -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(black_color);
    TPCLV_2 -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(black_color);
    TPC_1_layer_LV ->SetVisAttributes (G4VisAttributes::GetInvisible());
    TPC_2_layer_LV ->SetVisAttributes (G4VisAttributes::GetInvisible());
  
  return TPC_Construction_list;
}