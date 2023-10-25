// TPC
#include "DetectorConstruction.hh"
#include "Detector_Module/TPC_geometry.hh"
#include "Detector_Module/beampipe_geometry.hh"

#include "G4SDManager.hh"
#include "Sensitive_Detector/TPCSD.hh"

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
  G4double a; G4double z; G4double density;

  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elAr = nistManager->FindOrBuildElement("Ar");
  G4Element* elAl = nistManager->FindOrBuildElement("Al");
  // CO2
  G4Material* CO2 = new G4Material("CO2", density= 1.98*g/cm3, 2);
  CO2->AddElement(elO, 2); CO2->AddElement(elC, 1);
  // Ar/CO2 80/20
  G4Material* Gas = new G4Material("Gas", density=1.823*mg/cm3, 2);
  Gas->AddElement(elAr, .8); Gas->AddMaterial(CO2, .2);
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
    auto defaultMaterial = G4Material::GetMaterial("Galactic");
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
    G4double TPC_total_length = 2.0 * TPC_x_1 + TPC_x_2; // beware of this value, may affect the lead glass and scintillators
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
          new G4PVPlacement(0,TPC_layer_pos,TPC_1_layer_LV,"TPC_1_layer_PV",TPCLV_1,false,i,true);
        }

    std::cout << Beampipe_5_radius_2 <<"," <<TPC_x_2/2.<<" "<< TPC_y_2/2.<<" "<<TPC_z/2.<< "TPC_layer_x_2::"<<TPC_layer_x_2<<std::endl;
    //horizontal TPCs
    auto TPCS_2 = new G4Box("TPCS_2",TPC_x_2/2.,TPC_y_2/2.,TPC_z/2.); 
    auto TPCLV_2 = new G4LogicalVolume(TPCS_2,TPC_wall_Material,"TPCLV_2");
        auto TPC_2_layer_S = new G4Box("TPCS_2_layer",TPC_layer_x_2/2.,TPC_layer_y_2/2.,TPC_layer_z_2/2.);
        auto TPC_2_layer_LV = new G4LogicalVolume(TPC_2_layer_S,TPCMaterial,"TPC_2_layer_LV");
        for (int i=0;i<TPC_layer_2_number;i++){
          G4ThreeVector TPC_layer_pos = G4ThreeVector(0.,-TPC_drift_len/2.+(double(i)+0.5)*TPC_layer_thickness,0.); 
          new G4PVPlacement(0,TPC_layer_pos,TPC_2_layer_LV,"TPC_2_layer_PV",TPCLV_2,false,i,true);
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
    
    new G4PVPlacement(0,TPC_pos1,TPCLV_2,"TPCPV",mother,false,0,true);
    new G4PVPlacement(0,TPC_pos2,TPCLV_1,"TPCPV",mother,false,1,true);
    new G4PVPlacement(0,TPC_pos3,TPCLV_1,"TPCPV",mother,false,2,true);
    new G4PVPlacement(0,TPC_pos4,TPCLV_2,"TPCPV",mother,false,3,true);
    new G4PVPlacement(0,TPC_pos5,TPCLV_1,"TPCPV",mother,false,4,true);
    new G4PVPlacement(0,TPC_pos6,TPCLV_1,"TPCPV",mother,false,5,true);

    new G4PVPlacement(0,TPC_pos7,TPCLV_2,"TPCPV",mother,false,6,true);
    new G4PVPlacement(0,TPC_pos8,TPCLV_1,"TPCPV",mother,false,7,true);
    new G4PVPlacement(0,TPC_pos9,TPCLV_1,"TPCPV",mother,false,8,true);
    new G4PVPlacement(0,TPC_pos10,TPCLV_2,"TPCPV",mother,false,9,true);
    new G4PVPlacement(0,TPC_pos11,TPCLV_1,"TPCPV",mother,false,10,true);
    new G4PVPlacement(0,TPC_pos12,TPCLV_1,"TPCPV",mother,false,11,true);

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

    //Region for the TPC LiAr volume
    G4Region* TPC_region = new G4Region("TPC_region");
    TPC_region->AddRootLogicalVolume(TPC_1_layer_LV);
    TPC_region->AddRootLogicalVolume(TPC_2_layer_LV);

    G4Region* TPC_Region_ = G4RegionStore::GetInstance()->GetRegion("TPC_region");
    G4ProductionCuts* TPCcut = new G4ProductionCuts();
    TPCcut->SetProductionCut(0.0001*cm,"gamma");
    TPCcut->SetProductionCut(0.00001*nm,"e-");
    TPC_Region_->SetProductionCuts(TPCcut);


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