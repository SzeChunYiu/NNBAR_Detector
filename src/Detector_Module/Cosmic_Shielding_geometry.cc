// Cosmic_Shielding
#include "DetectorConstruction.hh"
#include "Detector_Module/Cosmic_Shielding_geometry.hh"
#include "Detector_Module/beampipe_geometry.hh"

#include "G4SDManager.hh"
//#include "Sensitive_Detector/Cosmic_ShieldingSD.hh"

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

CosmicShielding::CosmicShielding():G4VUserDetectorConstruction(){} //:G4VUserDetectorConstruction()
CosmicShielding::~CosmicShielding(){}

std::vector<G4LogicalVolume*> CosmicShielding_Construction_list;

extern G4double Beampipe_5_radius_2;
extern G4double Beampipe_4_radius_2;

/////////////////////////////////////////////////////////////////
//////////// Numbers for CosmicShielding construction ///////////
/////////////////////////////////////////////////////////////////
    
    // shield structure x=y cross section
    // ########%%%%%%%%%%%%%%%%%%%%%%%########
    // ########%%%%%%%%%%%%%%%%%%%%%%%########
    // ########***********************########
    // ########&                     &######## 
    // ########&                     &########
    // ########&                     &########
    // ########&                     &########
    // ########&                     &########
    // ########&                     &########
    // ########***********************########
    // ########%%%%%%%%%%%%%%%%%%%%%%%########
    // ########%%%%%%%%%%%%%%%%%%%%%%%########

    // distance of the shield from lead glass
    G4double shield_offset = 50.0*cm;
    G4double shield_thickness = 2.*m;

    // Distance from center to the side Lead Glass Modules
    G4double dist_center_lead_glass = 2.75*m;
    
    // Distance from center to the side Lead Glass Modules
    G4double dist_center_lead_glass_fb = 3.35*m;

    // Lead Shield on top surface
    G4double lead_top_x = 2.0*dist_center_lead_glass+2.0*shield_offset; 
    G4double lead_top_y = shield_thickness; 
    G4double lead_top_z = 2.0*dist_center_lead_glass_fb+2.0*shield_offset;

    // Lead shield Sides  
    G4double lead_side_x = shield_thickness; 
    G4double lead_side_y = 2.0*(dist_center_lead_glass+shield_offset+lead_top_y); 
    G4double lead_side_z = lead_top_z;

    //Lead Shield FB 
    G4double lead_fb_x = lead_top_x +2.0*lead_side_x;
    G4double lead_fb_y = 2.*(dist_center_lead_glass +shield_offset+lead_top_y);
    G4double lead_fb_z = shield_thickness;

    // Steel shield (inner layer)
    G4double steel_thickness = 30.0*cm;
    G4double steel_top_x = lead_top_x; 
    G4double steel_top_y = steel_thickness; 
    G4double steel_top_z = lead_top_z;

    // Steel shield side (inner layer)
    G4double steel_side_x = steel_thickness; 
    G4double steel_side_y = lead_side_y-2*lead_top_y-2*steel_top_y; 
    G4double steel_side_z = lead_side_z;

    // Steel shield fb (inner layer)
    G4double steel_fb_x = lead_fb_x-2.*(steel_thickness+shield_thickness); 
    G4double steel_fb_y = lead_fb_y-2.*(steel_thickness+shield_thickness); 
    G4double steel_fb_z = steel_thickness;

void CosmicShielding::DefineMaterials()
{
  auto nistManager = G4NistManager::Instance();
  G4double a; G4double z; G4double density;
  
  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elCa = nistManager->FindOrBuildElement("Ca");
  G4Element* elB = nistManager->FindOrBuildElement("B");
  G4Element* elTi = nistManager->FindOrBuildElement("Ti");
  G4Element* elAs = nistManager->FindOrBuildElement("As");
  G4Element* elPb = nistManager->FindOrBuildElement("Pb");
  G4Element* elO = nistManager->FindOrBuildElement("O");
  G4Element* elK = nistManager->FindOrBuildElement("K");
  G4Element* elSi = nistManager->FindOrBuildElement("Si");
  G4Element* elCl = nistManager->FindOrBuildElement("Cl");
  G4Element* elS = nistManager->FindOrBuildElement("S");
  G4Element* elN = nistManager->FindOrBuildElement("N");
  G4Element* elAr = nistManager->FindOrBuildElement("Ar");
  G4Element* elAl = nistManager->FindOrBuildElement("Al");
  G4Element* elMg = nistManager->FindOrBuildElement("Mg");
  G4Element* elMn = nistManager->FindOrBuildElement("Mn");
  G4Element* elNa = nistManager->FindOrBuildElement("Na");
  G4Element* elNi = nistManager->FindOrBuildElement("Ni");
  G4Element* elCr = nistManager->FindOrBuildElement("Cr");
  G4Element* elF = nistManager->FindOrBuildElement("F");
  G4Element* elFe = nistManager->FindOrBuildElement("Fe");
  G4Element* elP = nistManager->FindOrBuildElement("P");
  G4Element* elBe = nistManager->FindOrBuildElement("Be");
  G4Element* elCd = nistManager->FindOrBuildElement("Cd");
  G4Element* elLi = nistManager->FindOrBuildElement("Li");

  G4Material* Lead = new G4Material("Lead", z=82., a= 207.2*g/mole, density = 11.29*g/cm3);

  G4Material* PE_B4C_concrete = new G4Material("PE_B4C_concrete", density = 1.97*g/cm3,15);
  PE_B4C_concrete->AddElement(elO,0.4606); PE_B4C_concrete->AddElement(elCa,0.0805);
  PE_B4C_concrete->AddElement(elSi,0.2840); PE_B4C_concrete->AddElement(elAl,0.0234);
  PE_B4C_concrete->AddElement(elFe,0.00837); PE_B4C_concrete->AddElement(elMg,0.00195);
  PE_B4C_concrete->AddElement(elNa,0.00613); PE_B4C_concrete->AddElement(elK,0.0125);
  PE_B4C_concrete->AddElement(elS,0.00276); PE_B4C_concrete->AddElement(elCl,0.000353);
  PE_B4C_concrete->AddElement(elH,0.02362); PE_B4C_concrete->AddElement(elTi,0.00517);
  PE_B4C_concrete->AddElement(elP,0.00259); PE_B4C_concrete->AddElement(elC,0.0893);
  PE_B4C_concrete->AddElement(elB,0.0596);

  // Concrete material
  G4Material* MagnadenseHC = new G4Material("MagnadenseHC",   density = 3.8*g/cm3,10);
  MagnadenseHC->AddElement(elH, 0.0053);
  MagnadenseHC->AddElement(elO, 0.332);
  MagnadenseHC->AddElement(elNa, 0.0046);
  MagnadenseHC->AddElement(elAl, 0.0064);
  MagnadenseHC->AddElement(elSi, 0.0469);
  MagnadenseHC->AddElement(elP, 0.0044);
  MagnadenseHC->AddElement(elS, 0.0003);
  MagnadenseHC->AddElement(elK, 0.0015);
  MagnadenseHC->AddElement(elCa, 0.0198);
  MagnadenseHC->AddElement(elFe, 0.579);

  G4Material* StainlessSteel = new G4Material("StainlessSteel",   density = 8.02*g/cm3,5);
  StainlessSteel->AddElement(elMn, 0.02);
  StainlessSteel->AddElement(elSi, 0.01);
  StainlessSteel->AddElement(elCr, 0.19);
  StainlessSteel->AddElement(elNi, 0.10);
  StainlessSteel->AddElement(elFe, 0.68);

}

std::vector<G4LogicalVolume*> CosmicShielding::Construct_Volumes(G4LogicalVolume* mother)
{
  // OUR CosmicShielding GEOMETRY OUTPUT SENT BACK TO DETECTOR_CONSTRUCTION.CC //
  //std::vector<G4LogicalVolume*> CosmicShielding_Construction_list;
  //####################################################################//

  ////////////////////////////////////////////////////
  ///////////// Defining materials ///////////////////
  ////////////////////////////////////////////////////
    
    DefineMaterials();
    auto defaultMaterial = G4Material::GetMaterial("Galactic");
    auto CosmicShieldingMaterial = G4Material::GetMaterial("MagnadenseHC"); //PE_B4C_concrete
    auto SteelShieldingMaterial = G4Material::GetMaterial("StainlessSteel"); //PE_B4C_concrete

  ///////////////////////////////////////////////////////////////////////
  ///////////// CosmicShielding geometry construction ///////////////////
  ///////////////////////////////////////////////////////////////////////

    // Lead Shield top & bottom
    auto LeadShield_top_S = new G4Box("LeadShield_top_S",lead_top_x/2., lead_top_y/2., lead_top_z/2.);
    auto LeadShield_top_LV = new G4LogicalVolume(LeadShield_top_S,CosmicShieldingMaterial,"LeadShield_top_LV");
    auto LeadShield_top_y = dist_center_lead_glass+shield_offset+lead_top_y/2.;
    new G4PVPlacement(0,G4ThreeVector(0.,LeadShield_top_y,0.),LeadShield_top_LV,"LeadShield",mother,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0.,-LeadShield_top_y,0.),LeadShield_top_LV,"LeadShield",mother,false,0,true);

    auto SteelShield_top_S = new G4Box("SteelShield_top_S",steel_top_x/2., steel_top_y/2., steel_top_z/2.);
    auto SteelShield_top_LV = new G4LogicalVolume(SteelShield_top_S,SteelShieldingMaterial,"SteelShield_top_LV");
    auto SteelShield_top_y = LeadShield_top_y-lead_top_y/2.-steel_top_y/2.;
    new G4PVPlacement(0,G4ThreeVector(0.,SteelShield_top_y,0.),SteelShield_top_LV,"SteelShield",mother,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(0.,-SteelShield_top_y,0.),SteelShield_top_LV,"SteelShield",mother,false,0,true);

    // ===== 2 Sides ===== //
    auto LeadShield_side_S = new G4Box("LeadShield_side_S",lead_side_x/2., lead_side_y/2., lead_side_z/2.);
    auto LeadShield_side_LV = new G4LogicalVolume(LeadShield_side_S,CosmicShieldingMaterial,"LeadShield_side_LV");
    auto LeadShield_side_x = (lead_top_x/2.+lead_side_x/2.);
    new G4PVPlacement(0,G4ThreeVector(LeadShield_side_x,0.,0.),LeadShield_side_LV,"LeadShield_side",mother,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(-(LeadShield_side_x),0.,0.),LeadShield_side_LV,"LeadShield_side",mother,false,0,true);

    auto SteelShield_side_S = new G4Box("SteelShield_side_S",steel_side_x/2., steel_side_y/2., steel_side_z/2.);
    auto SteelShield_side_LV = new G4LogicalVolume(SteelShield_side_S,SteelShieldingMaterial,"SteelShield_side_LV");
    auto SteelShield_side_x = LeadShield_side_x - lead_side_x/2. - steel_side_x/2.;
    new G4PVPlacement(0,G4ThreeVector(SteelShield_side_x,0.,0.),SteelShield_side_LV,"SteelShield_side",mother,false,0,true);
    new G4PVPlacement(0,G4ThreeVector(-(SteelShield_side_x),0.,0.),SteelShield_side_LV,"SteelShield_side",mother,false,0,true);


    // ===== Front back surfaces
      // Define the shield board
        auto LeadShield_fb_temp_S = new G4Box("LeadShield_fb_temp_S",lead_fb_x/2., lead_fb_y/2., lead_fb_z/2.);
      // Need to make a hole on the surface
        auto virtual_beampipe_4_S = new G4Cons("virtual_Beampipe_4_S", 
                                                0.,Beampipe_4_radius_2,
                                                0.,Beampipe_4_radius_2,
                                                shield_thickness/2.,0.,360.*deg);

        auto LeadShield_fb_S = new G4SubtractionSolid("LeadShield_fb_S", LeadShield_fb_temp_S, virtual_beampipe_4_S,0, G4ThreeVector(0.,0.,0.));
        auto LeadShield_fb_LV = new G4LogicalVolume(LeadShield_fb_S,CosmicShieldingMaterial,"LeadShield_fb_LV");
        auto LeadShield_fb_z = (lead_top_z/2.+lead_fb_z/2.);
        new G4PVPlacement(0,G4ThreeVector(0.,0.,LeadShield_fb_z),LeadShield_fb_LV,"LeadShield_front",mother,false,0,true);
        new G4PVPlacement(0,G4ThreeVector(0.,0.,-LeadShield_fb_z),LeadShield_fb_LV,"LeadShield_back",mother,false,0,true);  
  
        auto SteelShield_fb_temp_S = new G4Box("SteelShield_fb_temp_S",steel_fb_x/2., steel_fb_y/2., steel_fb_z/2.);
      // Need to make a hole on the surface
        // auto virtual_beampipe_4_S = new G4Cons("virtual_Beampipe_4_S", 
        //                                         0.,Beampipe_4_radius_2,
        //                                         0.,Beampipe_4_radius_2,
        //                                         shield_thickness/2.,0.,360.*deg);

        auto SteelShield_fb_S = new G4SubtractionSolid("SteelShield_fb_S", SteelShield_fb_temp_S, virtual_beampipe_4_S,0, G4ThreeVector(0.,0.,0.));
        auto SteelShield_fb_LV = new G4LogicalVolume(SteelShield_fb_S,CosmicShieldingMaterial,"SteelShield_fb_LV");
        auto SteelShield_fb_z = LeadShield_fb_z - lead_fb_z/2. - steel_fb_z/2.; 
        new G4PVPlacement(0,G4ThreeVector(0.,0.,SteelShield_fb_z),SteelShield_fb_LV,"SteelShield_front",mother,false,0,true);
        new G4PVPlacement(0,G4ThreeVector(0.,0.,-SteelShield_fb_z),SteelShield_fb_LV,"SteelShield_back",mother,false,1,true);  

  //////////////////////////////////////////////////////////////////////
  //////////////////// Writing the outputs /////////////////////////////
  //////////////////////////////////////////////////////////////////////

    CosmicShielding_Construction_list.push_back(LeadShield_top_LV);
    CosmicShielding_Construction_list.push_back(LeadShield_side_LV);
    CosmicShielding_Construction_list.push_back(LeadShield_fb_LV);
    CosmicShielding_Construction_list.push_back(SteelShield_top_LV);
    CosmicShielding_Construction_list.push_back(SteelShield_side_LV);
    CosmicShielding_Construction_list.push_back(SteelShield_fb_LV);

  //////////////////////////////////////////////////////////////////////
  //////////////////// Region Settings /////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  
    // Region for cosmic shielding
    G4Region* shield_region = new G4Region("Shield_region");
    shield_region->AddRootLogicalVolume(LeadShield_top_LV);
    shield_region->AddRootLogicalVolume(LeadShield_side_LV);
    shield_region->AddRootLogicalVolume(LeadShield_fb_LV);
    shield_region->AddRootLogicalVolume(SteelShield_top_LV);
    shield_region->AddRootLogicalVolume(SteelShield_side_LV);
    shield_region->AddRootLogicalVolume(SteelShield_fb_LV);
      
      
    G4Region* shield_Region = G4RegionStore::GetInstance()->GetRegion("Shield_region");
    G4ProductionCuts* shieldcut = new G4ProductionCuts();
    shieldcut->SetProductionCut(5.0*cm,"gamma"); // 130cm -> 1.3 MeV   220cm -> 3.4 MeV 400*cm -> 50 MeV
    shieldcut->SetProductionCut(5.0*mm,"e-");
    shieldcut->SetProductionCut(5.0*mm,"e+");
    shieldcut->SetProductionCut(15.0*mm,"proton");
    shield_Region->SetProductionCuts(shieldcut);

  ////////////////////////////////////////////////////
  ///////////////// Color settings ///////////////////
  ////////////////////////////////////////////////////

    auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549));black_color->SetVisibility(true);
    auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true);
    auto green_color= new G4VisAttributes(G4Colour(0.788235,0.890196,0.741176)); green_color->SetVisibility(true);

    LeadShield_top_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(black_color);
    LeadShield_side_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(black_color);
    LeadShield_fb_LV -> SetVisAttributes (green_color);// SetVisAttributes(black_color);
    SteelShield_fb_LV-> SetVisAttributes (black_color);
    
  return CosmicShielding_Construction_list;
}