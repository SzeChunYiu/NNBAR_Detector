// Beampipe
// Beampipe geometry up to 2023.01.16 -> Yuri's drawing on the beampipe with different sections
// Beampipe geometry labelling and dimensions please find the supplementary document ** (Missing)
// Beampipe sections are coated with neutron absorber

#include "DetectorConstruction.hh"
#include "Detector_Module/beampipe_shielding_geometry.hh"
#include "Detector_Module/Cosmic_Shielding_geometry.hh"

#include "G4SDManager.hh"
#include "Sensitive_Detector/TubeSD.hh"

#include "G4NistManager.hh"
#include "G4Isotope.hh"
#include "G4Material.hh"
#include "G4Element.hh"

#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4AutoDelete.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

Beampipe_Shielding::Beampipe_Shielding(){}

Beampipe_Shielding::~Beampipe_Shielding(){}

extern G4double Beampipe_2_radius_2;
extern G4double Beampipe_2_len;
extern G4double Beampipe_2_len;
extern G4double Beampipe_2_pos_z;

extern G4double Beampipe_4_radius_2;
extern G4double Beampipe_4_len;
extern G4double Beampipe_4_pos_z;

extern G4double Beampipe_8_radius_2;
extern G4double Beampipe_8_len;
extern G4double Beampipe_8_pos_z;
extern G4double lead_top_z;
extern G4double lead_fb_z;

void Beampipe_Shielding::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  G4double a; G4double z; G4double density;
  G4Material* Lead = new G4Material("Lead_shield", z=82., a= 207.2*g/mole, density = 11.29*g/cm3);
}

std::vector<G4LogicalVolume*> Beampipe_Shielding::Construct_Volumes(G4LogicalVolume* mother)
{
  // OUR Beampipe Shielding GEOMETRY OUTPUT SENT BACK TO DETECTOR_CONSTRUCTION.CC //
  std::vector<G4LogicalVolume*> Beampipe_Shielding_Construction_list;
  //####################################################################//

  DefineMaterials();
  auto DefaultMaterial = G4Material::GetMaterial("Galactic");
  auto Beampipe_Shielding_Material = G4Material::GetMaterial("Lead_shield"); // Aluminum , el_Be
  

  //////////////////////////////////////////////////////////////////////
  //////////// Numbers for Beam shielding construction /////////////////
  //////////////////////////////////////////////////////////////////////
      
      G4double Beampipe_Shielding_2_thickness = 3.0*m;
      G4double Beampipe_Shielding_2_radius_1 = Beampipe_2_radius_2;
      G4double Beampipe_Shielding_2_radius_2 = Beampipe_2_radius_2 + Beampipe_Shielding_2_thickness;
      G4double Beampipe_Shielding_2_len = Beampipe_2_len;

      G4double Beampipe_Shielding_4_thickness = 3.5*m;
      G4double Beampipe_Shielding_4_radius_1 = Beampipe_4_radius_2;
      G4double Beampipe_Shielding_4_radius_2 = Beampipe_4_radius_2 + Beampipe_Shielding_4_thickness;
      G4double Beampipe_Shielding_4_len = (-(lead_top_z/2.+lead_fb_z))-(Beampipe_4_pos_z - Beampipe_4_len/2.);
      std::cout <<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<< std::endl;
      std::cout <<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<< std::endl;
      std::cout <<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<< std::endl;
      std::cout <<":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::"<< std::endl; 
      std::cout << Beampipe_Shielding_4_len << "::"<< lead_top_z << "," << lead_fb_z << "," << Beampipe_4_pos_z - Beampipe_4_len/2. << std::endl; 
      G4double Beampipe_Shielding_8_thickness = 3.5*m;
      G4double Beampipe_Shielding_8_radius_1 = Beampipe_8_radius_2; 
      G4double Beampipe_Shielding_8_radius_2 = Beampipe_8_radius_2+Beampipe_Shielding_8_thickness;

      G4double Beampipe_Shielding_8_len = (Beampipe_8_pos_z + Beampipe_8_len/2.)-(lead_top_z/2.+lead_fb_z);
      G4double Beampipe_Shielding_8_pos_z = Beampipe_8_pos_z + Beampipe_8_len/2. - Beampipe_Shielding_8_len/2.;
  
  
  //////////////////////////////////////////////////////////////////////
  ///////////// Shielding Geometry Construction ////////////////////////
  //////////////////////////////////////////////////////////////////////

      // Shielding for Beampipe 2

            auto Beampipe_Shielding_2_S = new G4Cons("Beampipe_Shielding_2_S", 
                                          Beampipe_Shielding_2_radius_1,Beampipe_Shielding_2_radius_2,
                                          Beampipe_Shielding_2_radius_1,Beampipe_Shielding_2_radius_2,
                                          Beampipe_Shielding_2_len/2.,0.,360.*deg);

            auto Beampipe_Shielding_2_LV = new G4LogicalVolume(Beampipe_Shielding_2_S,Beampipe_Shielding_Material,"Beampipe_Shielding_2_LV");
            new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_2_pos_z),Beampipe_Shielding_2_LV,"Beampipe_Shielding_2_PV",mother,false,0,true);


      // Shielding for Beampipe 4
  
            auto Beampipe_Shielding_4_S = new G4Cons("Beampipe_Shielding_4_S", 
                                          Beampipe_Shielding_4_radius_1,Beampipe_Shielding_4_radius_2,
                                          Beampipe_Shielding_4_radius_1,Beampipe_Shielding_4_radius_2,
                                          Beampipe_Shielding_4_len/2.,0.,360.*deg);

            auto Beampipe_Shielding_4_LV = new G4LogicalVolume(Beampipe_Shielding_4_S,Beampipe_Shielding_Material,"Beampipe_Shielding_4_LV");
            new G4PVPlacement(0,G4ThreeVector(0., 0.,(Beampipe_2_pos_z+Beampipe_2_len/2.+Beampipe_Shielding_4_len/2.)),Beampipe_Shielding_4_LV,"Beampipe_Shielding_4_PV",mother,false,0,true);
 

      // Shielding for Beampipe 8

            auto Beampipe_Shielding_8_S = new G4Cons("Beampipe_Shielding_8_S", 
                                          Beampipe_Shielding_8_radius_1,Beampipe_Shielding_8_radius_2,
                                          Beampipe_Shielding_8_radius_1,Beampipe_Shielding_8_radius_2,
                                          Beampipe_Shielding_8_len/2.,0.,360.*deg);

            auto Beampipe_Shielding_8_LV = new G4LogicalVolume(Beampipe_Shielding_8_S,Beampipe_Shielding_Material,"Beampipe_8_LV");
            new G4PVPlacement(0,G4ThreeVector(0., 0.,Beampipe_Shielding_8_pos_z),Beampipe_Shielding_8_LV,"Beampipe_Shielding_8_PV",mother,false,0,true);
 
  ///////////////////////////////////////////////////////////////////
  ////////////// Position of all the beampipe parts /////////////////
  ///////////////////////////////////////////////////////////////////

  auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549));black_color->SetVisibility(true);
  auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true);
  auto green_color= new G4VisAttributes(G4Colour(0.788235,0.890196,0.741176)); green_color->SetVisibility(true);
  auto orange_color= new G4VisAttributes(G4Colour(0.988235,0.686275,0.243137)); orange_color->SetVisibility(true);
  auto purple_color= new G4VisAttributes(G4Colour(0.48,0.27,0.833)); purple_color->SetVisibility(true);
  
  Beampipe_Shielding_2_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(green_color);
  Beampipe_Shielding_4_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());//SetVisAttributes(green_color);
  Beampipe_Shielding_8_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());//SetVisAttributes(green_color);


  return Beampipe_Shielding_Construction_list;
}