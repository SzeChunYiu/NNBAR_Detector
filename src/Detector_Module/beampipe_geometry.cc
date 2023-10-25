// Beampipe
// Beampipe geometry up to 2023.01.16 -> Yuri's drawing on the beampipe with different sections
// Beampiep geometry labelling and dimensions please find the supplementary document ** (Missing)
// Beampipe sections are coated with neutron absorber

#include "DetectorConstruction.hh"
#include "Detector_Module/beampipe_geometry.hh"

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


#include <G4ProductionCuts.hh>
#include "G4RegionStore.hh"
#include "G4Region.hh"

Beampipe::Beampipe(){}

Beampipe::~Beampipe(){}

//// All dimensions for beampipe ////
// Also shared by other detector module classes
// Dimensions for the beampipe
G4double Beampipe_thickness = 2.0*cm; G4double Beampipe_angle = 360. * deg;
G4double Beampipe_coating_thickness = 1.0*cm;

//Beampipe -- 1
G4double Beampipe_1_radius_1 = 1.0*m;
G4double Beampipe_1_radius_2 = 1.8*m;
G4double Beampipe_1_radius_11 = Beampipe_1_radius_1-Beampipe_thickness;
G4double Beampipe_1_radius_12 = Beampipe_1_radius_1;
G4double Beampipe_1_radius_21 = Beampipe_1_radius_2-Beampipe_thickness;
G4double Beampipe_1_radius_22 = Beampipe_1_radius_2;
G4double Beampipe_1_len = 14.4*m;

//Beampipe -- 2
G4double Beampipe_2_radius = 1.8*m;
G4double Beampipe_2_radius_1 = Beampipe_2_radius-Beampipe_thickness;
G4double Beampipe_2_radius_2 = Beampipe_2_radius;
G4double Beampipe_2_len = 171.*m; // original 171 m can change to 161m 

//Beampipe -- 4 (short Beampipe before the detector)
G4double Beampipe_4_radius_1 = 1.04*m; 
G4double Beampipe_4_radius_2 = Beampipe_4_radius_1+Beampipe_thickness;
G4double Beampipe_4_len = 3.5*m; //original 2.5 m can increase to 12.5 m

//Beampipe -- 3 (cap)
G4double Beampipe_3_radius_1 = Beampipe_4_radius_2; 
G4double Beampipe_3_radius_2 = Beampipe_2_radius_2;
G4double Beampipe_3_len = Beampipe_thickness;

    //Beampipe 3 Coating
    G4double Beampipe_3_coating_radius_1 = Beampipe_4_radius_1;
    G4double Beampipe_3_coating_radius_2 = Beampipe_2_radius_1-Beampipe_coating_thickness;
    G4double Beampipe_3_coating_thickness = 2.0*cm;

//Beampipe -- 5 (part attached to the detector)
G4double Beampipe_5_radius_1 = 1.12*m; 
G4double Beampipe_5_radius_2 = Beampipe_5_radius_1+Beampipe_thickness;
G4double Beampipe_5_len = 5.0*m;

//Beampipe -- 6 (Cover of the Beampipe 5)
G4double Beampipe_6_radius_1 = Beampipe_4_radius_2;
G4double Beampipe_6_radius_2 = Beampipe_5_radius_2;
G4double Beampipe_6_len = Beampipe_thickness;
    
      //Beampipe 6 Coating
      G4double Beampipe_6_coating_radius_1 = Beampipe_4_radius_1;
      G4double Beampipe_6_coating_radius_2 = Beampipe_5_radius_1-Beampipe_coating_thickness;
      G4double Beampipe_6_coating_thickness = 2.0*cm;

//Beampipe -- 7 (Cover of the Beampipe 5)
G4double Beampipe_7_radius_1 = Beampipe_4_radius_2;
G4double Beampipe_7_radius_2 = Beampipe_5_radius_2;
G4double Beampipe_7_len = Beampipe_thickness;

      //Beampipe 7 Coating
      G4double Beampipe_7_coating_radius_1 = Beampipe_4_radius_1;
      G4double Beampipe_7_coating_radius_2 = Beampipe_5_radius_1-Beampipe_coating_thickness;
      G4double Beampipe_7_coating_thickness = 2.0*cm;

//part 8 of the Beampipe (last part of the Beampipe)
G4double Beampipe_8_radius_1 = 1.02*m; 
G4double Beampipe_8_radius_2 = Beampipe_7_radius_1;
G4double Beampipe_8_len = 16.5*m;

//Beam stop
G4double BeamStop_radius_1 = 0.;
G4double BeamStop_radius_2 = Beampipe_8_radius_1-Beampipe_coating_thickness;
G4double BeamStop_Absorber_thickness = 30.0*cm;
G4double BeamStop_Metal_thickness = 3.0*m;
G4double BeamStop_len =  BeamStop_Absorber_thickness + BeamStop_Metal_thickness;


// Position of all the beampipe parts //
// Position all reference to Beampipe 5
G4double Beampipe_5_pos_z = 0.0*cm;

G4double Beampipe_6_pos_z = Beampipe_5_pos_z-Beampipe_5_len/2.-Beampipe_6_len/2.;
G4double Beampipe_6_coating_pos_z = Beampipe_5_pos_z-Beampipe_5_len/2.+Beampipe_6_coating_thickness/2.;

G4double Beampipe_7_pos_z = Beampipe_5_pos_z+Beampipe_5_len/2.+Beampipe_7_len/2.;
G4double Beampipe_7_coating_pos_z = Beampipe_5_pos_z+Beampipe_5_len/2.-Beampipe_7_coating_thickness/2.;

G4double Beampipe_8_pos_z = Beampipe_5_pos_z+Beampipe_5_len/2.+Beampipe_8_len/2.;

G4double BeamStop_pos_z = Beampipe_8_pos_z + Beampipe_8_len/2. - BeamStop_len/2.;

G4double Beampipe_4_pos_z = Beampipe_5_pos_z-Beampipe_5_len/2.-Beampipe_4_len/2.;

G4double Beampipe_3_pos_z = Beampipe_4_pos_z-Beampipe_4_len/2.+Beampipe_3_len/2.;

G4double Beampipe_3_coating_pos_z = Beampipe_4_pos_z-Beampipe_4_len/2.-Beampipe_3_coating_thickness/2.;

G4double Beampipe_2_pos_z = Beampipe_4_pos_z-Beampipe_4_len/2.-Beampipe_2_len/2.;

G4double Beampipe_1_pos_z = Beampipe_2_pos_z - Beampipe_2_len/2. - Beampipe_1_len/2.;  



void Beampipe::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  G4double a; G4double z; G4double density;
  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elB = nistManager->FindOrBuildElement("B");
  G4Element* elSi = nistManager->FindOrBuildElement("Si");
  G4Element* elMn = nistManager->FindOrBuildElement("Mn");
  G4Element* elNi = nistManager->FindOrBuildElement("Ni");
  G4Element* elCr = nistManager->FindOrBuildElement("Cr");
  G4Element* elF = nistManager->FindOrBuildElement("F");
  G4Element* elFe = nistManager->FindOrBuildElement("Fe");
  G4Element* elCd = nistManager->FindOrBuildElement("Cd");
  G4Element* elLi = nistManager->FindOrBuildElement("Li");

  std::cout<< "Working on beampipe materials 1 " << std::endl;
  // Beam stop copper
  new G4Material("Copper", z=29., a=63.5*g/mole, density=8.9*g/cm3); 
  
  // beamline coating materals
  // = = Li6
  G4Isotope* Li6_isotope = new G4Isotope("Li6_isotope",3,6,a=6.015*g/mole); 
  G4Element* el_Li6 = new G4Element("el_Li","li6",1);
  el_Li6->AddIsotope(Li6_isotope,1.0);
  G4Material* Li6 = new G4Material("el_Li6", 0.460*g/cm3, 1);
  Li6->AddElement(el_Li6, 1);
  
  // = = 6LiF 
  G4Material* Li6F = new G4Material("Li6F", 2.635*g/cm3, 3);
  Li6F->AddElement(el_Li6, 0.475); Li6F->AddElement(elF, 0.5);
  Li6F->AddElement(elLi, 0.025); 

  // = = LiF 
  G4Material* LiF = new G4Material("LiF", 2.635*g/cm3, 2);
  LiF->AddElement(elLi, 0.5); LiF->AddElement(elF, 0.5);

  // = = B4C
  G4Material *B4C = new G4Material("B4C", density = 2.52*g/cm3, 2);
  B4C->AddElement(elB, 0.8); B4C->AddElement(elC, 0.2);

  // = = Cadmium
  G4Material *Cadmium = new G4Material("el_Cd", density = 8.65*g/cm3, 1);
  Cadmium->AddElement(elCd, 1);

  // Beampipe
  // = = Aluminum
  new G4Material("Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);
  // = = stainless steel
  G4Material* StainlessSteel = new G4Material("StainlessSteel",   density = 8.02*g/cm3,5);
  StainlessSteel->AddElement(elMn, 0.02);
  StainlessSteel->AddElement(elSi, 0.01);
  StainlessSteel->AddElement(elCr, 0.19);
  StainlessSteel->AddElement(elNi, 0.10);
  StainlessSteel->AddElement(elFe, 0.68);

}

std::vector<G4LogicalVolume*> Beampipe::Construct_Volumes(G4LogicalVolume* mother)
{
  // OUR BEAMPIPE GEOMETRY OUTPUT SENT BACK TO DETECTOR_CONSTRUCTION.CC //
  std::vector<G4LogicalVolume*> Beampipe_Construction_list;
  //####################################################################//

  DefineMaterials();
  auto DefaultMaterial = G4Material::GetMaterial("Galactic");
  auto BeampipeMaterial = G4Material::GetMaterial("Aluminum"); // Aluminum , el_Be
  auto CopperMaterial = G4Material::GetMaterial("Copper");
  auto Li6FMaterial = G4Material::GetMaterial("Li6F"); // el_Li6,B4C,Cd,Galactic,Li6F,LiF
  auto B4CMaterial = G4Material::GetMaterial("B4C");

  //// The Structure goes like this
  //// Beampipe_S
  //// ---- Beampipe wall S
  //// ---- Beampipe Coating

  //// The Cap Structure goes like this (for 3,6,7)
  //// Beampipe_wall S

  auto Beampipe_1_S = new G4Cons("Beampipe_1_S", 
                                  Beampipe_1_radius_11-Beampipe_coating_thickness,Beampipe_1_radius_12,
                                  Beampipe_1_radius_21-Beampipe_coating_thickness,Beampipe_1_radius_22,
                                  Beampipe_1_len/2.,0.,Beampipe_angle);

  auto Beampipe_1_LV = new G4LogicalVolume(Beampipe_1_S,DefaultMaterial,"Beampipe_1_LV");
  
        auto Beampipe_1_wall_S = new G4Cons("Beampipe_1_wall_S", 
                                              Beampipe_1_radius_11,Beampipe_1_radius_12,
                                              Beampipe_1_radius_21,Beampipe_1_radius_22,
                                              Beampipe_1_len/2.,0.,Beampipe_angle);

        auto Beampipe_1_coating_S = new G4Cons("Beampipe_1_coating_S", 
                                              Beampipe_1_radius_11-Beampipe_coating_thickness,Beampipe_1_radius_11,
                                              Beampipe_1_radius_21-Beampipe_coating_thickness,Beampipe_1_radius_21,
                                              Beampipe_1_len/2.,0.,Beampipe_angle);

        auto Beampipe_1_wall_LV = new G4LogicalVolume(Beampipe_1_wall_S,BeampipeMaterial,"Beampipe_1_wall_LV");
        auto Beampipe_1_coating_LV = new G4LogicalVolume(Beampipe_1_coating_S,B4CMaterial,"Beampipe_1_coating_LV");

        auto Beampipe_1_wall_PV = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),Beampipe_1_wall_LV,"Beampipe_1_wall_PV",Beampipe_1_LV,false,0,true);
        auto Beampipe_1_coating_PV = new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),Beampipe_1_coating_LV,"Beampipe_1_coating_PV",Beampipe_1_LV,false,0,true);


  auto Beampipe_2_S = new G4Cons("Beampipe_2_S", 
                                  Beampipe_2_radius_1-Beampipe_coating_thickness,Beampipe_2_radius_2,
                                  Beampipe_2_radius_1-Beampipe_coating_thickness,Beampipe_2_radius_2,
                                  Beampipe_2_len/2.,0.,Beampipe_angle);

  auto Beampipe_2_LV = new G4LogicalVolume(Beampipe_2_S,DefaultMaterial,"Beampipe_2_LV");
  
        auto Beampipe_2_wall_S = new G4Cons("Beampipe_2_wall_S", 
                                            Beampipe_2_radius_1,Beampipe_2_radius_2,
                                            Beampipe_2_radius_1,Beampipe_2_radius_2,
                                            Beampipe_2_len/2.,0.,Beampipe_angle);

        auto Beampipe_2_coating_S = new G4Cons("Beampipe_2_coating_S", 
                                                Beampipe_2_radius_1-Beampipe_coating_thickness,Beampipe_2_radius_1,
                                                Beampipe_2_radius_1-Beampipe_coating_thickness,Beampipe_2_radius_1,
                                                Beampipe_2_len/2.,0.,Beampipe_angle);

        auto Beampipe_2_wall_LV = new G4LogicalVolume(Beampipe_2_S,BeampipeMaterial,"Beampipe_2_wall_LV");
        auto Beampipe_2_coating_LV = new G4LogicalVolume(Beampipe_2_S,B4CMaterial,"Beampipe_2_coating_LV");

        new G4PVPlacement(0,G4ThreeVector(),Beampipe_2_wall_LV,"Beampipe_2_wall_PV",Beampipe_2_LV,false,0,true);
        new G4PVPlacement(0,G4ThreeVector(),Beampipe_2_coating_LV,"Beampipe_2_coating_PV",Beampipe_2_LV,false,0,true);


  auto Beampipe_3_S = new G4Cons("Beampipe_3_S", 
                                  Beampipe_3_radius_1,Beampipe_3_radius_2,
                                  Beampipe_3_radius_1,Beampipe_3_radius_2,
                                  Beampipe_3_len/2.,0.,Beampipe_angle);
  
  auto Beampipe_3_LV = new G4LogicalVolume(Beampipe_3_S,BeampipeMaterial,"Beampipe_3_wall_LV");

        auto Beampipe_3_coating_S = new G4Cons("Beampipe_3_coating_S", 
                                                Beampipe_3_coating_radius_1,Beampipe_3_coating_radius_2,
                                                Beampipe_3_coating_radius_1,Beampipe_3_coating_radius_2,
                                                Beampipe_3_coating_thickness/2.,0.,Beampipe_angle);

        auto Beampipe_3_coating_LV = new G4LogicalVolume(Beampipe_3_coating_S,B4CMaterial,"Beampipe_3_coating_LV");

  
  auto Beampipe_4_S = new G4Cons("Beampipe_4_S", 
                                  Beampipe_4_radius_1-Beampipe_coating_thickness,Beampipe_4_radius_2,
                                  Beampipe_4_radius_1-Beampipe_coating_thickness,Beampipe_4_radius_2,
                                  Beampipe_4_len/2.,0.,Beampipe_angle);

  auto Beampipe_4_LV = new G4LogicalVolume(Beampipe_4_S,DefaultMaterial,"Beampipe_4_LV");
  
        auto Beampipe_4_wall_S = new G4Cons("Beampipe_4_wall_S", 
                                            Beampipe_4_radius_1,Beampipe_4_radius_2,
                                            Beampipe_4_radius_1,Beampipe_4_radius_2,
                                            Beampipe_4_len/2.,0.,Beampipe_angle);

        auto Beampipe_4_coating_S = new G4Cons("Beampipe_4_coating_S", 
                                                Beampipe_4_radius_1-Beampipe_coating_thickness,Beampipe_4_radius_1,
                                                Beampipe_4_radius_1-Beampipe_coating_thickness,Beampipe_4_radius_1,
                                                Beampipe_4_len/2.,0.,Beampipe_angle);

        auto Beampipe_4_wall_LV = new G4LogicalVolume(Beampipe_4_S,BeampipeMaterial,"Beampipe_4_wall_LV");
        auto Beampipe_4_coating_LV = new G4LogicalVolume(Beampipe_4_S,B4CMaterial,"Beampipe_4_coating_LV");

        new G4PVPlacement(0,G4ThreeVector(),Beampipe_4_wall_LV,"Beampipe_4_wall_PV",Beampipe_4_LV,false,0,true);
        new G4PVPlacement(0,G4ThreeVector(),Beampipe_4_coating_LV,"Beampipe_4_coating_PV",Beampipe_4_LV,false,0,true);

  auto Beampipe_5_S = new G4Cons("Beampipe_5_S", 
                                  Beampipe_5_radius_1-Beampipe_coating_thickness,Beampipe_5_radius_2,
                                  Beampipe_5_radius_1-Beampipe_coating_thickness,Beampipe_5_radius_2,
                                  Beampipe_5_len/2.,0.,Beampipe_angle);

  auto Beampipe_5_LV = new G4LogicalVolume(Beampipe_5_S,DefaultMaterial,"Beampipe_5_LV");
  
        auto Beampipe_5_wall_S = new G4Cons("Beampipe_5_wall_S", 
                                            Beampipe_5_radius_1,Beampipe_5_radius_2,
                                            Beampipe_5_radius_1,Beampipe_5_radius_2,
                                            Beampipe_5_len/2.,0.,Beampipe_angle);

        auto Beampipe_5_coating_S = new G4Cons("Beampipe_5_coating_S", 
                                                Beampipe_5_radius_1-Beampipe_coating_thickness,Beampipe_5_radius_1,
                                                Beampipe_5_radius_1-Beampipe_coating_thickness,Beampipe_5_radius_1,
                                                Beampipe_5_len/2.,0.,Beampipe_angle);

        auto Beampipe_5_wall_LV = new G4LogicalVolume(Beampipe_5_wall_S,BeampipeMaterial,"Beampipe_5_wall_LV");
        auto Beampipe_5_coating_LV = new G4LogicalVolume(Beampipe_5_coating_S,B4CMaterial,"Beampipe_5_coating_LV");

        new G4PVPlacement(0,G4ThreeVector(),Beampipe_5_wall_LV,"Beampipe_5_wall_PV",Beampipe_5_LV,false,0,true);
        new G4PVPlacement(0,G4ThreeVector(),Beampipe_5_coating_LV,"Beampipe_5_coating_PV",Beampipe_5_LV,false,0,true);

  
  auto Beampipe_6_S = new G4Cons("Beampipe_6_S", 
                                  Beampipe_6_radius_1,Beampipe_6_radius_2,
                                  Beampipe_6_radius_1,Beampipe_6_radius_2,
                                  Beampipe_6_len/2.,0.,Beampipe_angle);
  
  auto Beampipe_6_LV = new G4LogicalVolume(Beampipe_6_S,BeampipeMaterial,"Beampipe_6_wall_LV");

        auto Beampipe_6_coating_S = new G4Cons("Beampipe_6_coating_S", 
                                                Beampipe_6_coating_radius_1,Beampipe_6_coating_radius_2,
                                                Beampipe_6_coating_radius_1,Beampipe_6_coating_radius_2,
                                                Beampipe_6_coating_thickness/2.,0.,Beampipe_angle);
        auto Beampipe_6_coating_LV = new G4LogicalVolume(Beampipe_6_coating_S,B4CMaterial,"Beampipe_6_coating_LV");

    
  auto Beampipe_7_S = new G4Cons("Beampipe_7_S", 
                                  Beampipe_7_radius_1,Beampipe_7_radius_2,
                                  Beampipe_7_radius_1,Beampipe_7_radius_2,
                                  Beampipe_7_len/2.,0.,Beampipe_angle);
  
  auto Beampipe_7_LV = new G4LogicalVolume(Beampipe_7_S,BeampipeMaterial,"Beampipe_7_wall_LV");

        auto Beampipe_7_coating_S = new G4Cons("Beampipe_7_coating_S", 
                                                Beampipe_7_coating_radius_1,Beampipe_7_coating_radius_2,
                                                Beampipe_7_coating_radius_1,Beampipe_7_coating_radius_2,
                                                Beampipe_7_coating_thickness/2.,0.,Beampipe_angle);
        auto Beampipe_7_coating_LV = new G4LogicalVolume(Beampipe_7_coating_S,B4CMaterial,"Beampipe_7_coating_LV");


  auto Beampipe_8_S = new G4Cons("Beampipe_8_S", 
                                  Beampipe_8_radius_1-Beampipe_coating_thickness,Beampipe_8_radius_2,
                                  Beampipe_8_radius_1-Beampipe_coating_thickness,Beampipe_8_radius_2,
                                  Beampipe_8_len/2.,0.,Beampipe_angle);

  auto Beampipe_8_LV = new G4LogicalVolume(Beampipe_8_S,DefaultMaterial,"Beampipe_8_LV");
  
        auto Beampipe_8_wall_S = new G4Cons("Beampipe_8_wall_S", 
                                            Beampipe_8_radius_1,Beampipe_8_radius_2,
                                            Beampipe_8_radius_1,Beampipe_8_radius_2,
                                            Beampipe_8_len/2.,0.,Beampipe_angle);

        auto Beampipe_8_coating_S = new G4Cons("Beampipe_8_coating_S", 
                                                Beampipe_8_radius_1-Beampipe_coating_thickness,Beampipe_8_radius_1,
                                                Beampipe_8_radius_1-Beampipe_coating_thickness,Beampipe_8_radius_1,
                                                Beampipe_8_len/2.,0.,Beampipe_angle);

        auto Beampipe_8_wall_LV = new G4LogicalVolume(Beampipe_8_S,BeampipeMaterial,"Beampipe_8_wall_LV");
        auto Beampipe_8_coating_LV = new G4LogicalVolume(Beampipe_8_S,B4CMaterial,"Beampipe_8_coating_LV");

        new G4PVPlacement(0,G4ThreeVector(),Beampipe_8_wall_LV,"Beampipe_8_wall_PV",Beampipe_8_LV,false,0,true);
        new G4PVPlacement(0,G4ThreeVector(),Beampipe_8_coating_LV,"Beampipe_8_coating_PV",Beampipe_8_LV,false,0,true);


  auto BeamStop_S = new G4Cons("BeamStop_S", 
                                  BeamStop_radius_1,BeamStop_radius_2,
                                  BeamStop_radius_1,BeamStop_radius_2,
                                  BeamStop_len/2.,0.,Beampipe_angle);
  
  auto BeamStop_LV = new G4LogicalVolume(BeamStop_S,DefaultMaterial,"BeamStop_LV");

        auto BeamStop_Metal_S = new G4Cons("BeamStop_Metal_S", 
                                                BeamStop_radius_1,BeamStop_radius_2,
                                                BeamStop_radius_1,BeamStop_radius_2,
                                                BeamStop_Metal_thickness/2.,0.,Beampipe_angle);
        auto BeamStop_Metal_LV = new G4LogicalVolume(BeamStop_Metal_S,CopperMaterial,"BeamStop_Metal_LV");


        auto BeamStop_Absorber_S = new G4Cons("BeamStop_Absorber_S", 
                                          BeamStop_radius_1,BeamStop_radius_2,
                                          BeamStop_radius_1,BeamStop_radius_2,
                                          BeamStop_Absorber_thickness/2.,0.,Beampipe_angle);
        auto BeamStop_Absorber_LV = new G4LogicalVolume(BeamStop_Absorber_S,B4CMaterial,"BeamStop_Absorber_LV"); //CopperMaterial

  auto BeamStop_Absorber_PV = new G4PVPlacement(0,G4ThreeVector(0., 0.,-BeamStop_len/2.+BeamStop_Absorber_thickness/2.),BeamStop_Absorber_LV,"BeamStop_Absorber_PV",BeamStop_LV,false,0,true);
  auto BeamStop_Metal_PV = new G4PVPlacement(0,G4ThreeVector(0., 0.,-BeamStop_len/2.+BeamStop_Absorber_thickness+BeamStop_Metal_thickness/2.),BeamStop_Metal_LV,"BeamStop_Metal_PV",BeamStop_LV,false,0,true);

  auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549));black_color->SetVisibility(true);
  auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true);
  auto green_color= new G4VisAttributes(G4Colour(0.788235,0.890196,0.741176)); green_color->SetVisibility(true);
  auto orange_color= new G4VisAttributes(G4Colour(0.988235,0.686275,0.243137)); orange_color->SetVisibility(true);
  auto purple_color= new G4VisAttributes(G4Colour(0.48,0.27,0.833)); purple_color->SetVisibility(true);
  
  Beampipe_1_LV ->SetVisAttributes (G4VisAttributes::GetInvisible());// -> SetVisAttributes(black_color);
  Beampipe_1_coating_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(orange_color);
  Beampipe_1_wall_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(grey_color);

  Beampipe_2_LV ->SetVisAttributes (G4VisAttributes::GetInvisible());
  Beampipe_2_coating_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(orange_color);
  Beampipe_2_wall_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(grey_color);

  Beampipe_3_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());//SetVisAttributes (grey_color);
  Beampipe_3_coating_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());//SetVisAttributes (green_color);

  Beampipe_4_LV ->SetVisAttributes (G4VisAttributes::GetInvisible());
  Beampipe_4_coating_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(orange_color);
  Beampipe_4_wall_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(grey_color);

  Beampipe_5_LV ->SetVisAttributes (G4VisAttributes::GetInvisible());
  Beampipe_5_coating_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(orange_color);
  Beampipe_5_wall_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(grey_color);
  
  Beampipe_6_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(grey_color);
  Beampipe_6_coating_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(green_color);

  Beampipe_7_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(grey_color);
  Beampipe_7_coating_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(green_color);

  Beampipe_8_LV ->SetVisAttributes (G4VisAttributes::GetInvisible());
  Beampipe_8_coating_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(orange_color);
  Beampipe_8_wall_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(grey_color);

  BeamStop_LV ->SetVisAttributes (G4VisAttributes::GetInvisible());
  BeamStop_Absorber_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(purple_color);
  BeamStop_Metal_LV -> SetVisAttributes (G4VisAttributes::GetInvisible());// SetVisAttributes(purple_color);

  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_1_pos_z),Beampipe_1_LV,"Beampipe_1_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_2_pos_z),Beampipe_2_LV,"Beampipe_2_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_3_pos_z),Beampipe_3_LV,"Beampipe_3_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_3_coating_pos_z),Beampipe_3_coating_LV,"Beampipe_3_coating_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_4_pos_z),Beampipe_4_LV,"Beampipe_4_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_5_pos_z),Beampipe_5_LV,"Beampipe_5_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_6_pos_z),Beampipe_6_LV,"Beampipe_6_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_6_coating_pos_z),Beampipe_6_coating_LV,"Beampipe_6_coating_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_7_pos_z),Beampipe_7_LV,"Beampipe_7_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_7_coating_pos_z),Beampipe_7_coating_LV,"Beampipe_7_coating_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., Beampipe_8_pos_z),Beampipe_8_LV,"Beampipe_8_PV",mother,false,0,true);
  new G4PVPlacement(0,G4ThreeVector(0., 0., BeamStop_pos_z),BeamStop_LV,"BeamStop_PV",mother,false,0,true);
  
  //////////////////////////////////////////////////////////////////////
  //////////////////// Region Settings /////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  
    // Region for beampipe
    G4Region* Beampipe_region = new G4Region("Beampipe_region");
      Beampipe_region->AddRootLogicalVolume(Beampipe_1_wall_LV);
       Beampipe_region->AddRootLogicalVolume(Beampipe_1_LV);
      Beampipe_region->AddRootLogicalVolume(Beampipe_2_wall_LV);
       Beampipe_region->AddRootLogicalVolume(Beampipe_2_LV);
      Beampipe_region->AddRootLogicalVolume(Beampipe_3_LV);
      Beampipe_region->AddRootLogicalVolume(Beampipe_4_wall_LV);
       Beampipe_region->AddRootLogicalVolume(Beampipe_4_LV);
       Beampipe_region->AddRootLogicalVolume(Beampipe_5_LV);
      Beampipe_region->AddRootLogicalVolume(Beampipe_5_wall_LV);
      Beampipe_region->AddRootLogicalVolume(Beampipe_6_LV);
      Beampipe_region->AddRootLogicalVolume(Beampipe_7_LV);
      Beampipe_region->AddRootLogicalVolume(Beampipe_8_wall_LV);
       Beampipe_region->AddRootLogicalVolume(Beampipe_8_LV);

    G4Region* Beampipe_Region = G4RegionStore::GetInstance()->GetRegion("Beampipe_region");
    G4ProductionCuts* Beampipecut = new G4ProductionCuts();
    Beampipecut->SetProductionCut(1.0*cm,"gamma"); // 130cm -> 1.3 MeV   220cm -> 3.4 MeV 400*cm -> 50 MeV
    Beampipecut->SetProductionCut(1.0*mm,"e-");
    Beampipecut->SetProductionCut(1.0*mm,"e+");
    Beampipecut->SetProductionCut(1.0*mm,"proton");
    Beampipe_Region->SetProductionCuts(Beampipecut);

  Beampipe_Construction_list.push_back(Beampipe_1_LV);
  Beampipe_Construction_list.push_back(Beampipe_2_LV);
  Beampipe_Construction_list.push_back(Beampipe_3_LV);
  Beampipe_Construction_list.push_back(Beampipe_4_LV);
  Beampipe_Construction_list.push_back(Beampipe_5_LV);
  Beampipe_Construction_list.push_back(Beampipe_6_LV);
  Beampipe_Construction_list.push_back(Beampipe_7_LV);
  Beampipe_Construction_list.push_back(Beampipe_8_LV);
  Beampipe_Construction_list.push_back(Beampipe_3_coating_LV);
  Beampipe_Construction_list.push_back(Beampipe_7_coating_LV);
  Beampipe_Construction_list.push_back(Beampipe_6_coating_LV);
  Beampipe_Construction_list.push_back(BeamStop_LV);
  Beampipe_Construction_list.push_back(Beampipe_1_wall_LV);
  Beampipe_Construction_list.push_back(Beampipe_2_wall_LV);
  Beampipe_Construction_list.push_back(Beampipe_4_wall_LV);
  Beampipe_Construction_list.push_back(Beampipe_5_wall_LV);
  Beampipe_Construction_list.push_back(Beampipe_8_wall_LV);
  Beampipe_Construction_list.push_back(Beampipe_1_coating_LV);
  Beampipe_Construction_list.push_back(Beampipe_2_coating_LV);
  Beampipe_Construction_list.push_back(Beampipe_4_coating_LV);
  Beampipe_Construction_list.push_back(Beampipe_5_coating_LV);
  Beampipe_Construction_list.push_back(Beampipe_8_coating_LV);
  return Beampipe_Construction_list;
}