// Silicon
#include "DetectorConstruction.hh"
#include "Detector_Module/Silicon_geometry.hh"
#include "Detector_Module/Silicon_geometry.hh"

#include "G4SDManager.hh"
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


Silicon::Silicon():G4VUserDetectorConstruction(){} //:G4VUserDetectorConstruction()
Silicon::~Silicon(){}


// Beampipe dimensions
extern G4double Beampipe_5_len;
extern G4double Beampipe_5_radius_2;

void Silicon::DefineMaterials()
{
  auto nistManager = G4NistManager::Instance();
  G4double a; G4double z; G4double density;

  G4Element* elC = nistManager->FindOrBuildElement("C");
  G4Element* elB = nistManager->FindOrBuildElement("B");
  G4Element* elLi = nistManager->FindOrBuildElement("Li");
  G4Element* elF = nistManager->FindOrBuildElement("F");
  
  // Silicon
  G4Material* Silicon = new G4Material("Silicon", z=14., a= 28.0855*g/mole, density = 2.33*g/cm3);
  
  // = = B4C
  G4Material *B4C = new G4Material("B4C", density = 2.52*g/cm3, 2);
  B4C->AddElement(elB, 0.8); B4C->AddElement(elC, 0.2);

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

  
}

std::vector<G4LogicalVolume*> Silicon::Construct_Volumes(G4LogicalVolume* mother)
{
  // OUR Silicon GEOMETRY OUTPUT SENT BACK TO DETECTOR_CONSTRUCTION.CC //
  std::vector<G4LogicalVolume*> Silicon_Construction_list;
  //####################################################################//

  ////////////////////////////////////////////////////
  ///////////// Defining materials ///////////////////
  ////////////////////////////////////////////////////
    DefineMaterials();
    auto defaultMaterial = G4Material::GetMaterial("Galactic");
    auto SiliconMaterial = G4Material::GetMaterial("Silicon");
    auto beam_coating_Material = G4Material::GetMaterial("Li6F"); //Li6F

  /////////////////////////////////////////////////////
  //////////// Numbers for Silicon construction ///////////
  /////////////////////////////////////////////////////

  G4double silicon_radius_1 = 103.0*cm;//87.97*cm; 
  G4double silicon_radius_2 = 107.0*cm;//97.97*cm;
  G4double silicon_len = Beampipe_5_len - 10.*cm;
  G4double silicon_thickness = 0.2*cm; G4double silicon_angle = 360. * deg;
  G4double Silicon_coating_thickness = 0.1*cm;

    ///////////////////////////////////////
    //////////   Silicon coating  /////////
    ///////////////////////////////////////

      auto siliconS_1 = new G4Cons("siliconS_1", silicon_radius_1,silicon_radius_1+silicon_thickness,silicon_radius_1,silicon_radius_1+silicon_thickness,silicon_len/2.,0.,silicon_angle);
      auto siliconLV_1 = new G4LogicalVolume(siliconS_1,SiliconMaterial,"siliconLV_1");
      auto siliconPV_1 = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),siliconLV_1,"siliconPV_1",mother,false,0,true);
      
      auto siliconS_2 = new G4Cons("siliconS_2", silicon_radius_2,silicon_radius_2+silicon_thickness,silicon_radius_2,silicon_radius_2+silicon_thickness,silicon_len/2.,0.,silicon_angle);
      auto siliconLV_2 = new G4LogicalVolume(siliconS_2,SiliconMaterial,"siliconLV_2");
      auto siliconPV_2 = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),siliconLV_2,"siliconPV_2",mother,false,1,true);
  
    ///////////////////////////////////////
    //////////   Silicon coating  /////////
    ///////////////////////////////////////

      auto silicon_coating_1_S = new G4Cons("Silicon_Coating_1_S", 
                                            silicon_radius_1-Silicon_coating_thickness, silicon_radius_1,
                                            silicon_radius_1-Silicon_coating_thickness,silicon_radius_1,
                                            silicon_len/2.,0.,silicon_angle);    
      
      auto silicon_coating_1_LV = new G4LogicalVolume(silicon_coating_1_S,beam_coating_Material,"Silicon_Coating_LV");

      
      auto silicon_coating_2_S = new G4Cons("Silicon_Coating_2_S", silicon_radius_2-Silicon_coating_thickness, silicon_radius_2,
                                            silicon_radius_2-Silicon_coating_thickness,silicon_radius_2,
                                            silicon_len/2.,0.,silicon_angle);    
      
      auto silicon_coating_2_LV = new G4LogicalVolume(silicon_coating_2_S,beam_coating_Material,"Silicon_Coating_2_LV");

      auto silicon_coating_1_PV = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),silicon_coating_1_LV,"silicon_coating_PV",mother,false,0,true);
      auto silicon_coating_2_PV = new G4PVPlacement(0,G4ThreeVector(0., 0., 0.),silicon_coating_2_LV,"silicon_coating_2_PV",mother,false,1,true);
   

  //////////////////////////////////////////////////////////////////////
  //////////////////// Writing the outputs /////////////////////////////
  //////////////////////////////////////////////////////////////////////

    Silicon_Construction_list.push_back(siliconLV_1);
    Silicon_Construction_list.push_back(siliconLV_2);
  
  //////////////////////////////////////////////////////////////////////
  //////////////////// Region Settings /////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  
    G4Region* Silicon_region = new G4Region("Silicon_region");
      Silicon_region->AddRootLogicalVolume(siliconLV_1);
      Silicon_region->AddRootLogicalVolume(siliconLV_2);

    G4Region* Silicon_Region = G4RegionStore::GetInstance()->GetRegion("Silicon_region");
    G4ProductionCuts* Siliconcut = new G4ProductionCuts();
    Siliconcut->SetProductionCut(0.01*mm,"gamma"); // 130cm -> 1.3 MeV   220cm -> 3.4 MeV 400*cm -> 50 MeV
    Siliconcut->SetProductionCut(0.01*mm,"e-");
    Siliconcut->SetProductionCut(0.01*mm,"e+");
    Siliconcut->SetProductionCut(0.01*mm,"proton");
    Silicon_Region->SetProductionCuts(Siliconcut);



  ////////////////////////////////////////////////////
  ///////////////// Color settings ///////////////////
  ////////////////////////////////////////////////////

    auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549));black_color->SetVisibility(true);
    auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true);
    auto green_color= new G4VisAttributes(G4Colour(0.788235,0.890196,0.741176)); green_color->SetVisibility(true);
    auto orange_color= new G4VisAttributes(G4Colour(0.988235,0.686275,0.243137)); orange_color->SetVisibility(true);
    auto purple_color= new G4VisAttributes(G4Colour(0.48,0.27,0.833)); purple_color->SetVisibility(true);
    
    siliconLV_1 -> SetVisAttributes(orange_color);
    siliconLV_2 -> SetVisAttributes(orange_color);
    // silicon_coating_1_LV->SetVisAttributes (G4VisAttributes::GetInvisible());
    // silicon_coating_2_LV ->SetVisAttributes (G4VisAttributes::GetInvisible());
  
  return Silicon_Construction_list;
}