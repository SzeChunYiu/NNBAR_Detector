// Scintillator
#include "DetectorConstruction.hh"
#include "Detector_Module/TPC_geometry.hh"
#include "Detector_Module/beampipe_geometry.hh"
#include "Detector_Module/Scintillator_geometry.hh"
#include "G4SDManager.hh"

#include <G4ProductionCuts.hh>
#include "G4RegionStore.hh"
#include "G4Region.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"

#include "G4Box.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"

#include "G4AutoDelete.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

std::ofstream scintillator_module_pos_outFile;

Scintillator::Scintillator():G4VUserDetectorConstruction(){}
Scintillator::~Scintillator(){}

extern G4double Beampipe_5_len;
extern G4double Beampipe_6_len;
extern G4double Beampipe_5_radius_2; // need to calculate how many scintilator modules we put along z direction! 
extern G4double TPC_drift_len;
extern G4double TPC_wall_thickness;
extern G4double TPC_z; 

void Scintillator::DefineMaterials()
{
  auto nistManager = G4NistManager::Instance();
  G4double a; G4double z; G4double density;

  G4Element* elH = nistManager->FindOrBuildElement("H");
  G4Element* elC = nistManager->FindOrBuildElement("C");
  
  // BC-408 taken from datasheet
  G4Material* Scint = new G4Material("Scint", 1.023*g/cm3, 2);
  Scint->AddElement(elH, 0.524573); Scint->AddElement(elC, 1.-0.524573);

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
  scintMPT->AddConstProperty("SCINTILLATIONYIELD", 0./ MeV); //original 11136000. 17400
  scintMPT->AddConstProperty("RESOLUTIONSCALE", 1.0);
  scintMPT->AddConstProperty("FASTTIMECONSTANT", 1.0*ns); // org: 0.9
  scintMPT->AddConstProperty("SLOWTIMECONSTANT", 1.0*ns); // org: 2.1
  scintMPT->AddConstProperty("YIELDRATIO", 1.);
  Scint->SetMaterialPropertiesTable(scintMPT);
  scintMPT->DumpTable();

}

std::vector<G4LogicalVolume*> Scintillator::Construct_Volumes(G4LogicalVolume* mother)
{
  /////////////////////////////////////////////////////////////
  ///////////// The output is defined here: ///////////////////
  /////////////////////////////////////////////////////////////
        scintillator_module_pos_outFile.open("./output/Scintillator_Module_Position.txt");
        scintillator_module_pos_outFile << "index,x,y,z" << G4endl;
        std::vector<G4LogicalVolume*> Scintillator_output;

  /////////////////////////////////////////////////////////////
  /////////////////// Define Materials ////////////////////////
  /////////////////////////////////////////////////////////////

    DefineMaterials();
    auto defaultMaterial = G4Material::GetMaterial("Galactic");
    auto scintMaterial = G4Material::GetMaterial("Scint");

  ////////////////////////////////////////////////////////////////
  //////////// Numbers needed for scintillator construction //////
  ////////////////////////////////////////////////////////////////

    G4double scintillator_coating_thickness = 1.*mm;
    G4double scint_bar_x = 10.*cm;
    G4double scint_bar_y = 3.*cm; 
    G4double scint_bar_z = 40.*cm;
    G4double scint_layers = 10.; int no_of_Layers=10;
    G4double WLS_radius = 1.0*cm;
    G4double WLS_length = 40.0*cm;
      

    double n_bar_x = 4.; // how many bars a long x direction in a layer
    double n_bar_z = 1.; 

    G4double scint_module_x = n_bar_x*scint_bar_x; // 4 scint bar along x will make length 20 => a square! 
    G4double scint_module_y = scint_layers * scint_bar_y; 
    G4double scint_module_z = n_bar_z*scint_bar_z;
    
    double n_scint_module_x = 10.0; 
    double n_scint_module_z = 11.0; //number along z direction

    G4double dy= 20.*cm;

    G4double scint_base_dist = Beampipe_5_radius_2 + TPC_drift_len + 2.*TPC_wall_thickness;
    std::cout << " Scintillator distance from origin to scintillator top surface " << (scint_base_dist + dy + scint_module_y)/cm << std::endl;
    
    G4double total_width_x = 2.*Beampipe_5_radius_2+(2.*TPC_drift_len+4.*TPC_wall_thickness) + dy + scint_module_y;//width we have for placing scintillator along x
    G4double width_x = 2.*Beampipe_5_radius_2+(2.*TPC_drift_len+4.*TPC_wall_thickness); // width above the TPC
    std::cout<< "TPC_total_width " << width_x / cm << "drift len: " << TPC_drift_len << std::endl;
    G4double total_width_z = Beampipe_5_len + 2.*Beampipe_6_len;

    G4double dz = (total_width_z-n_scint_module_z*scint_module_z)/(n_scint_module_z-1.0);
    G4double dx = (total_width_x-n_scint_module_x*scint_module_x)/(n_scint_module_x-1.0);

    std::cout << "scintillator gap along z direction" << dz/cm << std::endl;
    std::cout << "for scint along x" <<  dx/cm << std::endl;

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////  Scintillator modules for Top, bottom, left right surfaces /////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////

  //The structrue goes like this 
  //// Scintillator module
  //// --- Scintillator layer
  //// ------ Scintillator bars
      int scint_index_count = 0;

      std::vector<G4VPhysicalVolume *> scint_PV_array;
      std::vector<G4VPhysicalVolume *> scint_Module_PV_array;
      std::vector<G4VPhysicalVolume *> scint_detector_PV_array; 

      auto scint_moduleS = new G4Box("Scint_moduleS",scint_module_x/2., scint_module_y/2., scint_module_z/2.);
      auto scint_moduleLV = new G4LogicalVolume(scint_moduleS,defaultMaterial,"Scint_moduleLV");

      // scintillator layers 
      // layer 1 and 2: (1 = horizontal, 2 = vertical bars)
      auto scint_layerH_S = new G4Box("ScintS",scint_module_x/2., scint_bar_y/2., scint_module_z/2.);
      auto scint_layerH_LV = new G4LogicalVolume(scint_layerH_S, defaultMaterial,"Scint_layerH_LV");
      auto scint_layerV_S = new G4Box("Scint2S",scint_module_x/2., scint_bar_y/2., scint_module_z/2.);
      auto scint_layerV_LV = new G4LogicalVolume(scint_layerV_S, defaultMaterial,"Scint_layerV_LV");

      // Two types of bars here, horizontal and vertical
      // --- horizontal bars (along x)
      auto scint_barH_x = scint_bar_x;
      auto scint_barH_z = scint_bar_z;
      auto scint_barH_S = new G4Box("Scint_barS",scint_barH_x/2., scint_bar_y/2., scint_barH_z/2.);
      auto scint_barH_LV = new G4LogicalVolume(scint_barH_S,scintMaterial,"Scint_barH_LV");
      
      G4cout<< "scint H bar: " << scint_barH_x << " z: " << scint_barH_z << G4endl;

      //--- vertical bars (divided along z)
      auto scint_barV_x = scint_bar_z;
      auto scint_barV_z = scint_bar_x;
      auto scint_barV_S = new G4Box("Scint_barS",scint_barV_x/2., scint_bar_y/2.,scint_barV_z/2.);
      auto scint_barV_LV = new G4LogicalVolume(scint_barV_S,scintMaterial,"Scint_barV_LV");

      G4cout<< "scint V bar: " << scint_barV_x << " z: " << scint_barV_z << G4endl;
      
      // // Variables to define dimensions
      // G4double start_angle = 0.0*deg;
      // G4double spanning_angle = 360.0*deg;

      // // Create the cylindrical shape (a "tub")
      // G4Tubs* WLS_S = new G4Tubs("WLS",0,WLS_radius,WLS_length/2.,start_angle,spanning_angle);
      // G4LogicalVolume* WLS_LV = new G4LogicalVolume(WLS_S,defaultMaterial,"WLS_LV");
      // new G4PVPlacement(0, G4ThreeVector(2.*cm,0.,0.),WLS_LV,"WLS_PV",scint_barH_LV,false,i,true)
      // new G4PVPlacement(0, G4ThreeVector(-2.*cm,0.,0.),WLS_LV,"WLS_PV",scint_barH_LV,false,i,true)

      // G4RotationMatrix* rotationMatrix = new G4RotationMatrix(); rotationMatrix->rotateY(90.*deg);
      // new G4PVPlacement(rotationMatrix, position, WLS_LV, "WLS_PV", scint_barV_LV, false, 0, true);
      // G4ThreeVector(0, 0, 2.*cm);

      for (int i=0; i<no_of_Layers; i++){
        if (i%2==0){new G4PVPlacement(0, G4ThreeVector(0.,-scint_module_y/2.+(2*i+1)/2.*scint_bar_y,0.),scint_layerH_LV,"Scint_layerPV",scint_moduleLV,false,i,true);}
        else{new G4PVPlacement(0, G4ThreeVector(0.,-scint_module_y/2.+(2*i+1)/2.*scint_bar_y,0.),scint_layerV_LV,"Scint_layerPV",scint_moduleLV,false,i,true);}
      }

      // placing the bars
      for (int i=0; i<n_bar_x; i++){
        new G4PVPlacement(0,G4ThreeVector(-scint_module_x/2.+(2*i+1)/2.*scint_barH_x,0.,0.),scint_barH_LV,"Scint_barPV_H",scint_layerH_LV,false,i,true);
        new G4PVPlacement(0,G4ThreeVector(0.,0.,-scint_module_z/2.+(2*i+1)/2.*scint_barV_z),scint_barV_LV,"Scint_barPV_V",scint_layerV_LV,false,i,true);
      }

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////// Position calculation part of the scintillator modules ////////////////////////
      /////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
      std::vector<G4ThreeVector> scint_module_pos_vector_horizontal; // for top and bottom
      std::vector<G4ThreeVector> scint_module_pos_vector_vertical; // for sides

      for (int i=0; i<int(n_scint_module_x); i++){
        for (int j=0; j<int(n_scint_module_z); j++){

          G4double x = -(n_scint_module_x*scint_module_x+(n_scint_module_x-1)*dx)/2. + 0.5*scint_module_x + double(i)*(scint_module_x+dx);
          G4double y = 0.*cm;
          G4double z = -(n_scint_module_z*scint_module_z+(n_scint_module_z-1)*dz)/2. + 0.5*scint_module_z + double(j)*(scint_module_z+dz);
          scint_module_pos_vector_horizontal.push_back(G4ThreeVector(x,y,z));
          scint_module_pos_vector_vertical.push_back(G4ThreeVector(y,x,z));
        }
      }

      //Rotation matrix for scintillators to be put of top, bottom, left and right sides
      
      std::vector<G4RotationMatrix *> scint_rot_array;
      G4RotationMatrix* zRot1 = new G4RotationMatrix; zRot1 -> rotateZ(0.*deg);
      G4RotationMatrix* zRot2 = new G4RotationMatrix; zRot2 -> rotateZ(270.*deg);
      G4RotationMatrix* zRot3 = new G4RotationMatrix; zRot3 -> rotateZ(180.*deg);
      G4RotationMatrix* zRot4 = new G4RotationMatrix; zRot4 -> rotateZ(90.*deg);

      auto scint_group_pos1 = G4ThreeVector((total_width_x-width_x)/2.,(scint_base_dist+scint_module_y/2.+dy),0.);
      auto scint_group_pos2 = G4ThreeVector(-(scint_base_dist+scint_module_y/2.+dy),(total_width_x-width_x)/2.,0.);
      auto scint_group_pos3 = G4ThreeVector(-(total_width_x-width_x)/2.,-(scint_base_dist+scint_module_y/2.+dy),0.);
      auto scint_group_pos4 = G4ThreeVector((scint_base_dist+scint_module_y/2.+dy),-(total_width_x-width_x)/2.,0.);

      for (int i; i<scint_module_pos_vector_horizontal.size(); i++){
        
        auto scint_mod_pos1 = scint_module_pos_vector_horizontal[i]+scint_group_pos1;
        auto scint_mod_pos2 = scint_module_pos_vector_vertical[i]+scint_group_pos2;
        auto scint_mod_pos3 = scint_module_pos_vector_horizontal[i]+scint_group_pos3;
        auto scint_mod_pos4 = scint_module_pos_vector_vertical[i]+scint_group_pos4;

        new G4PVPlacement(zRot1,scint_mod_pos1,scint_moduleLV,"ScintPV",mother,false,i,true);
        new G4PVPlacement(zRot2,scint_mod_pos2,scint_moduleLV,"ScintPV",mother,false,scint_module_pos_vector_horizontal.size()+i,true);
        new G4PVPlacement(zRot3,scint_mod_pos3,scint_moduleLV,"ScintPV",mother,false,2*scint_module_pos_vector_horizontal.size()+i,true);
        new G4PVPlacement(zRot4,scint_mod_pos4,scint_moduleLV,"ScintPV",mother,false,3*scint_module_pos_vector_horizontal.size()+i,true);

        scintillator_module_pos_outFile << i << "," << scint_mod_pos1[0]/cm << "," << scint_mod_pos1[1]/cm << "," << scint_mod_pos1[2]/cm << G4endl;
        scintillator_module_pos_outFile << scint_module_pos_vector_horizontal.size()+i << "," << scint_mod_pos2[0]/cm << "," << scint_mod_pos2[1]/cm << "," << scint_mod_pos2[2]/cm << G4endl;
        scintillator_module_pos_outFile << 2*scint_module_pos_vector_horizontal.size()+i << "," << scint_mod_pos3[0]/cm << "," << scint_mod_pos3[1]/cm << "," << scint_mod_pos3[2]/cm << G4endl;
        scintillator_module_pos_outFile << 3*scint_module_pos_vector_horizontal.size()+i << "," << scint_mod_pos4[0]/cm << "," << scint_mod_pos4[1]/cm << "," << scint_mod_pos4[2]/cm << G4endl;
        
        scint_index_count+=4;
      }

  ///////////////////////////////////////////////////////////////////////////////////////////
  ///////////////  Scintillator modules on front and back surfaces //////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////

      // - - - front and back scintillators
      G4double dz_fb = 20.*cm; //dist between the TPC and the front back surface
      G4double scint_bar_fb_x = 5.*cm;
      G4double scint_bar_fb_y = 30.*cm;
      G4double scint_bar_fb_z = scint_bar_y;
      G4double n_bar_x_fb = 6;

      G4double scint_module_fb_x = scint_bar_fb_x*n_bar_x_fb;
      G4double scint_module_fb_y = scint_bar_fb_y;
      G4double scint_module_fb_z = scint_bar_fb_z*scint_layers; 

      // Important number for lead glass position calculation:
      // Distance of the surface of the front scintillator from center
      G4double scint_fb_top_dist = Beampipe_5_len/2.+dz_fb+scint_module_fb_z;
      std::cout << "Scintillator :: Distance from center to front surface :: " << scint_fb_top_dist/cm << std::endl; 

            // Scint Fb group 1 
            G4double scint_fb_group_1_x = 2.*(Beampipe_5_radius_2+TPC_drift_len+2.*TPC_wall_thickness+scint_module_y+dy); // avaiable size in x
            G4double scint_fb_group_1_y = (TPC_drift_len+2.*TPC_wall_thickness+scint_module_y+dy); // avaiable size in y 
            G4double n_scint_module_fb_1_x = 15.;
            G4double n_scint_module_fb_1_y = 4.;

            G4double dx_fb_1 = (scint_fb_group_1_x-n_scint_module_fb_1_x*scint_module_fb_x)/(n_scint_module_fb_1_x-1);
            G4double dy_fb_1 = (scint_fb_group_1_y-n_scint_module_fb_1_y*scint_module_fb_y)/(n_scint_module_fb_1_y-1);

            std::cout << dx_fb_1/cm << "    "<< dy_fb_1/cm<< std::endl;

            // Scint Fb group 2 
            G4double dy_fb_12 = dy_fb_1/2.; // reserve a gap between group 1 and group 2, take the half value of y spacing in group 1
            G4double scint_fb_group_2_x = TPC_drift_len + 2.*TPC_wall_thickness+dy+scint_module_y;
            G4double scint_fb_group_2_y = 2.*Beampipe_5_radius_2-2.*dy_fb_1;
            G4double n_scint_module_fb_2_x = 4.;
            G4double n_scint_module_fb_2_y = 7.;

            G4double dx_fb_2 = (scint_fb_group_2_x-n_scint_module_fb_2_x*scint_module_fb_x)/(n_scint_module_fb_2_x-1);
            G4double dy_fb_2 = (scint_fb_group_2_y-n_scint_module_fb_2_y*scint_module_fb_y)/(n_scint_module_fb_2_y-1);

            std::cout << "group 2 " << dx_fb_2/cm << "  " << dy_fb_2/cm << std::endl;

      // virtual volume for the layers
      auto scint_module_fb_S = new G4Box("ScintS",scint_module_fb_x/2.,scint_module_fb_y/2.,scint_module_fb_z/2.);
      auto scint_module_fb_LV = new G4LogicalVolume(scint_module_fb_S,defaultMaterial,"Scint_fb_LV1");
      
          // defining the layers
          // type 1 layer 1 and 2: (1 = horizontal, 2 = vertical bars)
          auto scint_layer_fb1_S  = new G4Box("ScintS", scint_module_fb_x/2.,scint_module_fb_y/2., scint_bar_fb_z/2.);
          auto scint_layer_fb1_LV = new G4LogicalVolume(scint_layer_fb1_S,defaultMaterial,"Scint_bar_LV11");
          auto scint_layer_fb2_S  = new G4Box("ScintS", scint_module_fb_x/2.,scint_module_fb_y/2., scint_bar_fb_z/2.);
          auto scint_layer_fb2_LV = new G4LogicalVolume(scint_layer_fb2_S,defaultMaterial,"Scint_bar_LV12");

              // // bars for the front and back layers
              // bars for type 1 fb scint
              // --- vertical bars (divide x)
              G4double scint_bar_xV_fb1 = scint_bar_fb_x; 
              G4double scint_bar_yV_fb1 = scint_bar_fb_y;
              auto scint_bar_fb1V_S = new G4Box("ScintS", scint_bar_xV_fb1/2.,scint_bar_yV_fb1/2., scint_bar_fb_z/2.);
              auto scint_bar_fb1V_LV = new G4LogicalVolume(scint_bar_fb1V_S,scintMaterial,"Scint_fb_bar1V_LV");

              // --- horizontal bars (divide y)
              G4double scint_bar_xH_fb1 = scint_bar_fb_y;
              G4double scint_bar_yH_fb1 = scint_bar_fb_x;
              auto scint_bar_fb1H_S = new G4Box("ScintS", scint_bar_xH_fb1/2.,scint_bar_yH_fb1/2., scint_bar_fb_z/2.);
              auto scint_bar_fb1H_LV = new G4LogicalVolume(scint_bar_fb1H_S,scintMaterial,"Scint_fb_bar1H_LV");
              
              for (int i; i<no_of_Layers; i++){ //10
                if (i%2==0){new G4PVPlacement(0, G4ThreeVector(0.,0.,-scint_module_fb_z/2.+(2.*i+1.)/2.*scint_bar_fb_z),scint_layer_fb1_LV,"Scint_FB_layerPV",scint_module_fb_LV,false,i,true);}
                else{new G4PVPlacement(0, G4ThreeVector(0.,0.,-scint_module_fb_z/2.+(2.*i+1.)/2.*scint_bar_fb_z),scint_layer_fb2_LV,"Scint_FB_layerPV",scint_module_fb_LV,false,i,true);}
              }

              // placing the staves 
              for (int i; i<n_bar_x_fb; i++){
                // --- horizontal ones for odd number layers
                new G4PVPlacement(0,G4ThreeVector(0.,-scint_module_fb_y/2.+(2.0*i+1.0)/2.*scint_bar_yH_fb1,0.),scint_bar_fb1H_LV,"Scint_FB_barPV_H",scint_layer_fb1_LV,false,i,true);
                // --- vertical ones for even number layers
                new G4PVPlacement(0,G4ThreeVector(-scint_module_fb_x/2.+(2.0*i+1.0)/2.*scint_bar_xV_fb1,0.,0.),scint_bar_fb1V_LV,"Scint_FB_barPV_V",scint_layer_fb2_LV,false,i,true);
              }


  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////// Position calculation part of the scintillator modules ////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
    ////////////Scintillator fb group 1 Horizontal

    // rotation is introduced to make a better indexing of the layers
    G4RotationMatrix * scint_fb_front_rot =  new G4RotationMatrix; scint_fb_front_rot -> rotateX(0.0*deg);
    G4RotationMatrix * scint_fb_back_rot =  new G4RotationMatrix; scint_fb_back_rot -> rotateX(180.0*deg);

    std::vector<G4ThreeVector> scint_fb_group_1_pos_array; // position of each scintillator module
    for (int i=0; i<int(n_scint_module_fb_1_x); i++){
      for (int j=0; j<int(n_scint_module_fb_1_y); j++){

        G4double x = -scint_fb_group_1_x/2. + 0.5*scint_module_fb_x + double(i)*(scint_module_fb_x+dx_fb_1);
        G4double y = -scint_fb_group_1_y/2. + 0.5*scint_module_fb_y + double(j)*(scint_module_fb_y+dy_fb_1);
        G4double z = 0.*cm;
        scint_fb_group_1_pos_array.push_back(G4ThreeVector(x,y,z));
      }
    }

    // The group position of the modules (See the sketch for more information)
    auto scint_fb_group_1_pos_front_top = G4ThreeVector(0.,Beampipe_5_radius_2+scint_fb_group_1_y/2.,Beampipe_5_len/2.+dz_fb+scint_module_fb_z/2.);// front surface top
    auto scint_fb_group_1_pos_front_bottom = G4ThreeVector(0.,-(Beampipe_5_radius_2+scint_fb_group_1_y/2.),Beampipe_5_len/2.+dz_fb+scint_module_fb_z/2.);// front surface bottom
    auto scint_fb_group_1_pos_back_top = G4ThreeVector(0.,Beampipe_5_radius_2+scint_fb_group_1_y/2.,-(Beampipe_5_len/2.+dz_fb+scint_module_fb_z/2.));// back surface top
    auto scint_fb_group_1_pos_back_bottom = G4ThreeVector(0.,-(Beampipe_5_radius_2+scint_fb_group_1_y/2.),-(Beampipe_5_len/2.+dz_fb+scint_module_fb_z/2.));// back surface bottom
    
    for (int i; i<scint_fb_group_1_pos_array.size(); i++){

      int scint_fb_1_index_1 = scint_index_count + i;
      int scint_fb_1_index_2 = scint_index_count + scint_fb_group_1_pos_array.size()+i;
      int scint_fb_1_index_3 = scint_index_count + 2*scint_fb_group_1_pos_array.size()+i;
      int scint_fb_1_index_4 = scint_index_count + 3*scint_fb_group_1_pos_array.size()+i;
      
      auto scint_mod_pos_fb1 = scint_fb_group_1_pos_array[i]+scint_fb_group_1_pos_front_top;
      auto scint_mod_pos_fb2 = scint_fb_group_1_pos_array[i]+scint_fb_group_1_pos_front_bottom;
      auto scint_mod_pos_fb3 = scint_fb_group_1_pos_array[i]+scint_fb_group_1_pos_back_top;
      auto scint_mod_pos_fb4 = scint_fb_group_1_pos_array[i]+scint_fb_group_1_pos_back_bottom;

      scintillator_module_pos_outFile << scint_fb_1_index_1 << "," << scint_mod_pos_fb1[0]/cm << "," << scint_mod_pos_fb1[1]/cm << "," << scint_mod_pos_fb1[2]/cm << G4endl;
      scintillator_module_pos_outFile << scint_fb_1_index_2 << "," << scint_mod_pos_fb2[0]/cm << "," << scint_mod_pos_fb2[1]/cm << "," << scint_mod_pos_fb2[2]/cm << G4endl;
      scintillator_module_pos_outFile << scint_fb_1_index_3 << "," << scint_mod_pos_fb3[0]/cm << "," << scint_mod_pos_fb3[1]/cm << "," << scint_mod_pos_fb3[2]/cm << G4endl;
      scintillator_module_pos_outFile << scint_fb_1_index_4 << "," << scint_mod_pos_fb4[0]/cm << "," << scint_mod_pos_fb4[1]/cm << "," << scint_mod_pos_fb4[2]/cm << G4endl;

      new G4PVPlacement(scint_fb_front_rot,scint_fb_group_1_pos_array[i]+scint_fb_group_1_pos_front_top,scint_module_fb_LV,"Scint_FB_PV",mother,false,scint_fb_1_index_1,true);
      new G4PVPlacement(scint_fb_front_rot,scint_fb_group_1_pos_array[i]+scint_fb_group_1_pos_front_bottom,scint_module_fb_LV,"Scint_FB_PV",mother,false,scint_fb_1_index_2,true);
      new G4PVPlacement(scint_fb_back_rot,scint_fb_group_1_pos_array[i]+scint_fb_group_1_pos_back_top,scint_module_fb_LV,"Scint_FB_PV",mother,false,scint_fb_1_index_3,true);
      new G4PVPlacement(scint_fb_back_rot,scint_fb_group_1_pos_array[i]+scint_fb_group_1_pos_back_bottom,scint_module_fb_LV,"Scint_FB_PV",mother,false,scint_fb_1_index_4,true);

    }
    
    scint_index_count+=4*scint_fb_group_1_pos_array.size();

  ////////////Scintillator fb group 2 Vertical

    std::vector<G4ThreeVector> scint_fb_group_2_pos_array;
    for (int i=0; i<int(n_scint_module_fb_2_x); i++){
      for (int j=0; j<int(n_scint_module_fb_2_y); j++){

        G4double x = -scint_fb_group_2_x/2. + 0.5*scint_module_fb_x + double(i)*(scint_module_fb_x+dx_fb_2);
        G4double y = -scint_fb_group_2_y/2. + 0.5*scint_module_fb_y + double(j)*(scint_module_fb_y+dy_fb_2);
        G4double z = 0.*cm;
        scint_fb_group_2_pos_array.push_back(G4ThreeVector(x,y,z));
      }
    }

    auto scint_fb_group_2_pos_front_left  = G4ThreeVector(Beampipe_5_radius_2+scint_fb_group_2_x/2.,0.,Beampipe_5_len/2.+dz_fb+scint_module_fb_z/2.);// front surface top
    auto scint_fb_group_2_pos_front_right = G4ThreeVector(-(Beampipe_5_radius_2+scint_fb_group_2_x/2.),0.,Beampipe_5_len/2.+dz_fb+scint_module_fb_z/2.);// front surface bottom
    auto scint_fb_group_2_pos_back_left  = G4ThreeVector(Beampipe_5_radius_2+scint_fb_group_2_x/2.,0.,-(Beampipe_5_len/2.+dz_fb+scint_module_fb_z/2.));// front surface top
    auto scint_fb_group_2_pos_back_right = G4ThreeVector(-(Beampipe_5_radius_2+scint_fb_group_2_x/2.),0.,-(Beampipe_5_len/2.+dz_fb+scint_module_fb_z/2.));// front surface bottom
    
    for (int i; i<scint_fb_group_2_pos_array.size(); i++){

      int scint_fb_2_index_1 = scint_index_count + i;
      int scint_fb_2_index_2 = scint_index_count + scint_fb_group_2_pos_array.size()+i;
      int scint_fb_2_index_3 = scint_index_count + 2*scint_fb_group_2_pos_array.size()+i;
      int scint_fb_2_index_4 = scint_index_count + 3*scint_fb_group_2_pos_array.size()+i;
      
      auto scint_mod_pos_fb1 = scint_fb_group_2_pos_array[i]+scint_fb_group_2_pos_front_left;
      auto scint_mod_pos_fb2 = scint_fb_group_2_pos_array[i]+scint_fb_group_2_pos_front_right;
      auto scint_mod_pos_fb3 = scint_fb_group_2_pos_array[i]+scint_fb_group_2_pos_back_left;
      auto scint_mod_pos_fb4 = scint_fb_group_2_pos_array[i]+scint_fb_group_2_pos_back_right;

      scintillator_module_pos_outFile << scint_fb_2_index_1 << "," << scint_mod_pos_fb1[0]/cm << "," << scint_mod_pos_fb1[1]/cm << "," << scint_mod_pos_fb1[2]/cm << G4endl;
      scintillator_module_pos_outFile << scint_fb_2_index_2 << "," << scint_mod_pos_fb2[0]/cm << "," << scint_mod_pos_fb2[1]/cm << "," << scint_mod_pos_fb2[2]/cm << G4endl;
      scintillator_module_pos_outFile << scint_fb_2_index_3 << "," << scint_mod_pos_fb3[0]/cm << "," << scint_mod_pos_fb3[1]/cm << "," << scint_mod_pos_fb3[2]/cm << G4endl;
      scintillator_module_pos_outFile << scint_fb_2_index_4 << "," << scint_mod_pos_fb4[0]/cm << "," << scint_mod_pos_fb4[1]/cm << "," << scint_mod_pos_fb4[2]/cm << G4endl;

      new G4PVPlacement(scint_fb_front_rot,scint_mod_pos_fb1,scint_module_fb_LV,"Scint_FB_PV",mother,false,scint_fb_2_index_1,true);
      new G4PVPlacement(scint_fb_front_rot,scint_mod_pos_fb2,scint_module_fb_LV,"Scint_FB_PV",mother,false,scint_fb_2_index_2,true);
      new G4PVPlacement(scint_fb_back_rot,scint_mod_pos_fb3,scint_module_fb_LV,"Scint_FB_PV",mother,false,scint_fb_2_index_3,true);
      new G4PVPlacement(scint_fb_back_rot,scint_mod_pos_fb4,scint_module_fb_LV,"Scint_FB_PV",mother,false,scint_fb_2_index_4,true);
    }

  scintillator_module_pos_outFile.close();
  ////////////////////////////////////////////////////////////
  /////////// Optics settings of the scintillator bars //////
  //////////////////////////////////////////////////////////

    G4OpticalSurface* op_scintillator = new G4OpticalSurface("scintillator_optical_surface");
    op_scintillator->SetType(dielectric_metal);
    op_scintillator->SetFinish(polished);
    op_scintillator->SetModel(unified);

    G4double pp_scintillator[] = { 2.0 * eV, 3.5 * eV }; const G4int num_scintillator = sizeof(pp_scintillator) / sizeof(G4double);
    G4double reflectivity_scintillator[] = { 0.0, 0.0}; G4double efficiency_scintillator[] = { 0.0, 0.0 };
    G4MaterialPropertiesTable* Scintillator_Optical_Property = new G4MaterialPropertiesTable();
    Scintillator_Optical_Property->AddProperty("REFLECTIVITY", pp_scintillator, reflectivity_scintillator, num_scintillator);
    Scintillator_Optical_Property->AddProperty("EFFICIENCY", pp_scintillator, efficiency_scintillator, num_scintillator);
    op_scintillator ->SetMaterialPropertiesTable(Scintillator_Optical_Property);

    //scintillators
    new G4LogicalSkinSurface("name",scint_barH_LV,op_scintillator);
    new G4LogicalSkinSurface("name",scint_barV_LV,op_scintillator);
    new G4LogicalSkinSurface("name",scint_bar_fb1V_LV,op_scintillator);
    new G4LogicalSkinSurface("name",scint_bar_fb1H_LV,op_scintillator);

/////////////////////////////////////////////////////////
////////// Region settings for scintillator bars ///////
////////////////////////////////////////////////////////

    // cuts for scintillator
    G4Region* scint_region = new G4Region("Scint_region");
    scint_region->AddRootLogicalVolume(scint_barH_LV);
    scint_region->AddRootLogicalVolume(scint_barV_LV);
    scint_region->AddRootLogicalVolume(scint_bar_fb1V_LV);
    scint_region->AddRootLogicalVolume(scint_bar_fb1H_LV);

    G4Region* scint_Region = G4RegionStore::GetInstance()->GetRegion("Scint_region");
    G4ProductionCuts* scintcut = new G4ProductionCuts();
    scintcut->SetProductionCut(6.0*cm,"gamma");
    scintcut->SetProductionCut(2.0*mm,"e-");
    scintcut->SetProductionCut(2.0*mm,"e+");
    scintcut->SetProductionCut(2.0*mm,"proton");
    scint_Region->SetProductionCuts(scintcut);

/////////////////////////////////////////////////////////
///////////////// color settings ///////////////////////
/////////////////////////////////////////////////////////

  auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true);
  auto green_color= new G4VisAttributes(G4Colour(0.788235,0.890196,0.741176)); green_color->SetVisibility(true);
  // original green: 0.517647,0.772549,0.556863
  auto orange_color= new G4VisAttributes(G4Colour(0.988235,0.686275,0.243137)); orange_color->SetVisibility(true);
  auto red_color= new G4VisAttributes(G4Colour(0.956863,0.0901961,0.494118)); red_color->SetVisibility(true); 
  auto blue_color= new G4VisAttributes(G4Colour(0.447059,0.623529,0.811765)); blue_color->SetVisibility(true);
  auto dark_blue_color= new G4VisAttributes(G4Colour(0.116,0.43,0.895)); dark_blue_color->SetVisibility(true);
  auto pink_color= new G4VisAttributes(G4Colour(0.8,0.29,0.61)); pink_color->SetVisibility(true);
  auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549)); black_color->SetVisibility(true);

  scint_moduleLV->SetVisAttributes(blue_color);;
  scint_layerH_LV->SetVisAttributes(G4VisAttributes::GetInvisible());//
  scint_layerV_LV->SetVisAttributes(G4VisAttributes::GetInvisible());//
  scint_barH_LV->SetVisAttributes (G4VisAttributes::GetInvisible());
  scint_barV_LV->SetVisAttributes (G4VisAttributes::GetInvisible());

  scint_module_fb_LV->SetVisAttributes(dark_blue_color);
  scint_layer_fb1_LV->SetVisAttributes (G4VisAttributes::GetInvisible());
  scint_layer_fb2_LV->SetVisAttributes (G4VisAttributes::GetInvisible());
  scint_bar_fb1V_LV->SetVisAttributes (G4VisAttributes::GetInvisible());
  scint_bar_fb1H_LV->SetVisAttributes (G4VisAttributes::GetInvisible());


  Scintillator_output.push_back(scint_barH_LV);
  Scintillator_output.push_back(scint_barV_LV);
  Scintillator_output.push_back(scint_bar_fb1H_LV);
  Scintillator_output.push_back(scint_bar_fb1V_LV);
  
  return Scintillator_output;
  
}