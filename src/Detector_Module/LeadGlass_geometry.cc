// LeadGlass
#include "DetectorConstruction.hh"
#include "Detector_Module/LeadGlass_geometry.hh"
#include "Detector_Module/Scintillator_geometry.hh"
#include "Detector_Module/TPC_geometry.hh"

#include "G4SDManager.hh"

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
#include <boost/lexical_cast.hpp>

#include <G4ProductionCuts.hh>
#include "G4RegionStore.hh"
#include "G4Region.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"

using boost::lexical_cast;


extern G4double Beampipe_5_radius_2; // need to calculate how many scintilator modules we put along z direction! 
extern G4double TPC_drift_len;
extern G4double TPC_wall_thickness;

std::vector<G4LogicalVolume*> LeadGlass_Construction_list;

std::string filename_data_lead_glass_pos = "./lead_glass_position/lead_glass_position.csv";
std::string filename_data_lead_glass_pos_fb = "./lead_glass_position/lead_glass_position_fb.csv";

std::vector<std::vector<double>> data_lead_glass_pos;
std::vector<std::vector<double>> data_lead_glass_pos_fb;

void import_lead_glass_pos(std::string file_name, std::vector<std::vector<double> >& data) {
	
	std::string row;
	std::ifstream init_file(file_name.c_str());

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

LeadGlass::LeadGlass():G4VUserDetectorConstruction(){}

LeadGlass::~LeadGlass(){}

void LeadGlass::DefineMaterials()
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

  // Lead-glass (taken from PDG)
  G4Material* LeadGlass = new G4Material("LeadGlass", 6.22*g/cm3, 5);
  LeadGlass->AddElement(elO, 0.156453); LeadGlass->AddElement(elSi, 0.080866); LeadGlass->AddElement(elTi, 0.008092); LeadGlass->AddElement(elAs, .002651); LeadGlass->AddElement(elPb, 0.751938);

  // MgF2 coating of the lead glass
  G4Material *AlMgF2 = new G4Material("AlMgF2", density = 2.9007*g/cm3, 3);
  AlMgF2->AddElement(elAl, 0.331); AlMgF2->AddElement(elF, 0.408); AlMgF2->AddElement(elMg, 0.261);

  // PMT window quartz actually
  G4Material* PMT_window = new G4Material("PMT_window_mat",density= 2.200*g/cm3, 2);
  PMT_window->AddElement(elSi, 1); PMT_window->AddElement(elO , 2); 
  // -- optical properties of PMT
  G4double pmt_window_Energy[] = { 2.0 * eV, 7.0 * eV}; const G4int pmt_window_num = sizeof(pmt_window_Energy) / sizeof(G4double); G4double pmt_window_RIND[] = { 1.53,1.53 };
  G4MaterialPropertiesTable * pmt_window_mt = new G4MaterialPropertiesTable();
  pmt_window_mt->AddProperty("RINDEX", pmt_window_Energy, pmt_window_RIND, pmt_window_num); PMT_window->SetMaterialPropertiesTable(pmt_window_mt);
  
  ////////////////////////////////////////////////////
  //// Lead Glass Schott SF5 Optical properties //////
  ////////////////////////////////////////////////////

    G4double PhotonWavelength[] =
          { 2325.4*nm, 1970.1*nm, 1529.6*nm, 1060.0*nm,
            1014.0*nm, 852.10*nm, 706.50*nm, 656.30*nm,
            643.80*nm, 632.80*nm, 589.30*nm, 587.60*nm,
            546.10*nm, 486.10*nm, 480.00*nm, 435.80*nm,
            404.70*nm, 365.00*nm};

    const G4int nEntries = sizeof(PhotonWavelength)/sizeof(G4double);
    G4double PhotonEnergy[nEntries]; for (int i=0; i < nEntries; ++i) {PhotonEnergy[i] = (1240.*nm/PhotonWavelength[i])*eV;};

    G4double refractiveIndex[] =
            { 1.63289, 1.63785, 1.64359, 1.65104,
              1.65206, 1.65664, 1.66327, 1.66661,
              1.66756, 1.66846, 1.67252, 1.67270,
              1.67764, 1.68750, 1.68876, 1.69986,
              1.71069, 1.73056};

    G4MaterialPropertiesTable* LeadGlass_MPT = new G4MaterialPropertiesTable();
    LeadGlass_MPT->AddProperty("RINDEX", PhotonEnergy, refractiveIndex, nEntries)->SetSpline(true);
    LeadGlass->SetMaterialPropertiesTable(LeadGlass_MPT);
    LeadGlass_MPT->DumpTable();
}

std::vector<G4LogicalVolume*> LeadGlass::Construct_Volumes(G4LogicalVolume* mother)
{
  // OUR LeadGlass GEOMETRY OUTPUT SENT BACK TO DETECTOR_CONSTRUCTION.CC //
  //std::vector<G4LogicalVolume*> LeadGlass_Construction_list;
  //####################################################################//

  ////////////////////////////////////////////////////
  ///////////// Defining materials ///////////////////
  ////////////////////////////////////////////////////

    DefineMaterials();
    auto defaultMaterial = G4Material::GetMaterial("Galactic");
    auto LeadGlassMaterial = G4Material::GetMaterial("LeadGlass");
    auto CoatingMaterial = G4Material::GetMaterial("AlMgF2");
    auto PMTMaterial = G4Material::GetMaterial("PMT_window_mat");

  ///////////////////////////////////////////////////////////
  //////////// Numbers for LeadGlass construction ///////////
  ///////////////////////////////////////////////////////////

    G4double lead_glass_x = 8.*cm; 
    G4double lead_glass_y = 25.*cm;
    G4double lead_glass_z = 8.*cm;
    
    G4double coating_thickness = 0.01*mm;
    
    G4double PMT_thickness = 0.01*mm;
    G4double PMT_radius = 5.*cm;
    G4double offset_lead_glass = 0.0*mm;

    G4double lead_glass_y_level = Beampipe_5_radius_2 + 2.*TPC_wall_thickness + TPC_drift_len + 30.*cm + 1.0*cm; 
      
  ////////////////////////////////////////////////////////////
  ///////////// LeadGlass Construction ///////////////////////
  ////////////////////////////////////////////////////////////

    int lead_index = 0;

    import_lead_glass_pos(filename_data_lead_glass_pos, data_lead_glass_pos);
    import_lead_glass_pos(filename_data_lead_glass_pos_fb, data_lead_glass_pos_fb);
    // Sturcture goes like this:
    // LeadGlass_module
    // -- Lead Glass body
    // -- Coating (hallow)
    // -- PMT 
    // # Note: the coating is on the sides of the Lead glass, not the top or bottom
    // # Under the optical setting, the photons hitting the sides should be reflected
    // # while the photons hitting the bottom should be killed
 
    auto LeadGlass_module_S = new G4Box("LeadGlass_module_S",(lead_glass_x+2.*coating_thickness)/2., (lead_glass_y+PMT_thickness)/2., (lead_glass_z+2.*coating_thickness)/2.); // is the whole lead glass module including the coating
    auto LeadGlass_module_LV = new G4LogicalVolume(LeadGlass_module_S,defaultMaterial,"LeadGlass_module_LV");

        // Lead glass
        auto leadglassS = new G4Box("LeadGlassS",(lead_glass_x)/2., (lead_glass_y)/2., (lead_glass_z)/2.); // just the lead glass itself
        auto leadglassLV = new G4LogicalVolume(leadglassS,LeadGlassMaterial,"LeadGlassLV");
        new G4PVPlacement(0,G4ThreeVector(0., -(lead_glass_y+PMT_thickness)/2.+(lead_glass_y)/2., 0.),leadglassLV,"LeadGlassPV",LeadGlass_module_LV,false,0,true);

        // Coating of the lead glass 
        auto hallowS = new G4Box("HallowS",(lead_glass_x)/2., (lead_glass_y+PMT_thickness)/2., (lead_glass_z)/2.); // a copy of the lead glass for creating the hallow
  
        auto coatingS = new G4SubtractionSolid("CoatingS", LeadGlass_module_S, hallowS,0, G4ThreeVector(0.,0.,0.));
        auto coatingLV = new G4LogicalVolume(coatingS,defaultMaterial,"LeadGlass_CoatingLV");                              
        new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),coatingLV,"LeadGlass_CoatingPV",LeadGlass_module_LV,false,0,true);

        // virtual PMT 
        auto PMTS = new G4Box("PMTS",(lead_glass_x)/2., PMT_thickness/2., (lead_glass_z)/2.);
        auto PMTLV = new G4LogicalVolume(PMTS,PMTMaterial,"PMTLV");            
        new G4PVPlacement(0,G4ThreeVector(0., -(lead_glass_y+PMT_thickness)/2.+lead_glass_y+PMT_thickness/2., 0.),PMTLV,"PMTPV",LeadGlass_module_LV,false,0,true);
    
  //////////////////////////////////////////////////////////////////////
  //////////////////// Lead Glass Placement ////////////////////////////
  //////////////////////////////////////////////////////////////////////

    std::cout << " Constructing the lead glass blocks from data" << std::endl;
    std::vector<std::vector<G4RotationMatrix *>> rot_array_dir;
    std::vector<std::vector<G4VPhysicalVolume *>> lg_PV_array;
    G4double rot_angle = 90.0*deg;

        for (int i = 0; i < 4; i++){ 
          std::vector<G4RotationMatrix *>temp_rot_array_dir;
          std::vector<G4VPhysicalVolume*>temp_PV_array;
          for (int j = 0 ; j < data_lead_glass_pos.size();j++){temp_rot_array_dir.push_back(new G4RotationMatrix);}
          rot_array_dir.push_back(temp_rot_array_dir);
          lg_PV_array.push_back(temp_PV_array);
        }
      
        // Placing the lead glass blocks on the 4 surfaces 
        for (int i = 0 ; i<4;i++){
          for (int j = 0 ; j <data_lead_glass_pos.size();j++){ // data_lead_glass_pos.size()

            G4double lead_x = 0.; G4double lead_y=0.;
            G4double lead_x0 = 0. ; G4double lead_y0 = 0.*cm;
            lead_x0 = data_lead_glass_pos[j][0]*cm; lead_y0 = data_lead_glass_pos[j][1]*cm + 1.5*cm; 

            double i_angle = (double) i;
            lead_x = lead_x0*cos(i_angle*rot_angle) - (lead_y0)*sin(i_angle*rot_angle);
            lead_y = lead_x0*sin(i_angle*rot_angle) + (lead_y0)*cos(i_angle*rot_angle);
            rot_array_dir[i][j]->rotate(-data_lead_glass_pos[j][3]*deg, G4ThreeVector(cos(i_angle*rot_angle),sin(i_angle*rot_angle),0.));
            rot_array_dir[i][j]-> rotateZ(data_lead_glass_pos[j][5]*deg -i_angle*rot_angle); //data_lead_glass_pos[i][5]*deg 
            lg_PV_array[i].push_back(new G4PVPlacement(rot_array_dir[i][j],G4ThreeVector(lead_x, lead_y,data_lead_glass_pos[j][2]*cm) ,LeadGlass_module_LV,"LeadGlassPV",mother,false,i*data_lead_glass_pos.size()+j,false));
            
            //std::cout << i*data_lead_glass_pos.size()+j << std::endl;
            lead_index ++ ;
          }
        }
  
    // for the Front and back lead glass
    std::vector<G4RotationMatrix *> rot_array_dir111;
    std::vector<G4RotationMatrix *> rot_array_dir211;
    int lead_index_fb = 4*data_lead_glass_pos.size();

        for (int i = 0 ; i < data_lead_glass_pos_fb.size();i++){rot_array_dir111.push_back(new G4RotationMatrix);}
        for (int i = 0 ; i < data_lead_glass_pos_fb.size();i++){rot_array_dir211.push_back(new G4RotationMatrix);}
        
        for (int i = 0 ; i < data_lead_glass_pos_fb.size();i++){
              rot_array_dir111[i]-> rotateX(-data_lead_glass_pos_fb[i][3]*deg+90.*deg); rot_array_dir111[i]-> rotateZ(data_lead_glass_pos_fb[i][5]*deg); //rot_array_dir111[i]-> rotateY(-data_lead_glass_pos_fb[i][5]*deg);
              new G4PVPlacement(rot_array_dir111[i],
              G4ThreeVector(data_lead_glass_pos_fb[i][0]*cm,data_lead_glass_pos_fb[i][2]*cm,-(data_lead_glass_pos_fb[i][1]*cm))
              ,LeadGlass_module_LV,"LeadGlassPV",mother,false,lead_index_fb+i,false);    
        }
    
    lead_index_fb = 4*data_lead_glass_pos.size() + data_lead_glass_pos_fb.size();

        for (int i = 0 ; i < data_lead_glass_pos_fb.size();i++){
              rot_array_dir211[i]-> rotateX(data_lead_glass_pos_fb[i][3]*deg-90.*deg); rot_array_dir211[i]-> rotateZ(data_lead_glass_pos_fb[i][5]*deg); //rot_array_dir111[i]-> rotateY(-data_lead_glass_pos_fb[i][5]*deg);
              new G4PVPlacement(rot_array_dir211[i],
              G4ThreeVector(data_lead_glass_pos_fb[i][0]*cm,data_lead_glass_pos_fb[i][2]*cm,(data_lead_glass_pos_fb[i][1]*cm))
              ,LeadGlass_module_LV,"LeadGlassPV",mother,false,lead_index_fb+i,false);  
        }
        
  //////////////////////////////////////////////////////////////////////
  //////////////////// Optics settings /////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  
  G4OpticalSurface* op_coating= new G4OpticalSurface("coating");
  op_coating->SetType(dielectric_metal);
  op_coating->SetFinish(polished);
  op_coating->SetModel(unified);
  
  G4double pp_coating[] = { 2.0 * eV, 3.5 * eV }; const G4int num_coating = sizeof(pp_coating) / sizeof(G4double);
  G4double reflectivity_coating[] = { 0.95, 0.95}; G4double efficiency_coating[] = { 0.0, 0.0 };
  G4MaterialPropertiesTable* CoatingProperty = new G4MaterialPropertiesTable();
  CoatingProperty->AddProperty("REFLECTIVITY", pp_coating, reflectivity_coating, num_coating);
  CoatingProperty->AddProperty("EFFICIENCY", pp_coating, efficiency_coating, num_coating);
  op_coating ->SetMaterialPropertiesTable(CoatingProperty);
  
  new G4LogicalSkinSurface("name",coatingLV,op_coating);

  
  //////////////////////////////////////////////////////////////////////
  //////////////////// Region Settings /////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  
  G4Region* LeadGlass_region = new G4Region("LeadGlass_region");
  LeadGlass_region->AddRootLogicalVolume(leadglassLV); 
  G4Region* LeadGlassRegion = G4RegionStore::GetInstance()->GetRegion("LeadGlass_region");
  G4ProductionCuts* LeadGlasscut = new G4ProductionCuts();
  LeadGlasscut->SetProductionCut(5.0 *mm,"gamma");  // 15 cm 
  LeadGlasscut->SetProductionCut(5.0 *mm,"e-"); // 1 mm
  LeadGlasscut->SetProductionCut(5.0 *mm,"e+");
  LeadGlasscut->SetProductionCut(5.0 *mm,"proton");
  LeadGlassRegion->SetProductionCuts(LeadGlasscut);

  //////////////////////////////////////////////////////////////////////
  //////////////////// Writing the outputs /////////////////////////////
  //////////////////////////////////////////////////////////////////////

  LeadGlass_Construction_list.push_back(leadglassLV);
  LeadGlass_Construction_list.push_back(PMTLV);

  ////////////////////////////////////////////////////
  ///////////////// Color settings ///////////////////
  ////////////////////////////////////////////////////

    auto black_color = new G4VisAttributes(G4Colour(0.333333,0.341176,0.32549));black_color->SetVisibility(true);
    auto grey_color= new G4VisAttributes(G4Colour(0.533333,0.541176,0.521569)); grey_color->SetVisibility(true);
    auto green_color= new G4VisAttributes(G4Colour(0.788235,0.890196,0.741176)); green_color->SetVisibility(true);
    auto orange_color= new G4VisAttributes(G4Colour(0.988235,0.686275,0.243137)); orange_color->SetVisibility(true);
    auto purple_color= new G4VisAttributes(G4Colour(0.48,0.27,0.833)); purple_color->SetVisibility(true);
    
    leadglassLV->SetVisAttributes(green_color);
    coatingLV -> SetVisAttributes (G4VisAttributes::GetInvisible());
    PMTLV->SetVisAttributes (G4VisAttributes::GetInvisible());//->SetVisAttributes(grey_color);
    LeadGlass_module_LV->SetVisAttributes (G4VisAttributes::GetInvisible());

  return LeadGlass_Construction_list;
}