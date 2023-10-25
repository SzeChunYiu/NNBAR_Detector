#ifndef Silicon_geometry_h
#define Silicon_geometry_h 1

#include "globals.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "Detector_Module/beampipe_geometry.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

class Silicon : public G4VUserDetectorConstruction
{
  public:
    Silicon();
    virtual ~Silicon();

  public:
    std::vector<G4LogicalVolume*> Construct_Volumes(G4LogicalVolume* mother);

  private:
    void DefineMaterials();
    G4bool  fCheckOverlaps = true; // option to activate checking of volumes overlaps
};

#endif
