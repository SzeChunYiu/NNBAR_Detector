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
#include "detector/beampipe_geometry.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

class Silicon : public G4VUserDetectorConstruction
{
  public:
    Silicon();
    virtual ~Silicon();

  public:
    virtual G4VPhysicalVolume* Construct() override { return nullptr; }
    std::vector<G4LogicalVolume*> Construct_Volumes(G4LogicalVolume* mother);

  private:
    void DefineMaterials();
};

#endif
