#ifndef Beampipe_Shielding_geometry_h
#define Beampipe_Shielding_geometry_h 1

#include "globals.hh"
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Cons.hh"  
#include "G4SubtractionSolid.hh"  


class Beampipe_Shielding : public G4VUserDetectorConstruction
{
  public:
    Beampipe_Shielding();
    ~Beampipe_Shielding();

  public:
    std::vector<G4LogicalVolume*> Construct_Volumes(G4LogicalVolume* mother);

  private:
    void DefineMaterials();
    G4bool  fCheckOverlaps = true; // option to activate checking of volumes overlaps
};

#endif

