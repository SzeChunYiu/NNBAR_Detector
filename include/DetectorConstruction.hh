#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4GlobalMagFieldMessenger;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

  public:	
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();

  private:
    // methods
    void DefineMaterials();
    G4VPhysicalVolume* DefineVolumes();
    // data members
    G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
};

//....

#endif

