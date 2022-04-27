#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4DataVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "globals.hh"
//#include "RNGWrapper.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"

using namespace std;


class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
  PrimaryGeneratorAction();
  virtual ~PrimaryGeneratorAction();

  public:
	void GeneratePrimaries(G4Event* anEvent);

  void SetRandomFlag(G4bool value);

  private:
    G4ParticleGun*  fParticleGun; 

  private:
	G4GenericMessenger* fMessenger;
        G4ParticleTable* particleTable;
	G4ParticleGun* particleGun;
	PrimaryGeneratorMessenger* gunMessenger;
};

//....

#endif
