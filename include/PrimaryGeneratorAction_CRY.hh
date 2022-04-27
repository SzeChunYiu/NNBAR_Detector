#ifndef PrimaryGeneratorAction_CRY_h
#define PrimaryGeneratorAction_CRY_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "G4DataVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleGun.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "CRYGenerator.h"
#include "CRYParticle.h"
#include "CRYUtils.h"
#include "vector"
#include "RNGWrapper.hh"
#include "CRYSetup.h"
#include "PrimaryGeneratorMessenger_CRY.hh"
#include "G4ParticleDefinition.hh"
#include "G4GenericMessenger.hh"

using namespace std;


class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction_CRY : public G4VUserPrimaryGeneratorAction
{
  public:
  PrimaryGeneratorAction_CRY();
  virtual ~PrimaryGeneratorAction_CRY();

  public:
	void GeneratePrimaries(G4Event* anEvent);
	void InputCRY();
	void UpdateCRY(std::string* MessInput);
	void CRYFromFile(G4String newValue);
  
  // set methods
  void SetRandomFlag(G4bool value);

  private:
    G4ParticleGun*  fParticleGun; 

  private:
	G4GenericMessenger* fMessenger;
	std::vector<CRYParticle*>* vect; // vector of generated particles
	G4ParticleTable* particleTable;
	G4ParticleGun* particleGun;
	CRYGenerator* gen;
	G4int InputState;
	
	PrimaryGeneratorMessenger_CRY* gunMessenger;
};

//....

#endif
