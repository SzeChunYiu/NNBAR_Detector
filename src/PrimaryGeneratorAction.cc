#include "PrimaryGeneratorAction.hh"
#include <iomanip>
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4GenericMessenger.hh"
#include "Analysis.hh"


using namespace std;

extern std::ofstream Particle_outFile;
extern std::vector<std::vector<G4double>> particle_gun_record;
extern G4double event_number;
extern G4int run_number;

PrimaryGeneratorAction::PrimaryGeneratorAction():fParticleGun(nullptr){fParticleGun = new G4ParticleGun();}


PrimaryGeneratorAction::~PrimaryGeneratorAction(){}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event * anEvent)
{
	G4double x; G4double y; G4double z;
	G4double t; G4double px; G4double py; G4double pz;
	G4double KE;

	G4String particleName;

	std::vector<G4double> particle_gun_record_row;

	x = 0.0* m;
	y = 0.0*m;
	z = 0.0* m; // (*vect)[j]->z() * m
	KE = (50.0+ 50.0*std::floor(event_number/1000)) *MeV ; //250.0 * MeV;
	px = G4UniformRand();
	py = G4UniformRand();
	pz = G4UniformRand();
	t = 0.; 

	fParticleGun->SetParticleDefinition(particleTable->FindParticle("pi+"));
	//fParticleGun->SetParticleEnergy(KE);
	
	int n = event_number/100; // every energy, it has 50 entries
	fParticleGun->SetParticleEnergy(KE);	
	
	fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz)); //(*vect)[j]->w())
	fParticleGun->SetParticleTime(0.0);

	Particle_outFile <<  event_number << ",";
	Particle_outFile <<  event_number << ",";
	Particle_outFile <<  211 << ","; //PID
	Particle_outFile <<  particleTable -> FindParticle("pi+") -> GetPDGMass() << ","; //Mass
	Particle_outFile <<  0 << ","; //PID
	Particle_outFile <<  KE << ","; // 100+n*10 
	Particle_outFile <<  x << ",";
	Particle_outFile <<  y << ",";
	Particle_outFile <<  z << ",";
	Particle_outFile <<  t << ",";
	Particle_outFile <<  px << ",";
	Particle_outFile <<  py << ",";
	Particle_outFile <<  pz << G4endl;

	fParticleGun->GeneratePrimaryVertex(anEvent);
	

}

//....

