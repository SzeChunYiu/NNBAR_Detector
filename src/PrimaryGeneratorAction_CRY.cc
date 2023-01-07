#include "PrimaryGeneratorAction_CRY.hh"
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
#include <iostream>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include  "mcpl.h"

using namespace std;

extern std::ofstream Particle_outFile;
extern std::vector<std::vector<G4double>> particle_gun_record;
extern G4double event_number_global;
extern G4int run_number;

G4ThreadLocal G4int local_event_number_CRY;
namespace {G4Mutex PrimaryGeneratorCRYMutex = G4MUTEX_INITIALIZER;}
// extern mcpl_outfile_t f;
// extern mcpl_particle_t * p;

std::vector<std::vector<double>> Energy_range
{
	{0.1, 500.},
	{501., 1000.},
	{1001., 5000.},
	{5001.,10000.},
	{10001.,50000.},
	{50001.,100000.}
};

boost::random::mt19937 rng;

PrimaryGeneratorAction_CRY::PrimaryGeneratorAction_CRY():fParticleGun(nullptr)
{
	//G4AutoLock lock(&PrimaryGeneratorCRYMutex);
	const char* inputfile = "setup.file";
	fMessenger = new G4GenericMessenger(this, "/particle_generator/", "Name the particle for the file name");
	fParticleGun = new G4ParticleGun();

	std::ifstream inputFile;
	inputFile.open(inputfile, std::ios::in);

	char buffer[1000];

	if (inputFile.fail()) {
		if (*inputfile != 0)  //....only complain if a filename was given
			std::cerr << "PrimaryGeneratorAction: Failed to open CRY input file= " << inputfile << std::endl;
		InputState = -1;
	}

	else {
		std::string setupString("");
		while (!inputFile.getline(buffer, 1000).eof()) {
			setupString.append(buffer);
			setupString.append(" ");
		}

		CRYSetup* setup = new CRYSetup(setupString, "/home/billy/nnbar/cry_v1.7/data");

		gen = new CRYGenerator(setup);

		// set random number generator
		RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
		setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
		InputState = 0;

	}
	// create a vector to store the CRY particle properties
	vect = new std::vector<CRYParticle*>;
	particleTable = G4ParticleTable::GetParticleTable();
	gunMessenger = new PrimaryGeneratorMessenger_CRY(this);
	std::cerr << "Input state: " << InputState << std::endl;
	std::cerr << "particle table: " << particleTable << std::endl;
  
}

//....

PrimaryGeneratorAction_CRY::~PrimaryGeneratorAction_CRY()
{}

//....

void PrimaryGeneratorAction_CRY::InputCRY()
{
	InputState = 1;
}

//----------------------------------------------------------------------------//
void PrimaryGeneratorAction_CRY::UpdateCRY(std::string * MessInput)
{
	//G4AutoLock lock(&PrimaryGeneratorCRYMutex);
	CRYSetup* setup = new CRYSetup(*MessInput, "/home/billy/nnbar/cry_v1.7/data");
	gen = new CRYGenerator(setup);

	// set random number generator
	RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
	setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
	InputState = 0;

}

void PrimaryGeneratorAction_CRY::CRYFromFile(G4String newValue)
{
	//G4AutoLock lock(&PrimaryGeneratorCRYMutex);
	// Read the cry input file
	std::ifstream inputFile;
	inputFile.open(newValue, std::ios::in);
	char buffer[1000];

	if (inputFile.fail()) {
		std::cerr << "Failed to open input file " << newValue << std::endl;
		std::cerr << "Make sure to define the cry library on the command line" << std::endl;
		InputState = -1;
	}
	else {

		std::string setupString("");
		while (!inputFile.getline(buffer, 1000).eof()) {
			setupString.append(buffer);
			setupString.append(" ");
		}

		CRYSetup* setup = new CRYSetup(setupString, "/home/billy/nnbar/cry_v1.7/data");

		gen = new CRYGenerator(setup);

		// set random number generator
		RNGWrapper<CLHEP::HepRandomEngine>::set(CLHEP::HepRandom::getTheEngine(), &CLHEP::HepRandomEngine::flat);
		setup->setRandomFunction(RNGWrapper<CLHEP::HepRandomEngine>::rng);
		InputState = 0;
	}

	std::cerr << "Input state after CRYFromFile: " << InputState << std::endl;
 
}

void PrimaryGeneratorAction_CRY::GeneratePrimaries(G4Event * anEvent)
{	
	anEvent->SetEventID(event_number_global);
	local_event_number_CRY = anEvent->GetEventID(); //event_number; //_MCPL
	event_number_global++;

	if (InputState != 0) {
		G4String* str = new G4String("CRY library was not successfully initialized");
		std::cerr << "Error in Generate Primaries" << std::endl;
		G4Exception("PrimaryGeneratorAction", "1", RunMustBeAborted, *str);
	}

	G4double x; G4double y; G4double z;
	G4double t; G4double px; G4double py; G4double pz;
	G4double KE;

	G4AutoLock lock(&PrimaryGeneratorCRYMutex);
	G4String particleName;
	vect->clear();
	gen->genEvent(vect);

	int index_ = std::floor(run_number); ///100
	if (index_ > 5){index_=5;}

	std::vector<double> Energy_range_run = Energy_range[index_];
	boost::random::uniform_int_distribution<> KE_generator(Energy_range_run[0],Energy_range_run[1]); //std::floor(i/b)

	
	std::cout << "Worker:"<< G4Threading::G4GetThreadId() << " = = " << std::fixed << (int)local_event_number_CRY <<  " event size:" << vect->size() << std::endl; 
	//for (unsigned j = 0; j < vect->size(); j++) {std::cout << "Particle Time: " << (*vect)[j]->t() << std::endl;}

	for (unsigned  j = 0; j < vect->size(); j++) {
		

		std::vector<G4double> particle_gun_record_row;

		particleName = CRYUtils::partName((*vect)[j]->id());
		x = (*vect)[j]->x()* m; //(*vect)[j]->x()* m;
		y = 5.0*m;//5.0*m;
		z = (*vect)[j]->y()*m;// (*vect)[j]->y()* m; // (*vect)[j]->z() * m

		KE = KE_generator(rng)*MeV;//(*vect)[j]->ke()*MeV;//5000.0*MeV; // // here we need to customize the energy in order to get the desired energy
		px = (*vect)[j]->u();//(*vect)[j]->u();
		py = (*vect)[j]->w();//(*vect)[j]->v();
		pz = (*vect)[j]->v(); //(*vect)[j]->w();
		t =  (*vect)[j]->t()*s;


	    // p->position[0] = x; p->position[1] = y;p->position[2] = z;
        // p->direction[0] = px; p->direction[1] = py; p->direction[2] = pz;
        // p->ekin = KE;
        // p->time = t/ms;

        // p->pdgcode = (*vect)[j]->PDGid();
        // p->weight = 1;
        // p->userflags =  anEvent->GetEventID();
        // mcpl_add_particle(f,p);


		std::cout << 
			G4Threading::G4GetThreadId() << " = = "
			//<<std::fixed << (int)local_event_number_CRY << "  " << particleName << " ID: " << (*vect)[j]->PDGid() << " charge= " << (*vect)[j]->charge() << " "
			<< setprecision(4)
			<< " energy (MeV)=" << KE << " "
			//<< " pos (m)"
			//<< G4ThreeVector((*vect)[j]->x(), (*vect)[j]->y(), y/m)
			//<< " " << "direction cosines "
			//<< G4ThreeVector((*vect)[j]->u(), (*vect)[j]->w(), (*vect)[j]->v())
			<< " " << "Particle Time: " << (*vect)[j]->t() 
		<< std::endl;

		// fParticleGun->GeneratePrimaryVertex(anEvent);
		delete (*vect)[j];
	
	}
	
}

//....

