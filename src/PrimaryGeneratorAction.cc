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
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"
#include "G4Threading.hh"

using namespace std;

extern std::ofstream Particle_outFile;
extern std::vector<std::vector<G4double>> particle_gun_record;
extern G4double event_number_global;
extern G4int run_number;

extern int theta_bin_index;
extern int KE_bin_index;
extern G4int particle_name_input;

#include "G4AutoLock.hh"
namespace { G4Mutex PrimaryGeneratorMutex = G4MUTEX_INITIALIZER;}
G4ThreadLocal G4int local_event_number;
G4ThreadLocal G4int initial_local_event_number = 0; // used?
G4ThreadLocal G4int updated_local_event_number; // used? 
G4ThreadLocal G4int flag; // flaggin the event with event ID = 0 !! //used?


std::vector<G4String> particle_name{"pi+","e-","gamma"};

std::vector<std::vector<G4double>> Energy_range_
{	
	{0.5,0.5},
	{2.5,2.5},
	{100.,100.},
	{400.,400.},
	{50.,100.},
	{100.,200.},
	{200.,500.},
	{500.,700.}
	// {0.00001,0.0001}, //1
	// {0.0001,0.001}, //2
	// {0.001,0.01}, //3
	// {0.01,0.1}, //4
	// {0.1,1.}, //5 
	// {1.,  10.}, //6
	// {10., 100.}, //7
	// {100., 1000.}, //8
	// {1000.,10000.}, // 9, 0.001,0.01 MeV
	// {10000.,100000.},//10, 0.01 to 0.1 MeV
	// {100000.,1000000.},//11,0.01 to 1 MeV
	// {1000000.,5000000.}, //12,1 MeV to 5 MeV
	// {5000000.,10000000.}, //13,5 MeV to 10 MeV
	// {10000000.,50000000.}, //14,1 to 50 MeV
	// {50000000.,100000000.}, //15,50 to 100 MeV
	// {100000000.,200000000.}, // 16,100 MeV to 200 MeV
	// {200000000.,300000000.}, // 17,200 MeV to 300 MeV
	// {300000000.,500000000.}, // 18,300 MeV to 500 MeV
	// {500000000.,750000000.},  // 19,500 MeV to 750 MeV
	// {750000000.,1000000000.},  // 20,750 MeV to 1000 MeV
	// {1000000000.,1500000000.},  // 21,1000 MeV to 1500 MeV
	// {1500000000.,2000000000.}  // 22,1500 MeV to 2000 MeV
}; //total 22 bins, 0-21 for index (0,8)

std::vector<std::vector<G4double>> angle_range_
{	
	{0.,0.}
	// {-68.198,-65.8},
	// {-65.8,-30.},
	// {-30.,0.},
	// {0.,30.},
	// {30.,65.8},
	// {65.8,68.198}
	// {-90.,-30.}, //1
	// {-60.,-30.},
	// {-30.,-10.}, //2
	// {-10.,-5.}, //3
	// {-5.0,-0.2}, //4
	// {-0.2,0.2}, //directly hit the target
	// {0.2, 5.0}, //6
	// {5., 10.},  //7
	// {10.,30.}, //8
	// {30.,60.},  
	// {30.,90.} //9
//9 bins in total, 0-8 for index (0,2)
};

boost::random::mt19937 rng_;

PrimaryGeneratorAction::PrimaryGeneratorAction():fParticleGun(nullptr){

	fParticleGun = new G4ParticleGun();
	fMessenger = new G4GenericMessenger(this,"/Particle_control/","Primary generator controls");
}


PrimaryGeneratorAction::~PrimaryGeneratorAction(){}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event * anEvent)
{
	G4AutoLock lock(&PrimaryGeneratorMutex);
	anEvent->SetEventID(event_number_global);
	local_event_number = anEvent->GetEventID(); //event_number; //_MCPL
	event_number_global++;

	G4double x; G4double y; G4double z;
	G4double t; G4double px; G4double py; G4double pz;
	G4double KE;

	G4String particleName;

	std::vector<G4double> particle_gun_record_row;

	// int index_ = std::floor(run_number/100);
	// if (index_ > 7){index_=7;}

	std::vector<double> Energy_range_run = Energy_range_[KE_bin_index];
	//boost::random::uniform_int_distribution<> (Energy_range_run[0]+G4UniformRand()*(Energy_range_run[1]-Energy_range_run[0])); //std::floor(i/b) //
	
	std::vector<double> angle_run = angle_range_[theta_bin_index];
	//std::uniform_real_distribution<double> angle_dis(angle_run[0], angle_run[1]);
	
	double angle = (angle_run[0]+G4UniformRand()*(angle_run[1]-angle_run[0]));
	double sign_x; double sign_y;
	
	if (G4UniformRand()>0.5){sign_x = -1.0;}
	else{sign_x = 1.0;}

	if (G4UniformRand()>0.5){sign_y = -1.0;}
	else{sign_y = 1.0;}

	x = 0.*m; //2,19 m from matthias
	y = 1.2*m; //2,19 m from matthias
	z = 0.02*m; // (*vect)[j]->z() * m
	KE = (Energy_range_run[0]+G4UniformRand()*(Energy_range_run[1]-Energy_range_run[0]))*MeV;//KE_generator(rng_)*eV;///(1.0+ 50.0*std::floor(run_number/10)) *MeV ; //250.0 * MeV;
	px = 0.;
	py = 1.; //sin(angle/360.0*2.0*3.1415926);//angle);
	pz = 0.; //cos(angle/360.0*2.0*3.1415926);//angle/(360.0*2.0*3.1415926));
	double dt = 68.5;
	t = 0.0;//(std::floor(run_number/3.0)*(2.86+dt) + G4UniformRand()*2.86)*millisecond; 

	fParticleGun->SetParticleDefinition(particleTable->FindParticle(particle_name[particle_name_input]));
	
	//fParticleGun->SetParticleEnergy(KE);
	
	int n = event_number_global/100; // every energy, it has 50 entries
	fParticleGun->SetParticleEnergy(KE);	
	
	fParticleGun->SetParticlePosition(G4ThreeVector(x, y, z));
	fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px, py, pz)); //(*vect)[j]->w())
	fParticleGun->SetParticleTime(t);

	std::cout << particle_name[particle_name_input] << " :: "<< angle_run[0] << ","<< angle_run[1] << " :: " << angle << std::endl; 
	std::cout << "time" << t/millisecond << " ms :: " <<"KE " << KE << " angle " << angle << " py: " << py << " pz: " << pz<< std::endl;
	
	Particle_outFile <<  anEvent->GetEventID() << ",";
	Particle_outFile <<  event_number_global << ",";
	Particle_outFile <<  particleTable -> FindParticle(particle_name[particle_name_input])->GetPDGEncoding()  << ","; //PID
	Particle_outFile <<  particleTable -> FindParticle(particle_name[particle_name_input])-> GetPDGMass() << ","; //Mass
	Particle_outFile <<  0 << ","; //PID
	Particle_outFile <<  KE << ","; // 100+n*10 
	Particle_outFile <<  angle << ","; 
	Particle_outFile <<  x/10.0 << ",";
	Particle_outFile <<  y/10.0 << ",";
	Particle_outFile <<  z/10.0 << ",";
	Particle_outFile <<  std::setprecision(13) << t << std::setprecision(4) << ",";
	Particle_outFile <<  px << ",";
	Particle_outFile <<  py << ",";
	Particle_outFile <<  pz <<",";
	Particle_outFile <<  0 << G4endl;

	fParticleGun->GeneratePrimaryVertex(anEvent);
	

}

//....

