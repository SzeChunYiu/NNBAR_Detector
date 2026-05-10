#include "sensitive/ShieldSD.hh"

#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "G4ios.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4TrackVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4UserEventAction.hh"
#include "G4TransportationManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Scintillation.hh"
#include "G4Cerenkov.hh"

//.....
ShieldSD::ShieldSD(G4String name):
G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="ShieldHitCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
    
}

//.....
ShieldSD::~ShieldSD()
{}

//.....
void ShieldSD::Initialize(G4HCofThisEvent*)
{
    
    HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,collectionName[0]);
}

//.....
G4bool ShieldSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{

    //std::cout<< " shiled hit class SD " << std::endl;

    //if (aStep -> GetPreStepPoint() -> GetPhysicalVolume() -> GetName() != "Layer") return false;
    
    // Get Direction
    G4Track * theTrack = aStep  ->  GetTrack();
    
    // Get unique track_id (in an event)
    G4int trackID = theTrack -> GetTrackID();
   
    // Get Energy deposited
    G4double energyDeposit = aStep -> GetTotalEnergyDeposit();

    G4StepPoint* PreStep = aStep->GetPreStepPoint();
    G4StepPoint* PostStep = aStep->GetPostStepPoint();

    // Position
    G4ThreeVector pos1 = PreStep->GetPosition();
    G4ThreeVector pos2 = PostStep->GetPosition();
    G4double x = ((pos1+pos2)/2.).getX();
    G4double y = ((pos1+pos2)/2.).getY();
    G4double z = ((pos1+pos2)/2.).getZ();

    G4ThreeVector momentum = PreStep->GetMomentumDirection();
    G4double px = momentum.getX();
    G4double py = momentum.getY();
    G4double pz = momentum.getZ();

    G4ThreeVector vertex = theTrack->GetVertexPosition();
    G4double origin = vertex.getZ();
    G4double tracklength = z - origin;


    // Get Time
    G4double time = theTrack->GetGlobalTime() / CLHEP::ns;

    // Get Local Time
    G4double localTime = theTrack->GetLocalTime() / CLHEP::ns;

    // Get Name
    G4String name = theTrack->GetDynamicParticle()->GetParticleDefinition()->GetParticleName();

    G4TouchableHandle touchPreStep = PreStep->GetTouchableHandle();
    G4VPhysicalVolume* volumePre = touchPreStep->GetVolume();
    G4String namePre = volumePre->GetName();
   
    G4int parentID = 0;
    G4String proc = ""; 

	if (trackID > 1 && theTrack->GetOriginTouchable()->GetVolume()->GetName() != "World") {
		//std::cout << name << " ID : " << trackID << " step 1 " << theTrack->GetOriginTouchable()->GetVolume()->GetName() << std::endl;
		parentID = theTrack->GetParentID();
		if (parentID > 0) { proc = theTrack->GetCreatorProcess()->GetProcessName();}
		else { proc = "primary"; }
		//std::cout << " Scint hit : " << name << " ID : " << trackID  << theTrack->GetOriginTouchable()->GetVolume()->GetName() << "  " << proc << " the parent ID is : " << parentID << std::endl;
	}
	
	else {
        proc = "primary";
		parentID = 0;
    }

    // Get the post-step kinetic energy
    G4double eKinPost = aStep -> GetPostStepPoint() -> GetKineticEnergy();

    NNbarHit* detectorHit = new NNbarHit();

    // Make this kinetic energy and position
    detectorHit -> SetLocalTime(localTime);
    detectorHit -> SetParentID(parentID);
    detectorHit -> SetProcess(proc);
    detectorHit -> SetTime(time);
    detectorHit -> SetName(name);
    detectorHit -> SetVolName(theTrack -> GetVolume()-> GetName());
    detectorHit -> SetTrackID(trackID);
    //detectorHit -> SetXID(k);
    detectorHit -> SetPosZ(tracklength);
    detectorHit -> SetEDep(energyDeposit);
    detectorHit -> SetKinEn(eKinPost);
    detectorHit -> SetPosX(x);
    detectorHit -> SetPosY(y);
    detectorHit -> SetPosZ(z);
    detectorHit -> SetPX(px);
    detectorHit -> SetPY(py);
    detectorHit -> SetPZ(pz);
    HitsCollection -> insert(detectorHit);
    
    //}
    return true;
}

//......
void ShieldSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    
    static G4int HCID = -1;
    if(HCID < 0)
    {
        HCID = GetCollectionID(0);
    }
    
    HCE -> AddHitsCollection(HCID,HitsCollection);
}

