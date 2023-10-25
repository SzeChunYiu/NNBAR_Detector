#include "Sensitive_Detector/TPCSD.hh"

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
#include <random>
//.....
TPCSD::TPCSD(G4String name):G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="TPCHitCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
}

TPCSD::~TPCSD(){}

void TPCSD::Initialize(G4HCofThisEvent*)
{HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,collectionName[0]);}

G4bool TPCSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{

//Only interested in particles coming in and going out
    // [x,y,z,t,u,v,w,KE,step_length,e-]


    if (aStep->IsFirstStepInVolume() || (aStep->IsLastStepInVolume())){}
    else{return false;}

    G4Track * theTrack = aStep -> GetTrack();
   
    //Get particle name
    G4ParticleDefinition *particleDef = theTrack -> GetDefinition();
    G4String particleName =  particleDef -> GetParticleName();
    G4int trackID = theTrack -> GetTrackID();
    G4int parentID = theTrack->GetParentID();
    G4String proc = "primary"; 
    if (parentID > 0){proc = theTrack->GetCreatorProcess()->GetProcessName();} 

    //eDep of the step
    G4double energyDeposit = aStep -> GetTotalEnergyDeposit();

    // Position
    G4StepPoint* PreStep = aStep->GetPreStepPoint();
    G4StepPoint* PostStep = aStep->GetPostStepPoint();
    G4ThreeVector pos1 = PreStep->GetPosition();
    G4ThreeVector pos2 = PostStep->GetPosition();
    G4double x = ((pos1+pos2)/2.).getX();
    G4double y = ((pos1+pos2)/2.).getY();
    G4double z = ((pos1+pos2)/2.).getZ();
    G4double DX = aStep -> GetStepLength();

    // Get Time
    G4double time = theTrack->GetGlobalTime() / CLHEP::ns;

    // Get the pre-step kinetic energy
    G4double eKinPre = aStep -> GetPreStepPoint() -> GetKineticEnergy();
    // Get the post-step kinetic energy
    G4double eKinPost = aStep -> GetPostStepPoint() -> GetKineticEnergy();
    // Get the step average kinetic energy
    G4double eKinMean = (eKinPre + eKinPost) * 0.5;

    G4ThreeVector momentum = aStep -> GetPreStepPoint() ->GetMomentumDirection ();
    G4String origin_vol = theTrack->GetOriginTouchable()->GetVolume()->GetName();
    G4String current_vol = aStep->GetPreStepPoint()->GetTouchable()->GetVolume()->GetName();

    // Read TPC volume index
    const G4VTouchable* touchable = aStep->GetPreStepPoint()->GetTouchable();
    G4int k = touchable ->GetReplicaNumber(0); // which layer it is in
    G4int TPC_index = touchable -> GetReplicaNumber(1); // which TPC module it is in


    G4int electrons = 0;
    G4ParticleDefinition* particle;
    if (particleName != "opticalphoton") {
        
        std::default_random_engine generator(std::random_device{}());
        G4double pairs = energyDeposit/(23.6*eV);
        std::poisson_distribution<int> distribution(pairs);
        electrons = distribution(generator);

        //std::cout << "TPC e- " << electrons << "paris: " << pairs << " energyDeposit :" << energyDeposit << std::endl;
        // const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();
        // for (int j = 0; j < (*secondary).size(); j++) {
        //     particle = (*secondary)[j]->GetDefinition();
        //     if (particle->GetParticleName() == "e-") {
        //         auto process_name =  (*secondary)[j]->GetCreatorProcess()->GetProcessName ();
        //         if (process_name== "eIoni" || process_name!= "hIoni" || process_name!= "muIoni"){electrons++;}
        //     }
        // }

    }

    else
    {return false;}
    
    // [x,y,z,t,u,v,w,KE,step_length,e-,TPC layer,TPC module]

    NNbarHit* detectorHit = new NNbarHit();

    detectorHit -> SetName(particleName);
    detectorHit -> SetTrackID(trackID);
    detectorHit -> SetParentID(parentID);
    detectorHit -> SetProcess(proc);
    detectorHit -> SetOriginVolName(origin_vol);
    detectorHit -> SetVolName(current_vol);
    detectorHit -> SetPosX(x);
    detectorHit -> SetPosY(y);
    detectorHit -> SetPosZ(z);
    detectorHit -> SetTime(time);
    detectorHit -> SetKinEn(eKinMean);
    detectorHit -> SetPX(momentum.getX());
    detectorHit -> SetPY(momentum.getY());
    detectorHit -> SetPZ(momentum.getZ());
    detectorHit -> SetEDep(energyDeposit);

    detectorHit -> SetXID(k); // TPC layer
    detectorHit -> SetMod_ID(TPC_index); // TPC module
    
    //remarks!! here it is named set photons but actually it counts the electrons in the TPC, too lazy to add one more function
    detectorHit -> SetPhotons(electrons);
    detectorHit -> SetTrackLength(DX); // not actually trackLength but lets stay in this way..

    // check whether this is the first step in the volume
    if (origin_vol!="TPC_1_layer_PV" && origin_vol!="TPC_2_layer_PV"){ //make sure that it is from outside
        if (aStep->IsFirstStepInVolume()){detectorHit -> SetStepInfo(1);} // 1 means it is first step
        else{detectorHit -> SetStepInfo(0);} // not the first step
    }
    else{detectorHit -> SetStepInfo(999);} //For those coming from the inside, we make a special mark
    HitsCollection -> insert(detectorHit);
    return true;
}

//......
void TPCSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    static G4int HCID = -1;
    if(HCID < 0){HCID = GetCollectionID(0);}
    HCE -> AddHitsCollection(HCID,HitsCollection);
}

