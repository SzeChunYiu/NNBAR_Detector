#include "Sensitive_Detector/CarbonSD.hh"

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
CarbonSD::CarbonSD(G4String name):
G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="CarbonHitCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
    
}

//.....
CarbonSD::~CarbonSD()
{}

//.....
void CarbonSD::Initialize(G4HCofThisEvent*)
{HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,collectionName[0]);}

//.....
G4bool CarbonSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{
    //Only interested in particles coming in and going out
    if (aStep->IsFirstStepInVolume() || (aStep->IsLastStepInVolume())){}
    else{return false;}

    G4Track * theTrack = aStep  ->  GetTrack();
   
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

    NNbarHit* detectorHit = new NNbarHit();

        detectorHit -> SetName(particleName);
        detectorHit -> SetTrackID(trackID);
        detectorHit -> SetParentID(parentID);
        detectorHit -> SetProcess(proc);
        detectorHit -> SetOriginVolName(origin_vol);

        detectorHit -> SetPosX(x);
        detectorHit -> SetPosY(y);
        detectorHit -> SetPosZ(z);
        detectorHit -> SetTime(time);
        detectorHit -> SetKinEn(eKinMean);
        detectorHit -> SetPX(momentum.getX());
        detectorHit -> SetPY(momentum.getY());
        detectorHit -> SetPZ(momentum.getZ());
        detectorHit -> SetEDep(energyDeposit);
        
        if (aStep->IsFirstStepInVolume() && (aStep->IsLastStepInVolume())){
            detectorHit -> SetStepInfo(2);
            HitsCollection -> insert(detectorHit);
            return true;
        }
        
        else if (aStep->IsFirstStepInVolume()){
            detectorHit -> SetStepInfo(0);
            HitsCollection -> insert(detectorHit);
            return true;
        }
        
        else if (aStep->IsLastStepInVolume()){
            detectorHit -> SetStepInfo(1);
            HitsCollection -> insert(detectorHit);
            return true;
        }

        else {return false;}    
}

//......
void CarbonSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    static G4int HCID = -1;
    if(HCID < 0){HCID = GetCollectionID(0);}
    HCE -> AddHitsCollection(HCID,HitsCollection);
}

