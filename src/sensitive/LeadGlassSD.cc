#include "sensitive/LeadGlassSD.hh"

#include "hits/NNbarHit.hh"

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

#include "G4Cerenkov.hh"


//......
LeadGlassSD::LeadGlassSD(G4String name):
G4VSensitiveDetector(name)
{
	error_count = 0;
    G4String HCname;
    collectionName.insert(HCname="LeadGlassHitCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
}

//......
LeadGlassSD::~LeadGlassSD()
{}

//......
void LeadGlassSD::Initialize(G4HCofThisEvent*)
{
    
    HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,collectionName[0]);
}

//.....
G4bool LeadGlassSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{
    // Safety check for Celeritas GPU callbacks - some pointers may be null
    G4StepPoint* preStep = aStep->GetPreStepPoint();
    if (!preStep) return false;

    G4VPhysicalVolume* preVol = preStep->GetPhysicalVolume();
    if (!preVol || preVol->GetName() != "LeadGlassPV") return false;

    // Get Track - may be null when called from Celeritas with track=false
    G4Track* theTrack = aStep->GetTrack();

    // Get Energy deposited - always available
    G4double energyDeposit = aStep->GetTotalEnergyDeposit();
    if (energyDeposit <= 0) return false;  // Skip zero-energy deposits

    // Get Step Length
    G4double DX = aStep->GetStepLength();

    // Variables with defaults for when track is unavailable (Celeritas GPU mode)
    G4String particleName = "unknown";
    G4String name = "unknown";
    G4int trackID = -1;
    G4int parentID = 0;
    G4String proc = "GPU";
    G4double time = 0;
    G4double eKinPost = 0;
    G4int photons = 0;
    G4String origin_vol = "GPU";
    G4int k = 0;
    G4ThreeVector lead_pos(0, 0, 0);

    // If track is available (normal Geant4 mode), get full info
    if (theTrack) {
        G4ParticleDefinition* particleDef = theTrack->GetDefinition();
        if (particleDef) {
            particleName = particleDef->GetParticleName();
            if (particleName == "opticalphoton") return false;
        }

        trackID = theTrack->GetTrackID();
        time = theTrack->GetGlobalTime() / CLHEP::ns;

        G4DynamicParticle* dynPart = (G4DynamicParticle*)theTrack->GetDynamicParticle();
        if (dynPart && dynPart->GetParticleDefinition()) {
            name = dynPart->GetParticleDefinition()->GetParticleName();
        }

        // Get optical photons (only if track available)
        if (particleName != "opticalphoton") {
            const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();
            if (secondary) {
                for (size_t j = 0; j < secondary->size(); j++) {
                    G4ParticleDefinition* particle = (*secondary)[j]->GetDefinition();
                    if (particle && particle->GetParticleName() == "opticalphoton") {
                        const G4VProcess* creatorProc = (*secondary)[j]->GetCreatorProcess();
                        if (creatorProc && creatorProc->GetProcessName() == "Cerenkov") {
                            photons++;
                        }
                    }
                }
            }
        }

        // Get Process and parentID safely
        if (trackID > 1) {
            G4VTouchable* originTouch = (G4VTouchable*)theTrack->GetOriginTouchable();
            if (originTouch) {
                G4VPhysicalVolume* originVol = originTouch->GetVolume();
                if (originVol && originVol->GetName() != "World") {
                    parentID = theTrack->GetParentID();
                    if (parentID != 0) {
                        const G4VProcess* creatorProc = theTrack->GetCreatorProcess();
                        if (creatorProc) proc = creatorProc->GetProcessName();
                        else proc = "unknown";
                    } else {
                        proc = "primary";
                    }
                    origin_vol = originVol->GetName();
                }
            }
        }

        // Get post-step kinetic energy
        G4StepPoint* postStep = aStep->GetPostStepPoint();
        if (postStep) {
            eKinPost = postStep->GetKineticEnergy();
        }
    }

    // Get touchable info if available (may be null in Celeritas mode)
    const G4VTouchable* touchable = preStep->GetTouchable();
    if (touchable) {
        k = touchable->GetReplicaNumber(1);
        lead_pos = touchable->GetTranslation(0);
    } else {
        // Use pre-step position as fallback
        lead_pos = preStep->GetPosition();
    }

    // Create and fill hit
    NNbarHit* detectorHit = new NNbarHit();
    detectorHit->SetParentID(parentID);
    detectorHit->SetProcess(proc);
    detectorHit->SetTime(time);
    detectorHit->SetName(name);
    detectorHit->SetTrackID(trackID);
    detectorHit->SetXID(k);
    detectorHit->SetPosX(lead_pos.getX());
    detectorHit->SetPosY(lead_pos.getY());
    detectorHit->SetPosZ(lead_pos.getZ());
    detectorHit->SetTrackLength(DX);
    detectorHit->SetEDep(energyDeposit);
    detectorHit->SetKinEn(eKinPost);
    detectorHit->SetPhotons(photons);
    detectorHit->SetOriginVolName(origin_vol);
    if (aStep->IsFirstStepInVolume()) { detectorHit->SetStepInfo(1); }
    else { detectorHit->SetStepInfo(0); }
    HitsCollection->insert(detectorHit);
    return true;
}

void LeadGlassSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    static G4int HCID = -1;
    if(HCID < 0)
    {HCID = GetCollectionID(0);}
    HCE -> AddHitsCollection(HCID,HitsCollection);
}

