#include "sensitive/ScintillatorSD.hh"

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
ScintillatorSD::ScintillatorSD(G4String name):
G4VSensitiveDetector(name)
{
    G4String HCname;
    collectionName.insert(HCname="ScintillatorHitCollection");
    HitsCollection = NULL;
    sensitiveDetectorName = name;
    
}

//.....
ScintillatorSD::~ScintillatorSD()
{}

//.....
void ScintillatorSD::Initialize(G4HCofThisEvent*)
{
    
    HitsCollection = new NNbarHitsCollection(sensitiveDetectorName,
                                                             collectionName[0]);
}

//.....
G4bool ScintillatorSD::ProcessHits(G4Step* aStep, G4TouchableHistory* )
{
    // Safety check for Celeritas GPU callbacks - some pointers may be null
    G4StepPoint* PreStep = aStep->GetPreStepPoint();
    if (!PreStep) return false;

    // Get Energy deposited - always available
    G4double energyDeposit = aStep->GetTotalEnergyDeposit();
    if (energyDeposit <= 0) return false;

    // Position from pre-step (always available)
    G4ThreeVector pos = PreStep->GetPosition();
    G4double x = pos.getX();
    G4double y = pos.getY();
    G4double z = pos.getZ();

    // Get Track - may be null when called from Celeritas with track=false
    G4Track* theTrack = aStep->GetTrack();

    // Variables with defaults for when track is unavailable (Celeritas GPU mode)
    G4String particleName = "unknown";
    G4String name = "unknown";
    G4int trackID = -1;
    G4int parentID = 0;
    G4String proc = "GPU";
    G4double time = 0;
    G4double localTime = 0;
    G4double eKinPost = 0;
    G4String origin_vol = "GPU";
    G4double tracklength = 0;
    G4int stave_ID = 0, layer_ID = 0, module_ID = 0;
    G4String current_vol = "unknown";
    G4ThreeVector scint_pos(x, y, z);
    G4double local_x = 0, local_y = 0, local_z = 0;

    // If track is available (normal Geant4 mode), get full info
    if (theTrack) {
        G4ParticleDefinition* particleDef = theTrack->GetDefinition();
        if (particleDef) {
            particleName = particleDef->GetParticleName();
            if (particleName == "opticalphoton") return false;
        }

        trackID = theTrack->GetTrackID();
        time = theTrack->GetGlobalTime();
        localTime = theTrack->GetLocalTime();

        G4ThreeVector vertex = theTrack->GetVertexPosition();
        tracklength = z - vertex.getZ();

        G4DynamicParticle* dynPart = (G4DynamicParticle*)theTrack->GetDynamicParticle();
        if (dynPart && dynPart->GetParticleDefinition()) {
            name = dynPart->GetParticleDefinition()->GetParticleName();
        }

        // Get Process and parentID safely
        if (trackID > 1) {
            G4VTouchable* originTouch = (G4VTouchable*)theTrack->GetOriginTouchable();
            if (originTouch) {
                G4VPhysicalVolume* originVol = originTouch->GetVolume();
                if (originVol && originVol->GetName() != "World") {
                    parentID = theTrack->GetParentID();
                    if (parentID > 0) {
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
    const G4VTouchable* touchable = PreStep->GetTouchable();
    if (touchable) {
        stave_ID = touchable->GetReplicaNumber(0);
        layer_ID = touchable->GetReplicaNumber(1);
        module_ID = touchable->GetReplicaNumber(2);
        G4VPhysicalVolume* vol = touchable->GetVolume();
        if (vol) current_vol = vol->GetName();
        scint_pos = touchable->GetTranslation(0);

        G4NavigationHistory* history = (G4NavigationHistory*)touchable->GetHistory();
        if (history) {
            G4ThreeVector localPosition = history->GetTopTransform().TransformPoint(pos);
            local_x = localPosition.getX();
            local_y = localPosition.getY();
            local_z = localPosition.getZ();
        }
    }

    G4int photons = static_cast<G4int>(energyDeposit * 11136.);

    NNbarHit* detectorHit = new NNbarHit();
    detectorHit->SetLocalTime(localTime);
    detectorHit->SetParentID(parentID);
    detectorHit->SetProcess(proc);
    detectorHit->SetTime(time);
    detectorHit->SetName(name);
    detectorHit->SetTrackID(trackID);
    detectorHit->SetStave_ID(stave_ID);
    detectorHit->SetXID(layer_ID);
    detectorHit->SetMod_ID(module_ID);
    detectorHit->SetPosX(scint_pos.getX());
    detectorHit->SetPosY(scint_pos.getY());
    detectorHit->SetPosZ(scint_pos.getZ());
    detectorHit->SetLocalPosX(local_x);
    detectorHit->SetLocalPosY(local_y);
    detectorHit->SetLocalPosZ(local_z);
    detectorHit->SetPosX_particle(x);
    detectorHit->SetPosY_particle(y);
    detectorHit->SetPosZ_particle(z);
    detectorHit->SetVolName(current_vol);
    detectorHit->SetOriginVolName(origin_vol);
    if (aStep->IsFirstStepInVolume()) { detectorHit->SetStepInfo(1); }
    else { detectorHit->SetStepInfo(0); }
    detectorHit->SetTrackLength(tracklength);
    detectorHit->SetEDep(energyDeposit);
    detectorHit->SetKinEn(eKinPost);
    detectorHit->SetPhotons(photons);
    HitsCollection->insert(detectorHit);

    return true;
}

//......
void ScintillatorSD::EndOfEvent(G4HCofThisEvent* HCE)
{
    
    static G4int HCID = -1;
    if(HCID < 0)
    {
        HCID = GetCollectionID(0);
    }
    
    HCE -> AddHitsCollection(HCID,HitsCollection);
}

