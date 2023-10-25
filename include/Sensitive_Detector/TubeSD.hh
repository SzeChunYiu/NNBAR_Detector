#ifndef TubeSD_h
#define TubeSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

#include "NNbarHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class TubeSD : public G4VSensitiveDetector
{
public:
    TubeSD(G4String name);
    ~TubeSD();
    
    
    std::ofstream ofs;
    void Initialize(G4HCofThisEvent*);
    
    G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
    void EndOfEvent(G4HCofThisEvent*HCE);
    
private:
    NNbarHitsCollection *HitsCollection;
    G4String sensitiveDetectorName;
};
#endif


