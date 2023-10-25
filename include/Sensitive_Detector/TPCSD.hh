#ifndef TPCSD_h
#define TPCSD_h 1

#include "NNbarHit.hh"
#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class TPCSD : public G4VSensitiveDetector{
    
    public:
    TPCSD(G4String name);
    ~TPCSD();

    std::ofstream ofs;
    void Initialize(G4HCofThisEvent*);
    
    G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
    
    void EndOfEvent(G4HCofThisEvent*HCE);
    
private:
    NNbarHitsCollection *HitsCollection;
    G4String sensitiveDetectorName;
	int error_count;
};
#endif


