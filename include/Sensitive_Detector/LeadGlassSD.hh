#ifndef LeadGlassSD_h
#define LeadGlassSD_h 1

#include "NNbarHit.hh"
//#include "Sensitive_Detector/LeadGlassHit.hh"


#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class LeadGlassSD : public G4VSensitiveDetector
{
public:
    LeadGlassSD(G4String name);
    ~LeadGlassSD();
    
    
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


