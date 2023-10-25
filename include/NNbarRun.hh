#include "G4Run.hh"

class NNbarRun : public G4Run
{
	public:
		NNbarRun();
		virtual ~NNbarRun();
		virtual void RecordEvent(const G4Event*);
		virtual void Merge(const G4Run*);
	private:
  // data members                   
  G4int  fAbsoEdepHCID;
  G4int  fGapEdepHCID;
  G4int  fAbsoTrackLengthHCID;
  G4int  fCerenkovHCID;
  G4int  fScintTrackLengthHCID;

  G4int  SiliconHitsCollectionID;

  G4int  CarbonHitsCollectionID;
  G4int  BeampipeHitsCollectionID;
  G4int  TPCHitsCollectionID;
  
  G4int  scintHitsCollectionID;
  G4int  LeadGlassHitsCollectionID;
  G4int  PMTHitsCollectionID;
  G4int  ShieldHitsCollectionID;
  
  G4int b = 1;
  G4int ltime     = 0.;
  G4int parentID  = 0;
  G4String proc   = "";
  G4String name   = "";
  G4double time   = 0.;
  G4int trID      = 0;
  G4int i         = 0;
  G4double kinEn  = 0.;
  G4double eDep   = 0.;
  G4double trackl = 0.;	
  G4int hitCount  = 0;
  G4int org_replica = 99;
  G4int group_ID = 999;
  G4int module_ID = 999;
  G4double x=0;
  G4double y=0;
  G4double z=0;
  
  // Book vector to keep track of Edep in each Scintillator Sheet
  G4double EdepPerSheet[10] = {0., 0., 0., 0., 0.,0., 0., 0., 0., 0.};
  G4int ScintPerSheet[10] = { 0,0,0,0,0,0,0,0,0,0 };
  G4int ScintPerSheet_new[10] = { 0,0,0,0,0,0,0,0,0,0 };
  G4double totEdep   = 0.;     
  G4double eDepScint = 0.;
  G4double eDepAbs   = 0.;
  G4double eDepTube  = 0.; 
  G4double extraEdep = 0.;
  G4double eDepCompt = 0.;
  G4double eDepInelastic= 0.;
  G4double eDephIoni = 0.;
  G4double eDepHadElas = 0.;
  G4double eDepPrimary = 0.;
  G4double eDepOther = 0.;
  G4int cerenkovCounter = 0;
  G4int scint_photons = 0;
  G4int scint_photons_check = 0;
  G4int PMT_photons = 0;
  
};
