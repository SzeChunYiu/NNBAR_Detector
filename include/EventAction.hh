//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"

#include "G4THitsMap.hh"
#include "globals.hh"

class EventAction : public G4UserEventAction
{
public:
  EventAction();
  virtual ~EventAction();

  virtual void  BeginOfEventAction(const G4Event* event);
  virtual void    EndOfEventAction(const G4Event* event);
    
private:
  // methods
  G4THitsMap<G4double>* GetHitsCollection(G4int hcID,
                                          const G4Event* event) const;
  G4double GetSum(G4THitsMap<G4double>* hitsMap) const;
  void PrintEventStatistics(G4double absoEdep, G4double absoTrackLength,
                            G4double gapEdep, G4double scintTrackLength) const; //, G4double gapTrackLength) const;
  
  // data members                   
  G4int  fAbsoEdepHCID;
  G4int  fGapEdepHCID;
  G4int  fAbsoTrackLengthHCID;
  G4int  fCerenkovHCID;
  G4int  fScintTrackLengthHCID;
  G4int  scintHitsCollectionID;
  G4int  absHitsCollectionID;
  G4int  tubeHitsCollectionID;
  G4int  TPCHitsCollectionID;
  G4int  PMTHitsCollectionID;
  G4int  ShieldHitsCollectionID;
  G4int  SiliconHitsCollectionID;
  
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
                     
//....

#endif

    
