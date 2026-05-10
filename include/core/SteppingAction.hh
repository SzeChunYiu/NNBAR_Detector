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

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <map>
#include <string>


class SteppingAction : public G4UserSteppingAction
{
  public:
    SteppingAction();
    virtual ~SteppingAction();

    virtual void UserSteppingAction(const G4Step*);

    // Energy deposition tracking per detector
    void ResetEdep();
    G4double GetEdepTPC() const { return fEdepTPC; }
    G4double GetEdepScintillator() const { return fEdepScintillator; }
    G4double GetEdepLeadGlass() const { return fEdepLeadGlass; }
    G4double GetEdepOther() const { return fEdepOther; }
    G4int GetPhotonsTPC() const { return fPhotonsTPC; }
    G4int GetPhotonsScintillator() const { return fPhotonsScintillator; }
    G4int GetPhotonsLeadGlass() const { return fPhotonsLeadGlass; }
    G4int GetScintPhotonCount() const { return fScintillationCounter; }
    G4int GetCerenkovPhotonCount() const { return fCerenkovCounter; }

  private:
    G4int fScintillationCounter;
    G4int fCerenkovCounter;
    G4int fEventNumber;

    // Energy deposition accumulators
    G4double fEdepTPC;
    G4double fEdepScintillator;
    G4double fEdepLeadGlass;
    G4double fEdepOther;

    // Optical photon counters per detector
    G4int fPhotonsTPC;
    G4int fPhotonsScintillator;
    G4int fPhotonsLeadGlass;
};

//....

#endif
