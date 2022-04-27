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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#ifndef NNbarHit_h
#define NNbarHit_h 1


#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"


class NNbarHit : public G4VHit
{
public:
    NNbarHit();
    NNbarHit(const NNbarHit&);
    virtual ~NNbarHit();
 
 
    const NNbarHit& operator=(const NNbarHit&);
 
    G4bool operator==(const NNbarHit&) const;

//******************************MT
inline void* operator new(size_t);
inline void operator delete(void*);
//******************************MT

private:
    G4double localTime;
    G4int parentID;
    G4String process;
    G4String name;
    G4double time;
    G4int trackID; 
    G4int xHitID; // Hit x voxel 
    G4double posX;
    G4double posY;
    G4double posZ;
    G4double px;
    G4double py;
    G4double pz;
    G4double posX_particle;
    G4double posY_particle;
    G4double posZ_particle;
    G4int stave_ID_;
    G4int group_ID_;
    G4int module_ID_;
    G4int origin_rp;
    G4double energyDeposit;
    G4double kinEnergy;
    G4int photons;
    G4double TrackLength;
    G4String vol_name;
    
    
public:
    inline G4double GetLocalTime(){return localTime;}
    inline G4int GetParentID(){return parentID;}
    inline G4String GetProcess(){return process;}
    inline G4String GetName(){return name;}
    inline G4double GetTime(){return time;}
    inline G4int GetTrackID(){return trackID;}

    inline G4int GetStave_ID(){return stave_ID_;}
    inline G4int GetXID(){return xHitID;}
    inline G4int GetGroup_ID(){return group_ID_;}
    inline G4int GetMod_ID(){return module_ID_;}
    
    inline G4double GetPosX(){return posX;}
    inline G4double GetPosY(){return posY;}
    inline G4double GetPosZ(){return posZ;}
    inline G4double GetPX(){return px;}
    inline G4double GetPY(){return py;}
    inline G4double GetPZ(){return pz;}

    inline G4double GetPosX_particle(){return posX_particle;}
    inline G4double GetPosY_particle(){return posY_particle;}
    inline G4double GetPosZ_particle(){return posZ_particle;}

    inline G4double GetTrackLength(){return TrackLength;}
    inline G4double GetEdep(){return energyDeposit;}
    inline G4double GetKinEn(){return kinEnergy;}
    inline G4double GetOrigin(){return origin_rp;} // get which layer this hit particle is from
    inline G4int GetPhotons(){return photons;}
    inline G4String GetVolName(){return vol_name;}

    inline void SetLocalTime(G4double ltime){localTime = ltime;}
    inline void SetParentID(G4int parent){parentID = parent;}
    inline void SetProcess(G4String p){process = p;}
    inline void SetName(G4String n){name = n;}
    inline void SetTime(G4double t){time = t;}
    inline void SetTrackID(G4int track){trackID = track;}

    inline void SetStave_ID(G4int stave_ID){stave_ID_ = stave_ID;}
    inline void SetXID(G4int xID){xHitID = xID;}
    inline void SetGroup_ID(G4int groupID){group_ID_=groupID;}
    inline void SetMod_ID(G4int ModID){module_ID_=ModID;}
    
    inline void SetPosX(G4double x){posX = x;}
    inline void SetPosY(G4double y){posY = y;}
    inline void SetPosZ(G4double z){posZ = z;}

    inline void SetPX(G4double pX){px = pX;}
    inline void SetPY(G4double pY){py = pY;}
    inline void SetPZ(G4double pZ){pz = pZ;}
    
    inline void SetPosX_particle(G4double x_particle){posX_particle = x_particle;}
    inline void SetPosY_particle(G4double y_particle){posY_particle = y_particle;}
    inline void SetPosZ_particle(G4double z_particle){posZ_particle = z_particle;}

    inline void SetTrackLength(G4double tracklength) {TrackLength = tracklength;}
    inline void SetEDep(G4double eDep){energyDeposit = eDep;}
    inline void SetKinEn(G4double kinEn){kinEnergy = kinEn;}
    inline void SetOrigin(G4int origin) {origin_rp = origin;}
    inline void SetPhotons(G4int photon) {photons = photon;}
    inline void SetVolName(G4String name){vol_name=name;}

};

typedef G4THitsCollection<NNbarHit> NNbarHitsCollection;
//******************************MT
extern G4ThreadLocal G4Allocator<NNbarHit>* NNbarHitAllocator;
//******************************MT

inline void* NNbarHit::operator new(size_t){ 
  
    if(!NNbarHitAllocator) {
        NNbarHitAllocator= new G4Allocator<NNbarHit>;
    }

    void *aHit;

    aHit = (void *) NNbarHitAllocator->MallocSingle();
    return aHit;

}

inline void NNbarHit::operator delete(void *aHit) {
    if(!NNbarHitAllocator){ 
        NNbarHitAllocator= new G4Allocator<NNbarHit>;
    }

    NNbarHitAllocator->FreeSingle((NNbarHit*) aHit);
}

#endif
