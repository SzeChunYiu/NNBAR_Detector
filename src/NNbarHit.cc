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

#include "NNbarHit.hh"


//**********************MT
G4ThreadLocal G4Allocator<NNbarHit>* NNbarHitAllocator=0;
//**********************MT

NNbarHit::NNbarHit()
: G4VHit()
{
    localTime = 0.;
    parentID = 0;
    process = "";
    name = "";
    time = 0.;
    trackID = 0;
    xHitID = 0;
    energyDeposit = 0.;
    kinEnergy = 0.;
    posZ = 0.;
    origin_rp = 99;
    photons = 0;
    stave_ID_ = 0;
    group_ID_ = 0;
    module_ID_=0;
}

NNbarHit::~NNbarHit()
{}

NNbarHit::NNbarHit(const NNbarHit& right)
  : G4VHit()
{
    localTime = right.localTime;
    parentID = right.parentID;
    process = right.process;
    name = right.name;
    time = right.time;
    trackID = right.trackID;
    xHitID = right.xHitID;
    origin_rp = right.origin_rp;
    stave_ID_ = right.stave_ID_;
    group_ID_ = right.group_ID_;
    module_ID_= right.module_ID_;
    vol_name = right.vol_name;
    posX = right.posX;
    posY = right.posY;
    posZ = right.posZ;
    posX_particle = right.posX_particle;
    posY_particle = right.posY_particle;
    posZ_particle = right.posZ_particle;
    px = right.px;
    py = right.py;
    pz = right.pz;



    TrackLength = right.TrackLength;
    energyDeposit = right.energyDeposit;
    kinEnergy = right.kinEnergy;
    photons = right.photons;
}

const NNbarHit& NNbarHit::operator=(const NNbarHit& right)
{
    localTime = right.localTime;
    parentID = right.parentID;
    process = right.process;
    name = right.name;
    time = right.time;
    trackID = right.trackID;
    xHitID = right.xHitID;
    origin_rp = right.origin_rp;
    stave_ID_ = right.stave_ID_;
    group_ID_ = right.group_ID_;
    module_ID_= right.module_ID_;
    posX = right.posX;
    posY = right.posY;
    posZ = right.posZ;
    px = right.px;
    py = right.py;
    pz = right.pz;
    posX_particle = right.posX_particle;
    posY_particle = right.posY_particle;
    posZ_particle = right.posZ_particle;
    TrackLength = right.TrackLength;
    vol_name = right.vol_name;
    energyDeposit = right.energyDeposit;
    kinEnergy = right.kinEnergy;
    photons = right.photons;

    return *this;
}

G4bool NNbarHit::operator==(const NNbarHit& right) const
{
return(xHitID==right.xHitID);
       	//return((xHitID==right.xHitID)&&(zHitID==right.zHitID)&&(yHitID==right.yHitID));
}
