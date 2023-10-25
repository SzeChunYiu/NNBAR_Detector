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

NNbarHit::NNbarHit(const NNbarHit& right): G4VHit()
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
    origin_vol_name = right.origin_vol_name;
    posX = right.posX;
    posY = right.posY;
    posZ = right.posZ;
    posX_local = right.posX_local;
    posY_local = right.posY_local;
    posZ_local = right.posZ_local;
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
    step_info = right.step_info;
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
    posX_local = right.posX_local;
    posY_local = right.posY_local;
    posZ_local = right.posZ_local;
    px = right.px;
    py = right.py;
    pz = right.pz;
    posX_particle = right.posX_particle;
    posY_particle = right.posY_particle;
    posZ_particle = right.posZ_particle;
    TrackLength = right.TrackLength;
    vol_name = right.vol_name;
    origin_vol_name = right.origin_vol_name;
    energyDeposit = right.energyDeposit;
    kinEnergy = right.kinEnergy;
    photons = right.photons;
    step_info = right.step_info;

    return *this;
}

G4bool NNbarHit::operator==(const NNbarHit& right) const
{
return(xHitID==right.xHitID);
       	//return((xHitID==right.xHitID)&&(zHitID==right.zHitID)&&(yHitID==right.yHitID));
}
