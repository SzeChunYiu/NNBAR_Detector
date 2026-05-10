#include "core/ActionInitialization.hh"
#include "config.h"

#include "core/PrimaryGeneratorAction.hh"

#include "generator/G4MCPLGenerator.hh"
#include "core/RunAction.hh"
#include "core/EventAction.hh"
#include "core/SteppingAction.hh"

#if WITH_GARFIELD_GPU
#include "physics/TPCDriftManager.hh"
#endif

//....

ActionInitialization::ActionInitialization()
 : G4VUserActionInitialization()
{}

//....

ActionInitialization::~ActionInitialization()
{;}

//....

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
}

//....

// Global flag to force particle gun mode (set via pre-init messenger command)
// This allows runtime switching even when compiled with MCPL_BUILD=1
G4bool g_useParticleGun = false;

void ActionInitialization::Build() const
{
  // Create single PrimaryGeneratorAction instance per thread
  // Calibration settings are shared via static members
  #if MCPL_BUILD==1
    // Check if particle gun mode is forced via /generator/use_particle_gun command
    if (g_useParticleGun) {
      std::cout << " Particle Gun mode is activated (forced)" << std::endl;
      SetUserAction(new PrimaryGeneratorAction());
    } else {
      std::cout << " MCPL generator is activated" << std::endl;
      SetUserAction(new G4MCPLGenerator("./mcpl_files/NNBAR_rwag_signal_GBL_jbar_100k_9009.mcpl"));
    }
  #else
    SetUserAction(new PrimaryGeneratorAction());
  #endif


  SetUserAction(new RunAction);

  // Create SteppingAction and pass to EventAction for energy deposition tracking
  SteppingAction* steppingAction = new SteppingAction();
  SetUserAction(steppingAction);

  EventAction* eventAction = new EventAction(steppingAction);
  SetUserAction(eventAction);

#if WITH_GARFIELD_GPU
  // Initialize TPC drift manager with Ar/CO2 80/20 gas properties
  // TPC Geometry from TPC_geometry.cc:
  //   - Beampipe inner radius: 1120 mm (112 cm)
  //   - TPC drift length: 850 mm (85 cm)
  //   - TPC z-extent: ~2504 mm (250 cm half-length)
  //   - 12 rectangular modules arranged around beampipe
  //   - Drift field ~250 V/cm along radial direction (x or y depending on module)
  auto* driftMgr = nnbar::TPCDriftManager::Instance();
  driftMgr->Initialize();
  // Set TPC geometry (in cm): inner radius=112, outer~197 (112+85), halfZ=250, drift=85
  driftMgr->SetTPCGeometry(112.0f, 197.0f, 250.0f, 85.0f);
  // Set Ar/CO2 80/20 gas at 760 Torr (1 atm), 293.15 K (20°C)
  driftMgr->SetGasProperties(0.80f, 0.20f, 760.0f, 293.15f);
#endif

}  

//....
