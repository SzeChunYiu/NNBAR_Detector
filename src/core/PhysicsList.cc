#include <iomanip>
#include "config.h"

#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4RadioactiveDecayPhysics.hh"

#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4StepLimiterPhysics.hh"

#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "core/PhysicsList.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

#include "config.h"

// Celeritas GPU acceleration
#if WITH_CELERITAS
#include "physics/CeleritasInterface.hh"
#endif


PhysicsList::PhysicsList(): G4VModularPhysicsList()
{
  defaultCutValue = 0.7*CLHEP::mm;
  fConfig = G4LossTableManager::Instance()->EmConfigurator();
  defaultCutValue = 1.0*mm;
  SetVerboseLevel(0);

  ConstructParticle();

 // EM Physics
#if WITH_CELERITAS
  if (nnbar::CeleritasInterface::IsEnabled()) {
    // Use G4EmStandardPhysics for Celeritas compatibility
    // G4EmStandardPhysics_option4 has processes that Celeritas can't handle
    RegisterPhysics(new G4EmStandardPhysics());
  } else {
    // Celeritas disabled - use higher accuracy Geant4 physics
    RegisterPhysics(new G4EmStandardPhysics_option4());
  }
#else
  RegisterPhysics(new G4EmStandardPhysics_option4());
#endif
  //RegisterPhysics(new G4EmExtraPhysics());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4HadronElasticPhysics());
  RegisterPhysics(new G4HadronPhysicsFTFP_BERT()); //_HP will slow down a lot 
  RegisterPhysics(new G4StoppingPhysics() );
  RegisterPhysics(new G4IonPhysics());
  RegisterPhysics(new G4NeutronTrackingCut());
  RegisterPhysics(new G4RadioactiveDecayPhysics());
  // Enable user limits (step limits, min kinetic energy, etc.)
  RegisterPhysics(new G4StepLimiterPhysics());
  G4cout << "PhysicsList: G4StepLimiterPhysics enabled for user limits" << G4endl;
  //AddPAIModel("pai");
  //

  // Optical physics - only enabled when NOT using Celeritas
  // TrackingManagerConstructor is now registered from main.cc to match
  // the official Celeritas example pattern
#if WITH_CELERITAS
  if (!nnbar::CeleritasInterface::IsEnabled()) {
    // Optical Physics for Geant4 11.x (when Celeritas disabled)
    G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
    auto opticalParams = G4OpticalParameters::Instance();
    opticalParams->SetWLSTimeProfile("delta");
    opticalParams->SetCerenkovMaxPhotonsPerStep(2000);
    opticalParams->SetCerenkovMaxBetaChange(100.0);
    opticalParams->SetCerenkovTrackSecondariesFirst(true);
    opticalParams->SetScintTrackSecondariesFirst(true);
    RegisterPhysics(opticalPhysics);
  }
  // Note: Celeritas TrackingManagerConstructor is registered from main.cc
#else
  // No Celeritas - check if optical physics is enabled via WITH_SCINTILLATION
  #if WITH_SCINTILLATION
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  auto opticalParams = G4OpticalParameters::Instance();
  opticalParams->SetWLSTimeProfile("delta");
  opticalParams->SetCerenkovMaxPhotonsPerStep(2000);
  opticalParams->SetCerenkovMaxBetaChange(100.0);
  opticalParams->SetCerenkovTrackSecondariesFirst(true);
  opticalParams->SetScintTrackSecondariesFirst(true);
  RegisterPhysics(opticalPhysics);
  G4cout << "PhysicsList: Optical physics ENABLED (scintillation + Cerenkov)" << G4endl;
  #else
  G4cout << "PhysicsList: Optical physics DISABLED (fast mode)" << G4endl;
  #endif
#endif
}

PhysicsList::~PhysicsList()
{}

void PhysicsList::ConstructParticle()
{
  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();

  // Only construct optical photon when NOT using Celeritas
  // Celeritas imports all particles and the mismatch between
  // ParticlesParams (4: gamma, opticalphoton, e-, e+) and
  // offload list (3: e-, e+, gamma) causes tracking manager attachment to fail
#if WITH_CELERITAS
  if (!nnbar::CeleritasInterface::IsEnabled()) {
    G4OpticalPhoton::OpticalPhotonDefinition();
  }
#else
  G4OpticalPhoton::OpticalPhotonDefinition();
#endif

  // CRITICAL: Call base class method to invoke ConstructParticle() on all
  // registered physics constructors (including TrackingManagerConstructor
  // for Celeritas). Without this, Celeritas's offload_particles_ is never
  // populated and tracking managers are not attached to particles.
  G4VModularPhysicsList::ConstructParticle();
}



void PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma

  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
  SetCutValue(defaultCutValue, "proton");

  if (verboseLevel>0) DumpCutValuesTable();
}

#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"

void PhysicsList::AddPAIModel(const G4String& modname)
{
  G4ParticleTable::G4PTblDicIterator* theParticleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
  theParticleIterator->reset();
  
  while ((*theParticleIterator)())
  { 
     //std::cout<< "PAI ..... " << std::endl;
     G4ParticleDefinition* particle = theParticleIterator->value();

     std::cout<< particle->GetParticleName() << std::endl;

     G4String partname = particle->GetParticleName();
     if(partname == "e-" || partname == "e+") {NewPAIModel(particle, modname, "eIoni");} 
     else if(partname == "mu-" || partname == "mu+") {NewPAIModel(particle, modname, "muIoni");} 
     else if(partname == "proton" ||partname == "pi+" || partname == "pi-") {NewPAIModel(particle, modname, "hIoni");}
   }
 }
 
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void PhysicsList::NewPAIModel(const G4ParticleDefinition* part, 
                               const G4String& modname,
                               const G4String& procname){

  //auto fEmPhysicsList = new G4EmLivermorePhysics();

  G4String partname = part->GetParticleName();
   if(modname == "pai") {
     std::cout<< "Adding pai model" << std::endl;
     G4PAIModel* pai = new G4PAIModel(part,"PAIModel");
     fConfig->SetExtraEmModel(partname,procname,pai,"TPC_region",0.0,100.*GeV,pai);
     fConfig->SetExtraEmModel(partname,procname,pai,"Silicon_region",0.0,100.*GeV,pai);
     //fConfig->SetExtraEmModel(partname,procname,fEmPhysicsList,"TPC_region",0.0,100.*TeV,fEmPhysicsList);
   } 
   else if(modname == "pai_photon") {
     G4PAIPhotModel* pai = new G4PAIPhotModel(part,"PAIPhotModel");
     fConfig->SetExtraEmModel(partname,procname,pai,"TPC_region",0.0,100.*TeV,pai);
   }
}