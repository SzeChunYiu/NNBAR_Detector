#include "ActionInitialization.hh"
#include "config.h"

#include "PrimaryGeneratorAction.hh"

#include "G4MCPLGenerator.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

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

void ActionInitialization::Build() const
{ 

  
  SetUserAction(new PrimaryGeneratorAction);  

  #if MCPL_BUILD==1
    std::cout << " MCPL generator is activated" << std::endl;
    SetUserAction(new G4MCPLGenerator("./mcpl_files/cosmic_muon_0.mcpl"));
  
  #else
    SetUserAction(new PrimaryGeneratorAction());
  #endif
  
  
  SetUserAction(new RunAction);
  EventAction* eventAction = new EventAction();
  SetUserAction(eventAction);
  //SetUserAction(new SteppingAction());
  
}  

//....
