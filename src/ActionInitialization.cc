#include "ActionInitialization.hh"
#include "config.h"

#if CRY_BUILD==0
  #include "PrimaryGeneratorAction.hh"
#else 
  #include "PrimaryGeneratorAction_CRY.hh"   //CRY is not specified => we will use the normal particle gun :)
#endif

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

  #if CRY_BUILD==1
  std::cout << "check .... CRY build is activated "<<std::endl;
    SetUserAction(new PrimaryGeneratorAction_CRY);
  #else
  std::cout << "check .... CRY build is activated "<<std::endl;
    SetUserAction(new PrimaryGeneratorAction);  
  #endif

  SetUserAction(new RunAction);
  EventAction* eventAction = new EventAction();
  SetUserAction(eventAction);
  SetUserAction(new SteppingAction());
  
}  

//....
