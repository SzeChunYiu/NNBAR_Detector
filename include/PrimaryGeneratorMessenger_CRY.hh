
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorMessenger_CRY_h
#define PrimaryGeneratorMessenger_CRY_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction_CRY;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorMessenger_CRY: public G4UImessenger
{
  public:
    PrimaryGeneratorMessenger_CRY(PrimaryGeneratorAction_CRY*);
   ~PrimaryGeneratorMessenger_CRY();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    PrimaryGeneratorAction_CRY*      Action;
    G4UIdirectory*               CRYDir;
    G4UIcmdWithAString*          FileCmd; 
    G4UIcmdWithAString*          InputCmd;
    G4UIcmdWithoutParameter*     UpdateCmd;
    std::string* MessInput;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

