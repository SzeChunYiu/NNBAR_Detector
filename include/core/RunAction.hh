// ============================================================================
// RunAction.hh
// Run-level action class for NNBAR detector simulation
// ============================================================================

#ifndef RUNACTION_HH
#define RUNACTION_HH

#include "config.h"
#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class G4GenericMessenger;

class RunAction : public G4UserRunAction {
public:
    RunAction();
    virtual ~RunAction();

    // G4UserRunAction interface
    virtual G4Run* GenerateRun() override;
    virtual void BeginOfRunAction(const G4Run* aRun) override;
    virtual void EndOfRunAction(const G4Run* aRun) override;

    // Messenger commands
    void SetOutputFolder(G4String& folderName);

#if MCPL_BUILD
    void SetMCPLFile(G4String& filename);
#endif

private:
    G4GenericMessenger* fMessenger = nullptr;
};

#endif // RUNACTION_HH
