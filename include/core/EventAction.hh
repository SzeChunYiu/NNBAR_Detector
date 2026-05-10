// ============================================================================
// EventAction.hh
// Event-level action for NNBAR detector simulation
// Collects track/hit data and sends to EventDisplay
// ============================================================================

#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"
#include "config.h"

#if WITH_CELERITAS
#include "physics/CeleritasCalorimeter.hh"
#endif

class G4Event;
class SteppingAction;

class EventAction : public G4UserEventAction {
public:
    EventAction(SteppingAction* steppingAction = nullptr);
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event) override;
    virtual void EndOfEventAction(const G4Event* event) override;

    void SetSteppingAction(SteppingAction* sa) { fSteppingAction = sa; }

private:
    void SendEventToDisplay(const G4Event* event);
    SteppingAction* fSteppingAction;

#if WITH_CELERITAS
    // Stores GPU energy from current event for use in SendEventToDisplay
    nnbar::GPUEnergyData fCurrentGPUEnergy;
#endif
};

#endif // EVENTACTION_HH
