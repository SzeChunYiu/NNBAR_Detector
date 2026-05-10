// ============================================================================
// NNbarRun.hh
// Custom G4Run class for NNBAR detector simulation
// ============================================================================

#ifndef NNBARRUN_HH
#define NNBARRUN_HH

#include "G4Run.hh"
#include "G4Event.hh"
#include "hits/NNbarHit.hh"
#include <map>

// Forward declarations
class G4HCofThisEvent;

namespace nnbar {
    class ParquetOutputManager;
}

class NNbarRun : public G4Run {
public:
    NNbarRun();
    virtual ~NNbarRun();

    // G4Run interface
    virtual void RecordEvent(const G4Event* event) override;
    virtual void Merge(const G4Run* aRun) override;

private:
    // Hit collection processing methods
    void ProcessCarbonHits(G4HCofThisEvent* hce, int eventId,
                           nnbar::ParquetOutputManager& output);
    void ProcessBeampipeHits(G4HCofThisEvent* hce, int eventId,
                             nnbar::ParquetOutputManager& output);
    void ProcessSiliconHits(G4HCofThisEvent* hce, int eventId,
                            nnbar::ParquetOutputManager& output);
    void ProcessTPCHits(G4HCofThisEvent* hce, int eventId,
                        nnbar::ParquetOutputManager& output);
    void ProcessScintillatorHits(G4HCofThisEvent* hce, int eventId,
                                 nnbar::ParquetOutputManager& output);
    void ProcessLeadGlassHits(G4HCofThisEvent* hce, int eventId,
                              nnbar::ParquetOutputManager& output);
    void ProcessPMTHits(G4HCofThisEvent* hce, int eventId,
                        nnbar::ParquetOutputManager& output);

    // Hit collection IDs
    G4int fCarbonHCID;
    G4int fSiliconHCID;
    G4int fBeampipeHCID;
    G4int fTPCHCID;
    G4int fScintHCID;
    G4int fLeadGlassHCID;
    G4int fPMTHCID;
};

#endif // NNBARRUN_HH
