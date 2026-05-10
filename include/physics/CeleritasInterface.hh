// ============================================================================
// CeleritasInterface.hh
// Interface for Celeritas GPU-accelerated EM physics
// ============================================================================
// Celeritas offloads electromagnetic physics (e-, e+, gamma) to GPU,
// providing 10-100x speedup for EM shower-dominated simulations.
//
// Requirements:
//   - Geant4 >= 11.0
//   - CUDA Toolkit
//   - Celeritas library (https://github.com/celeritas-project/celeritas)
//
// Usage:
//   1. In PhysicsList constructor (after other physics):
//      #ifdef WITH_CELERITAS
//      RegisterPhysics(nnbar::CeleritasInterface::CreatePhysicsConstructor());
//      #endif
//
//   2. In RunAction::BeginOfRunAction:
//      #ifdef WITH_CELERITAS
//      nnbar::CeleritasInterface::BeginRun(run);
//      #endif
//
//   3. In RunAction::EndOfRunAction:
//      #ifdef WITH_CELERITAS
//      nnbar::CeleritasInterface::EndRun(run);
//      #endif
// ============================================================================

#ifndef CELERITAS_INTERFACE_HH
#define CELERITAS_INTERFACE_HH

#include "config.h"

#if WITH_CELERITAS

#include <memory>

class G4Run;
class G4VPhysicsConstructor;

namespace nnbar {

// ============================================================================
// CeleritasInterface - Static interface for GPU EM physics offloading
// ============================================================================
// This provides a simplified interface to Celeritas TrackingManagerIntegration.
// Use this for simple integration with NNBAR simulation.
// ============================================================================
class CeleritasInterface {
public:
    // Configuration - call before G4RunManager::Initialize()
    static void SetMaxNumTracks(unsigned int n);
    static void SetMaxSteps(unsigned int n);
    static void SetOutputFile(const char* filename);

    // Apply configuration options - call before G4RunManager::Initialize()
    static void Configure();

    // Create physics constructor to register with physics list
    // Call this in PhysicsList constructor
    static G4VPhysicsConstructor* CreatePhysicsConstructor();

    // Run action hooks - call from UserRunAction
    static void BeginRun(const G4Run* run);
    static void EndRun(const G4Run* run);

    // Check if Celeritas is enabled (via CELER_DISABLE env var)
    static bool IsEnabled();

private:
    // Configuration values (applied in Configure())
    static unsigned int s_maxNumTracks;
    static unsigned int s_maxSteps;
    static const char* s_outputFile;
    static bool s_configured;
};

} // namespace nnbar

#endif // WITH_CELERITAS (value check)

#endif // CELERITAS_INTERFACE_HH
