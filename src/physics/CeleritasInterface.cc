// ============================================================================
// CeleritasInterface.cc
// Implementation of Celeritas GPU-accelerated EM physics interface
// ============================================================================

#include "config.h"

#if WITH_CELERITAS

#include "physics/CeleritasInterface.hh"

// Celeritas headers (modern accel API)
#include <accel/SetupOptions.hh>
#include <accel/TrackingManagerIntegration.hh>
#include <accel/TrackingManagerConstructor.hh>
#include <accel/AlongStepFactory.hh>

// Geant4 headers
#include <G4Run.hh>
#include <G4ios.hh>
#include <G4Electron.hh>
#include <G4Positron.hh>
#include <G4Gamma.hh>

#include <cstdlib>

namespace nnbar {

// ============================================================================
// Static member initialization
// ============================================================================
unsigned int CeleritasInterface::s_maxNumTracks = 262144;  // 256K tracks
unsigned int CeleritasInterface::s_maxSteps = 1000;
const char* CeleritasInterface::s_outputFile = "celeritas.out.json";
bool CeleritasInterface::s_configured = false;

// ============================================================================
// Configuration
// ============================================================================
void CeleritasInterface::SetMaxNumTracks(unsigned int n) {
    s_maxNumTracks = n;
}

void CeleritasInterface::SetMaxSteps(unsigned int n) {
    s_maxSteps = n;
}

void CeleritasInterface::SetOutputFile(const char* filename) {
    s_outputFile = filename;
}

// ============================================================================
// Apply configuration
// ============================================================================
void CeleritasInterface::Configure() {
    if (s_configured) {
        return;
    }

    // Check if disabled via environment
    if (!IsEnabled()) {
        G4cout << "[Celeritas] Disabled via CELER_DISABLE environment variable" << G4endl;
        s_configured = true;
        return;
    }

    G4cout << "\n";
    G4cout << "======================================================================\n";
    G4cout << "              Celeritas GPU Acceleration Configuration               \n";
    G4cout << "======================================================================\n";

    // Create setup options with proper along-step factory
    celeritas::SetupOptions options;
    options.max_num_tracks = s_maxNumTracks;
    options.max_steps = s_maxSteps;
    options.initializer_capacity = s_maxNumTracks * 2;
    options.secondary_stack_factor = 2.0;
    options.output_file = s_outputFile;

    // Use UniformAlongStepFactory with no field (linear propagation)
    // NNBAR detector has no magnetic field
    options.make_along_step = celeritas::UniformAlongStepFactory();

    // SD callbacks - Try enabling for hit-level data recording
    // Previously disabled due to ORANGE geometry incompatibility, but let's test
    // Check CELER_SD_DISABLE env var to allow runtime control
    const char* sd_disable = std::getenv("CELER_SD_DISABLE");
    if (sd_disable && sd_disable[0] != '\0') {
        options.sd.enabled = false;
        G4cout << "[Celeritas] SD callbacks DISABLED via CELER_SD_DISABLE" << G4endl;
    } else {
        options.sd.enabled = true;
        G4cout << "[Celeritas] SD callbacks ENABLED for hit-level recording" << G4endl;
    }

    // Skip processes that Celeritas cannot handle with complex geometries
    // CoulombScat causes validation errors with inconsistent cross-section limits
    options.ignore_processes = {"CoulombScat"};

    // Do NOT set offload_particles - let Celeritas use its defaults
    // Setting offload_particles = {} explicitly clears it, which breaks things
    // The default behavior handles e-, e+, gamma automatically
    G4cout << "[Celeritas] Particle offload: using defaults (e-, e+, gamma)" << G4endl;

    // Note: Optical physics is disabled when Celeritas is enabled (in PhysicsList)

    G4cout << "[Celeritas] Max tracks:     " << options.max_num_tracks << G4endl;
    G4cout << "[Celeritas] Max steps:      " << options.max_steps << G4endl;
    G4cout << "[Celeritas] Initializers:   " << options.initializer_capacity << G4endl;
    G4cout << "[Celeritas] Field:          None (linear propagation)" << G4endl;
    G4cout << "[Celeritas] Output file:    " << options.output_file << G4endl;
    G4cout << "======================================================================\n";
    G4cout << G4endl;

    // Apply options to the integration singleton
    G4cout << "[Celeritas] Setting options on TrackingManagerIntegration..." << G4endl;
    celeritas::TrackingManagerIntegration::Instance().SetOptions(std::move(options));
    G4cout << "[Celeritas] Options set successfully" << G4endl;

    s_configured = true;
}

// ============================================================================
// Create physics constructor
// ============================================================================
G4VPhysicsConstructor* CeleritasInterface::CreatePhysicsConstructor() {
    if (!IsEnabled()) {
        G4cout << "[Celeritas] Not creating physics constructor (disabled)" << G4endl;
        return nullptr;
    }

    // NOTE: Configure() should be called from main.cc AFTER the physics list
    // is registered with the run manager. This ensures proper TrackingManager attachment.

    G4cout << "[Celeritas] Creating TrackingManagerConstructor" << G4endl;
    return new celeritas::TrackingManagerConstructor{
        &celeritas::TrackingManagerIntegration::Instance()
    };
}

// ============================================================================
// Run action hooks
// ============================================================================
void CeleritasInterface::BeginRun(const G4Run* run) {
    if (!IsEnabled()) return;

    G4cout << "[Celeritas] BeginOfRunAction - initializing GPU transport" << G4endl;
    celeritas::TrackingManagerIntegration::Instance().BeginOfRunAction(run);
}

void CeleritasInterface::EndRun(const G4Run* run) {
    if (!IsEnabled()) return;

    G4cout << "[Celeritas] EndOfRunAction - finalizing GPU transport" << G4endl;
    celeritas::TrackingManagerIntegration::Instance().EndOfRunAction(run);
}

// ============================================================================
// Check if enabled
// ============================================================================
bool CeleritasInterface::IsEnabled() {
    // Check CELER_DISABLE environment variable (same as Celeritas's SharedParams::GetMode())
    const char* disable_env = std::getenv("CELER_DISABLE");
    if (disable_env && disable_env[0] != '\0') {
        return false;
    }
    // Enabled by default when compiled with Celeritas support
    return true;
}

} // namespace nnbar

#endif // WITH_CELERITAS
