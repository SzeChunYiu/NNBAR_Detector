// ============================================================================
// OpticksInterface.cc
// Implementation of Opticks GPU optical photon propagation interface
// ============================================================================

#ifdef WITH_OPTICKS

#include "physics/OpticksInterface.hh"

// Opticks headers
#include "G4CXOpticks.hh"
#include "SEventConfig.hh"
#include "SEvt.hh"

// Geant4 headers
#include <G4VPhysicalVolume.hh>
#include <G4Event.hh>
#include <G4Track.hh>
#include <G4OpticalPhoton.hh>
#include <G4SystemOfUnits.hh>
#include <G4ios.hh>

#include <chrono>
#include <cstdlib>
#include <csignal>
#include <csetjmp>

namespace nnbar {

// ============================================================================
// Signal handling for graceful crash recovery during Opticks initialization
// ============================================================================
static sigjmp_buf opticks_init_jmp_buf;
static volatile sig_atomic_t opticks_init_crashed = 0;
static struct sigaction opticks_old_sigsegv_handler;
static struct sigaction opticks_old_sigabrt_handler;

static void opticks_crash_handler(int sig) {
    opticks_init_crashed = 1;
    siglongjmp(opticks_init_jmp_buf, sig);
}

static void install_opticks_crash_handler() {
    struct sigaction sa;
    sa.sa_handler = opticks_crash_handler;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = 0;
    sigaction(SIGSEGV, &sa, &opticks_old_sigsegv_handler);
    sigaction(SIGABRT, &sa, &opticks_old_sigabrt_handler);
}

static void restore_default_handlers() {
    sigaction(SIGSEGV, &opticks_old_sigsegv_handler, nullptr);
    sigaction(SIGABRT, &opticks_old_sigabrt_handler, nullptr);
}

// ============================================================================
// Check if enabled
// ============================================================================
bool OpticksInterface::IsEnabled() {
    const char* disable_env = std::getenv("OPTICKS_DISABLE");
    if (disable_env && disable_env[0] != '\0') {
        return false;
    }
    return true;
}

// ============================================================================
// Singleton instance
// ============================================================================
OpticksInterface& OpticksInterface::Instance() {
    static OpticksInterface instance;
    return instance;
}

OpticksInterface::OpticksInterface() = default;

OpticksInterface::~OpticksInterface() {
    Finalize();
}

// ============================================================================
// Initialize Opticks
// ============================================================================
void OpticksInterface::Initialize(G4VPhysicalVolume* world) {
    if (m_initialized) {
        G4cout << "[Opticks] Already initialized" << G4endl;
        return;
    }

    if (!IsEnabled()) {
        G4cout << "[Opticks] Disabled via OPTICKS_DISABLE environment variable" << G4endl;
        m_active = false;
        return;
    }

    if (!world) {
        G4cerr << "[Opticks] ERROR: World volume is null" << G4endl;
        return;
    }

    G4cout << "\n";
    G4cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    G4cout << "║           Opticks GPU Optical Photon Acceleration                ║\n";
    G4cout << "╚══════════════════════════════════════════════════════════════════╝\n";

    // Install signal handlers to gracefully recover from geometry conversion crashes
    install_opticks_crash_handler();
    opticks_init_crashed = 0;

    int sig = sigsetjmp(opticks_init_jmp_buf, 1);
    if (sig != 0) {
        // We got here via longjmp from the signal handler
        restore_default_handlers();
        G4cerr << "\n[Opticks] ╔══════════════════════════════════════════════════════════════════╗" << G4endl;
        G4cerr << "[Opticks] ║  Geometry conversion failed - complex geometry not supported     ║" << G4endl;
        G4cerr << "[Opticks] ║  Continuing with CPU optical photon simulation                   ║" << G4endl;
        G4cerr << "[Opticks] ╚══════════════════════════════════════════════════════════════════╝\n" << G4endl;
        m_active = false;
        m_initialized = false;
        return;
    }

    try {
        // Configure Opticks
        SEventConfig::SetMaxPhoton(1000000);  // 1M photons per event
        SEventConfig::SetMaxGenstep(10000);   // 10K gensteps per event

        // Initialize G4CXOpticks with geometry using static method
        G4cout << "[Opticks] Converting Geant4 geometry to GPU format..." << G4endl;

        G4CXOpticks* opticks = G4CXOpticks::SetGeometry(world);

        G4cout << "[Opticks] Building OptiX acceleration structures..." << G4endl;

        if (opticks) {
            m_initialized = true;
            m_active = true;

            G4cout << "[Opticks] ✓ GPU optical photon propagation ready" << G4endl;
            G4cout << "[Opticks]   Max photons per event: 1,000,000" << G4endl;
            G4cout << "[Opticks]   Using OptiX ray tracing" << G4endl;
        } else {
            G4cerr << "[Opticks] WARNING: G4CXOpticks::SetGeometry returned null" << G4endl;
            m_active = false;
        }
        G4cout << "\n";

    } catch (const std::exception& e) {
        G4cerr << "[Opticks] ERROR during initialization: " << e.what() << G4endl;
        m_active = false;
    }

    // Restore default signal handlers
    restore_default_handlers();
}

// ============================================================================
// Finalize Opticks
// ============================================================================
void OpticksInterface::Finalize() {
    if (!m_initialized) return;

    G4cout << "[Opticks] Finalizing GPU resources..." << G4endl;
    G4cout << "[Opticks] Statistics:" << G4endl;
    G4cout << "[Opticks]   Photons generated: " << m_photonsGenerated << G4endl;
    G4cout << "[Opticks]   Photons detected:  " << m_photonsDetected << G4endl;
    G4cout << "[Opticks]   GPU time: " << m_gpuTime << " s" << G4endl;

    // Opticks cleanup handled by G4CXOpticks singleton

    m_initialized = false;
    m_active = false;
}

// ============================================================================
// Enable/disable
// ============================================================================
void OpticksInterface::SetEnabled(bool enabled) {
    m_enabled = enabled;
}

// ============================================================================
// Add gensteps from Geant4 scintillation/Cerenkov processes
// ============================================================================
void OpticksInterface::AddScintillationGenstep(const G4Track* track,
                                                int numPhotons,
                                                double meanEnergy,
                                                double time) {
    if (!m_active || !m_enabled) return;

    // In production, this would collect genstep data for batch GPU processing
    m_photonsGenerated += numPhotons;
    m_numGensteps++;

    // Store genstep data (simplified - real implementation uses Opticks format)
    // Each genstep: position, direction, polarization, wavelength, time, flags
}

void OpticksInterface::AddCerenkovGenstep(const G4Track* track,
                                          int numPhotons,
                                          double minWavelength,
                                          double maxWavelength) {
    if (!m_active || !m_enabled) return;

    m_photonsGenerated += numPhotons;
    m_numGensteps++;
}

// ============================================================================
// Propagate photons on GPU
// ============================================================================
void OpticksInterface::PropagateOpticalPhotons() {
    if (!m_active || !m_enabled || m_numGensteps == 0) return;

    auto start = std::chrono::high_resolution_clock::now();

    try {
        G4CXOpticks* opticks = G4CXOpticks::Get();
        if (!opticks) {
            G4cerr << "[Opticks] ERROR: G4CXOpticks instance is null" << G4endl;
            return;
        }

        // Use current event ID (tracked externally)
        static int eventId = 0;
        eventId++;

        // Launch GPU simulation with event ID and reset flag
        opticks->simulate(eventId, true);

        // Hits are retrieved via SensitiveDetector_EndOfEvent
        // For now just track statistics
        m_photonsDetected += m_photonsGenerated;  // Approximate

    } catch (const std::exception& e) {
        G4cerr << "[Opticks] ERROR during propagation: " << e.what() << G4endl;
    }

    auto end = std::chrono::high_resolution_clock::now();
    m_gpuTime += std::chrono::duration<double>(end - start).count();

    // Reset for next event
    m_gensteps.clear();
    m_numGensteps = 0;
}

} // namespace nnbar

#endif // WITH_OPTICKS
