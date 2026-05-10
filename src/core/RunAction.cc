// ============================================================================
// RunAction.cc
// Manages run-level operations for NNBAR detector simulation
// ============================================================================

#include "core/RunAction.hh"
#include "hits/NNbarRun.hh"
#include "output/ParquetOutputManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4EventManager.hh"
#include "G4TrackingManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"
#include "G4Allocator.hh"
#include "hits/NNbarHit.hh"

#include <iostream>
#include <iomanip>
#include <filesystem>
#include <string>

#include "config.h"

#if MCPL_BUILD
#include "generator/G4MCPLGenerator.hh"
#include "mcpl.h"
extern mcpl_file_t m_mcplfile;
extern const mcpl_particle_t* m_p;
#endif

// Celeritas GPU acceleration
#if WITH_CELERITAS
#include "physics/CeleritasInterface.hh"
#include "physics/CeleritasCalorimeter.hh"
#endif

// TPC GPU drift simulation
#if WITH_GARFIELD_GPU
#include "physics/TPCDriftManager.hh"
#endif

// Opticks GPU optical photon simulation
#if WITH_OPTICKS
#include "physics/OpticksInterface.hh"
#include "G4TransportationManager.hh"
#endif

// GPU computation engines
#include "gpu/GPUManager.hh"
#include "gpu/OpticalPhotonGPU.hh"
#include "gpu/TPCDriftGPU.hh"

// ============================================================================
// Global simulation settings
// ============================================================================
extern G4int run_number;
extern G4double event_number_global;

// Internal state
static G4String s_folderName = "";
static std::string s_particleNames[7] = {
    "", "_neutron", "_proton", "_gamma", "_electron", "_muon", "_pion"
};

extern G4ThreadLocal G4Allocator<NNbarHit>* NNbarHitAllocator;

// ============================================================================
// Constructor
// ============================================================================

RunAction::RunAction() : G4UserRunAction() {
    // Create messenger for run control
    fMessenger = new G4GenericMessenger(this, "/particle_generator/",
                                        "Particle generator control");

    // Run number control
    auto& runCmd = fMessenger->DeclareProperty("set_run_number", run_number,
                                                "Run number of the simulation");
    runCmd.SetParameterName("run_number", true);
    runCmd.SetRange("run_number>=0");
    runCmd.SetDefaultValue("0");

    // Event number control
    auto& eventCmd = fMessenger->DeclareProperty("set_event_number", event_number_global,
                                                  "Global event number");
    eventCmd.SetParameterName("Event_number", true);
    eventCmd.SetRange("Event_number>=0");
    eventCmd.SetDefaultValue("0.");

    // Output folder control
    fMessenger->DeclareMethod("set_folder_name", &RunAction::SetOutputFolder,
                              "Set output folder name");

#if MCPL_BUILD
    // MCPL file control
    fMessenger->DeclareMethod("set_mcpl_file", &RunAction::SetMCPLFile,
                              "Set and read MCPL input file");
#endif

    G4RunManager::GetRunManager()->SetPrintProgress(1);
}

// ============================================================================
// Destructor
// ============================================================================

RunAction::~RunAction() {
    delete fMessenger;
}

// ============================================================================
// Generate Run
// ============================================================================

G4Run* RunAction::GenerateRun() {
    return new NNbarRun;
}

// ============================================================================
// Begin of Run Action
// ============================================================================

void RunAction::BeginOfRunAction(const G4Run* aRun) {
#if !WITH_CELERITAS
    (void)aRun;
#endif

    // ========================================================================
    // GPU Manager Initialization (Master thread only, once)
    // ========================================================================
    static bool gpuManagerInitialized = false;
    if (!gpuManagerInitialized && IsMaster()) {
        auto& gpuMgr = nnbar::gpu::GPUManager::Instance();
        gpuMgr.SetComputeMode(nnbar::gpu::ComputeMode::GPU_PREFERRED);
        gpuMgr.Initialize();
        gpuManagerInitialized = true;
    }

    // ========================================================================
    // GPU Optical Photon Engine Initialization
    // ========================================================================
    static bool opticalGPUInitialized = false;
    if (!opticalGPUInitialized && IsMaster()) {
        auto& opticalGPU = nnbar::gpu::OpticalPhotonGPU::Instance();
        opticalGPU.Initialize();
        opticalGPUInitialized = true;
    }

    // ========================================================================
    // GPU TPC Drift Engine Initialization
    // ========================================================================
    static bool tpcGPUInitialized = false;
    if (!tpcGPUInitialized && IsMaster()) {
        auto& tpcGPU = nnbar::gpu::TPCDriftGPU::Instance();
        tpcGPU.Initialize();
        // Set drift parameters
        tpcGPU.SetDriftVelocity(5.0f);      // mm/us
        tpcGPU.SetDiffusionL(0.02f);        // mm/sqrt(cm)
        tpcGPU.SetDiffusionT(0.02f);        // mm/sqrt(cm)
        tpcGPU.SetElectricField(500.0f);    // V/cm
        tpcGPU.SetGasGain(1000.0f);
        tpcGPUInitialized = true;
    }

    // ========================================================================
    // Celeritas GPU EM Physics (must be called from all threads)
    // ========================================================================
#if WITH_CELERITAS
    nnbar::CeleritasInterface::BeginRun(aRun);
    // Create thread-local calorimeters for GPU energy recording
    nnbar::CeleritasCalorimeter::Instance().CreateThreadLocalSD();
#endif

    // Initialize Opticks GPU optical photon simulation
#if WITH_OPTICKS
    static bool opticksInitialized = false;
    if (!opticksInitialized && IsMaster()) {
        G4VPhysicalVolume* world = G4TransportationManager::GetTransportationManager()
            ->GetNavigatorForTracking()->GetWorldVolume();
        if (world) {
            nnbar::OpticksInterface::Instance().Initialize(world);
            opticksInitialized = true;
        }
    }
#endif

    // Initialize TPC drift manager (GPU/OpenMP acceleration)
#if WITH_GARFIELD_GPU
    static bool tpcDriftInitialized = false;
    if (!tpcDriftInitialized) {
        auto* driftMgr = nnbar::TPCDriftManager::Instance();
        if (driftMgr->Initialize()) {
            // Set TPC geometry (inner=114cm, outer=200cm, half-length=300cm, drift=85cm)
            driftMgr->SetTPCGeometry(114.0f, 200.0f, 300.0f, 85.0f);
            // Set gas properties (Ar/CO2 80/20 at 1 atm, room temp)
            driftMgr->SetGasProperties(0.8f, 0.2f, 760.0f, 293.15f);
            G4cout << "\n╔══════════════════════════════════════════════════════════════════╗" << G4endl;
            G4cout << "║  TPC DRIFT SIMULATION INITIALIZED                                ║" << G4endl;
            G4cout << "║    Device: " << std::left << std::setw(54) << driftMgr->GetDeviceName() << " ║" << G4endl;
            G4cout << "║    GPU Available: " << (driftMgr->IsGPUAvailable() ? "YES" : "NO (using CPU)")
                   << std::setw(driftMgr->IsGPUAvailable() ? 42 : 33) << " ║" << G4endl;
            G4cout << "╚══════════════════════════════════════════════════════════════════╝\n" << G4endl;
            tpcDriftInitialized = true;
        } else {
            G4cerr << "WARNING: TPCDriftManager initialization failed!" << G4endl;
        }
    }
#endif

    // Enable trajectory storage for EventDisplay visualization
#if WITH_DASHBOARD
    G4UImanager::GetUIpointer()->ApplyCommand("/tracking/storeTrajectory 1");
#endif

    // Worker thread: Reset hit allocator if needed
    if (!IsMaster()) {
        if (NNbarHitAllocator != nullptr) {
#if DEBUG_VERBOSE
            G4cout << "Resetting hit allocator, size: "
                   << NNbarHitAllocator->GetAllocatedSize() << G4endl;
#endif
            NNbarHitAllocator->ResetStorage();
        }
        return;
    }

    // Master thread: Initialize output
    std::string outputDir = "./output/" + std::string(s_folderName);

    // Create output directory
    std::error_code ec;
    std::filesystem::create_directories(outputDir, ec);
    if (ec) {
        G4cerr << "Warning: Failed to create output directory: " << outputDir
               << " - " << ec.message() << G4endl;
    }

    G4cout << "========================================" << G4endl;
    G4cout << "Run " << run_number << " starting" << G4endl;
    G4cout << "Output directory: " << outputDir << G4endl;
    G4cout << "========================================" << G4endl;

    // Initialize Parquet output manager
    nnbar::ParquetOutputManager::Instance().Initialize(
        outputDir, run_number, "./config/output_layouts");
}

// ============================================================================
// End of Run Action
// ============================================================================

void RunAction::EndOfRunAction(const G4Run* run) {
#if !WITH_CELERITAS
    (void)run;
#endif

    // Celeritas: must be called from all threads
#if WITH_CELERITAS
    nnbar::CeleritasInterface::EndRun(run);
#endif

    if (!IsMaster()) {
        return;
    }

    // ========================================================================
    // GPU Statistics and Shutdown
    // ========================================================================
    auto& gpuMgr = nnbar::gpu::GPUManager::Instance();
    if (gpuMgr.IsInitialized()) {
        gpuMgr.PrintStatistics();
    }

    // Print optical photon GPU stats
    auto& opticalGPU = nnbar::gpu::OpticalPhotonGPU::Instance();
    if (opticalGPU.IsInitialized()) {
        opticalGPU.PrintStatistics();
    }

    // Print TPC drift GPU stats
    auto& tpcGPU = nnbar::gpu::TPCDriftGPU::Instance();
    if (tpcGPU.IsInitialized()) {
        tpcGPU.PrintStatistics();
    }

    // Finalize Opticks (report statistics)
#if WITH_OPTICKS
    // Note: Full finalization happens at program exit via destructor
    // This just prints statistics for this run
#endif

    // Finalize Parquet output
    nnbar::ParquetOutputManager::Instance().Finalize();

    // Increment run number for next run
    run_number++;

    G4cout << "========================================" << G4endl;
    G4cout << "Run completed. Next run number: " << run_number << G4endl;
    G4cout << "========================================" << G4endl;
}

// ============================================================================
// Set Output Folder
// ============================================================================

void RunAction::SetOutputFolder(G4String& folderName) {
    if (!IsMaster()) {
        return;
    }

    s_folderName = folderName;

    std::string outputDir = "./output/" + std::string(s_folderName);
    std::error_code ec;
    std::filesystem::create_directories(outputDir, ec);

    if (ec) {
        G4cerr << "Warning: Failed to create directory: " << outputDir << G4endl;
    } else {
        G4cout << "Output folder set to: " << outputDir << G4endl;
    }
}

// ============================================================================
// Set MCPL File
// ============================================================================

#if MCPL_BUILD
void RunAction::SetMCPLFile(G4String& filename) {
    G4MCPLGenerator::SetInputFile(filename);
}
#endif
