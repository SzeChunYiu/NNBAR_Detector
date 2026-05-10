#include "core/DetectorConstruction.hh"
#include "core/ActionInitialization.hh"
#include "generator/G4MCPLGenerator.hh"
#include "core/Analysis.hh"
#include "G4Types.hh"
#include "G4OpticalPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "core/PhysicsList.hh"

#include "G4RunManagerFactory.hh"
#include <G4ProductionCuts.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>
#include "G4UImanager.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4SystemOfUnits.hh"
#include "generator/G4MCPLGenerator.hh"
#include "generator/G4MCPLWriter.hh"
#include  "mcpl.h"

#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"

#include "config.h"

#include <filesystem>
#include <fstream>
#include <sstream>

#if WITH_CELERITAS
#include "physics/CeleritasInterface.hh"
#include <accel/TrackingManagerConstructor.hh>
#include <accel/TrackingManagerIntegration.hh>
#include <accel/SetupOptions.hh>
#include <accel/AlongStepFactory.hh>
#endif

#include "arrow/io/file.h"
#include "parquet/stream_writer.h"
#include "parquet_writer.h"
#include "stdint.h"
#include "nlohmann/json.hpp"

// Dashboard must be included AFTER Arrow/Parquet to avoid Qt signal conflicts
#if WITH_DASHBOARD
#include "gui/DashboardWindow.hh"
#include "G4UIQt.hh"
#endif

// Global simulation state
G4double event_number_global = 1;
G4int run_number = 0;

extern mcpl_outfile_t f;

int theta_bin_index = 0;
int KE_bin_index = 0;
G4int particle_name_input = 0;

// Global flag for generator selection (defined in ActionInitialization.cc)
extern G4bool g_useParticleGun;

// Pre-init messenger for generator selection (must be set BEFORE /run/initialize)
class GeneratorPreInitMessenger : public G4UImessenger {
public:
  GeneratorPreInitMessenger() {
    fDir = new G4UIdirectory("/generator/");
    fDir->SetGuidance("Generator pre-initialization commands");

    fUseParticleGunCmd = new G4UIcmdWithABool("/generator/use_particle_gun", this);
    fUseParticleGunCmd->SetGuidance("Use particle gun instead of MCPL (set BEFORE /run/initialize)");
    fUseParticleGunCmd->SetParameterName("use_gun", false);
    fUseParticleGunCmd->AvailableForStates(G4State_PreInit);
  }
  ~GeneratorPreInitMessenger() {
    delete fUseParticleGunCmd;
    delete fDir;
  }
  void SetNewValue(G4UIcommand* cmd, G4String newValue) override {
    if (cmd == fUseParticleGunCmd) {
      g_useParticleGun = fUseParticleGunCmd->GetNewBoolValue(newValue);
      G4cout << "Generator mode set to: "
             << (g_useParticleGun ? "Particle Gun" : "MCPL") << G4endl;
    }
  }
private:
  G4UIdirectory* fDir;
  G4UIcmdWithABool* fUseParticleGunCmd;
};

namespace {
  void PrintUsage() {
    G4cerr << " Usage: " << G4endl;
    G4cerr << " nnbar-calo-sim [-m macro] [-u UIsession] [-t nThreads] [-g]" << G4endl;
    G4cerr << "   -m macro     : batch mode with macro file" << G4endl;
    G4cerr << "   -u session   : UI session type (tcsh, csh, etc.)" << G4endl;
    G4cerr << "   -t nThreads  : number of threads (1=sequential)" << G4endl;
    G4cerr << "   -g / --gun   : use particle gun instead of MCPL input" << G4endl;
  }

  bool RequireOptionValue(int argc, int index, const char* option) {
    if (index + 1 < argc) {
      return true;
    }
    G4cerr << "Missing value for option " << option << G4endl;
    PrintUsage();
    return false;
  }

  bool ConfigureMCPLForBatchMacro(const G4String& macro) {
#if MCPL_BUILD
    if (g_useParticleGun || macro.empty()) {
      return true;
    }

    std::ifstream input(macro.c_str());
    if (!input) {
      G4cerr << "[MCPL] Cannot open macro for MCPL preflight: " << macro << G4endl;
      return false;
    }

    G4String mcplFile;
    std::string line;
    while (std::getline(input, line)) {
      const auto comment = line.find('#');
      if (comment != std::string::npos) {
        line = line.substr(0, comment);
      }
      std::istringstream iss(line);
      std::string command;
      iss >> command;
      if (command == "/particle_generator/set_mcpl_file") {
        std::string path;
        iss >> path;
        if (!path.empty()) {
          mcplFile = path;
        }
      }
    }

    if (mcplFile.empty()) {
      G4cerr << "[MCPL] MCPL build selected, but macro " << macro
             << " does not set /particle_generator/set_mcpl_file." << G4endl;
      G4cerr << "[MCPL] Use a signal/cosmic MCPL macro, or run particle-gun macros with -g/--gun." << G4endl;
      return false;
    }

    if (!std::filesystem::exists(std::filesystem::path(mcplFile.c_str()))) {
      G4cerr << "[MCPL] Configured MCPL input does not exist: " << mcplFile << G4endl;
      G4cerr << "[MCPL] The path is resolved relative to the current run directory." << G4endl;
      return false;
    }

    G4MCPLGenerator::SetInputFile(mcplFile);
#else
    (void)macro;
#endif
    return true;
  }
}

int main(int argc, char** argv)
{  
  G4String macro;
  G4String session;
  G4int nThreads = 0;

  for (G4int i = 1; i < argc; ++i) {
    const G4String option(argv[i]);

    if (option == "-m") {
      if (!RequireOptionValue(argc, i, "-m")) {
        return 1;
      }
      macro = argv[++i];
    } else if (option == "-u") {
      if (!RequireOptionValue(argc, i, "-u")) {
        return 1;
      }
      session = argv[++i];
    } else if (option == "-t") {
      if (!RequireOptionValue(argc, i, "-t")) {
        return 1;
      }
      nThreads = G4UIcommand::ConvertToInt(argv[++i]);
    } else if (option == "-g" || option == "--gun") {
      g_useParticleGun = true;
    } else {
      PrintUsage();
      return 1;
    }
  }  

  if (!ConfigureMCPLForBatchMacro(macro)) {
    return 2;
  }
    
  // Detect interactive mode (if no macro provided) and define UI session
  G4UIExecutive* ui = nullptr;
#if WITH_DASHBOARD
  // Use Qt mode to get G4UIQt with proper OpenGL visualization
  // G4UIExecutive will create the QApplication internally
  ui = new G4UIExecutive(argc, argv, "qt");
#else
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv, session);
  }
#endif

  // Optionally: choose a different Random engine...
  //
  // G4Random::setTheEngine(new CLHEP::MTwistEngine);
  
  // Construct the run manager using factory (supports both MT and sequential)
  // Celeritas supports MT mode - each worker thread gets its own LocalTransporter
  G4RunManagerType rmType = G4RunManagerType::Default;

  // Determine number of threads: -t N sets N threads, -t 1 means sequential mode
  int numThreads = nThreads > 0 ? nThreads : G4Threading::G4GetNumberOfCores();

  // Use sequential mode if explicitly requested with -t 1
  if (numThreads == 1) {
    rmType = G4RunManagerType::Serial;
    G4cout << "[Main] Sequential mode (single thread)" << G4endl;
  }

  auto runManager = G4RunManagerFactory::CreateRunManager(rmType);

  // Set number of threads for MT mode (Celeritas will auto-detect and create streams)
  if (numThreads > 1) {
    runManager->SetNumberOfThreads(numThreads);
    G4cout << "[Main] Multi-threaded mode with " << numThreads << " threads" << G4endl;
  }

  // Set mandatory initialization classes

  // Create physics list
  G4VModularPhysicsList* physicsList = new PhysicsList();
  // Set energy range for production cuts
  // Minimum raised from 30 eV to 1 keV to avoid tracking sub-keV particles
  // which have no meaningful physics and cause simulation slowdown
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(1.0*keV, 10.0*TeV);

#if WITH_CELERITAS
  // Initialize Celeritas following the EXACT pattern from the official example:
  // https://github.com/celeritas-project/celeritas/blob/develop/example/offload-template/run-offload.cc
  if (nnbar::CeleritasInterface::IsEnabled()) {
    // 1. Get TMI instance FIRST (exactly like official example)
    auto& tmi = celeritas::TrackingManagerIntegration::Instance();
    G4cout << "[Main] Got TrackingManagerIntegration instance" << G4endl;

    // 2. Register TrackingManagerConstructor with physics list
    G4cout << "[Main] Registering Celeritas TrackingManagerConstructor" << G4endl;
    physicsList->RegisterPhysics(new celeritas::TrackingManagerConstructor(&tmi));

    // 3. Register physics list with run manager
    runManager->SetUserInitialization(physicsList);
    G4cout << "[Main] Physics list registered" << G4endl;

    // 4. SetOptions with optimized parameters for RTX A3000 (6GB VRAM)
    celeritas::SetupOptions opts;

    // GPU track capacity optimized for RTX A3000
    opts.max_num_tracks = 1024 * 64;      // 64K tracks per batch
    opts.initializer_capacity = 1024 * 256 * 4;  // 1M initializer slots
    opts.secondary_stack_factor = 3.0;    // More secondary capacity

    // Batching for efficiency
    opts.auto_flush = 1024 * 8;           // Flush at 8K tracks

    // Skip problematic processes
    opts.ignore_processes = {"CoulombScat"};

    // Use efficient along-step (no field in NNBAR)
    opts.make_along_step = celeritas::UniformAlongStepFactory();

    // SD callbacks crash with NNBAR geometry due to ORANGE incompatibility
    // (core dump when hits need to be reconstructed)
    // Using CELER_DISABLE=1 environment variable for full hit tracking instead
    opts.sd.enabled = false;

    // Output files
    opts.output_file = "celeritas.out.json";

    G4cout << "[Main] Calling tmi.SetOptions()..." << G4endl;
    tmi.SetOptions(std::move(opts));
    G4cout << "[Main] Celeritas configured for GPU acceleration" << G4endl;
  } else {
    // Celeritas disabled - just register physics list
    runManager->SetUserInitialization(physicsList);
  }
#else
  // Register physics list with run manager (no Celeritas)
  runManager->SetUserInitialization(physicsList);
#endif

  auto detConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detConstruction);

  auto actionInitialization = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);

  runManager->Initialize();
  
  // Initialize visualization
  auto visManager = new G4VisExecutive;
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  auto UImanager = G4UImanager::GetUIpointer();

#if WITH_DASHBOARD
  // Get G4UIQt session for accessing the Qt main window
  G4UIQt* g4uiqt = dynamic_cast<G4UIQt*>(ui->GetSession());
  if (g4uiqt) {
    G4cout << "[Main] G4UIQt session obtained successfully" << G4endl;

    // Initialize the OpenGL viewer using vis.mac
    UImanager->ApplyCommand("/control/execute vis.mac");
    G4cout << "[Main] OpenGL visualization initialized" << G4endl;

    // Get G4UIQt's main window and add our dashboard panels as dock widgets
    QMainWindow* g4MainWindow = g4uiqt->GetMainWindow();
    if (g4MainWindow) {
      // Initialize our dashboard in embedded mode (adds dock widgets to G4UIQt window)
      nnbar::DashboardWindow::Instance().InitializeEmbedded(g4MainWindow);
      G4cout << "[Main] Dashboard panels added to G4UIQt window" << G4endl;
    }
  } else {
    G4cerr << "[Main] Warning: Could not get G4UIQt session" << G4endl;
    // Fallback to standalone dashboard
    nnbar::DashboardWindow::Instance().Initialize();
    nnbar::DashboardWindow::Instance().Show();
  }

  // Execute macro if provided
  if (macro.size()) {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }

  // Start G4UIQt session (runs Qt event loop)
  if (ui) {
    ui->SessionStart();
    delete ui;
  }
#else
  if (macro.size()) {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
  else  {
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }
#endif

  // Job termination
  delete visManager;
  delete runManager;
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
