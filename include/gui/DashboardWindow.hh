// ============================================================================
// DashboardWindow.hh
// Professional Qt-based monitoring dashboard for NNBAR detector simulation
// Commercial-grade interface for physics experts
// ============================================================================

#ifndef DASHBOARD_WINDOW_HH
#define DASHBOARD_WINDOW_HH

#include "config.h"

#include <string>
#include <vector>
#include <cstdint>
#include <array>

// Data structures - always available
namespace nnbar {

struct ParticleInfo {
    int32_t event_id;
    int32_t track_id;
    std::string name;
    double ke;
    double x, y, z;
    double px, py, pz;
    std::string volume;
    std::string process;
};

struct EventStatistics {
    int32_t event_id = 0;
    int32_t n_primaries = 0;
    int32_t n_secondaries = 0;
    int32_t n_scintillator_hits = 0;
    int32_t n_tpc_hits = 0;
    int32_t n_leadglass_hits = 0;
    double total_edep = 0.0;
    double cpu_time_ms = 0.0;
    int32_t n_optical_photons = 0;
};

struct DetectorHitMap {
    std::vector<double> scintillator_edep;
    std::vector<double> leadglass_edep;
    std::vector<int32_t> tpc_hits;
};

} // namespace nnbar

// ============================================================================
// Qt Dashboard - only compiled when WITH_DASHBOARD is enabled
// ============================================================================
#if WITH_DASHBOARD

#include <QMainWindow>
#include <QDockWidget>
#include <QTableWidget>
#include <QLabel>
#include <QTimer>
#include <QProgressBar>
#include <QThread>
#include <QApplication>
#include <atomic>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QSplitter>
#include <QPushButton>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QListWidget>
#include <QLineEdit>
#include <QTextEdit>
#include <QToolButton>
#include <QTabWidget>
#include <QScrollArea>
#include <QStackedWidget>
#include <QToolBar>
#include <QButtonGroup>

#include "G4Types.hh"
#include <mutex>

// Forward declarations
namespace nnbar {
class EventDisplay;
class G4CommandPanel;
class DashboardPanel;
class KPICard;
class MiniHistogram;
class TPCPadDisplay;
}

namespace nnbar {

// Worker class for running simulation in separate thread
class SimulationWorker : public QObject {
    Q_OBJECT
public:
    explicit SimulationWorker(QObject* parent = nullptr) : QObject(parent) {}

public slots:
    void runEvents(int nEvents);
    void runSingleEvent();

signals:
    void started();
    void finished();
    void eventCompleted(int eventNum);
    void error(const QString& message);
};

class DashboardWindow : public QMainWindow {
    Q_OBJECT

public:
    static DashboardWindow& Instance();
    void Initialize();
    void InitializeEmbedded(QMainWindow* g4MainWindow);  // Embed in G4UIQt window
    void Show();
    void Hide();

    // Data update methods (thread-safe)
    void AddParticle(const ParticleInfo& info);
    void UpdateEventStats(const EventStatistics& stats);
    void UpdateHitMap(const DetectorHitMap& hits);
    void SetProgress(int current, int total);
    void UpdateEventCount(int count);

    // Event display access
    EventDisplay* GetEventDisplay() const { return m_eventDisplay; }

    // Track and vertex data for event display (thread-safe)
    void AddTrackToDisplay(int32_t trackId, int32_t parentId, int32_t pdg,
                           float energy, const std::string& name,
                           const std::string& creatorProcess,
                           const std::vector<std::array<float, 3>>& points);
    void AddVertexToDisplay(float x, float y, float z, float t,
                            const std::string& process,
                            int32_t parentTrackId,
                            const std::vector<int32_t>& daughterPDGs,
                            float totalEnergy);
    void AddCaloHitToDisplay(float x, float y, float z, float energy, int detType);
    void ClearEventDisplay();
    void FinalizeEvent(int32_t eventId, int32_t runId,
                       const std::string& primaryParticle, double primaryEnergy);

    DashboardWindow(const DashboardWindow&) = delete;
    DashboardWindow& operator=(const DashboardWindow&) = delete;

signals:
    void particleAdded();
    void statsUpdated();
    void hitMapUpdated();
    void progressChanged(int current, int total);
    void startSimulation(int nEvents);
    void startSingleEvent();

private slots:
    void onParticleAdded();
    void onStatsUpdated();
    void onHitMapUpdated();
    void onProgressChanged(int current, int total);
    void onRefreshTimer();

    // Simulation worker slots
    void onSimulationStarted();
    void onSimulationFinished();
    void onEventCompleted(int eventNum);
    void onSimulationError(const QString& msg);

    // Control slots
    void onRunEvents();
    void onSingleEvent();
    void onStopSimulation();
    void onPauseSimulation();
    void onResetSimulation();
    void onClearData();
    void onApplySource();
    void onExportData();
    void onExportCSV();
    void onSaveConfig();
    void onOpenOutputFolder();
    void onInitVis();
    void onToggleG4Viewer();

    // Command console slots
    void onExecuteCommand();
    void onCommandHistoryClicked(QListWidgetItem* item);
    void onQuickCommand(const QString& cmd);
    void onLoadMacro();
    void onClearConsole();

    // Event display slots
    void onViewTransverse();
    void onViewLongitudinal();
    void onView3D();
    void onResetView();
    void onPageChanged(int index);

    // Timer-based simulation slot
    void onSimTimerTick();

private:
    DashboardWindow();
    ~DashboardWindow();

    void CreateMenus();
    void CreateToolBar();
    void CreateMainLayout();
    void CreateEventDisplayPanel();
    void CreateControlPage();
    void CreateCommandPage();
    void CreateAnalysisPage();
    void CreateMaterialBudgetPage();
    void CreateStatsPanel();
    void CreateTPCPadPanel();
    void CreateStatusBar();
    void AddToCommandHistory(const QString& cmd);

    // Main layout components
    QSplitter* m_mainSplitter = nullptr;
    EventDisplay* m_eventDisplay = nullptr;
    QStackedWidget* m_pageStack = nullptr;
    QToolBar* m_pageToolBar = nullptr;
    QButtonGroup* m_pageButtonGroup = nullptr;

    // View control buttons
    QPushButton* m_viewXYButton = nullptr;
    QPushButton* m_viewRZButton = nullptr;
    QPushButton* m_view3DButton = nullptr;

    // Dock widgets (for stats)
    QDockWidget* m_statsDock = nullptr;

    // Tables and displays
    QTableWidget* m_particleTable = nullptr;

    // Statistics labels
    QLabel* m_eventLabel = nullptr;
    QLabel* m_primaryLabel = nullptr;
    QLabel* m_secondaryLabel = nullptr;
    QLabel* m_edepLabel = nullptr;
    QLabel* m_cpuLabel = nullptr;
    QLabel* m_scintHitsLabel = nullptr;
    QLabel* m_tpcHitsLabel = nullptr;
    QLabel* m_lgHitsLabel = nullptr;
    QLabel* m_photonLabel = nullptr;
    QLabel* m_evtRateLabel = nullptr;
    QLabel* m_gpuStatusLabel = nullptr;
    QLabel* m_statusLabel = nullptr;

    // Control widgets
    QPushButton* m_runButton = nullptr;
    QPushButton* m_singleEventButton = nullptr;
    QPushButton* m_stopButton = nullptr;
    QPushButton* m_pauseButton = nullptr;
    QPushButton* m_resetButton = nullptr;
    QSpinBox* m_eventsSpinBox = nullptr;
    QComboBox* m_particleCombo = nullptr;
    QDoubleSpinBox* m_energySpinBox = nullptr;
    QDoubleSpinBox* m_posXSpin = nullptr;
    QDoubleSpinBox* m_posYSpin = nullptr;
    QDoubleSpinBox* m_posZSpin = nullptr;
    QCheckBox* m_autoSaveCheck = nullptr;
    QProgressBar* m_progressBar = nullptr;

    // Command console widgets
    QLineEdit* m_commandInput = nullptr;
    QListWidget* m_commandHistory = nullptr;
    QTextEdit* m_commandOutput = nullptr;
    QTabWidget* m_commandTabs = nullptr;

    // New PowerBI-style components
    G4CommandPanel* m_g4CommandPanel = nullptr;
    DashboardPanel* m_eventDisplayPanel = nullptr;
    DashboardPanel* m_controlPanel = nullptr;
    DashboardPanel* m_statisticsPanel = nullptr;
    DashboardPanel* m_commandPanel = nullptr;

    // KPI Cards
    KPICard* m_eventCountKPI = nullptr;
    KPICard* m_edepKPI = nullptr;
    KPICard* m_hitRateKPI = nullptr;
    KPICard* m_eventRateKPI = nullptr;

    // Mini histograms for live data
    MiniHistogram* m_edepHistogram = nullptr;
    MiniHistogram* m_hitHistogram = nullptr;

    // TPC Pad Display
    TPCPadDisplay* m_tpcPadDisplay = nullptr;
    DashboardPanel* m_tpcPadPanel = nullptr;

    // Display filter widgets
    QCheckBox* m_showTracksCheck = nullptr;
    QCheckBox* m_showVerticesCheck = nullptr;
    QCheckBox* m_showHitsCheck = nullptr;
    QCheckBox* m_showArrowsCheck = nullptr;
    QCheckBox* m_showElectronsCheck = nullptr;
    QCheckBox* m_showGammasCheck = nullptr;
    QCheckBox* m_showNeutronsCheck = nullptr;
    QCheckBox* m_showProtonsCheck = nullptr;
    QCheckBox* m_showPionsCheck = nullptr;
    QCheckBox* m_showMuonsCheck = nullptr;
    QDoubleSpinBox* m_minEnergySpin = nullptr;

    // Energy deposition table (embedded mode)
    QTableWidget* m_edepTable = nullptr;

    // TPC drift stats labels (embedded mode)
    QLabel* m_tpcClustersLabel = nullptr;
    QLabel* m_tpcElectronsInLabel = nullptr;
    QLabel* m_tpcCollectedLabel = nullptr;
    QLabel* m_tpcDriftTimeLabel = nullptr;

    // Data buffers
    std::mutex m_mutex;
    std::vector<ParticleInfo> m_particles;
    EventStatistics m_currentStats;
    DetectorHitMap m_hitMap;

    QTimer* m_refreshTimer = nullptr;
    bool m_initialized = false;
    int m_maxParticles = 1000;
    QMainWindow* m_g4MainWindow = nullptr;  // G4UIQt's main window when embedded
    bool m_embeddedMode = false;

    // Threading for non-blocking simulation
    QThread* m_simThread = nullptr;
    SimulationWorker* m_simWorker = nullptr;
    std::atomic<bool> m_isSimulating{false};

    // Timer-based simulation (main thread, non-blocking)
    QTimer* m_simTimer = nullptr;
    int m_pendingEvents = 0;
    int m_completedEvents = 0;
};

// Global helper functions
void RecordParticleForDashboard(
    int32_t event_id, int32_t track_id,
    const std::string& name,
    double ke, double x, double y, double z,
    double px, double py, double pz,
    const std::string& volume,
    const std::string& process
);

void UpdateDashboardEventStats(const EventStatistics& stats);
void UpdateDashboardHitMap(const DetectorHitMap& hits);

} // namespace nnbar

#else // !WITH_DASHBOARD

namespace nnbar {

class EventDisplay;  // Forward declaration

class DashboardWindow {
public:
    static DashboardWindow& Instance() {
        static DashboardWindow instance;
        return instance;
    }
    void Initialize() {}
    void Show() {}
    void Hide() {}
    void UpdateEventCount(int) {}
    EventDisplay* GetEventDisplay() const { return nullptr; }
};

inline void RecordParticleForDashboard(
    int32_t, int32_t, const std::string&,
    double, double, double, double,
    double, double, double,
    const std::string&, const std::string&) {}

inline void UpdateDashboardEventStats(const EventStatistics&) {}
inline void UpdateDashboardHitMap(const DetectorHitMap&) {}

} // namespace nnbar

#endif // WITH_DASHBOARD

#endif // DASHBOARD_WINDOW_HH
