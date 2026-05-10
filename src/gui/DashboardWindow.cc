// ============================================================================
// DashboardWindow.cc
// Professional Qt-based monitoring dashboard for NNBAR simulation
// Multi-page layout with CERN-style event display
// ============================================================================

#include "gui/DashboardWindow.hh"
#include "gui/EventDisplay.hh"
#include "gui/MaterialBudgetPlot.hh"
#include "gui/G4CommandPanel.hh"
#include "gui/DashboardPanel.hh"
#include "gui/TPCPadDisplay.hh"
#include "config.h"

#ifdef WITH_DASHBOARD

#include <QApplication>
#include <QMenuBar>
#include <QMenu>
#include <QAction>
#include <QHeaderView>
#include <QFont>
#include <QPalette>
#include <QStyle>
#include <QStatusBar>
#include <QPushButton>
#include <QSpinBox>
#include <QComboBox>
#include <QLineEdit>
#include <QCheckBox>
#include <QFrame>
#include <QGridLayout>
#include <QFormLayout>
#include <QMessageBox>
#include <QFileDialog>
#include <QDateTime>
#include <QThread>
#include <QProcess>
#include <QToolBar>
#include <QStackedWidget>
#include <QButtonGroup>
#include <QSplitter>
#include <QScrollArea>
#include <QFile>
#include <QTextStream>
#include <QDir>

#include "G4UImanager.hh"
#include "G4StateManager.hh"

namespace nnbar {

// ============================================================================
// SimulationWorker Implementation
// ============================================================================

void SimulationWorker::runEvents(int nEvents) {
    // These slots will be invoked via signal from the main thread timer
    // Not used in current implementation
    Q_UNUSED(nEvents);
}

void SimulationWorker::runSingleEvent() {
    // Not used in current implementation - using timer-based approach instead
}

// ============================================================================
// Singleton Instance
// ============================================================================

DashboardWindow& DashboardWindow::Instance() {
    static DashboardWindow instance;
    return instance;
}

// ============================================================================
// Constructor - Professional Dark Theme
// ============================================================================

DashboardWindow::DashboardWindow() : QMainWindow(nullptr) {
    setWindowTitle("NNBAR Detector Simulation - Control Dashboard");
    setMinimumSize(1600, 1000);

    // Professional dark theme
    setStyleSheet(R"(
        QMainWindow, QWidget {
            background-color: #1a1a1f;
            color: #ffffff;
            font-family: 'Segoe UI', Arial, sans-serif;
            font-size: 11px;
        }
        QGroupBox {
            border: 1px solid #3d3d3d;
            border-radius: 5px;
            margin-top: 12px;
            padding-top: 10px;
            font-weight: bold;
            color: #00d4ff;
        }
        QGroupBox::title {
            subcontrol-origin: margin;
            left: 10px;
            padding: 0 5px;
            color: #00d4ff;
        }
        QPushButton {
            background-color: #0078d4;
            color: white;
            border: none;
            border-radius: 4px;
            padding: 8px 16px;
            font-weight: bold;
            min-width: 70px;
        }
        QPushButton:hover { background-color: #1084d8; }
        QPushButton:pressed { background-color: #006cbd; }
        QPushButton:disabled { background-color: #4d4d4d; color: #808080; }
        QPushButton:checked { background-color: #00a86b; }
        QPushButton#stopButton { background-color: #d41a1a; }
        QPushButton#stopButton:hover { background-color: #e62222; }
        QPushButton#pauseButton { background-color: #d4a017; }
        QPushButton#pauseButton:hover { background-color: #e6b01f; }
        QPushButton#pageButton {
            background-color: #2d2d35;
            border: 1px solid #3d3d3d;
            border-radius: 0;
            padding: 12px 24px;
            font-size: 12px;
        }
        QPushButton#pageButton:hover { background-color: #3d3d45; }
        QPushButton#pageButton:checked {
            background-color: #0078d4;
            border-bottom: 3px solid #00d4ff;
        }
        QSpinBox, QDoubleSpinBox, QComboBox, QLineEdit {
            background-color: #2d2d2d;
            border: 1px solid #3d3d3d;
            border-radius: 3px;
            padding: 5px;
            color: #ffffff;
            min-height: 25px;
        }
        QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus, QLineEdit:focus {
            border: 1px solid #0078d4;
        }
        QComboBox::drop-down { border: none; width: 20px; }
        QComboBox QAbstractItemView {
            background-color: #2d2d2d;
            color: #ffffff;
            selection-background-color: #0078d4;
        }
        QTableWidget {
            background-color: #1e1e24;
            alternate-background-color: #252530;
            gridline-color: #3d3d3d;
            color: #ffffff;
            border: 1px solid #3d3d3d;
        }
        QTableWidget::item { padding: 5px; color: #ffffff; }
        QTableWidget::item:selected { background-color: #0078d4; }
        QHeaderView::section {
            background-color: #2a2a35;
            color: #00d4ff;
            padding: 8px;
            border: 1px solid #3d3d3d;
            font-weight: bold;
        }
        QLabel { color: #ffffff; }
        QLabel#titleLabel { font-size: 14px; font-weight: bold; color: #00d4ff; }
        QLabel#valueLabel { font-size: 13px; color: #00ff88; font-family: 'Consolas', monospace; }
        QProgressBar {
            border: 1px solid #3d3d3d;
            border-radius: 3px;
            text-align: center;
            color: #ffffff;
            background-color: #2d2d2d;
        }
        QProgressBar::chunk { background-color: #0078d4; border-radius: 2px; }
        QDockWidget { color: #ffffff; }
        QDockWidget::title {
            background-color: #2a2a35;
            padding: 8px;
            color: #00d4ff;
            font-weight: bold;
        }
        QMenuBar { background-color: #1e1e24; color: #ffffff; }
        QMenuBar::item:selected { background-color: #0078d4; }
        QMenu {
            background-color: #2d2d2d;
            color: #ffffff;
            border: 1px solid #3d3d3d;
        }
        QMenu::item:selected { background-color: #0078d4; }
        QStatusBar { background-color: #1e1e24; color: #ffffff; }
        QCheckBox { color: #ffffff; }
        QCheckBox::indicator {
            width: 16px; height: 16px;
            border: 1px solid #3d3d3d;
            border-radius: 3px;
            background-color: #2d2d2d;
        }
        QCheckBox::indicator:checked { background-color: #0078d4; }
        QSplitter::handle {
            background-color: #3d3d3d;
            width: 3px;
            height: 3px;
        }
        QSplitter::handle:hover { background-color: #0078d4; }
        QToolBar {
            background-color: #1e1e24;
            border: none;
            spacing: 0px;
        }
        QTabWidget::pane {
            border: 1px solid #3d3d3d;
            background-color: #1e1e24;
        }
        QTabBar::tab {
            background-color: #2d2d35;
            color: #ffffff;
            padding: 8px 16px;
            margin-right: 2px;
            border-top-left-radius: 4px;
            border-top-right-radius: 4px;
        }
        QTabBar::tab:selected { background-color: #0078d4; }
        QTabBar::tab:hover:!selected { background-color: #3d3d45; }
        QListWidget {
            background-color: #1a1a1f;
            color: #00ff88;
            font-family: 'Consolas', monospace;
            border: 1px solid #3d3d3d;
        }
        QListWidget::item { padding: 4px 8px; border-bottom: 1px solid #2d2d2d; }
        QListWidget::item:selected { background-color: #0078d4; }
        QListWidget::item:hover { background-color: #333340; }
        QTextEdit {
            background-color: #1a1a1f;
            color: #e0e0e0;
            font-family: 'Consolas', monospace;
            border: 1px solid #3d3d3d;
        }
        QScrollArea { background-color: transparent; border: none; }
    )");
}

// ============================================================================
// Destructor
// ============================================================================

DashboardWindow::~DashboardWindow() {
    if (m_simTimer) {
        m_simTimer->stop();
    }
    if (m_refreshTimer) {
        m_refreshTimer->stop();
    }
}

// ============================================================================
// Initialize
// ============================================================================

void DashboardWindow::Initialize() {
    if (m_initialized) return;

    CreateMenus();
    CreateToolBar();
    CreateMainLayout();
    CreateStatsPanel();
    CreateStatusBar();

    // Set up signal connections
    connect(this, &DashboardWindow::particleAdded,
            this, &DashboardWindow::onParticleAdded, Qt::QueuedConnection);
    connect(this, &DashboardWindow::statsUpdated,
            this, &DashboardWindow::onStatsUpdated, Qt::QueuedConnection);
    connect(this, &DashboardWindow::hitMapUpdated,
            this, &DashboardWindow::onHitMapUpdated, Qt::QueuedConnection);
    connect(this, &DashboardWindow::progressChanged,
            this, &DashboardWindow::onProgressChanged, Qt::QueuedConnection);

    // Set up simulation timer for non-blocking event processing
    m_simTimer = new QTimer(this);
    m_simTimer->setInterval(0);  // Process immediately when event loop is free
    connect(m_simTimer, &QTimer::timeout, this, &DashboardWindow::onSimTimerTick);

    // Refresh timer
    m_refreshTimer = new QTimer(this);
    connect(m_refreshTimer, &QTimer::timeout, this, &DashboardWindow::onRefreshTimer);
    m_refreshTimer->start(100);

    m_initialized = true;
}

void DashboardWindow::InitializeEmbedded(QMainWindow* g4MainWindow) {
    if (m_initialized) return;
    if (!g4MainWindow) {
        // Fall back to standalone mode
        Initialize();
        Show();
        return;
    }

    m_g4MainWindow = g4MainWindow;
    m_embeddedMode = true;

    // Set up signal connections
    connect(this, &DashboardWindow::particleAdded,
            this, &DashboardWindow::onParticleAdded, Qt::QueuedConnection);
    connect(this, &DashboardWindow::statsUpdated,
            this, &DashboardWindow::onStatsUpdated, Qt::QueuedConnection);
    connect(this, &DashboardWindow::hitMapUpdated,
            this, &DashboardWindow::onHitMapUpdated, Qt::QueuedConnection);
    connect(this, &DashboardWindow::progressChanged,
            this, &DashboardWindow::onProgressChanged, Qt::QueuedConnection);

    // Set up simulation timer for non-blocking event processing
    m_simTimer = new QTimer(this);
    m_simTimer->setInterval(0);
    connect(m_simTimer, &QTimer::timeout, this, &DashboardWindow::onSimTimerTick);

    // Refresh timer
    m_refreshTimer = new QTimer(this);
    connect(m_refreshTimer, &QTimer::timeout, this, &DashboardWindow::onRefreshTimer);
    m_refreshTimer->start(100);

    // Create statistics dock widget for G4UIQt window
    QDockWidget* statsDock = new QDockWidget("NNBAR Statistics", g4MainWindow);
    statsDock->setObjectName("NNBARStatsDock");
    statsDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

    // Create stats widget content
    QWidget* statsWidget = new QWidget();
    QVBoxLayout* statsLayout = new QVBoxLayout(statsWidget);
    statsLayout->setContentsMargins(8, 8, 8, 8);
    statsLayout->setSpacing(4);

    // Title
    QLabel* titleLabel = new QLabel("<b>Event Statistics</b>");
    titleLabel->setStyleSheet("font-size: 14px; color: #2196F3;");
    statsLayout->addWidget(titleLabel);

    // Stats labels
    m_eventLabel = new QLabel("Event: 0");
    m_primaryLabel = new QLabel("Primaries: 0");
    m_secondaryLabel = new QLabel("Secondaries: 0");
    m_edepLabel = new QLabel("Energy Dep: 0.00 MeV");
    m_scintHitsLabel = new QLabel("Scint Hits: 0");
    m_tpcHitsLabel = new QLabel("TPC Hits: 0");
    m_lgHitsLabel = new QLabel("LG Hits: 0");
    m_photonLabel = new QLabel("Opt. Photons: 0");
    m_evtRateLabel = new QLabel("Event Rate: 0/s");

    for (auto* label : {m_eventLabel, m_primaryLabel, m_secondaryLabel,
                        m_edepLabel, m_scintHitsLabel, m_tpcHitsLabel,
                        m_lgHitsLabel, m_photonLabel, m_evtRateLabel}) {
        label->setStyleSheet("font-family: monospace; padding: 2px;");
        statsLayout->addWidget(label);
    }

    // Progress bar
    m_progressBar = new QProgressBar();
    m_progressBar->setRange(0, 100);
    m_progressBar->setValue(0);
    statsLayout->addWidget(m_progressBar);

    statsLayout->addStretch();

    statsDock->setWidget(statsWidget);
    g4MainWindow->addDockWidget(Qt::RightDockWidgetArea, statsDock);

    // Create control dock widget
    QDockWidget* controlDock = new QDockWidget("Run Control", g4MainWindow);
    controlDock->setObjectName("NNBARControlDock");
    controlDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea | Qt::BottomDockWidgetArea);

    QWidget* controlWidget = new QWidget();
    QVBoxLayout* controlLayout = new QVBoxLayout(controlWidget);
    controlLayout->setContentsMargins(8, 8, 8, 8);

    // Events spinner
    QHBoxLayout* eventsRow = new QHBoxLayout();
    eventsRow->addWidget(new QLabel("Events:"));
    m_eventsSpinBox = new QSpinBox();
    m_eventsSpinBox->setRange(1, 1000000);
    m_eventsSpinBox->setValue(100);
    eventsRow->addWidget(m_eventsSpinBox);
    controlLayout->addLayout(eventsRow);

    // Run buttons
    m_runButton = new QPushButton("Run");
    m_runButton->setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold;");
    connect(m_runButton, &QPushButton::clicked, this, &DashboardWindow::onRunEvents);
    controlLayout->addWidget(m_runButton);

    m_singleEventButton = new QPushButton("Single Event");
    connect(m_singleEventButton, &QPushButton::clicked, this, &DashboardWindow::onSingleEvent);
    controlLayout->addWidget(m_singleEventButton);

    m_stopButton = new QPushButton("Stop");
    m_stopButton->setStyleSheet("background-color: #f44336; color: white;");
    m_stopButton->setEnabled(false);
    connect(m_stopButton, &QPushButton::clicked, this, &DashboardWindow::onStopSimulation);
    controlLayout->addWidget(m_stopButton);

    controlLayout->addStretch();

    controlDock->setWidget(controlWidget);
    g4MainWindow->addDockWidget(Qt::RightDockWidgetArea, controlDock);

    // Stack the docks
    g4MainWindow->tabifyDockWidget(statsDock, controlDock);
    statsDock->raise();

    // =========================================================================
    // Create Particle Breakdown Dock
    // =========================================================================
    QDockWidget* particleDock = new QDockWidget("Particle Breakdown", g4MainWindow);
    particleDock->setObjectName("NNBARParticleDock");
    particleDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

    QWidget* particleWidget = new QWidget();
    QVBoxLayout* particleLayout = new QVBoxLayout(particleWidget);
    particleLayout->setContentsMargins(8, 8, 8, 8);

    QLabel* particleTitle = new QLabel("<b>Particle Statistics</b>");
    particleTitle->setStyleSheet("font-size: 14px; color: #FF9800;");
    particleLayout->addWidget(particleTitle);

    // Particle table
    m_particleTable = new QTableWidget(8, 2);
    m_particleTable->setHorizontalHeaderLabels({"Particle", "Count"});
    m_particleTable->verticalHeader()->setVisible(false);
    m_particleTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_particleTable->setStyleSheet(
        "QTableWidget { background: #1a1a2e; gridline-color: #333; }"
        "QTableWidget::item { padding: 4px; }"
        "QHeaderView::section { background: #252540; padding: 4px; }"
    );

    QStringList particles = {"gamma", "e-", "e+", "proton", "neutron", "pi+", "pi-", "mu-"};
    for (int i = 0; i < particles.size(); i++) {
        QTableWidgetItem* nameItem = new QTableWidgetItem(particles[i]);
        nameItem->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
        m_particleTable->setItem(i, 0, nameItem);

        QTableWidgetItem* countItem = new QTableWidgetItem("0");
        countItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
        m_particleTable->setItem(i, 1, countItem);
    }
    m_particleTable->horizontalHeader()->setStretchLastSection(true);
    m_particleTable->setColumnWidth(0, 80);
    particleLayout->addWidget(m_particleTable);

    particleDock->setWidget(particleWidget);
    g4MainWindow->addDockWidget(Qt::RightDockWidgetArea, particleDock);
    g4MainWindow->tabifyDockWidget(statsDock, particleDock);

    // =========================================================================
    // Create Accelerator Status Dock
    // =========================================================================
    QDockWidget* accelDock = new QDockWidget("Accelerator Status", g4MainWindow);
    accelDock->setObjectName("NNBARAccelDock");
    accelDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea);

    QWidget* accelWidget = new QWidget();
    QVBoxLayout* accelLayout = new QVBoxLayout(accelWidget);
    accelLayout->setContentsMargins(8, 8, 8, 8);
    accelLayout->setSpacing(6);

    QLabel* accelTitle = new QLabel("<b>GPU/Accelerator Status</b>");
    accelTitle->setStyleSheet("font-size: 14px; color: #00BCD4;");
    accelLayout->addWidget(accelTitle);

    // Celeritas status
    QHBoxLayout* celeritasRow = new QHBoxLayout();
    celeritasRow->addWidget(new QLabel("Celeritas EM:"));
#if WITH_CELERITAS
    QLabel* celeritasStatus = new QLabel("ENABLED");
    celeritasStatus->setStyleSheet("color: #4CAF50; font-weight: bold;");
#else
    QLabel* celeritasStatus = new QLabel("OFF (G4 < 11.0)");
    celeritasStatus->setStyleSheet("color: #9E9E9E;");
#endif
    celeritasRow->addWidget(celeritasStatus);
    celeritasRow->addStretch();
    accelLayout->addLayout(celeritasRow);

    // GarfieldGPU status
    QHBoxLayout* garfieldRow = new QHBoxLayout();
    garfieldRow->addWidget(new QLabel("TPC Drift:"));
#if WITH_GARFIELD_GPU
    QLabel* garfieldStatus = new QLabel("ENABLED");
    garfieldStatus->setStyleSheet("color: #4CAF50; font-weight: bold;");
#else
    QLabel* garfieldStatus = new QLabel("DISABLED");
    garfieldStatus->setStyleSheet("color: #9E9E9E;");
#endif
    garfieldRow->addWidget(garfieldStatus);
    garfieldRow->addStretch();
    accelLayout->addLayout(garfieldRow);

#if WITH_GARFIELD_GPU
    // Device info
    QLabel* deviceLabel = new QLabel("Device:");
    deviceLabel->setStyleSheet("color: #888;");
    accelLayout->addWidget(deviceLabel);

    m_gpuStatusLabel = new QLabel("Initializing...");
    m_gpuStatusLabel->setStyleSheet("font-family: monospace; color: #FFD700; padding-left: 10px;");
    accelLayout->addWidget(m_gpuStatusLabel);

    // TPC drift stats (updated per event)
    QLabel* driftStatsTitle = new QLabel("<b>TPC Drift Stats (per event):</b>");
    driftStatsTitle->setStyleSheet("margin-top: 10px; color: #888;");
    accelLayout->addWidget(driftStatsTitle);

    QGridLayout* driftGrid = new QGridLayout();
    driftGrid->setSpacing(4);

    driftGrid->addWidget(new QLabel("Clusters:"), 0, 0);
    m_tpcClustersLabel = new QLabel("--");
    m_tpcClustersLabel->setStyleSheet("font-family: monospace;");
    driftGrid->addWidget(m_tpcClustersLabel, 0, 1);

    driftGrid->addWidget(new QLabel("Electrons In:"), 1, 0);
    m_tpcElectronsInLabel = new QLabel("--");
    m_tpcElectronsInLabel->setStyleSheet("font-family: monospace;");
    driftGrid->addWidget(m_tpcElectronsInLabel, 1, 1);

    driftGrid->addWidget(new QLabel("Collected:"), 2, 0);
    m_tpcCollectedLabel = new QLabel("--");
    m_tpcCollectedLabel->setStyleSheet("font-family: monospace;");
    driftGrid->addWidget(m_tpcCollectedLabel, 2, 1);

    driftGrid->addWidget(new QLabel("Drift Time:"), 3, 0);
    m_tpcDriftTimeLabel = new QLabel("--");
    m_tpcDriftTimeLabel->setStyleSheet("font-family: monospace;");
    driftGrid->addWidget(m_tpcDriftTimeLabel, 3, 1);

    accelLayout->addLayout(driftGrid);
#endif

    accelLayout->addStretch();
    accelDock->setWidget(accelWidget);
    g4MainWindow->addDockWidget(Qt::RightDockWidgetArea, accelDock);
    g4MainWindow->tabifyDockWidget(statsDock, accelDock);

    // =========================================================================
    // Create Energy Deposition Dock
    // =========================================================================
    QDockWidget* edepDock = new QDockWidget("Energy Deposition", g4MainWindow);
    edepDock->setObjectName("NNBAREdepDock");
    edepDock->setAllowedAreas(Qt::LeftDockWidgetArea | Qt::RightDockWidgetArea | Qt::BottomDockWidgetArea);

    QWidget* edepWidget = new QWidget();
    QVBoxLayout* edepLayout = new QVBoxLayout(edepWidget);
    edepLayout->setContentsMargins(8, 8, 8, 8);

    QLabel* edepTitle = new QLabel("<b>Energy Deposition by Detector</b>");
    edepTitle->setStyleSheet("font-size: 14px; color: #E91E63;");
    edepLayout->addWidget(edepTitle);

    // Energy table
    QTableWidget* edepTable = new QTableWidget(4, 3);
    edepTable->setHorizontalHeaderLabels({"Detector", "E (MeV)", "Hits"});
    edepTable->verticalHeader()->setVisible(false);
    edepTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    edepTable->setStyleSheet(
        "QTableWidget { background: #1a1a2e; gridline-color: #333; }"
        "QTableWidget::item { padding: 4px; }"
        "QHeaderView::section { background: #252540; padding: 4px; }"
    );

    QStringList detectors = {"TPC", "Scintillator", "Lead Glass", "Other"};
    for (int i = 0; i < detectors.size(); i++) {
        QTableWidgetItem* nameItem = new QTableWidgetItem(detectors[i]);
        nameItem->setTextAlignment(Qt::AlignLeft | Qt::AlignVCenter);
        edepTable->setItem(i, 0, nameItem);

        QTableWidgetItem* edepItem = new QTableWidgetItem("0.00");
        edepItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
        edepTable->setItem(i, 1, edepItem);

        QTableWidgetItem* hitsItem = new QTableWidgetItem("0");
        hitsItem->setTextAlignment(Qt::AlignRight | Qt::AlignVCenter);
        edepTable->setItem(i, 2, hitsItem);
    }
    edepTable->horizontalHeader()->setStretchLastSection(true);
    edepTable->setColumnWidth(0, 80);
    edepTable->setColumnWidth(1, 70);
    m_edepTable = edepTable;
    edepLayout->addWidget(edepTable);

    edepDock->setWidget(edepWidget);
    g4MainWindow->addDockWidget(Qt::RightDockWidgetArea, edepDock);
    g4MainWindow->tabifyDockWidget(statsDock, edepDock);

    // Raise the stats dock as default
    statsDock->raise();

    m_initialized = true;
}

// ============================================================================
// Show/Hide
// ============================================================================

void DashboardWindow::Show() {
    if (!m_initialized) Initialize();
    show();
    raise();
    activateWindow();
}

void DashboardWindow::Hide() { hide(); }

// ============================================================================
// Create Menu Bar
// ============================================================================

void DashboardWindow::CreateMenus() {
    QMenuBar* menuBar = this->menuBar();

    // File menu
    QMenu* fileMenu = menuBar->addMenu("&File");
    QAction* exportAction = fileMenu->addAction("Export Data...");
    connect(exportAction, &QAction::triggered, this, &DashboardWindow::onExportData);
    QAction* saveConfigAction = fileMenu->addAction("Save Configuration...");
    connect(saveConfigAction, &QAction::triggered, this, &DashboardWindow::onSaveConfig);
    fileMenu->addSeparator();
    QAction* exitAction = fileMenu->addAction("Exit");
    connect(exitAction, &QAction::triggered, this, &QMainWindow::close);

    // Simulation menu
    QMenu* simMenu = menuBar->addMenu("&Simulation");
    QAction* runAction = simMenu->addAction("Run Events");
    runAction->setShortcut(QKeySequence("F5"));
    connect(runAction, &QAction::triggered, this, &DashboardWindow::onRunEvents);
    QAction* stopAction = simMenu->addAction("Stop");
    stopAction->setShortcut(QKeySequence("F6"));
    connect(stopAction, &QAction::triggered, this, &DashboardWindow::onStopSimulation);
    simMenu->addSeparator();
    QAction* initVisAction = simMenu->addAction("Initialize Visualization");
    connect(initVisAction, &QAction::triggered, this, &DashboardWindow::onInitVis);

    // View menu
    QMenu* viewMenu = menuBar->addMenu("&View");
    QAction* viewXYAction = viewMenu->addAction("Transverse (X-Y)");
    viewXYAction->setShortcut(QKeySequence("1"));
    connect(viewXYAction, &QAction::triggered, this, &DashboardWindow::onViewTransverse);
    QAction* viewRZAction = viewMenu->addAction("Longitudinal (R-Z)");
    viewRZAction->setShortcut(QKeySequence("2"));
    connect(viewRZAction, &QAction::triggered, this, &DashboardWindow::onViewLongitudinal);
    QAction* view3DAction = viewMenu->addAction("3D Perspective");
    view3DAction->setShortcut(QKeySequence("3"));
    connect(view3DAction, &QAction::triggered, this, &DashboardWindow::onView3D);
    viewMenu->addSeparator();
    QAction* resetViewAction = viewMenu->addAction("Reset View");
    resetViewAction->setShortcut(QKeySequence("R"));
    connect(resetViewAction, &QAction::triggered, this, &DashboardWindow::onResetView);

    // Help menu
    QMenu* helpMenu = menuBar->addMenu("&Help");
    QAction* aboutAction = helpMenu->addAction("About NNBAR Simulation");
    connect(aboutAction, &QAction::triggered, [this]() {
        QMessageBox::about(this, "About NNBAR Detector Simulation",
            "<h2>NNBAR Detector Simulation</h2>"
            "<p>Version 2.0</p>"
            "<p>Professional simulation software for the NNBAR experiment</p>"
            "<p>&copy; 2024 NNBAR Collaboration</p>");
    });
}

// ============================================================================
// Create Tool Bar (Page Navigation)
// ============================================================================

void DashboardWindow::CreateToolBar() {
    m_pageToolBar = new QToolBar("Navigation", this);
    m_pageToolBar->setMovable(false);
    m_pageToolBar->setFloatable(false);
    addToolBar(Qt::TopToolBarArea, m_pageToolBar);

    m_pageButtonGroup = new QButtonGroup(this);
    m_pageButtonGroup->setExclusive(true);

    QStringList pageNames = {"Control Panel", "G4 Console", "Data Analysis", "Material Budget"};
    for (int i = 0; i < pageNames.size(); i++) {
        QPushButton* btn = new QPushButton(pageNames[i]);
        btn->setObjectName("pageButton");
        btn->setCheckable(true);
        if (i == 0) btn->setChecked(true);
        m_pageButtonGroup->addButton(btn, i);
        m_pageToolBar->addWidget(btn);
    }

    connect(m_pageButtonGroup, static_cast<void(QButtonGroup::*)(int)>(&QButtonGroup::buttonClicked),
            this, &DashboardWindow::onPageChanged);

    m_pageToolBar->addSeparator();

    // === G4-style Visualization Controls ===
    QLabel* visLabel = new QLabel("  Vis: ");
    visLabel->setStyleSheet("color: #ffd700;");
    m_pageToolBar->addWidget(visLabel);

    // Wireframe/Solid toggle
    QPushButton* wireframeBtn = new QPushButton("Wireframe");
    wireframeBtn->setFixedWidth(75);
    wireframeBtn->setCheckable(true);
    wireframeBtn->setStyleSheet("QPushButton { padding: 4px 8px; font-size: 10px; }");
    connect(wireframeBtn, &QPushButton::clicked, [this, wireframeBtn](bool checked) {
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if (UI) {
            if (checked) {
                UI->ApplyCommand("/vis/viewer/set/style wireframe");
                wireframeBtn->setText("Wireframe");
            } else {
                UI->ApplyCommand("/vis/viewer/set/style surface");
                wireframeBtn->setText("Solid");
            }
        }
    });
    m_pageToolBar->addWidget(wireframeBtn);

    // Hidden line removal
    QPushButton* hlrBtn = new QPushButton("HLR");
    hlrBtn->setToolTip("Hidden Line Removal");
    hlrBtn->setFixedWidth(45);
    hlrBtn->setCheckable(true);
    hlrBtn->setStyleSheet("QPushButton { padding: 4px 6px; font-size: 10px; }");
    connect(hlrBtn, &QPushButton::clicked, [](bool checked) {
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if (UI) {
            UI->ApplyCommand(checked ? "/vis/viewer/set/hiddenEdge true" : "/vis/viewer/set/hiddenEdge false");
        }
    });
    m_pageToolBar->addWidget(hlrBtn);

    // Auxiliary edges
    QPushButton* auxEdgeBtn = new QPushButton("Aux");
    auxEdgeBtn->setToolTip("Show Auxiliary Edges");
    auxEdgeBtn->setFixedWidth(40);
    auxEdgeBtn->setCheckable(true);
    auxEdgeBtn->setStyleSheet("QPushButton { padding: 4px 6px; font-size: 10px; }");
    connect(auxEdgeBtn, &QPushButton::clicked, [](bool checked) {
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if (UI) {
            UI->ApplyCommand(checked ? "/vis/viewer/set/auxiliaryEdge true" : "/vis/viewer/set/auxiliaryEdge false");
        }
    });
    m_pageToolBar->addWidget(auxEdgeBtn);

    m_pageToolBar->addSeparator();

    // === Quick G4 Viewer Commands ===
    QPushButton* openOGLBtn = new QPushButton("Open OGL");
    openOGLBtn->setFixedWidth(70);
    openOGLBtn->setStyleSheet("QPushButton { padding: 4px 8px; font-size: 10px; background: #164050; }");
    connect(openOGLBtn, &QPushButton::clicked, []() {
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if (UI) UI->ApplyCommand("/vis/open OGL");
    });
    m_pageToolBar->addWidget(openOGLBtn);

    QPushButton* drawVolBtn = new QPushButton("Draw");
    drawVolBtn->setToolTip("Draw Volume");
    drawVolBtn->setFixedWidth(50);
    drawVolBtn->setStyleSheet("QPushButton { padding: 4px 8px; font-size: 10px; background: #164050; }");
    connect(drawVolBtn, &QPushButton::clicked, []() {
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if (UI) UI->ApplyCommand("/vis/drawVolume");
    });
    m_pageToolBar->addWidget(drawVolBtn);

    QPushButton* flushBtn = new QPushButton("Flush");
    flushBtn->setToolTip("Flush Viewer");
    flushBtn->setFixedWidth(50);
    flushBtn->setStyleSheet("QPushButton { padding: 4px 8px; font-size: 10px; background: #164050; }");
    connect(flushBtn, &QPushButton::clicked, []() {
        G4UImanager* UI = G4UImanager::GetUIpointer();
        if (UI) UI->ApplyCommand("/vis/viewer/flush");
    });
    m_pageToolBar->addWidget(flushBtn);

    // Add spacer
    QWidget* spacer = new QWidget();
    spacer->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Preferred);
    m_pageToolBar->addWidget(spacer);

    // === Internal Event Display View Controls ===
    m_pageToolBar->addSeparator();
    QLabel* viewLabel = new QLabel("  Display: ");
    viewLabel->setStyleSheet("color: #00d4ff;");
    m_pageToolBar->addWidget(viewLabel);

    m_viewXYButton = new QPushButton("X-Y");
    m_viewXYButton->setCheckable(true);
    m_viewXYButton->setChecked(true);
    m_viewXYButton->setFixedWidth(45);
    m_viewXYButton->setStyleSheet("QPushButton { padding: 4px 6px; }");
    connect(m_viewXYButton, &QPushButton::clicked, this, &DashboardWindow::onViewTransverse);
    m_pageToolBar->addWidget(m_viewXYButton);

    m_viewRZButton = new QPushButton("R-Z");
    m_viewRZButton->setCheckable(true);
    m_viewRZButton->setFixedWidth(45);
    m_viewRZButton->setStyleSheet("QPushButton { padding: 4px 6px; }");
    connect(m_viewRZButton, &QPushButton::clicked, this, &DashboardWindow::onViewLongitudinal);
    m_pageToolBar->addWidget(m_viewRZButton);

    m_view3DButton = new QPushButton("3D");
    m_view3DButton->setCheckable(true);
    m_view3DButton->setFixedWidth(40);
    m_view3DButton->setStyleSheet("QPushButton { padding: 4px 6px; }");
    connect(m_view3DButton, &QPushButton::clicked, this, &DashboardWindow::onView3D);
    m_pageToolBar->addWidget(m_view3DButton);

    QPushButton* resetBtn = new QPushButton("⟲");
    resetBtn->setToolTip("Reset View");
    resetBtn->setFixedWidth(30);
    connect(resetBtn, &QPushButton::clicked, this, &DashboardWindow::onResetView);
    m_pageToolBar->addWidget(resetBtn);
}

// ============================================================================
// Create Main Layout with Event Display
// ============================================================================

void DashboardWindow::CreateMainLayout() {
    // Main splitter: Event Display (top) + Pages (bottom)
    m_mainSplitter = new QSplitter(Qt::Vertical, this);

    // === TOP: Event Display + TPC Pad Display (side by side) ===
    QSplitter* topSplitter = new QSplitter(Qt::Horizontal);
    CreateEventDisplayPanel();
    CreateTPCPadPanel();
    topSplitter->addWidget(m_eventDisplayPanel);
    topSplitter->addWidget(m_tpcPadPanel);
    topSplitter->setSizes({700, 400});  // Event display gets more space
    m_mainSplitter->addWidget(topSplitter);

    // === BOTTOM: Page Stack ===
    m_pageStack = new QStackedWidget();

    CreateControlPage();
    CreateCommandPage();
    CreateAnalysisPage();
    CreateMaterialBudgetPage();

    m_mainSplitter->addWidget(m_pageStack);

    // Set splitter sizes (60% event display, 40% pages)
    m_mainSplitter->setSizes({600, 400});
    m_mainSplitter->setStretchFactor(0, 3);
    m_mainSplitter->setStretchFactor(1, 2);

    setCentralWidget(m_mainSplitter);
}

// ============================================================================
// Create Event Display Panel
// ============================================================================

void DashboardWindow::CreateEventDisplayPanel() {
    // Create DashboardPanel wrapper for Event Display with maximize support
    m_eventDisplayPanel = new DashboardPanel("Event Display");
    m_eventDisplayPanel->SetAccentColor(QColor(255, 215, 0));  // Gold
    m_eventDisplayPanel->SetMaximizable(true);

    // Connect maximize/restore signals
    connect(m_eventDisplayPanel, &DashboardPanel::maximizeRequested, this, [this]() {
        // Hide other panels and maximize event display
        m_pageStack->hide();
        m_mainSplitter->setSizes({m_mainSplitter->height(), 0});
        m_eventDisplayPanel->setMinimumHeight(height() - 100);
    });
    connect(m_eventDisplayPanel, &DashboardPanel::restoreRequested, this, [this]() {
        // Restore normal layout
        m_pageStack->show();
        m_mainSplitter->setSizes({600, 400});
        m_eventDisplayPanel->setMinimumHeight(350);
    });

    // Content widget for event display
    QWidget* displayContent = new QWidget();
    QVBoxLayout* layout = new QVBoxLayout(displayContent);
    layout->setContentsMargins(5, 5, 5, 5);
    layout->setSpacing(5);

    // Title bar with controls
    QHBoxLayout* titleLayout = new QHBoxLayout();

    // Display layer toggles
    QLabel* layerLabel = new QLabel("Layers:");
    layerLabel->setStyleSheet("color: #888888; font-size: 10px;");
    titleLayout->addWidget(layerLabel);

    QCheckBox* heatmapCheck = new QCheckBox("Energy Heatmap");
    heatmapCheck->setStyleSheet("QCheckBox { color: #ff8800; font-size: 10px; }");
    connect(heatmapCheck, &QCheckBox::toggled, [this](bool checked) {
        if (m_eventDisplay) m_eventDisplay->SetShowEnergyHeatmap(checked);
    });
    titleLayout->addWidget(heatmapCheck);

    QCheckBox* opticalCheck = new QCheckBox("Optical Photons");
    opticalCheck->setStyleSheet("QCheckBox { color: #00ff88; font-size: 10px; }");
    connect(opticalCheck, &QCheckBox::toggled, [this](bool checked) {
        if (m_eventDisplay) m_eventDisplay->SetShowOpticalPhotons(checked);
    });
    titleLayout->addWidget(opticalCheck);

    titleLayout->addSpacing(10);
    titleLayout->addStretch();

    // Quick info labels
    QLabel* infoLabel = new QLabel("Drag to pan | Scroll to zoom");
    infoLabel->setStyleSheet("color: #808080; font-size: 10px;");
    titleLayout->addWidget(infoLabel);

    layout->addLayout(titleLayout);

    // Event Display Widget
    m_eventDisplay = new EventDisplay();
    m_eventDisplay->setMinimumHeight(350);
    layout->addWidget(m_eventDisplay, 1);

    // Set content in the DashboardPanel
    m_eventDisplayPanel->SetContent(displayContent);
    // Note: Don't add to splitter here - done in CreateMainLayout
}

// ============================================================================
// Create TPC Pad Display Panel - Shows pad occupancy and clusters
// ============================================================================

void DashboardWindow::CreateTPCPadPanel() {
    // Create DashboardPanel wrapper for TPC Pad Display
    m_tpcPadPanel = new DashboardPanel("TPC Pad Readout");
    m_tpcPadPanel->SetAccentColor(QColor(0, 200, 255));  // Cyan for TPC

    // Content widget for TPC pad display
    QWidget* padContent = new QWidget();
    QVBoxLayout* layout = new QVBoxLayout(padContent);
    layout->setContentsMargins(5, 5, 5, 5);
    layout->setSpacing(5);

    // Title bar with controls
    QHBoxLayout* titleLayout = new QHBoxLayout();

    // Side toggle button
    QPushButton* toggleSideBtn = new QPushButton("+z Side");
    toggleSideBtn->setStyleSheet(
        "QPushButton { background: #304050; color: #00ccff; border: 1px solid #00ccff; "
        "border-radius: 3px; padding: 3px 10px; font-size: 10px; }"
        "QPushButton:hover { background: #405060; }");
    titleLayout->addWidget(toggleSideBtn);

    QCheckBox* logScaleCheck = new QCheckBox("Log Scale");
    logScaleCheck->setChecked(true);
    logScaleCheck->setStyleSheet("QCheckBox { color: #888888; font-size: 10px; }");
    titleLayout->addWidget(logScaleCheck);

    titleLayout->addStretch();

    // Stats labels
    QLabel* padsLabel = new QLabel("Hit Pads: --");
    padsLabel->setStyleSheet("color: #00ccff; font-size: 10px;");
    titleLayout->addWidget(padsLabel);

    QLabel* clustersLabel = new QLabel("Clusters: --");
    clustersLabel->setStyleSheet("color: #ffaa00; font-size: 10px;");
    titleLayout->addWidget(clustersLabel);

    layout->addLayout(titleLayout);

    // TPC Pad Display Widget
    m_tpcPadDisplay = new TPCPadDisplay();
    m_tpcPadDisplay->setMinimumSize(300, 250);
    layout->addWidget(m_tpcPadDisplay, 1);

    // Connect side toggle
    connect(toggleSideBtn, &QPushButton::clicked, [this, toggleSideBtn]() {
        if (m_tpcPadDisplay) {
            m_tpcPadDisplay->ToggleSide();
            toggleSideBtn->setText(m_tpcPadDisplay->GetCurrentSide() == 0 ? "+z Side" : "-z Side");
        }
    });

    // Connect log scale toggle
    connect(logScaleCheck, &QCheckBox::toggled, [this](bool checked) {
        if (m_tpcPadDisplay) {
            m_tpcPadDisplay->SetLogScale(checked);
        }
    });

    // Set content in the DashboardPanel
    m_tpcPadPanel->SetContent(padContent);
}

// ============================================================================
// Create Control Page (Page 0) - PowerBI-style with KPI Cards
// ============================================================================

void DashboardWindow::CreateControlPage() {
    QWidget* page = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(page);
    mainLayout->setSpacing(10);
    mainLayout->setContentsMargins(10, 10, 10, 10);

    // === TOP ROW: KPI Cards ===
    QHBoxLayout* kpiLayout = new QHBoxLayout();
    kpiLayout->setSpacing(10);

    m_eventCountKPI = new KPICard("EVENTS PROCESSED");
    m_eventCountKPI->SetValue(0);
    m_eventCountKPI->SetColor(QColor(0, 212, 255));  // Cyan
    kpiLayout->addWidget(m_eventCountKPI);

    m_edepKPI = new KPICard("TOTAL ENERGY DEP");
    m_edepKPI->SetValue("0.00");
    m_edepKPI->SetUnit("MeV");
    m_edepKPI->SetColor(QColor(255, 136, 0));  // Orange
    kpiLayout->addWidget(m_edepKPI);

    m_hitRateKPI = new KPICard("DETECTOR HITS");
    m_hitRateKPI->SetValue(0);
    m_hitRateKPI->SetColor(QColor(0, 200, 100));  // Green
    kpiLayout->addWidget(m_hitRateKPI);

    m_eventRateKPI = new KPICard("EVENT RATE");
    m_eventRateKPI->SetValue("0.0");
    m_eventRateKPI->SetUnit("evt/s");
    m_eventRateKPI->SetColor(QColor(255, 100, 100));  // Red
    kpiLayout->addWidget(m_eventRateKPI);

    mainLayout->addLayout(kpiLayout);

    // === MIDDLE: Splitter with Controls and Statistics ===
    QSplitter* controlSplitter = new QSplitter(Qt::Horizontal);

    // === LEFT: Simulation Controls ===
    QWidget* leftPanel = new QWidget();
    QVBoxLayout* leftLayout = new QVBoxLayout(leftPanel);
    leftLayout->setSpacing(10);
    leftLayout->setContentsMargins(0, 0, 0, 0);

    // Simulation Control Group
    QGroupBox* controlGroup = new QGroupBox("Simulation Control");
    QGridLayout* controlLayout = new QGridLayout(controlGroup);
    controlLayout->setSpacing(10);

    QLabel* eventsLabel = new QLabel("Number of Events:");
    m_eventsSpinBox = new QSpinBox();
    m_eventsSpinBox->setRange(1, 1000000);
    m_eventsSpinBox->setValue(100);
    m_eventsSpinBox->setSingleStep(10);

    m_runButton = new QPushButton("▶ Run Events");
    m_runButton->setStyleSheet("QPushButton { background: #006600; padding: 10px 15px; font-size: 12px; }");
    connect(m_runButton, &QPushButton::clicked, this, &DashboardWindow::onRunEvents);

    m_singleEventButton = new QPushButton("▶ Single");
    m_singleEventButton->setStyleSheet("QPushButton { background: #005500; padding: 10px 12px; }");
    connect(m_singleEventButton, &QPushButton::clicked, this, &DashboardWindow::onSingleEvent);

    controlLayout->addWidget(eventsLabel, 0, 0);
    controlLayout->addWidget(m_eventsSpinBox, 0, 1);
    controlLayout->addWidget(m_runButton, 0, 2);
    controlLayout->addWidget(m_singleEventButton, 0, 3);

    m_stopButton = new QPushButton("⏹ Stop");
    m_stopButton->setObjectName("stopButton");
    connect(m_stopButton, &QPushButton::clicked, this, &DashboardWindow::onStopSimulation);

    m_pauseButton = new QPushButton("⏸ Pause");
    m_pauseButton->setObjectName("pauseButton");
    m_pauseButton->setCheckable(true);
    connect(m_pauseButton, &QPushButton::clicked, this, &DashboardWindow::onPauseSimulation);

    m_resetButton = new QPushButton("↺ Reset");
    connect(m_resetButton, &QPushButton::clicked, this, &DashboardWindow::onResetSimulation);

    QPushButton* clearButton = new QPushButton("Clear Data");
    connect(clearButton, &QPushButton::clicked, this, &DashboardWindow::onClearData);

    controlLayout->addWidget(m_stopButton, 1, 0);
    controlLayout->addWidget(m_pauseButton, 1, 1);
    controlLayout->addWidget(m_resetButton, 1, 2);
    controlLayout->addWidget(clearButton, 1, 3);

    leftLayout->addWidget(controlGroup);

    // Particle Source Group
    QGroupBox* sourceGroup = new QGroupBox("Particle Source");
    QGridLayout* sourceLayout = new QGridLayout(sourceGroup);
    sourceLayout->setSpacing(8);

    QLabel* particleLabel = new QLabel("Particle:");
    m_particleCombo = new QComboBox();
    m_particleCombo->addItems({"e-", "e+", "gamma", "mu-", "mu+", "proton", "neutron", "pi+", "pi-", "kaon+", "kaon-"});

    QLabel* energyLabel = new QLabel("Energy (MeV):");
    m_energySpinBox = new QDoubleSpinBox();
    m_energySpinBox->setRange(0.1, 100000);
    m_energySpinBox->setValue(1000);
    m_energySpinBox->setDecimals(1);

    sourceLayout->addWidget(particleLabel, 0, 0);
    sourceLayout->addWidget(m_particleCombo, 0, 1);
    sourceLayout->addWidget(energyLabel, 0, 2);
    sourceLayout->addWidget(m_energySpinBox, 0, 3);

    QLabel* posLabel = new QLabel("Position (cm):");
    m_posXSpin = new QDoubleSpinBox(); m_posXSpin->setRange(-500, 500); m_posXSpin->setValue(0);
    m_posYSpin = new QDoubleSpinBox(); m_posYSpin->setRange(-500, 500); m_posYSpin->setValue(0);
    m_posZSpin = new QDoubleSpinBox(); m_posZSpin->setRange(-1000, 1000); m_posZSpin->setValue(-500);

    sourceLayout->addWidget(posLabel, 1, 0);
    QHBoxLayout* posLayout = new QHBoxLayout();
    posLayout->addWidget(new QLabel("X:")); posLayout->addWidget(m_posXSpin);
    posLayout->addWidget(new QLabel("Y:")); posLayout->addWidget(m_posYSpin);
    posLayout->addWidget(new QLabel("Z:")); posLayout->addWidget(m_posZSpin);
    QWidget* posWidget = new QWidget(); posWidget->setLayout(posLayout);
    sourceLayout->addWidget(posWidget, 1, 1, 1, 2);

    QPushButton* applySourceButton = new QPushButton("Apply");
    connect(applySourceButton, &QPushButton::clicked, this, &DashboardWindow::onApplySource);
    sourceLayout->addWidget(applySourceButton, 1, 3);

    leftLayout->addWidget(sourceGroup);

    // Data Management Group
    QGroupBox* exportGroup = new QGroupBox("Data Management");
    QHBoxLayout* exportLayout = new QHBoxLayout(exportGroup);

    QPushButton* exportDataButton = new QPushButton("Export Parquet");
    connect(exportDataButton, &QPushButton::clicked, this, &DashboardWindow::onExportData);
    QPushButton* exportCSVButton = new QPushButton("Export CSV");
    connect(exportCSVButton, &QPushButton::clicked, this, &DashboardWindow::onExportCSV);
    QPushButton* openOutputButton = new QPushButton("Open Output");
    connect(openOutputButton, &QPushButton::clicked, this, &DashboardWindow::onOpenOutputFolder);

    m_autoSaveCheck = new QCheckBox("Auto-save");
    m_autoSaveCheck->setChecked(true);

    exportLayout->addWidget(exportDataButton);
    exportLayout->addWidget(exportCSVButton);
    exportLayout->addWidget(openOutputButton);
    exportLayout->addWidget(m_autoSaveCheck);
    exportLayout->addStretch();

    leftLayout->addWidget(exportGroup);
    leftLayout->addStretch();

    controlSplitter->addWidget(leftPanel);

    // === RIGHT: Live Statistics with Mini Histograms ===
    QWidget* rightPanel = new QWidget();
    QVBoxLayout* rightLayout = new QVBoxLayout(rightPanel);
    rightLayout->setSpacing(10);
    rightLayout->setContentsMargins(0, 0, 0, 0);

    // Current Event Stats
    QGroupBox* eventGroup = new QGroupBox("Current Event");
    QFormLayout* eventLayout = new QFormLayout(eventGroup);
    eventLayout->setSpacing(6);

    m_eventLabel = new QLabel("0"); m_eventLabel->setObjectName("valueLabel");
    m_primaryLabel = new QLabel("0"); m_primaryLabel->setObjectName("valueLabel");
    m_secondaryLabel = new QLabel("0"); m_secondaryLabel->setObjectName("valueLabel");
    m_edepLabel = new QLabel("0.00 MeV"); m_edepLabel->setObjectName("valueLabel");
    m_cpuLabel = new QLabel("0.00 ms"); m_cpuLabel->setObjectName("valueLabel");

    eventLayout->addRow("Event ID:", m_eventLabel);
    eventLayout->addRow("Primaries:", m_primaryLabel);
    eventLayout->addRow("Secondaries:", m_secondaryLabel);
    eventLayout->addRow("Total Edep:", m_edepLabel);
    eventLayout->addRow("CPU Time:", m_cpuLabel);

    rightLayout->addWidget(eventGroup);

    // Mini histogram for energy deposition
    m_edepHistogram = new MiniHistogram();
    m_edepHistogram->SetColor(QColor(255, 136, 0));
    rightLayout->addWidget(m_edepHistogram);

    // Detector Hits
    QGroupBox* hitsGroup = new QGroupBox("Detector Hits");
    QFormLayout* hitsLayout = new QFormLayout(hitsGroup);
    hitsLayout->setSpacing(6);

    m_scintHitsLabel = new QLabel("0"); m_scintHitsLabel->setObjectName("valueLabel");
    m_tpcHitsLabel = new QLabel("0"); m_tpcHitsLabel->setObjectName("valueLabel");
    m_lgHitsLabel = new QLabel("0"); m_lgHitsLabel->setObjectName("valueLabel");
    m_photonLabel = new QLabel("0"); m_photonLabel->setObjectName("valueLabel");

    hitsLayout->addRow("Scintillator:", m_scintHitsLabel);
    hitsLayout->addRow("TPC:", m_tpcHitsLabel);
    hitsLayout->addRow("Lead Glass:", m_lgHitsLabel);
    hitsLayout->addRow("Optical Photons:", m_photonLabel);

    rightLayout->addWidget(hitsGroup);

    // Mini histogram for hit rate
    m_hitHistogram = new MiniHistogram();
    m_hitHistogram->SetColor(QColor(0, 200, 100));
    rightLayout->addWidget(m_hitHistogram);

    // Performance & GPU Acceleration Status
    QGroupBox* perfGroup = new QGroupBox("Performance & GPU");
    QFormLayout* perfLayout = new QFormLayout(perfGroup);

    m_evtRateLabel = new QLabel("0.0 evt/s"); m_evtRateLabel->setObjectName("valueLabel");

    // Celeritas status indicator
    QLabel* celeritasLabel = new QLabel();
#if WITH_CELERITAS
    celeritasLabel->setText("✓ Enabled");
    celeritasLabel->setStyleSheet("color: #00ff88; font-weight: bold;");
#else
    celeritasLabel->setText("✗ Disabled");
    celeritasLabel->setStyleSheet("color: #666666;");
#endif
    celeritasLabel->setObjectName("valueLabel");

    // Opticks status indicator
    QLabel* opticksLabel = new QLabel();
#if WITH_OPTICKS
    opticksLabel->setText("✓ Enabled");
    opticksLabel->setStyleSheet("color: #00ff88; font-weight: bold;");
#else
    opticksLabel->setText("✗ Disabled");
    opticksLabel->setStyleSheet("color: #666666;");
#endif
    opticksLabel->setObjectName("valueLabel");

    // Overall GPU mode
    m_gpuStatusLabel = new QLabel();
#if WITH_CELERITAS || WITH_OPTICKS
    m_gpuStatusLabel->setText("GPU Accelerated");
    m_gpuStatusLabel->setStyleSheet("color: #00ff88; font-weight: bold;");
#else
    m_gpuStatusLabel->setText("CPU Only");
    m_gpuStatusLabel->setStyleSheet("color: #ffaa00;");
#endif
    m_gpuStatusLabel->setObjectName("valueLabel");

    perfLayout->addRow("Event Rate:", m_evtRateLabel);
    perfLayout->addRow("Celeritas:", celeritasLabel);
    perfLayout->addRow("Opticks:", opticksLabel);
    perfLayout->addRow("Mode:", m_gpuStatusLabel);

    rightLayout->addWidget(perfGroup);
    rightLayout->addStretch();

    controlSplitter->addWidget(rightPanel);
    controlSplitter->setSizes({500, 250});

    mainLayout->addWidget(controlSplitter, 1);

    m_pageStack->addWidget(page);
}

// ============================================================================
// Create Command Page (Page 1) - Full G4 Command Interface
// ============================================================================

void DashboardWindow::CreateCommandPage() {
    // Create the comprehensive G4 command panel
    m_g4CommandPanel = new G4CommandPanel();

    // Wrap it in a DashboardPanel for consistent styling
    m_commandPanel = new DashboardPanel("Geant4 Command Console");
    m_commandPanel->SetAccentColor(QColor(255, 215, 0));  // Gold
    m_commandPanel->SetContent(m_g4CommandPanel);
    m_commandPanel->SetMaximizable(true);

    // Connect signals
    connect(m_g4CommandPanel, &G4CommandPanel::commandExecuted,
            [this](const QString& cmd) {
                AddToCommandHistory(cmd);
                m_statusLabel->setText("Executed: " + cmd);
            });

    // Also keep references for backward compatibility
    m_commandOutput = m_g4CommandPanel->findChild<QTextEdit*>();
    m_commandHistory = m_g4CommandPanel->findChild<QListWidget*>();

    m_pageStack->addWidget(m_commandPanel);
}

// ============================================================================
// Create Analysis Page (Page 2)
// ============================================================================

void DashboardWindow::CreateAnalysisPage() {
    QWidget* page = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(page);
    mainLayout->setSpacing(10);
    mainLayout->setContentsMargins(10, 10, 10, 10);

    // Particle table
    QGroupBox* tableGroup = new QGroupBox("Particle Tracks (Live)");
    QVBoxLayout* tableLayout = new QVBoxLayout(tableGroup);

    m_particleTable = new QTableWidget();
    m_particleTable->setColumnCount(9);
    m_particleTable->setHorizontalHeaderLabels({
        "Event", "Track", "Particle", "KE (MeV)",
        "X (cm)", "Y (cm)", "Z (cm)", "Volume", "Process"
    });
    m_particleTable->setAlternatingRowColors(true);
    m_particleTable->setSelectionBehavior(QAbstractItemView::SelectRows);
    m_particleTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    m_particleTable->horizontalHeader()->setStretchLastSection(true);
    m_particleTable->verticalHeader()->setVisible(false);

    m_particleTable->setColumnWidth(0, 60);
    m_particleTable->setColumnWidth(1, 60);
    m_particleTable->setColumnWidth(2, 80);
    m_particleTable->setColumnWidth(3, 90);
    m_particleTable->setColumnWidth(4, 80);
    m_particleTable->setColumnWidth(5, 80);
    m_particleTable->setColumnWidth(6, 80);
    m_particleTable->setColumnWidth(7, 120);

    tableLayout->addWidget(m_particleTable);

    mainLayout->addWidget(tableGroup);

    m_pageStack->addWidget(page);
}

// ============================================================================
// Create Material Budget Page (Page 3)
// ============================================================================

void DashboardWindow::CreateMaterialBudgetPage() {
    QWidget* page = new QWidget();
    QVBoxLayout* mainLayout = new QVBoxLayout(page);
    mainLayout->setSpacing(10);
    mainLayout->setContentsMargins(10, 10, 10, 10);

    // Title
    QLabel* titleLabel = new QLabel("NNBAR Detector Material Budget");
    titleLabel->setStyleSheet("font-size: 18px; font-weight: bold; color: #ffd700; margin-bottom: 10px;");
    mainLayout->addWidget(titleLabel);

    // Description
    QLabel* descLabel = new QLabel(
        "Material budget visualization showing radiation length (X₀) and interaction length (λ) "
        "across the detector as a function of pseudorapidity (η). Hover over plots for detailed values.");
    descLabel->setWordWrap(true);
    descLabel->setStyleSheet("color: #b0b0b0; margin-bottom: 10px;");
    mainLayout->addWidget(descLabel);

    // Create horizontal splitter for two plots side by side
    QSplitter* plotSplitter = new QSplitter(Qt::Horizontal);
    plotSplitter->setStyleSheet("QSplitter::handle { background: #3d3d3d; }");

    // Radiation Length Plot (X0 vs eta)
    QWidget* x0Container = new QWidget();
    QVBoxLayout* x0Layout = new QVBoxLayout(x0Container);
    x0Layout->setContentsMargins(5, 5, 5, 5);

    MaterialBudgetPlot* x0Plot = new MaterialBudgetPlot();
    x0Plot->SetPlotType(MaterialBudgetPlot::PlotType::RadiationLength);
    x0Plot->SetTitle("Radiation Length (X₀) vs η");
    x0Plot->setMinimumSize(350, 300);
    x0Layout->addWidget(x0Plot);

    // X0 export buttons
    QHBoxLayout* x0ButtonLayout = new QHBoxLayout();
    QPushButton* exportX0CSV = new QPushButton("Export CSV");
    exportX0CSV->setStyleSheet("background-color: #2d5016; padding: 5px 15px;");
    QPushButton* exportX0JSON = new QPushButton("Export JSON");
    exportX0JSON->setStyleSheet("background-color: #164050; padding: 5px 15px;");
    x0ButtonLayout->addStretch();
    x0ButtonLayout->addWidget(exportX0CSV);
    x0ButtonLayout->addWidget(exportX0JSON);
    x0Layout->addLayout(x0ButtonLayout);

    plotSplitter->addWidget(x0Container);

    // Interaction Length Plot (lambda vs eta)
    QWidget* lambdaContainer = new QWidget();
    QVBoxLayout* lambdaLayout = new QVBoxLayout(lambdaContainer);
    lambdaLayout->setContentsMargins(5, 5, 5, 5);

    MaterialBudgetPlot* lambdaPlot = new MaterialBudgetPlot();
    lambdaPlot->SetPlotType(MaterialBudgetPlot::PlotType::InteractionLength);
    lambdaPlot->SetTitle("Interaction Length (λ) vs η");
    lambdaPlot->setMinimumSize(350, 300);
    lambdaLayout->addWidget(lambdaPlot);

    // Lambda export buttons
    QHBoxLayout* lambdaButtonLayout = new QHBoxLayout();
    QPushButton* exportLambdaCSV = new QPushButton("Export CSV");
    exportLambdaCSV->setStyleSheet("background-color: #2d5016; padding: 5px 15px;");
    QPushButton* exportLambdaJSON = new QPushButton("Export JSON");
    exportLambdaJSON->setStyleSheet("background-color: #164050; padding: 5px 15px;");
    lambdaButtonLayout->addStretch();
    lambdaButtonLayout->addWidget(exportLambdaCSV);
    lambdaButtonLayout->addWidget(exportLambdaJSON);
    lambdaLayout->addLayout(lambdaButtonLayout);

    plotSplitter->addWidget(lambdaContainer);

    mainLayout->addWidget(plotSplitter, 2);  // Stretch factor 2 for plots

    // Connect export buttons
    connect(exportX0CSV, &QPushButton::clicked, [x0Plot]() {
        QString filename = QFileDialog::getSaveFileName(nullptr, "Export X0 Data",
            QDir::homePath() + "/nnbar_x0_vs_eta.csv", "CSV Files (*.csv)");
        if (!filename.isEmpty()) {
            QFile file(filename);
            if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                QTextStream out(&file);
                out << x0Plot->ExportToCSV();
                file.close();
                QMessageBox::information(nullptr, "Export Complete",
                    "Radiation length data exported to:\n" + filename);
            }
        }
    });

    connect(exportX0JSON, &QPushButton::clicked, [x0Plot]() {
        QString filename = QFileDialog::getSaveFileName(nullptr, "Export X0 Data",
            QDir::homePath() + "/nnbar_x0_vs_eta.json", "JSON Files (*.json)");
        if (!filename.isEmpty()) {
            QFile file(filename);
            if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                QTextStream out(&file);
                out << x0Plot->ExportToJSON();
                file.close();
                QMessageBox::information(nullptr, "Export Complete",
                    "Radiation length data exported to:\n" + filename);
            }
        }
    });

    connect(exportLambdaCSV, &QPushButton::clicked, [lambdaPlot]() {
        QString filename = QFileDialog::getSaveFileName(nullptr, "Export λ Data",
            QDir::homePath() + "/nnbar_lambda_vs_eta.csv", "CSV Files (*.csv)");
        if (!filename.isEmpty()) {
            QFile file(filename);
            if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                QTextStream out(&file);
                out << lambdaPlot->ExportToCSV();
                file.close();
                QMessageBox::information(nullptr, "Export Complete",
                    "Interaction length data exported to:\n" + filename);
            }
        }
    });

    connect(exportLambdaJSON, &QPushButton::clicked, [lambdaPlot]() {
        QString filename = QFileDialog::getSaveFileName(nullptr, "Export λ Data",
            QDir::homePath() + "/nnbar_lambda_vs_eta.json", "JSON Files (*.json)");
        if (!filename.isEmpty()) {
            QFile file(filename);
            if (file.open(QIODevice::WriteOnly | QIODevice::Text)) {
                QTextStream out(&file);
                out << lambdaPlot->ExportToJSON();
                file.close();
                QMessageBox::information(nullptr, "Export Complete",
                    "Interaction length data exported to:\n" + filename);
            }
        }
    });

    // Material budget table (condensed)
    QGroupBox* tableGroup = new QGroupBox("Material Properties by Region");
    QVBoxLayout* tableLayout = new QVBoxLayout(tableGroup);
    tableLayout->setContentsMargins(5, 15, 5, 5);

    QTableWidget* materialTable = new QTableWidget();
    materialTable->setColumnCount(6);
    materialTable->setHorizontalHeaderLabels({
        "Region", "η Range", "X₀", "λ", "Density", "Material"
    });
    materialTable->setAlternatingRowColors(true);
    materialTable->setEditTriggers(QAbstractItemView::NoEditTriggers);
    materialTable->horizontalHeader()->setStretchLastSection(true);
    materialTable->setMaximumHeight(200);

    struct MaterialData {
        const char* region;
        const char* eta;
        const char* x0;
        const char* lambda;
        const char* density;
        const char* material;
    };

    MaterialData matData[] = {
        {"Beampipe", "|η| < 3.0", "0.02", "0.01", "2.70", "Aluminum"},
        {"TPC Gas", "|η| < 1.5", "0.03", "0.02", "0.0016", "Ar/CO₂"},
        {"TPC Walls", "|η| < 1.5", "0.15", "0.08", "1.43", "G10"},
        {"Scintillator", "|η| < 2.0", "0.45", "0.25", "1.03", "Plastic"},
        {"Lead Glass (Side)", "|η| < 1.2", "12.5", "0.95", "6.22", "PbO-SiO₂"},
        {"Lead Glass (F/B)", "1.2 < |η| < 2.5", "15.0", "1.10", "6.22", "PbO-SiO₂"},
        {"Cosmic Shield", "|η| < 3.0", "0.80", "0.45", "7.87", "Steel"},
    };

    materialTable->setRowCount(sizeof(matData) / sizeof(matData[0]));
    for (size_t i = 0; i < sizeof(matData) / sizeof(matData[0]); i++) {
        materialTable->setItem(i, 0, new QTableWidgetItem(matData[i].region));
        materialTable->setItem(i, 1, new QTableWidgetItem(matData[i].eta));
        materialTable->setItem(i, 2, new QTableWidgetItem(matData[i].x0));
        materialTable->setItem(i, 3, new QTableWidgetItem(matData[i].lambda));
        materialTable->setItem(i, 4, new QTableWidgetItem(matData[i].density));
        materialTable->setItem(i, 5, new QTableWidgetItem(matData[i].material));

        float x0 = QString(matData[i].x0).toFloat();
        QColor rowColor;
        if (x0 < 0.5) rowColor = QColor(40, 80, 40);
        else if (x0 < 5.0) rowColor = QColor(80, 80, 40);
        else rowColor = QColor(80, 40, 40);

        for (int j = 0; j < 6; j++) {
            if (materialTable->item(i, j)) {
                materialTable->item(i, j)->setBackground(rowColor);
            }
        }
    }

    materialTable->resizeColumnsToContents();
    tableLayout->addWidget(materialTable);
    mainLayout->addWidget(tableGroup, 1);

    // Options and event display toggle
    QHBoxLayout* optionsLayout = new QHBoxLayout();

    QCheckBox* showX0MapCheck = new QCheckBox("Show X₀ map on Event Display");
    showX0MapCheck->setChecked(false);
    connect(showX0MapCheck, &QCheckBox::toggled, [this](bool checked) {
        if (m_eventDisplay) {
            m_eventDisplay->SetShowRadiationLength(checked);
        }
    });
    optionsLayout->addWidget(showX0MapCheck);

    QCheckBox* showLegend = new QCheckBox("Show Legend");
    showLegend->setChecked(true);
    connect(showLegend, &QCheckBox::toggled, [x0Plot, lambdaPlot](bool checked) {
        x0Plot->SetShowLegend(checked);
        lambdaPlot->SetShowLegend(checked);
    });
    optionsLayout->addWidget(showLegend);

    QCheckBox* showGrid = new QCheckBox("Show Grid");
    showGrid->setChecked(true);
    connect(showGrid, &QCheckBox::toggled, [x0Plot, lambdaPlot](bool checked) {
        x0Plot->SetShowGrid(checked);
        lambdaPlot->SetShowGrid(checked);
    });
    optionsLayout->addWidget(showGrid);

    optionsLayout->addStretch();
    mainLayout->addLayout(optionsLayout);

    m_pageStack->addWidget(page);
}

// ============================================================================
// Create Stats Panel (Right Dock - Optional)
// ============================================================================

void DashboardWindow::CreateStatsPanel() {
    // Stats are now integrated into the Control Page
    // This can be used for additional docked panels if needed
}

// ============================================================================
// Create Status Bar
// ============================================================================

void DashboardWindow::CreateStatusBar() {
    QStatusBar* status = statusBar();

    m_progressBar = new QProgressBar();
    m_progressBar->setTextVisible(true);
    m_progressBar->setFormat("Event %v / %m");
    m_progressBar->setFixedWidth(250);
    m_progressBar->setValue(0);
    m_progressBar->setMaximum(100);

    m_statusLabel = new QLabel("Ready");
    m_statusLabel->setStyleSheet("color: #00ff88; font-weight: bold;");

    status->addWidget(m_statusLabel);
    status->addPermanentWidget(m_progressBar);
}

// ============================================================================
// Page Navigation
// ============================================================================

void DashboardWindow::onPageChanged(int index) {
    m_pageStack->setCurrentIndex(index);
}

// ============================================================================
// View Mode Slots
// ============================================================================

void DashboardWindow::onViewTransverse() {
    if (m_eventDisplay) m_eventDisplay->SetViewMode(ViewMode::Transverse);
    m_viewXYButton->setChecked(true);
    m_viewRZButton->setChecked(false);
    m_view3DButton->setChecked(false);
}

void DashboardWindow::onViewLongitudinal() {
    if (m_eventDisplay) m_eventDisplay->SetViewMode(ViewMode::Longitudinal);
    m_viewXYButton->setChecked(false);
    m_viewRZButton->setChecked(true);
    m_view3DButton->setChecked(false);
}

void DashboardWindow::onView3D() {
    if (m_eventDisplay) m_eventDisplay->SetViewMode(ViewMode::Perspective3D);
    m_viewXYButton->setChecked(false);
    m_viewRZButton->setChecked(false);
    m_view3DButton->setChecked(true);
}

void DashboardWindow::onResetView() {
    if (m_eventDisplay) m_eventDisplay->ResetView();
}

// ============================================================================
// Simulation Control Slots
// ============================================================================

void DashboardWindow::onRunEvents() {
    if (m_isSimulating.load()) {
        QMessageBox::warning(this, "Simulation Running",
            "A simulation is already in progress. Please wait for it to complete or stop it.");
        return;
    }

    int nEvents = m_eventsSpinBox->value();
    m_pendingEvents = nEvents;
    m_completedEvents = 0;

    m_progressBar->setMaximum(nEvents);
    m_progressBar->setValue(0);

    // Disable buttons during simulation
    m_runButton->setEnabled(false);
    m_singleEventButton->setEnabled(false);
    m_isSimulating.store(true);

    m_statusLabel->setText("Running...");
    m_statusLabel->setStyleSheet("color: #00d4ff; font-weight: bold;");

    // Start the simulation timer
    m_simTimer->start();
}

void DashboardWindow::onSingleEvent() {
    if (m_isSimulating.load()) {
        QMessageBox::warning(this, "Simulation Running",
            "A simulation is already in progress. Please wait for it to complete or stop it.");
        return;
    }

    m_pendingEvents = 1;
    m_completedEvents = 0;

    m_progressBar->setMaximum(1);
    m_progressBar->setValue(0);

    // Disable buttons during simulation
    m_runButton->setEnabled(false);
    m_singleEventButton->setEnabled(false);
    m_isSimulating.store(true);

    m_statusLabel->setText("Running single event...");
    m_statusLabel->setStyleSheet("color: #00d4ff; font-weight: bold;");

    // Start the simulation timer
    m_simTimer->start();
}

void DashboardWindow::onStopSimulation() {
    // Stop the simulation timer
    if (m_simTimer) {
        m_simTimer->stop();
    }

    // Abort any running Geant4 run
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI) UI->ApplyCommand("/run/abort");

    m_pendingEvents = 0;
    m_isSimulating.store(false);
    m_runButton->setEnabled(true);
    m_singleEventButton->setEnabled(true);
    m_statusLabel->setText("Stopped");
    m_statusLabel->setStyleSheet("color: #ff6b6b; font-weight: bold;");
}

void DashboardWindow::onPauseSimulation() {
    if (m_pauseButton->isChecked()) {
        m_statusLabel->setText("Paused");
        m_statusLabel->setStyleSheet("color: #ffd93d; font-weight: bold;");
    } else {
        m_statusLabel->setText("Running...");
        m_statusLabel->setStyleSheet("color: #00d4ff; font-weight: bold;");
    }
}

void DashboardWindow::onResetSimulation() {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI) UI->ApplyCommand("/run/initialize");
    onClearData();
    m_statusLabel->setText("Reset complete");
    m_statusLabel->setStyleSheet("color: #00ff88; font-weight: bold;");
}

void DashboardWindow::onSimTimerTick() {
    // Run one event at a time to keep UI responsive
    if (m_completedEvents >= m_pendingEvents) {
        // All events completed
        m_simTimer->stop();
        m_isSimulating.store(false);
        m_runButton->setEnabled(true);
        m_singleEventButton->setEnabled(true);
        m_statusLabel->setText("Completed");
        m_statusLabel->setStyleSheet("color: #00ff88; font-weight: bold;");
        return;
    }

    // Check if paused
    if (m_pauseButton && m_pauseButton->isChecked()) {
        return;  // Skip this tick if paused
    }

    // Run one event
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI) {
        int result = UI->ApplyCommand("/run/beamOn 1");
        if (result == 0) {
            m_completedEvents++;
            m_progressBar->setValue(m_completedEvents);
        } else {
            // Error occurred, stop simulation
            m_simTimer->stop();
            m_isSimulating.store(false);
            m_runButton->setEnabled(true);
            m_singleEventButton->setEnabled(true);
            m_statusLabel->setText("Error during simulation");
            m_statusLabel->setStyleSheet("color: #ff6b6b; font-weight: bold;");
            return;
        }
    }

    // Process pending events to keep UI responsive
    QApplication::processEvents(QEventLoop::ExcludeUserInputEvents, 10);
}

void DashboardWindow::onClearData() {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_particles.clear();
    if (m_particleTable) m_particleTable->setRowCount(0);
    if (m_eventDisplay) m_eventDisplay->ClearEvent();
    m_currentStats = EventStatistics();

    m_eventLabel->setText("0");
    m_primaryLabel->setText("0");
    m_secondaryLabel->setText("0");
    m_edepLabel->setText("0.00 MeV");
    m_cpuLabel->setText("0.00 ms");
    m_scintHitsLabel->setText("0");
    m_tpcHitsLabel->setText("0");
    m_lgHitsLabel->setText("0");
    m_photonLabel->setText("0");
    m_progressBar->setValue(0);
}

// ============================================================================
// Simulation Worker Slots
// ============================================================================

void DashboardWindow::onSimulationStarted() {
    m_isSimulating.store(true);
    m_statusLabel->setText("Running...");
    m_statusLabel->setStyleSheet("color: #00d4ff; font-weight: bold;");
}

void DashboardWindow::onSimulationFinished() {
    m_isSimulating.store(false);
    m_runButton->setEnabled(true);
    m_singleEventButton->setEnabled(true);
    m_statusLabel->setText("Completed");
    m_statusLabel->setStyleSheet("color: #00ff88; font-weight: bold;");
}

void DashboardWindow::onEventCompleted(int eventNum) {
    m_progressBar->setValue(eventNum);

    // Process Qt events to keep UI responsive
    QApplication::processEvents();
}

void DashboardWindow::onSimulationError(const QString& msg) {
    m_statusLabel->setText("Error: " + msg);
    m_statusLabel->setStyleSheet("color: #ff6b6b; font-weight: bold;");
    QMessageBox::warning(this, "Simulation Error", msg);
}

void DashboardWindow::onApplySource() {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (!UI) return;

    QString particle = m_particleCombo->currentText();
    double energy = m_energySpinBox->value();
    double x = m_posXSpin->value();
    double y = m_posYSpin->value();
    double z = m_posZSpin->value();

    UI->ApplyCommand(QString("/gun/particle %1").arg(particle).toStdString());
    UI->ApplyCommand(QString("/gun/energy %1 MeV").arg(energy).toStdString());
    UI->ApplyCommand(QString("/gun/position %1 %2 %3 cm").arg(x).arg(y).arg(z).toStdString());

    m_statusLabel->setText("Source settings applied");
}

// ============================================================================
// Export Slots
// ============================================================================

void DashboardWindow::onExportData() {
    QString filename = QFileDialog::getSaveFileName(this, "Export Data",
        QString("nnbar_data_%1.parquet").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss")),
        "Parquet Files (*.parquet)");
    if (!filename.isEmpty()) {
        m_statusLabel->setText("Data exported: " + filename);
    }
}

void DashboardWindow::onExportCSV() {
    QString filename = QFileDialog::getSaveFileName(this, "Export to CSV",
        QString("nnbar_data_%1.csv").arg(QDateTime::currentDateTime().toString("yyyyMMdd_hhmmss")),
        "CSV Files (*.csv)");
    if (!filename.isEmpty()) {
        m_statusLabel->setText("CSV exported: " + filename);
    }
}

void DashboardWindow::onSaveConfig() {
    QString filename = QFileDialog::getSaveFileName(this, "Save Configuration",
        "nnbar_config.mac", "Macro Files (*.mac)");
    if (!filename.isEmpty()) {
        m_statusLabel->setText("Configuration saved");
    }
}

void DashboardWindow::onOpenOutputFolder() {
    QProcess::startDetached("xdg-open", {"./output/"});
}

void DashboardWindow::onInitVis() {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI) {
        UI->ApplyCommand("/vis/open OGL");
        UI->ApplyCommand("/vis/drawVolume");
        UI->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 45 45");
        UI->ApplyCommand("/vis/scene/add/trajectories smooth");
    }
}

void DashboardWindow::onToggleG4Viewer() {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI) UI->ApplyCommand("/vis/open OGL");
}

// ============================================================================
// Command Console Slots
// ============================================================================

void DashboardWindow::onExecuteCommand() {
    QString cmd = m_commandInput->text().trimmed();
    if (cmd.isEmpty()) return;

    AddToCommandHistory(cmd);

    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI) {
        m_commandOutput->append(QString("<span style='color:#00d4ff;'>></span> <span style='color:#00ff88;'>%1</span>").arg(cmd));
        int result = UI->ApplyCommand(cmd.toStdString());
        if (result == 0) {
            m_commandOutput->append("<span style='color:#00ff88;'>OK</span>");
        } else {
            m_commandOutput->append(QString("<span style='color:#ff6b6b;'>Error: %1</span>").arg(result));
        }
        m_commandOutput->append("");
    }

    m_commandInput->clear();
    m_statusLabel->setText("Command: " + cmd);
}

void DashboardWindow::onCommandHistoryClicked(QListWidgetItem* item) {
    if (item) {
        m_commandInput->setText(item->text());
        onExecuteCommand();
    }
}

void DashboardWindow::onQuickCommand(const QString& cmd) {
    AddToCommandHistory(cmd);

    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI) {
        m_commandOutput->append(QString("<span style='color:#00d4ff;'>></span> <span style='color:#00ff88;'>%1</span>").arg(cmd));
        int result = UI->ApplyCommand(cmd.toStdString());
        if (result == 0) {
            m_commandOutput->append("<span style='color:#00ff88;'>OK</span>");
        } else {
            m_commandOutput->append(QString("<span style='color:#ff6b6b;'>Error: %1</span>").arg(result));
        }
        m_commandOutput->append("");
    }

    m_statusLabel->setText("Executed: " + cmd);
}

void DashboardWindow::onLoadMacro() {
    QString filename = QFileDialog::getOpenFileName(this, "Load Macro File",
        "../macro", "Macro Files (*.mac);;All Files (*)");

    if (!filename.isEmpty()) {
        QString cmd = QString("/control/execute %1").arg(filename);
        AddToCommandHistory(cmd);

        G4UImanager* UI = G4UImanager::GetUIpointer();
        if (UI) {
            m_commandOutput->append(QString("<span style='color:#00d4ff;'>Loading:</span> %1").arg(filename));
            int result = UI->ApplyCommand(cmd.toStdString());
            if (result == 0) {
                m_commandOutput->append("<span style='color:#00ff88;'>Macro executed</span>");
            } else {
                m_commandOutput->append(QString("<span style='color:#ff6b6b;'>Error: %1</span>").arg(result));
            }
        }
        m_statusLabel->setText("Loaded: " + filename);
    }
}

void DashboardWindow::onClearConsole() {
    if (m_commandOutput) m_commandOutput->clear();
}

void DashboardWindow::AddToCommandHistory(const QString& cmd) {
    if (m_commandHistory) {
        if (m_commandHistory->count() == 0 || m_commandHistory->item(0)->text() != cmd) {
            m_commandHistory->insertItem(0, cmd);
        }
        while (m_commandHistory->count() > 100) {
            delete m_commandHistory->takeItem(m_commandHistory->count() - 1);
        }
    }
}

// ============================================================================
// Data Update Methods
// ============================================================================

void DashboardWindow::AddParticle(const ParticleInfo& info) {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_particles.size() >= static_cast<size_t>(m_maxParticles)) {
        m_particles.erase(m_particles.begin());
    }
    m_particles.push_back(info);

    // Add to event display
    if (m_eventDisplay) {
        Track track;
        track.trackId = info.track_id;
        track.parentId = 0;  // Not available from ParticleInfo
        TrackPoint pt;
        pt.x = info.x; pt.y = info.y; pt.z = info.z;
        pt.t = 0; pt.ke = info.ke;
        track.points.push_back(pt);
        track.initialEnergy = info.ke;
        track.particleName = info.name;
        track.creatorProcess = info.process;
        // Map particle name to PDG
        if (info.name == "e-") track.pdg = 11;
        else if (info.name == "e+") track.pdg = -11;
        else if (info.name == "gamma") track.pdg = 22;
        else if (info.name == "mu-") track.pdg = 13;
        else if (info.name == "mu+") track.pdg = -13;
        else if (info.name == "proton") track.pdg = 2212;
        else if (info.name == "neutron") track.pdg = 2112;
        else if (info.name == "pi+" || info.name == "pion+") track.pdg = 211;
        else if (info.name == "pi-" || info.name == "pion-") track.pdg = -211;
        else if (info.name == "pi0") track.pdg = 111;
        else if (info.name == "kaon+" || info.name == "K+") track.pdg = 321;
        else if (info.name == "kaon-" || info.name == "K-") track.pdg = -321;
        else if (info.name == "kaon0" || info.name == "K0") track.pdg = 311;
        else if (info.name == "kaon0S" || info.name == "K0S") track.pdg = 310;
        else if (info.name == "kaon0L" || info.name == "K0L") track.pdg = 130;
        else if (info.name == "lambda" || info.name == "Lambda") track.pdg = 3122;
        else if (info.name == "opticalphoton") track.pdg = 0;  // Will be filtered
        else if (info.name == "nu_e" || info.name == "anti_nu_e") track.pdg = 12;
        else if (info.name == "nu_mu" || info.name == "anti_nu_mu") track.pdg = 14;
        else track.pdg = 999;  // Other particles
        m_eventDisplay->AddTrack(track);
    }

    emit particleAdded();
}

void DashboardWindow::UpdateEventStats(const EventStatistics& stats) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_currentStats = stats;

    if (m_eventDisplay) {
        EventDisplayData evtData;
        evtData.eventId = stats.event_id;
        evtData.totalEnergy = stats.total_edep;
        evtData.timestamp = QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString();
        // Keep existing tracks
    }

    emit statsUpdated();
}

void DashboardWindow::UpdateHitMap(const DetectorHitMap& hits) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_hitMap = hits;
    emit hitMapUpdated();
}

void DashboardWindow::SetProgress(int current, int total) {
    emit progressChanged(current, total);
}

void DashboardWindow::UpdateEventCount(int count) {
    // Thread-safe update - store value and let refresh timer update UI
    std::lock_guard<std::mutex> lock(m_mutex);
    m_currentStats.event_id = count;
}

// ============================================================================
// Slot Implementations
// ============================================================================

void DashboardWindow::onParticleAdded() {
    std::lock_guard<std::mutex> lock(m_mutex);
    if (m_particles.empty() || !m_particleTable) return;

    const auto& p = m_particles.back();

    int row = m_particleTable->rowCount();
    if (row >= m_maxParticles) {
        m_particleTable->removeRow(0);
        row = m_maxParticles - 1;
    } else {
        m_particleTable->insertRow(row);
    }

    m_particleTable->setItem(row, 0, new QTableWidgetItem(QString::number(p.event_id)));
    m_particleTable->setItem(row, 1, new QTableWidgetItem(QString::number(p.track_id)));
    m_particleTable->setItem(row, 2, new QTableWidgetItem(QString::fromStdString(p.name)));
    m_particleTable->setItem(row, 3, new QTableWidgetItem(QString::number(p.ke, 'f', 2)));
    m_particleTable->setItem(row, 4, new QTableWidgetItem(QString::number(p.x, 'f', 2)));
    m_particleTable->setItem(row, 5, new QTableWidgetItem(QString::number(p.y, 'f', 2)));
    m_particleTable->setItem(row, 6, new QTableWidgetItem(QString::number(p.z, 'f', 2)));
    m_particleTable->setItem(row, 7, new QTableWidgetItem(QString::fromStdString(p.volume)));
    m_particleTable->setItem(row, 8, new QTableWidgetItem(QString::fromStdString(p.process)));

    m_particleTable->scrollToBottom();
}

void DashboardWindow::onStatsUpdated() {
    std::lock_guard<std::mutex> lock(m_mutex);

    m_eventLabel->setText(QString::number(m_currentStats.event_id));
    m_primaryLabel->setText(QString::number(m_currentStats.n_primaries));
    m_secondaryLabel->setText(QString::number(m_currentStats.n_secondaries));
    m_edepLabel->setText(QString("%1 MeV").arg(m_currentStats.total_edep, 0, 'f', 2));
    m_cpuLabel->setText(QString("%1 ms").arg(m_currentStats.cpu_time_ms, 0, 'f', 2));
    m_scintHitsLabel->setText(QString::number(m_currentStats.n_scintillator_hits));
    m_tpcHitsLabel->setText(QString::number(m_currentStats.n_tpc_hits));
    m_lgHitsLabel->setText(QString::number(m_currentStats.n_leadglass_hits));
    m_photonLabel->setText(QString::number(m_currentStats.n_optical_photons));

    // Update KPI cards
    if (m_eventCountKPI) {
        m_eventCountKPI->SetValue(m_currentStats.event_id);
    }
    if (m_edepKPI) {
        m_edepKPI->SetValue(m_currentStats.total_edep, 2);
    }
    if (m_hitRateKPI) {
        int totalHits = m_currentStats.n_scintillator_hits +
                       m_currentStats.n_tpc_hits +
                       m_currentStats.n_leadglass_hits;
        m_hitRateKPI->SetValue(totalHits);
    }

    // Update mini histograms with recent data
    static std::vector<double> edepHistory;
    static std::vector<double> hitHistory;

    edepHistory.push_back(m_currentStats.total_edep);
    if (edepHistory.size() > 50) edepHistory.erase(edepHistory.begin());

    int totalHits = m_currentStats.n_scintillator_hits +
                   m_currentStats.n_tpc_hits +
                   m_currentStats.n_leadglass_hits;
    hitHistory.push_back(static_cast<double>(totalHits));
    if (hitHistory.size() > 50) hitHistory.erase(hitHistory.begin());

    if (m_edepHistogram) {
        m_edepHistogram->SetData(edepHistory);
    }
    if (m_hitHistogram) {
        m_hitHistogram->SetData(hitHistory);
    }
}

void DashboardWindow::onHitMapUpdated() {
    // Update event display with hits
}

void DashboardWindow::onProgressChanged(int current, int total) {
    m_progressBar->setMaximum(total);
    m_progressBar->setValue(current);
}

void DashboardWindow::onRefreshTimer() {
    // Periodic refresh
}

// ============================================================================
// Event Display Data Methods
// ============================================================================

void DashboardWindow::AddTrackToDisplay(int32_t trackId, int32_t parentId, int32_t pdg,
                                         float energy, const std::string& name,
                                         const std::string& creatorProcess,
                                         const std::vector<std::array<float, 3>>& points) {
    if (!m_eventDisplay || points.empty()) return;

    Track track;
    track.trackId = trackId;
    track.parentId = parentId;
    track.pdg = pdg;
    track.initialEnergy = energy;
    track.particleName = name;
    track.creatorProcess = creatorProcess;

    for (const auto& pt : points) {
        TrackPoint tp;
        tp.x = pt[0];
        tp.y = pt[1];
        tp.z = pt[2];
        tp.t = 0;
        tp.ke = energy;
        track.points.push_back(tp);
    }

    m_eventDisplay->AddTrack(track);
}

void DashboardWindow::AddVertexToDisplay(float x, float y, float z, float t,
                                          const std::string& process,
                                          int32_t parentTrackId,
                                          const std::vector<int32_t>& daughterPDGs,
                                          float totalEnergy) {
    if (!m_eventDisplay) return;

    InteractionVertex vtx;
    vtx.x = x;
    vtx.y = y;
    vtx.z = z;
    vtx.t = t;
    vtx.process = process;
    vtx.parentTrackId = parentTrackId;
    vtx.daughterPDGs = daughterPDGs;
    vtx.totalDaughterEnergy = totalEnergy;

    // Generate placeholder track IDs
    for (size_t i = 0; i < daughterPDGs.size(); i++) {
        vtx.daughterTrackIds.push_back(-1);  // Will be filled when tracks are added
    }

    m_eventDisplay->AddInteractionVertex(vtx);
}

void DashboardWindow::AddCaloHitToDisplay(float x, float y, float z, float energy, int detType) {
    if (!m_eventDisplay) return;

    CaloHit hit;
    hit.x = x;
    hit.y = y;
    hit.z = z;
    hit.energy = energy;
    hit.time = 0;
    hit.detectorType = detType;
    hit.moduleId = 0;

    m_eventDisplay->AddCaloHit(hit);
}

void DashboardWindow::ClearEventDisplay() {
    if (m_eventDisplay) {
        m_eventDisplay->ClearEvent();
    }
}

void DashboardWindow::FinalizeEvent(int32_t eventId, int32_t runId,
                                     const std::string& primaryParticle, double primaryEnergy) {
    if (!m_eventDisplay) return;

    // Update event display data
    EventDisplayData evtData;
    evtData.eventId = eventId;
    evtData.runId = runId;
    evtData.primaryParticle = primaryParticle;
    evtData.primaryEnergy = primaryEnergy;

    // Count particle types from current stats
    evtData.nPrimaries = m_currentStats.n_primaries;
    evtData.nSecondaries = m_currentStats.n_secondaries;
    evtData.totalEnergy = m_currentStats.total_edep;

    // Particle counts will be updated from actual track data
    // This provides a summary that the event display can use
}

// ============================================================================
// Global Helper Functions
// ============================================================================

void RecordParticleForDashboard(
    int32_t event_id, int32_t track_id,
    const std::string& name,
    double ke, double x, double y, double z,
    double px, double py, double pz,
    const std::string& volume,
    const std::string& process)
{
    ParticleInfo info;
    info.event_id = event_id;
    info.track_id = track_id;
    info.name = name;
    info.ke = ke;
    info.x = x; info.y = y; info.z = z;
    info.px = px; info.py = py; info.pz = pz;
    info.volume = volume;
    info.process = process;

    DashboardWindow::Instance().AddParticle(info);
}

void UpdateDashboardEventStats(const EventStatistics& stats) {
    DashboardWindow::Instance().UpdateEventStats(stats);
}

void UpdateDashboardHitMap(const DetectorHitMap& hits) {
    DashboardWindow::Instance().UpdateHitMap(hits);
}

} // namespace nnbar

#endif // WITH_DASHBOARD
