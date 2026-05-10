// ============================================================================
// G4CommandPanel.cc
// Geant4 command tree and session interface panel
// ============================================================================

#include "gui/G4CommandPanel.hh"

#ifdef WITH_DASHBOARD

#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QGroupBox>
#include <QLabel>
#include <QHeaderView>
#include <QScrollBar>
#include <QDateTime>
#include <QKeyEvent>

#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIcommand.hh"

namespace nnbar {

G4CommandPanel::G4CommandPanel(QWidget* parent)
    : QWidget(parent)
{
    SetupUI();
    BuildCommandTree();
    UpdateCompleter();
}

void G4CommandPanel::SetupUI() {
    QVBoxLayout* mainLayout = new QVBoxLayout(this);
    mainLayout->setContentsMargins(5, 5, 5, 5);
    mainLayout->setSpacing(5);

    m_mainSplitter = new QSplitter(Qt::Horizontal);

    // ========== LEFT PANEL: Command Tree ==========
    QWidget* treePanel = new QWidget();
    QVBoxLayout* treeLayout = new QVBoxLayout(treePanel);
    treeLayout->setContentsMargins(0, 0, 0, 0);

    QLabel* treeLabel = new QLabel("Command Tree");
    treeLabel->setStyleSheet("font-weight: bold; color: #ffd700;");
    treeLayout->addWidget(treeLabel);

    m_treeSearch = new QLineEdit();
    m_treeSearch->setPlaceholderText("Search commands...");
    m_treeSearch->setStyleSheet("background: #2d2d2d; border: 1px solid #444; padding: 3px;");
    connect(m_treeSearch, &QLineEdit::textChanged, this, &G4CommandPanel::onFilterChanged);
    treeLayout->addWidget(m_treeSearch);

    m_commandTree = new QTreeWidget();
    m_commandTree->setHeaderLabel("G4 Commands");
    m_commandTree->setStyleSheet(
        "QTreeWidget { background: #1e1e1e; border: 1px solid #333; }"
        "QTreeWidget::item { padding: 2px; }"
        "QTreeWidget::item:selected { background: #0066aa; }"
    );
    m_commandTree->setAlternatingRowColors(true);
    connect(m_commandTree, &QTreeWidget::itemClicked, this, &G4CommandPanel::onTreeItemClicked);
    connect(m_commandTree, &QTreeWidget::itemDoubleClicked, this, &G4CommandPanel::onTreeItemDoubleClicked);
    treeLayout->addWidget(m_commandTree);

    m_mainSplitter->addWidget(treePanel);

    // ========== RIGHT PANEL: Session ==========
    QWidget* sessionPanel = new QWidget();
    QVBoxLayout* sessionLayout = new QVBoxLayout(sessionPanel);
    sessionLayout->setContentsMargins(0, 0, 0, 0);

    // Quick command buttons
    QHBoxLayout* quickLayout = new QHBoxLayout();
    m_runBtn = new QPushButton("/run/beamOn 1");
    m_runBtn->setStyleSheet("background: #2d5016; padding: 5px;");
    connect(m_runBtn, &QPushButton::clicked, this, &G4CommandPanel::onQuickCommand);

    m_visBtn = new QPushButton("/vis/drawVolume");
    m_visBtn->setStyleSheet("background: #164050; padding: 5px;");
    connect(m_visBtn, &QPushButton::clicked, this, &G4CommandPanel::onQuickCommand);

    m_gunBtn = new QPushButton("/gun/particle");
    m_gunBtn->setStyleSheet("background: #503016; padding: 5px;");
    connect(m_gunBtn, &QPushButton::clicked, this, &G4CommandPanel::onQuickCommand);

    quickLayout->addWidget(m_runBtn);
    quickLayout->addWidget(m_visBtn);
    quickLayout->addWidget(m_gunBtn);
    sessionLayout->addLayout(quickLayout);

    // Output area with filter
    QHBoxLayout* filterLayout = new QHBoxLayout();
    QLabel* outputLabel = new QLabel("Output:");
    outputLabel->setStyleSheet("font-weight: bold; color: #00d4ff;");
    filterLayout->addWidget(outputLabel);

    m_filterCombo = new QComboBox();
    m_filterCombo->addItems({"All", "Info", "Warning", "Error"});
    m_filterCombo->setStyleSheet("background: #2d2d2d;");
    filterLayout->addWidget(m_filterCombo);

    QPushButton* clearBtn = new QPushButton("Clear");
    clearBtn->setStyleSheet("background: #501616; padding: 3px 10px;");
    connect(clearBtn, &QPushButton::clicked, this, &G4CommandPanel::onClearOutput);
    filterLayout->addWidget(clearBtn);
    filterLayout->addStretch();
    sessionLayout->addLayout(filterLayout);

    m_outputArea = new QTextEdit();
    m_outputArea->setReadOnly(true);
    m_outputArea->setStyleSheet(
        "QTextEdit { background: #0a0a0a; color: #00ff00; "
        "font-family: 'Courier New', monospace; font-size: 10px; "
        "border: 1px solid #333; }"
    );
    sessionLayout->addWidget(m_outputArea, 3);

    // Command history
    QLabel* histLabel = new QLabel("History:");
    histLabel->setStyleSheet("color: #888;");
    sessionLayout->addWidget(histLabel);

    m_historyList = new QListWidget();
    m_historyList->setMaximumHeight(80);
    m_historyList->setStyleSheet("background: #1a1a1a; border: 1px solid #333;");
    connect(m_historyList, &QListWidget::itemClicked, this, &G4CommandPanel::onHistoryItemClicked);
    sessionLayout->addWidget(m_historyList);

    // Command input
    QHBoxLayout* inputLayout = new QHBoxLayout();
    QLabel* promptLabel = new QLabel("Session:");
    promptLabel->setStyleSheet("color: #ffd700; font-weight: bold;");
    inputLayout->addWidget(promptLabel);

    m_commandInput = new QLineEdit();
    m_commandInput->setPlaceholderText("Enter G4 command (e.g., /run/beamOn 10)");
    m_commandInput->setStyleSheet(
        "background: #1a1a1a; color: #00ff00; "
        "font-family: 'Courier New', monospace; "
        "border: 1px solid #444; padding: 5px;"
    );

    m_completerModel = new QStringListModel(this);
    m_completer = new QCompleter(m_completerModel, this);
    m_completer->setCaseSensitivity(Qt::CaseInsensitive);
    m_commandInput->setCompleter(m_completer);

    connect(m_commandInput, &QLineEdit::returnPressed, this, &G4CommandPanel::onCommandEntered);
    inputLayout->addWidget(m_commandInput);

    QPushButton* execBtn = new QPushButton("Execute");
    execBtn->setStyleSheet("background: #0066aa; padding: 5px 15px;");
    connect(execBtn, &QPushButton::clicked, this, &G4CommandPanel::onCommandEntered);
    inputLayout->addWidget(execBtn);

    sessionLayout->addLayout(inputLayout);

    m_mainSplitter->addWidget(sessionPanel);
    m_mainSplitter->setSizes({250, 450});

    mainLayout->addWidget(m_mainSplitter);
}

void G4CommandPanel::BuildCommandTree() {
    m_commandTree->clear();
    BuildTreeFromG4();
}

void G4CommandPanel::BuildTreeFromG4() {
    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (!UI) return;

    G4UIcommandTree* tree = UI->GetTree();
    if (!tree) return;

    // Build tree recursively
    std::function<void(G4UIcommandTree*, QTreeWidgetItem*)> buildBranch;
    buildBranch = [&](G4UIcommandTree* cmdTree, QTreeWidgetItem* parent) {
        // Add subdirectories
        for (G4int i = 1; i <= cmdTree->GetTreeEntry(); i++) {
            G4UIcommandTree* subTree = cmdTree->GetTree(i);
            if (subTree) {
                QString path = QString::fromStdString(subTree->GetPathName());
                QString name = path.split('/').filter(QString()).last();

                QTreeWidgetItem* item;
                if (parent) {
                    item = new QTreeWidgetItem(parent);
                } else {
                    item = new QTreeWidgetItem(m_commandTree);
                }
                item->setText(0, name);
                item->setData(0, Qt::UserRole, path);
                item->setIcon(0, style()->standardIcon(QStyle::SP_DirIcon));

                buildBranch(subTree, item);
            }
        }

        // Add commands
        for (G4int i = 1; i <= cmdTree->GetCommandEntry(); i++) {
            G4UIcommand* cmd = cmdTree->GetCommand(i);
            if (cmd) {
                QString path = QString::fromStdString(cmd->GetCommandPath());
                QString name = path.split('/').last();

                QTreeWidgetItem* item;
                if (parent) {
                    item = new QTreeWidgetItem(parent);
                } else {
                    item = new QTreeWidgetItem(m_commandTree);
                }
                item->setText(0, name);
                item->setData(0, Qt::UserRole, path);
                item->setToolTip(0, QString::fromStdString(cmd->GetTitle()));
            }
        }
    };

    buildBranch(tree, nullptr);
}

void G4CommandPanel::onCommandEntered() {
    QString cmd = m_commandInput->text().trimmed();
    if (cmd.isEmpty()) return;

    ExecuteCommand(cmd);
    m_commandInput->clear();
}

void G4CommandPanel::ExecuteCommand(const QString& cmd) {
    AddToHistory(cmd);

    AppendOutput("> " + cmd, "command");

    G4UImanager* UI = G4UImanager::GetUIpointer();
    if (UI) {
        G4int result = UI->ApplyCommand(cmd.toStdString());
        if (result != 0) {
            AppendOutput("Command failed with code: " + QString::number(result), "error");
        }
    }

    emit commandExecuted(cmd);
}

void G4CommandPanel::AddToHistory(const QString& cmd) {
    // Remove if already exists
    auto it = std::find(m_commandHistory.begin(), m_commandHistory.end(), cmd);
    if (it != m_commandHistory.end()) {
        m_commandHistory.erase(it);
    }

    m_commandHistory.push_front(cmd);
    if (m_commandHistory.size() > MAX_HISTORY) {
        m_commandHistory.pop_back();
    }

    // Update history list widget
    m_historyList->clear();
    for (const auto& c : m_commandHistory) {
        m_historyList->addItem(c);
    }

    m_historyIndex = -1;
}

void G4CommandPanel::onTreeItemClicked(QTreeWidgetItem* item, int) {
    QString path = item->data(0, Qt::UserRole).toString();
    emit commandTreeItemSelected(path);
}

void G4CommandPanel::onTreeItemDoubleClicked(QTreeWidgetItem* item, int) {
    QString path = item->data(0, Qt::UserRole).toString();
    if (!path.endsWith('/')) {
        // It's a command, put it in the input
        m_commandInput->setText(path);
        m_commandInput->setFocus();
    }
}

void G4CommandPanel::onHistoryItemClicked(QListWidgetItem* item) {
    m_commandInput->setText(item->text());
}

void G4CommandPanel::onQuickCommand() {
    QPushButton* btn = qobject_cast<QPushButton*>(sender());
    if (btn) {
        QString cmd = btn->text();
        if (cmd.startsWith("/gun/particle")) {
            m_commandInput->setText("/gun/particle ");
            m_commandInput->setFocus();
        } else {
            ExecuteCommand(cmd);
        }
    }
}

void G4CommandPanel::onFilterChanged(const QString& text) {
    // Filter tree items
    std::function<bool(QTreeWidgetItem*)> filterItem;
    filterItem = [&](QTreeWidgetItem* item) -> bool {
        bool match = item->text(0).contains(text, Qt::CaseInsensitive);
        bool childMatch = false;

        for (int i = 0; i < item->childCount(); i++) {
            if (filterItem(item->child(i))) {
                childMatch = true;
            }
        }

        item->setHidden(!match && !childMatch && !text.isEmpty());
        if (childMatch && !text.isEmpty()) {
            item->setExpanded(true);
        }

        return match || childMatch;
    };

    for (int i = 0; i < m_commandTree->topLevelItemCount(); i++) {
        filterItem(m_commandTree->topLevelItem(i));
    }
}

void G4CommandPanel::onClearOutput() {
    m_outputArea->clear();
}

void G4CommandPanel::AppendOutput(const QString& text, const QString& type) {
    std::lock_guard<std::mutex> lock(m_outputMutex);

    QString timestamp = QDateTime::currentDateTime().toString("hh:mm:ss");
    QString color = "#00ff00";  // Green default

    if (type == "error") color = "#ff4444";
    else if (type == "warning") color = "#ffaa00";
    else if (type == "command") color = "#00d4ff";

    QString formatted = QString("<span style='color: #666;'>[%1]</span> "
                                "<span style='color: %2;'>%3</span>")
                        .arg(timestamp, color, text.toHtmlEscaped());

    m_outputArea->append(formatted);

    // Auto-scroll to bottom
    QScrollBar* sb = m_outputArea->verticalScrollBar();
    sb->setValue(sb->maximum());
}

void G4CommandPanel::AppendG4cout(const QString& text) {
    AppendOutput(text, "info");
}

void G4CommandPanel::AppendG4cerr(const QString& text) {
    AppendOutput(text, "error");
}

void G4CommandPanel::ClearOutput() {
    m_outputArea->clear();
}

void G4CommandPanel::UpdateCompleter() {
    QStringList commands;

    // Common G4 commands
    commands << "/run/beamOn" << "/run/initialize" << "/run/verbose"
             << "/gun/particle" << "/gun/energy" << "/gun/position" << "/gun/direction"
             << "/vis/drawVolume" << "/vis/viewer/flush" << "/vis/viewer/refresh"
             << "/tracking/verbose" << "/process/list"
             << "/control/execute" << "/control/verbose";

    m_completerModel->setStringList(commands);
}

QString G4CommandPanel::GetCommandPath(QTreeWidgetItem* item) {
    return item->data(0, Qt::UserRole).toString();
}

} // namespace nnbar

#endif // WITH_DASHBOARD
