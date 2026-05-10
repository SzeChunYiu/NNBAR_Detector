// ============================================================================
// G4CommandPanel.hh
// Geant4 command tree and session interface panel
// ============================================================================

#ifndef G4_COMMAND_PANEL_HH
#define G4_COMMAND_PANEL_HH

#include "config.h"

#ifdef WITH_DASHBOARD

#include <QWidget>
#include <QTreeWidget>
#include <QLineEdit>
#include <QTextEdit>
#include <QListWidget>
#include <QComboBox>
#include <QSplitter>
#include <QPushButton>
#include <QCompleter>
#include <QStringListModel>
#include <mutex>
#include <deque>

namespace nnbar {

class G4CommandPanel : public QWidget {
    Q_OBJECT

public:
    explicit G4CommandPanel(QWidget* parent = nullptr);
    ~G4CommandPanel() = default;

    // Add output to console
    void AppendOutput(const QString& text, const QString& type = "info");
    void AppendG4cout(const QString& text);
    void AppendG4cerr(const QString& text);

    // Execute a G4 command
    void ExecuteCommand(const QString& cmd);

    // Build command tree from G4UIcommandTree
    void BuildCommandTree();

    // Clear output
    void ClearOutput();

signals:
    void commandExecuted(const QString& cmd);
    void commandTreeItemSelected(const QString& path);

private slots:
    void onCommandEntered();
    void onTreeItemClicked(QTreeWidgetItem* item, int column);
    void onTreeItemDoubleClicked(QTreeWidgetItem* item, int column);
    void onHistoryItemClicked(QListWidgetItem* item);
    void onQuickCommand();
    void onFilterChanged(const QString& text);
    void onClearOutput();

private:
    void SetupUI();
    void BuildTreeFromG4();
    void AddToHistory(const QString& cmd);
    void UpdateCompleter();
    QString GetCommandPath(QTreeWidgetItem* item);

    // UI Components
    QSplitter* m_mainSplitter;

    // Left side: Command tree
    QTreeWidget* m_commandTree;
    QLineEdit* m_treeSearch;

    // Right side: Session
    QTextEdit* m_outputArea;
    QLineEdit* m_commandInput;
    QListWidget* m_historyList;
    QComboBox* m_filterCombo;
    QCompleter* m_completer;
    QStringListModel* m_completerModel;

    // Quick command buttons
    QPushButton* m_runBtn;
    QPushButton* m_visBtn;
    QPushButton* m_gunBtn;

    // Command history
    std::deque<QString> m_commandHistory;
    int m_historyIndex = -1;
    static const int MAX_HISTORY = 100;

    // Output buffer
    std::mutex m_outputMutex;
};

} // namespace nnbar

#endif // WITH_DASHBOARD
#endif // G4_COMMAND_PANEL_HH
