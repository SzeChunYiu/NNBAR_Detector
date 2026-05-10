// ============================================================================
// DashboardPanel.hh
// Reusable panel widget with maximize/restore functionality
// PowerBI-style design for physics monitoring
// ============================================================================

#ifndef DASHBOARD_PANEL_HH
#define DASHBOARD_PANEL_HH

#include "config.h"

#ifdef WITH_DASHBOARD

#include <QWidget>
#include <QFrame>
#include <QLabel>
#include <QPushButton>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <functional>

namespace nnbar {

class DashboardPanel : public QFrame {
    Q_OBJECT

public:
    explicit DashboardPanel(const QString& title, QWidget* parent = nullptr);
    ~DashboardPanel() = default;

    // Set the content widget
    void SetContent(QWidget* content);
    QWidget* GetContent() const { return m_content; }

    // Title
    void SetTitle(const QString& title);
    QString GetTitle() const { return m_titleLabel->text(); }

    // Subtitle/info
    void SetSubtitle(const QString& subtitle);

    // Icon color (accent color for the panel)
    void SetAccentColor(const QColor& color);

    // KPI value display (for summary panels)
    void SetKPIValue(const QString& value, const QString& unit = "");

    // Enable/disable maximize button
    void SetMaximizable(bool enabled);

    // Check if currently maximized
    bool IsMaximized() const { return m_isMaximized; }

signals:
    void maximizeRequested();
    void restoreRequested();
    void refreshRequested();

public slots:
    void Maximize();
    void Restore();

protected:
    void paintEvent(QPaintEvent* event) override;

private slots:
    void onMaximizeClicked();
    void onRefreshClicked();

private:
    void SetupUI();
    void ApplyStyle();

    // UI Components
    QVBoxLayout* m_mainLayout;
    QWidget* m_headerWidget;
    QHBoxLayout* m_headerLayout;
    QLabel* m_titleLabel;
    QLabel* m_subtitleLabel;
    QLabel* m_kpiLabel;
    QPushButton* m_maximizeBtn;
    QPushButton* m_refreshBtn;
    QWidget* m_content;
    QFrame* m_contentFrame;

    // State
    bool m_isMaximized = false;
    bool m_maximizable = true;
    QColor m_accentColor = QColor(0, 212, 255);  // Cyan default

    // For restoring position
    QWidget* m_originalParent = nullptr;
    int m_originalIndex = -1;
};

// ============================================================================
// KPI Card - Small summary card showing a single metric
// ============================================================================

class KPICard : public QFrame {
    Q_OBJECT

public:
    explicit KPICard(const QString& label, QWidget* parent = nullptr);

    void SetValue(const QString& value);
    void SetValue(double value, int precision = 2);
    void SetValue(int value);
    void SetUnit(const QString& unit);
    void SetColor(const QColor& color);
    void SetTrend(double percentChange);  // Shows up/down arrow with %

private:
    void SetupUI();

    QLabel* m_labelText;
    QLabel* m_valueText;
    QLabel* m_unitText;
    QLabel* m_trendText;
    QColor m_color = QColor(0, 212, 255);
};

// ============================================================================
// MiniHistogram - Small inline histogram for dashboard cards
// ============================================================================

class MiniHistogram : public QWidget {
    Q_OBJECT

public:
    explicit MiniHistogram(QWidget* parent = nullptr);

    void SetData(const std::vector<double>& data);
    void SetColor(const QColor& color);
    void SetRange(double min, double max);

protected:
    void paintEvent(QPaintEvent* event) override;

private:
    std::vector<double> m_data;
    QColor m_color = QColor(0, 212, 255);
    double m_min = 0.0;
    double m_max = 0.0;
    bool m_autoRange = true;
};

} // namespace nnbar

#endif // WITH_DASHBOARD
#endif // DASHBOARD_PANEL_HH
