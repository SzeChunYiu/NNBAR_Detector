// ============================================================================
// DashboardPanel.cc
// PowerBI-style reusable panel widget implementation
// ============================================================================

#include "gui/DashboardPanel.hh"

#ifdef WITH_DASHBOARD

#include <QPainter>
#include <QPainterPath>
#include <QGraphicsDropShadowEffect>
#include <QPropertyAnimation>

namespace nnbar {

// ============================================================================
// DashboardPanel Implementation
// ============================================================================

DashboardPanel::DashboardPanel(const QString& title, QWidget* parent)
    : QFrame(parent), m_content(nullptr)
{
    SetupUI();
    SetTitle(title);
    ApplyStyle();
}

void DashboardPanel::SetupUI() {
    m_mainLayout = new QVBoxLayout(this);
    m_mainLayout->setContentsMargins(0, 0, 0, 0);
    m_mainLayout->setSpacing(0);

    // Header
    m_headerWidget = new QWidget();
    m_headerWidget->setFixedHeight(36);
    m_headerLayout = new QHBoxLayout(m_headerWidget);
    m_headerLayout->setContentsMargins(12, 6, 8, 6);
    m_headerLayout->setSpacing(8);

    m_titleLabel = new QLabel();
    m_titleLabel->setStyleSheet("font-size: 13px; font-weight: 600; color: #ffffff;");

    m_subtitleLabel = new QLabel();
    m_subtitleLabel->setStyleSheet("font-size: 10px; color: #b0b0b0;");  // Improved contrast
    m_subtitleLabel->hide();

    m_kpiLabel = new QLabel();
    m_kpiLabel->setStyleSheet("font-size: 18px; font-weight: bold; color: #00d4ff;");
    m_kpiLabel->hide();

    m_refreshBtn = new QPushButton("↻");
    m_refreshBtn->setFixedSize(24, 24);
    m_refreshBtn->setStyleSheet(
        "QPushButton { background: transparent; border: none; color: #999; font-size: 14px; }"
        "QPushButton:hover { color: #00d4ff; }"
    );
    m_refreshBtn->setToolTip("Refresh");
    connect(m_refreshBtn, &QPushButton::clicked, this, &DashboardPanel::onRefreshClicked);

    m_maximizeBtn = new QPushButton("⬜");
    m_maximizeBtn->setFixedSize(24, 24);
    m_maximizeBtn->setStyleSheet(
        "QPushButton { background: transparent; border: none; color: #999; font-size: 12px; }"
        "QPushButton:hover { color: #00d4ff; }"
    );
    m_maximizeBtn->setToolTip("Maximize");
    connect(m_maximizeBtn, &QPushButton::clicked, this, &DashboardPanel::onMaximizeClicked);

    m_headerLayout->addWidget(m_titleLabel);
    m_headerLayout->addWidget(m_subtitleLabel);
    m_headerLayout->addStretch();
    m_headerLayout->addWidget(m_kpiLabel);
    m_headerLayout->addWidget(m_refreshBtn);
    m_headerLayout->addWidget(m_maximizeBtn);

    m_mainLayout->addWidget(m_headerWidget);

    // Content frame
    m_contentFrame = new QFrame();
    m_contentFrame->setStyleSheet("background: transparent;");
    QVBoxLayout* contentLayout = new QVBoxLayout(m_contentFrame);
    contentLayout->setContentsMargins(8, 4, 8, 8);
    m_mainLayout->addWidget(m_contentFrame, 1);

    // Add subtle shadow effect
    QGraphicsDropShadowEffect* shadow = new QGraphicsDropShadowEffect(this);
    shadow->setBlurRadius(15);
    shadow->setColor(QColor(0, 0, 0, 80));
    shadow->setOffset(0, 2);
    setGraphicsEffect(shadow);
}

void DashboardPanel::ApplyStyle() {
    setStyleSheet(QString(
        "DashboardPanel {"
        "   background: qlineargradient(x1:0, y1:0, x2:0, y2:1,"
        "       stop:0 #2a2a2f, stop:1 #1e1e22);"
        "   border: 1px solid #3a3a3f;"
        "   border-radius: 8px;"
        "}"
    ));

    m_headerWidget->setStyleSheet(QString(
        "background: transparent;"
        "border-bottom: 1px solid %1;"
    ).arg(m_accentColor.darker(200).name()));
}

void DashboardPanel::paintEvent(QPaintEvent* event) {
    QFrame::paintEvent(event);

    // Draw accent line at top
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    QPainterPath path;
    path.moveTo(8, 0);
    path.lineTo(width() - 8, 0);

    QPen pen(m_accentColor);
    pen.setWidth(2);
    p.setPen(pen);
    p.drawPath(path);
}

void DashboardPanel::SetContent(QWidget* content) {
    if (m_content) {
        m_contentFrame->layout()->removeWidget(m_content);
    }
    m_content = content;
    if (m_content) {
        m_contentFrame->layout()->addWidget(m_content);
    }
}

void DashboardPanel::SetTitle(const QString& title) {
    m_titleLabel->setText(title);
}

void DashboardPanel::SetSubtitle(const QString& subtitle) {
    m_subtitleLabel->setText(subtitle);
    m_subtitleLabel->setVisible(!subtitle.isEmpty());
}

void DashboardPanel::SetAccentColor(const QColor& color) {
    m_accentColor = color;
    ApplyStyle();
    update();
}

void DashboardPanel::SetKPIValue(const QString& value, const QString& unit) {
    QString text = value;
    if (!unit.isEmpty()) {
        text += " <span style='font-size: 10px; color: #a0a0a0;'>" + unit + "</span>";
    }
    m_kpiLabel->setText(text);
    m_kpiLabel->show();
}

void DashboardPanel::SetMaximizable(bool enabled) {
    m_maximizable = enabled;
    m_maximizeBtn->setVisible(enabled);
}

void DashboardPanel::onMaximizeClicked() {
    if (m_isMaximized) {
        Restore();
    } else {
        Maximize();
    }
}

void DashboardPanel::Maximize() {
    if (!m_maximizable || m_isMaximized) return;

    m_isMaximized = true;
    m_maximizeBtn->setText("❐");
    m_maximizeBtn->setToolTip("Restore");

    emit maximizeRequested();
}

void DashboardPanel::Restore() {
    if (!m_isMaximized) return;

    m_isMaximized = false;
    m_maximizeBtn->setText("⬜");
    m_maximizeBtn->setToolTip("Maximize");

    emit restoreRequested();
}

void DashboardPanel::onRefreshClicked() {
    emit refreshRequested();
}

// ============================================================================
// KPICard Implementation
// ============================================================================

KPICard::KPICard(const QString& label, QWidget* parent)
    : QFrame(parent)
{
    SetupUI();
    m_labelText->setText(label);
}

void KPICard::SetupUI() {
    setFixedHeight(80);
    setMinimumWidth(120);

    setStyleSheet(
        "KPICard {"
        "   background: qlineargradient(x1:0, y1:0, x2:0, y2:1,"
        "       stop:0 #2a2a2f, stop:1 #1e1e22);"
        "   border: 1px solid #3a3a3f;"
        "   border-radius: 6px;"
        "}"
    );

    QVBoxLayout* layout = new QVBoxLayout(this);
    layout->setContentsMargins(12, 8, 12, 8);
    layout->setSpacing(2);

    m_labelText = new QLabel();
    m_labelText->setStyleSheet("font-size: 10px; color: #b0b0b0; font-weight: 500;");  // Improved contrast

    QHBoxLayout* valueLayout = new QHBoxLayout();
    valueLayout->setSpacing(4);

    m_valueText = new QLabel("--");
    m_valueText->setStyleSheet("font-size: 24px; font-weight: bold; color: #00d4ff;");

    m_unitText = new QLabel();
    m_unitText->setStyleSheet("font-size: 11px; color: #a0a0a0;");  // Improved contrast

    valueLayout->addWidget(m_valueText);
    valueLayout->addWidget(m_unitText);
    valueLayout->addStretch();

    m_trendText = new QLabel();
    m_trendText->setStyleSheet("font-size: 10px;");
    m_trendText->hide();

    layout->addWidget(m_labelText);
    layout->addLayout(valueLayout);
    layout->addWidget(m_trendText);
    layout->addStretch();
}

void KPICard::SetValue(const QString& value) {
    m_valueText->setText(value);
}

void KPICard::SetValue(double value, int precision) {
    m_valueText->setText(QString::number(value, 'f', precision));
}

void KPICard::SetValue(int value) {
    m_valueText->setText(QString::number(value));
}

void KPICard::SetUnit(const QString& unit) {
    m_unitText->setText(unit);
}

void KPICard::SetColor(const QColor& color) {
    m_color = color;
    m_valueText->setStyleSheet(QString("font-size: 24px; font-weight: bold; color: %1;").arg(color.name()));
}

void KPICard::SetTrend(double percentChange) {
    QString arrow = percentChange >= 0 ? "▲" : "▼";
    QString color = percentChange >= 0 ? "#00cc66" : "#ff4444";
    m_trendText->setText(QString("<span style='color: %1;'>%2 %3%</span>")
                         .arg(color, arrow, QString::number(std::abs(percentChange), 'f', 1)));
    m_trendText->show();
}

// ============================================================================
// MiniHistogram Implementation
// ============================================================================

MiniHistogram::MiniHistogram(QWidget* parent)
    : QWidget(parent)
{
    setMinimumHeight(40);
    setMaximumHeight(60);
}

void MiniHistogram::SetData(const std::vector<double>& data) {
    m_data = data;

    if (m_autoRange && !m_data.empty()) {
        m_min = *std::min_element(m_data.begin(), m_data.end());
        m_max = *std::max_element(m_data.begin(), m_data.end());
        if (m_max == m_min) m_max = m_min + 1.0;
    }

    update();
}

void MiniHistogram::SetColor(const QColor& color) {
    m_color = color;
    update();
}

void MiniHistogram::SetRange(double min, double max) {
    m_min = min;
    m_max = max;
    m_autoRange = false;
    update();
}

void MiniHistogram::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    // Background
    p.fillRect(rect(), QColor(20, 20, 25));

    if (m_data.empty()) return;

    int n = m_data.size();
    double barWidth = (double)width() / n;
    double range = m_max - m_min;

    // Draw bars
    for (int i = 0; i < n; i++) {
        double normalized = (m_data[i] - m_min) / range;
        int barHeight = normalized * (height() - 4);

        QRectF bar(i * barWidth + 1, height() - barHeight - 2, barWidth - 2, barHeight);

        QLinearGradient grad(bar.topLeft(), bar.bottomLeft());
        grad.setColorAt(0, m_color);
        grad.setColorAt(1, m_color.darker(150));

        p.fillRect(bar, grad);
    }

    // Border
    p.setPen(QPen(QColor(60, 60, 65), 1));
    p.drawRect(rect().adjusted(0, 0, -1, -1));
}

} // namespace nnbar

#endif // WITH_DASHBOARD
