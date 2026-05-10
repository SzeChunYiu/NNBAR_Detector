// ============================================================================
// MaterialBudgetPlot.cc
// Custom Qt widget for plotting material budget (X0, λ) vs pseudorapidity (η)
// ============================================================================

#include "gui/MaterialBudgetPlot.hh"

#ifdef WITH_DASHBOARD

#include <QPainter>
#include <QPainterPath>
#include <QToolTip>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iomanip>

namespace nnbar {

// Region colors (ATLAS/CMS style palette)
const std::vector<std::pair<std::string, QColor>> MaterialBudgetPlot::s_regionColors = {
    {"Beampipe",        QColor(100, 100, 100)},    // Gray
    {"Silicon",         QColor(255, 200, 100)},    // Orange
    {"TPC Gas",         QColor(150, 220, 255)},    // Light blue
    {"TPC Walls",       QColor(100, 180, 220)},    // Blue
    {"Scintillator",    QColor(120, 200, 120)},    // Green
    {"Lead Glass Side", QColor(220, 80, 80)},      // Red
    {"Lead Glass F/B",  QColor(180, 60, 60)},      // Dark red
    {"Cosmic Shield",   QColor(80, 80, 120)},      // Dark blue
    {"Support",         QColor(180, 180, 180)},    // Light gray
    {"Other",           QColor(200, 200, 100)}     // Yellow
};

MaterialBudgetPlot::MaterialBudgetPlot(QWidget* parent)
    : QWidget(parent)
{
    setMouseTracking(true);
    setMinimumSize(400, 300);

    // Default title
    m_title = "Radiation Length vs η";

    // Initialize with NNBAR detector data
    InitializeNNBARData();
}

void MaterialBudgetPlot::InitializeNNBARData() {
    // NNBAR detector material budget data
    // Based on actual detector geometry
    // SORTED BY GEOMETRIC RADIUS: innermost (bottom) to outermost (top)
    m_data.clear();

    // Create eta bins from -3 to +3
    const int nBins = 60;
    const double etaStep = 6.0 / nBins;

    for (int i = 0; i < nBins; i++) {
        EtaBin bin;
        bin.etaMin = -3.0 + i * etaStep;
        bin.etaMax = bin.etaMin + etaStep;
        double etaCenter = (bin.etaMin + bin.etaMax) / 2.0;
        double absEta = std::abs(etaCenter);

        // ========== SORTED BY RADIUS (innermost first, for bottom of stack) ==========

        // 1. Beampipe - innermost, r ~ 112-114 cm
        bin.contributions.push_back({"Beampipe", 0.02});

        // 2. TPC Gas - inside TPC modules, r ~ 114-200 cm
        if (absEta < 1.5) {
            bin.contributions.push_back({"TPC Gas", 0.03});
        }

        // 3. TPC Walls - TPC structure
        if (absEta < 1.5) {
            bin.contributions.push_back({"TPC Walls", 0.15 * (1.0 - absEta/2.0)});
        }

        // 4. Scintillator - r ~ 197 cm
        if (absEta < 2.0) {
            double scintX0 = 0.45 * std::exp(-absEta * 0.5);
            bin.contributions.push_back({"Scintillator", scintX0});
        }

        // 5. Lead Glass Side - r ~ 275 cm (central region)
        if (absEta < 1.2) {
            double lgSideX0 = 12.5 * (1.0 - absEta/1.5);
            bin.contributions.push_back({"Lead Glass Side", lgSideX0});
        }

        // 6. Lead Glass Front/Back - z ~ 335 cm (forward regions)
        if (absEta > 1.0 && absEta < 2.5) {
            double lgFBX0 = 15.0 * std::min(1.0, (absEta - 0.8) / 0.5);
            if (absEta > 2.0) lgFBX0 *= (2.5 - absEta) / 0.5;
            bin.contributions.push_back({"Lead Glass F/B", lgFBX0});
        }

        // 7. Support structures - distributed
        bin.contributions.push_back({"Support", 0.1});

        // 8. Cosmic Shield - outermost, r ~ 375 cm
        if (absEta < 2.5) {
            bin.contributions.push_back({"Cosmic Shield", 0.8 * (1.0 - absEta/3.0)});
        }

        m_data.push_back(bin);
    }

    // Auto-scale Y axis
    UpdateYScale();
}

void MaterialBudgetPlot::UpdateYScale() {
    m_valueMax = 1.0;  // Minimum

    for (const auto& bin : m_data) {
        double total = 0.0;
        for (const auto& contrib : bin.contributions) {
            if (m_plotType == PlotType::RadiationLength) {
                total += contrib.second;
            } else {
                // Convert X0 to lambda (rough approximation: lambda ≈ X0 / 10 for lead glass)
                total += contrib.second * GetLambdaFactor(contrib.first);
            }
        }
        m_valueMax = std::max(m_valueMax, total * 1.2);
    }
}

double MaterialBudgetPlot::GetLambdaFactor(const std::string& region) const {
    // Ratio of interaction length to radiation length for each material
    if (region.find("Lead Glass") != std::string::npos) return 0.075;  // λ/X0 for PbO-SiO2
    if (region == "Scintillator") return 0.55;
    if (region == "TPC Gas") return 0.67;
    if (region == "TPC Walls") return 0.53;
    if (region == "Beampipe") return 0.5;
    if (region == "Silicon") return 0.6;
    if (region == "Cosmic Shield") return 0.56;
    return 0.5;  // Default
}

void MaterialBudgetPlot::SetPlotType(PlotType type) {
    m_plotType = type;
    m_title = (type == PlotType::RadiationLength)
        ? "Radiation Length (X₀) vs η"
        : "Interaction Length (λ) vs η";
    UpdateYScale();
    update();
}

void MaterialBudgetPlot::SetData(const std::vector<EtaBin>& bins) {
    m_data = bins;
    UpdateYScale();
    update();
}

void MaterialBudgetPlot::ClearData() {
    m_data.clear();
    update();
}

void MaterialBudgetPlot::paintEvent(QPaintEvent*) {
    QPainter p(this);
    p.setRenderHint(QPainter::Antialiasing);

    // Background
    p.fillRect(rect(), QColor(30, 30, 35));

    // Plot area background
    QRect plotArea(m_marginLeft, m_marginTop,
                   width() - m_marginLeft - m_marginRight,
                   height() - m_marginTop - m_marginBottom);
    p.fillRect(plotArea, QColor(20, 20, 25));

    if (m_showGrid) DrawGrid(p);
    DrawAxes(p);
    DrawData(p);
    DrawTitle(p);
    if (m_showLegend) DrawLegend(p);
    if (m_hoverX >= 0) DrawHoverInfo(p);
}

void MaterialBudgetPlot::DrawGrid(QPainter& p) {
    p.setPen(QPen(QColor(50, 50, 55), 1, Qt::DotLine));

    int plotWidth = width() - m_marginLeft - m_marginRight;
    int plotHeight = height() - m_marginTop - m_marginBottom;

    // Vertical lines (eta)
    for (double eta = -3.0; eta <= 3.0; eta += 0.5) {
        int x = m_marginLeft + (eta - m_etaMin) / (m_etaMax - m_etaMin) * plotWidth;
        p.drawLine(x, m_marginTop, x, height() - m_marginBottom);
    }

    // Horizontal lines (value)
    int nLines = 5;
    for (int i = 0; i <= nLines; i++) {
        double value = m_valueMin + (m_valueMax - m_valueMin) * i / nLines;
        int y = height() - m_marginBottom - (value - m_valueMin) / (m_valueMax - m_valueMin) * plotHeight;
        p.drawLine(m_marginLeft, y, width() - m_marginRight, y);
    }
}

void MaterialBudgetPlot::DrawAxes(QPainter& p) {
    p.setPen(QPen(QColor(200, 200, 200), 2));

    int plotWidth = width() - m_marginLeft - m_marginRight;
    int plotHeight = height() - m_marginTop - m_marginBottom;

    // X axis
    p.drawLine(m_marginLeft, height() - m_marginBottom,
               width() - m_marginRight, height() - m_marginBottom);

    // Y axis
    p.drawLine(m_marginLeft, m_marginTop,
               m_marginLeft, height() - m_marginBottom);

    // X axis labels
    p.setFont(QFont("Sans", 9));
    for (double eta = -3.0; eta <= 3.0; eta += 1.0) {
        int x = m_marginLeft + (eta - m_etaMin) / (m_etaMax - m_etaMin) * plotWidth;
        QString label = QString::number(eta, 'f', 0);
        if (eta > 0) label = "+" + label;
        p.drawText(x - 15, height() - m_marginBottom + 20, 30, 20,
                   Qt::AlignCenter, label);
    }

    // X axis title
    p.setFont(QFont("Sans", 10, QFont::Bold));
    p.drawText(m_marginLeft, height() - 15, plotWidth, 20,
               Qt::AlignCenter, "Pseudorapidity η");

    // Y axis labels
    p.setFont(QFont("Sans", 9));
    int nLabels = 5;
    for (int i = 0; i <= nLabels; i++) {
        double value = m_valueMin + (m_valueMax - m_valueMin) * i / nLabels;
        int y = height() - m_marginBottom - (value - m_valueMin) / (m_valueMax - m_valueMin) * plotHeight;
        QString label = QString::number(value, 'f', 1);
        p.drawText(5, y - 10, m_marginLeft - 10, 20,
                   Qt::AlignRight | Qt::AlignVCenter, label);
    }

    // Y axis title (rotated)
    p.save();
    p.translate(15, height() / 2);
    p.rotate(-90);
    p.setFont(QFont("Sans", 10, QFont::Bold));
    QString yTitle = (m_plotType == PlotType::RadiationLength) ? "X₀" : "λ";
    p.drawText(-plotHeight/2, 0, plotHeight, 20, Qt::AlignCenter, yTitle);
    p.restore();
}

void MaterialBudgetPlot::DrawData(QPainter& p) {
    if (m_data.empty()) return;

    int plotWidth = width() - m_marginLeft - m_marginRight;
    int plotHeight = height() - m_marginTop - m_marginBottom;

    // Collect all unique regions
    std::vector<std::string> regions;
    for (const auto& bin : m_data) {
        for (const auto& contrib : bin.contributions) {
            if (std::find(regions.begin(), regions.end(), contrib.first) == regions.end()) {
                regions.push_back(contrib.first);
            }
        }
    }

    // Draw stacked histogram
    for (size_t binIdx = 0; binIdx < m_data.size(); binIdx++) {
        const auto& bin = m_data[binIdx];

        int x1 = m_marginLeft + (bin.etaMin - m_etaMin) / (m_etaMax - m_etaMin) * plotWidth;
        int x2 = m_marginLeft + (bin.etaMax - m_etaMin) / (m_etaMax - m_etaMin) * plotWidth;
        int barWidth = x2 - x1;

        double cumulative = 0.0;

        for (const auto& region : regions) {
            // Find contribution for this region in this bin
            double value = 0.0;
            for (const auto& contrib : bin.contributions) {
                if (contrib.first == region) {
                    value = contrib.second;
                    if (m_plotType == PlotType::InteractionLength) {
                        value *= GetLambdaFactor(region);
                    }
                    break;
                }
            }

            if (value <= 0) continue;

            double y1 = cumulative;
            double y2 = cumulative + value;
            cumulative = y2;

            int screenY1 = height() - m_marginBottom - y1 / m_valueMax * plotHeight;
            int screenY2 = height() - m_marginBottom - y2 / m_valueMax * plotHeight;

            QRect bar(x1, screenY2, barWidth, screenY1 - screenY2);

            QColor color = ColorForRegion(region);
            p.fillRect(bar, color);

            // Border
            p.setPen(QPen(color.darker(120), 1));
            p.drawRect(bar);
        }
    }
}

void MaterialBudgetPlot::DrawLegend(QPainter& p) {
    // Collect unique regions from data
    std::vector<std::string> regions;
    for (const auto& bin : m_data) {
        for (const auto& contrib : bin.contributions) {
            if (std::find(regions.begin(), regions.end(), contrib.first) == regions.end()) {
                regions.push_back(contrib.first);
            }
        }
    }

    int legendX = width() - m_marginRight - 150;
    int legendY = m_marginTop + 10;
    int boxSize = 12;
    int lineHeight = 16;
    int legendHeight = regions.size() * lineHeight + 10;

    // Legend background
    p.fillRect(legendX - 5, legendY - 5, 155, legendHeight, QColor(40, 40, 45, 200));
    p.setPen(QPen(QColor(80, 80, 80), 1));
    p.drawRect(legendX - 5, legendY - 5, 155, legendHeight);

    p.setFont(QFont("Sans", 8));

    int y = legendY;
    for (const auto& region : regions) {
        QColor color = ColorForRegion(region);
        p.fillRect(legendX, y, boxSize, boxSize, color);
        p.setPen(QPen(color.darker(120), 1));
        p.drawRect(legendX, y, boxSize, boxSize);

        p.setPen(QColor(200, 200, 200));
        p.drawText(legendX + boxSize + 5, y, 130, lineHeight,
                   Qt::AlignLeft | Qt::AlignVCenter, QString::fromStdString(region));
        y += lineHeight;
    }
}

void MaterialBudgetPlot::DrawTitle(QPainter& p) {
    p.setPen(QColor(255, 215, 0));  // Gold
    p.setFont(QFont("Sans", 12, QFont::Bold));
    p.drawText(m_marginLeft, 5, width() - m_marginLeft - m_marginRight, 30,
               Qt::AlignCenter, m_title);
}

void MaterialBudgetPlot::DrawHoverInfo(QPainter& p) {
    QString info = GetDataAtEta(m_hoverEta);
    if (info.isEmpty()) return;

    // Info box
    QFont font("Sans", 9);
    p.setFont(font);
    QFontMetrics fm(font);

    QStringList lines = info.split('\n');
    int boxWidth = 0;
    for (const auto& line : lines) {
#if QT_VERSION >= QT_VERSION_CHECK(5, 11, 0)
        boxWidth = std::max(boxWidth, fm.horizontalAdvance(line));
#else
        boxWidth = std::max(boxWidth, fm.width(line));
#endif
    }
    int boxHeight = lines.size() * fm.height() + 10;
    boxWidth += 20;

    int boxX = std::min(m_hoverX + 10, width() - boxWidth - 10);
    int boxY = std::max(m_hoverY - boxHeight - 10, 10);

    p.fillRect(boxX, boxY, boxWidth, boxHeight, QColor(50, 50, 55, 230));
    p.setPen(QPen(QColor(100, 100, 100), 1));
    p.drawRect(boxX, boxY, boxWidth, boxHeight);

    p.setPen(QColor(200, 200, 200));
    int y = boxY + 5;
    for (const auto& line : lines) {
        p.drawText(boxX + 5, y, boxWidth - 10, fm.height(),
                   Qt::AlignLeft | Qt::AlignVCenter, line);
        y += fm.height();
    }

    // Vertical line at hover position
    p.setPen(QPen(QColor(255, 255, 0, 150), 1, Qt::DashLine));
    p.drawLine(m_hoverX, m_marginTop, m_hoverX, height() - m_marginBottom);
}

void MaterialBudgetPlot::mouseMoveEvent(QMouseEvent* event) {
    int x = event->pos().x();
    int y = event->pos().y();

    // Check if within plot area
    if (x >= m_marginLeft && x <= width() - m_marginRight &&
        y >= m_marginTop && y <= height() - m_marginBottom) {
        m_hoverX = x;
        m_hoverY = y;
        m_hoverEta = ScreenToEta(x);

        QString info = GetDataAtEta(m_hoverEta);
        emit etaHovered(m_hoverEta, info);
    } else {
        m_hoverX = -1;
    }

    update();
}

void MaterialBudgetPlot::mousePressEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        double eta = ScreenToEta(event->pos().x());
        emit plotClicked(eta);
    }
}

void MaterialBudgetPlot::resizeEvent(QResizeEvent*) {
    update();
}

QPointF MaterialBudgetPlot::DataToScreen(double eta, double value) {
    int plotWidth = width() - m_marginLeft - m_marginRight;
    int plotHeight = height() - m_marginTop - m_marginBottom;

    double x = m_marginLeft + (eta - m_etaMin) / (m_etaMax - m_etaMin) * plotWidth;
    double y = height() - m_marginBottom - (value - m_valueMin) / (m_valueMax - m_valueMin) * plotHeight;

    return QPointF(x, y);
}

double MaterialBudgetPlot::ScreenToEta(int x) {
    int plotWidth = width() - m_marginLeft - m_marginRight;
    return m_etaMin + (x - m_marginLeft) / (double)plotWidth * (m_etaMax - m_etaMin);
}

double MaterialBudgetPlot::ScreenToValue(int y) {
    int plotHeight = height() - m_marginTop - m_marginBottom;
    return m_valueMin + (height() - m_marginBottom - y) / (double)plotHeight * (m_valueMax - m_valueMin);
}

QColor MaterialBudgetPlot::ColorForRegion(const std::string& region) {
    for (const auto& rc : s_regionColors) {
        if (region.find(rc.first) != std::string::npos || rc.first.find(region) != std::string::npos) {
            return rc.second;
        }
    }
    return QColor(150, 150, 150);  // Default gray
}

QString MaterialBudgetPlot::GetDataAtEta(double eta) const {
    // Find bin containing this eta
    for (const auto& bin : m_data) {
        if (eta >= bin.etaMin && eta < bin.etaMax) {
            std::ostringstream oss;
            oss << std::fixed << std::setprecision(2);
            oss << "η = " << eta << "\n";

            double total = 0.0;
            for (const auto& contrib : bin.contributions) {
                double value = contrib.second;
                if (m_plotType == PlotType::InteractionLength) {
                    value *= GetLambdaFactor(contrib.first);
                }
                if (value > 0.001) {
                    oss << contrib.first << ": " << value << "\n";
                    total += value;
                }
            }
            oss << "Total: " << total;

            return QString::fromStdString(oss.str());
        }
    }
    return QString();
}

QString MaterialBudgetPlot::ExportToCSV() const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4);

    // Collect all regions
    std::vector<std::string> regions;
    for (const auto& bin : m_data) {
        for (const auto& contrib : bin.contributions) {
            if (std::find(regions.begin(), regions.end(), contrib.first) == regions.end()) {
                regions.push_back(contrib.first);
            }
        }
    }

    // Header
    oss << "eta_min,eta_max,eta_center";
    for (const auto& region : regions) {
        oss << "," << region;
    }
    oss << ",Total\n";

    // Data
    for (const auto& bin : m_data) {
        double etaCenter = (bin.etaMin + bin.etaMax) / 2.0;
        oss << bin.etaMin << "," << bin.etaMax << "," << etaCenter;

        double total = 0.0;
        for (const auto& region : regions) {
            double value = 0.0;
            for (const auto& contrib : bin.contributions) {
                if (contrib.first == region) {
                    value = contrib.second;
                    if (m_plotType == PlotType::InteractionLength) {
                        value *= GetLambdaFactor(region);
                    }
                    break;
                }
            }
            oss << "," << value;
            total += value;
        }
        oss << "," << total << "\n";
    }

    return QString::fromStdString(oss.str());
}

QString MaterialBudgetPlot::ExportToJSON() const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4);

    oss << "{\n";
    oss << "  \"plot_type\": \"" << (m_plotType == PlotType::RadiationLength ? "radiation_length" : "interaction_length") << "\",\n";
    oss << "  \"unit\": \"" << (m_plotType == PlotType::RadiationLength ? "X0" : "lambda") << "\",\n";
    oss << "  \"eta_range\": [" << m_etaMin << ", " << m_etaMax << "],\n";
    oss << "  \"bins\": [\n";

    for (size_t i = 0; i < m_data.size(); i++) {
        const auto& bin = m_data[i];
        oss << "    {\n";
        oss << "      \"eta_min\": " << bin.etaMin << ",\n";
        oss << "      \"eta_max\": " << bin.etaMax << ",\n";
        oss << "      \"contributions\": {\n";

        double total = 0.0;
        for (size_t j = 0; j < bin.contributions.size(); j++) {
            const auto& contrib = bin.contributions[j];
            double value = contrib.second;
            if (m_plotType == PlotType::InteractionLength) {
                value *= GetLambdaFactor(contrib.first);
            }
            oss << "        \"" << contrib.first << "\": " << value;
            if (j < bin.contributions.size() - 1) oss << ",";
            oss << "\n";
            total += value;
        }

        oss << "      },\n";
        oss << "      \"total\": " << total << "\n";
        oss << "    }";
        if (i < m_data.size() - 1) oss << ",";
        oss << "\n";
    }

    oss << "  ]\n";
    oss << "}\n";

    return QString::fromStdString(oss.str());
}

} // namespace nnbar

#endif // WITH_DASHBOARD
