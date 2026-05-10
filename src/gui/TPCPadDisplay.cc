// ============================================================================
// TPCPadDisplay.cc
// Qt widget for TPC pad plane readout visualization
// ============================================================================

#include "config.h"

#ifdef WITH_DASHBOARD

#include "gui/TPCPadDisplay.hh"
#include <QPainter>
#include <QPaintEvent>
#include <QMouseEvent>
#include <QWheelEvent>
#include <cmath>
#include <algorithm>

namespace nnbar {

// ============================================================================
// Constructor
// ============================================================================

TPCPadDisplay::TPCPadDisplay(QWidget* parent) : QWidget(parent) {
    setMinimumSize(400, 300);
    setMouseTracking(true);
    setAutoFillBackground(true);

    QPalette pal = palette();
    pal.setColor(QPalette::Window, QColor(20, 20, 25));
    setPalette(pal);

    // Initialize pad data
    SetGeometry(m_nRows, m_nPads);
}

// ============================================================================
// Configuration
// ============================================================================

void TPCPadDisplay::SetGeometry(int nRows, int nPads) {
    m_nRows = nRows;
    m_nPads = nPads;

    // Initialize data arrays for both sides
    m_padData.resize(2);
    for (int side = 0; side < 2; side++) {
        m_padData[side].resize(nRows);
        for (int row = 0; row < nRows; row++) {
            m_padData[side][row].resize(nPads, 0.0f);
        }
    }

    m_cacheValid = false;
    update();
}

void TPCPadDisplay::SetPadData(const std::vector<std::tuple<int, int, float>>& data, int side) {
    if (side < 0 || side > 1) return;

    // Clear existing data for this side
    for (auto& row : m_padData[side]) {
        std::fill(row.begin(), row.end(), 0.0f);
    }

    // Set new data
    for (const auto& [row, pad, charge] : data) {
        if (row >= 0 && row < m_nRows && pad >= 0 && pad < m_nPads) {
            m_padData[side][row][pad] = charge;
        }
    }

    m_cacheValid = false;
    update();
}

void TPCPadDisplay::AddCluster(int row, int pad, float charge) {
    m_clusters.emplace_back(row, pad, charge);
    m_cacheValid = false;
    update();
}

void TPCPadDisplay::Clear() {
    for (int side = 0; side < 2; side++) {
        for (auto& row : m_padData[side]) {
            std::fill(row.begin(), row.end(), 0.0f);
        }
    }
    m_clusters.clear();
    m_cacheValid = false;
    update();
}

void TPCPadDisplay::SetColorRange(float minCharge, float maxCharge) {
    m_minCharge = minCharge;
    m_maxCharge = maxCharge;
    m_cacheValid = false;
    update();
}

void TPCPadDisplay::SetSide(int side) {
    if (side != m_currentSide && side >= 0 && side < 2) {
        m_currentSide = side;
        m_cacheValid = false;
        update();
    }
}

void TPCPadDisplay::ToggleSide() {
    SetSide(1 - m_currentSide);
}

// ============================================================================
// Color Mapping
// ============================================================================

QColor TPCPadDisplay::ChargeToColor(float charge) const {
    if (charge < m_minCharge) {
        return QColor(20, 20, 40);  // Dark background for empty pads
    }

    float value;
    if (m_useLogScale) {
        float logMin = std::log10(m_minCharge);
        float logMax = std::log10(m_maxCharge);
        float logVal = std::log10(std::max(charge, m_minCharge));
        value = (logVal - logMin) / (logMax - logMin);
    } else {
        value = (charge - m_minCharge) / (m_maxCharge - m_minCharge);
    }

    value = std::max(0.0f, std::min(1.0f, value));

    // Blue -> Cyan -> Green -> Yellow -> Red colormap
    int r, g, b;
    if (value < 0.25f) {
        float t = value / 0.25f;
        r = 0;
        g = static_cast<int>(255 * t);
        b = 255;
    } else if (value < 0.5f) {
        float t = (value - 0.25f) / 0.25f;
        r = 0;
        g = 255;
        b = static_cast<int>(255 * (1 - t));
    } else if (value < 0.75f) {
        float t = (value - 0.5f) / 0.25f;
        r = static_cast<int>(255 * t);
        g = 255;
        b = 0;
    } else {
        float t = (value - 0.75f) / 0.25f;
        r = 255;
        g = static_cast<int>(255 * (1 - t));
        b = 0;
    }

    return QColor(r, g, b);
}

// ============================================================================
// Coordinate Conversion
// ============================================================================

QPointF TPCPadDisplay::PadToScreen(int row, int pad) const {
    // Map pad coordinates to screen coordinates
    float margin = 60.0f;
    float plotW = width() - 2 * margin - 60;  // Leave room for colorbar
    float plotH = height() - 2 * margin;

    float x = margin + (pad + 0.5f) / m_nPads * plotW;
    float y = margin + (m_nRows - row - 0.5f) / m_nRows * plotH;

    // Apply zoom and pan
    x = width() / 2.0f + (x - width() / 2.0f) * m_zoom + m_panX;
    y = height() / 2.0f + (y - height() / 2.0f) * m_zoom + m_panY;

    return QPointF(x, y);
}

void TPCPadDisplay::ScreenToPad(const QPoint& screen, int& row, int& pad) const {
    // Inverse of PadToScreen
    float margin = 60.0f;
    float plotW = width() - 2 * margin - 60;
    float plotH = height() - 2 * margin;

    float x = (screen.x() - m_panX - width() / 2.0f) / m_zoom + width() / 2.0f;
    float y = (screen.y() - m_panY - height() / 2.0f) / m_zoom + height() / 2.0f;

    pad = static_cast<int>((x - margin) / plotW * m_nPads);
    row = static_cast<int>(m_nRows - (y - margin) / plotH * m_nRows);
}

// ============================================================================
// Rendering
// ============================================================================

void TPCPadDisplay::RenderToCache() {
    m_renderCache = QPixmap(size());
    m_renderCache.fill(QColor(20, 20, 25));

    QPainter p(&m_renderCache);
    p.setRenderHint(QPainter::Antialiasing);

    // Draw title
    DrawTitle(p);

    // Calculate pad dimensions on screen
    float margin = 60.0f;
    float plotW = width() - 2 * margin - 60;
    float plotH = height() - 2 * margin;
    float padW = plotW / m_nPads * m_zoom;
    float padH = plotH / m_nRows * m_zoom;

    // Draw pads as colored rectangles
    const auto& sideData = m_padData[m_currentSide];
    for (int row = 0; row < m_nRows; row++) {
        for (int pad = 0; pad < m_nPads; pad++) {
            float charge = sideData[row][pad];
            if (charge < m_minCharge * 0.5f) continue;  // Skip empty pads

            QPointF center = PadToScreen(row, pad);
            QRectF rect(center.x() - padW / 2, center.y() - padH / 2, padW, padH);

            if (rect.right() < 0 || rect.left() > width() ||
                rect.bottom() < 0 || rect.top() > height()) {
                continue;  // Skip off-screen pads
            }

            p.fillRect(rect, ChargeToColor(charge));
        }
    }

    // Draw cluster markers
    p.setPen(QPen(Qt::white, 2));
    for (const auto& [row, pad, charge] : m_clusters) {
        QPointF center = PadToScreen(row, pad);
        float radius = std::min(padW, padH) * 2;
        p.drawEllipse(center, radius, radius);
    }

    // Draw grid (if zoomed in enough)
    if (m_zoom > 2.0f) {
        DrawGrid(p);
    }

    // Draw axes labels
    p.setPen(Qt::lightGray);
    p.setFont(QFont("Arial", 10));

    // X-axis label
    p.drawText(QRect(0, height() - 25, width() - 60, 20),
               Qt::AlignCenter, "Pad Number");

    // Y-axis label (rotated)
    p.save();
    p.translate(15, height() / 2);
    p.rotate(-90);
    p.drawText(QRect(-100, 0, 200, 20), Qt::AlignCenter, "Pad Row");
    p.restore();

    // Draw colorbar
    DrawColorBar(p);

    m_cacheValid = true;
}

void TPCPadDisplay::DrawTitle(QPainter& p) const {
    QString sideStr = (m_currentSide == 0) ? "+z Side" : "-z Side";
    QString title = QString("TPC Pad Readout (%1)").arg(sideStr);

    p.setPen(Qt::white);
    p.setFont(QFont("Arial", 12, QFont::Bold));
    p.drawText(QRect(0, 5, width(), 25), Qt::AlignCenter, title);
}

void TPCPadDisplay::DrawGrid(QPainter& p) const {
    p.setPen(QPen(QColor(60, 60, 80), 1, Qt::DotLine));

    // Draw row lines
    for (int row = 0; row <= m_nRows; row += 10) {
        QPointF start = PadToScreen(row, 0);
        QPointF end = PadToScreen(row, m_nPads - 1);
        p.drawLine(start, end);
    }

    // Draw pad lines
    for (int pad = 0; pad <= m_nPads; pad += 25) {
        QPointF start = PadToScreen(0, pad);
        QPointF end = PadToScreen(m_nRows - 1, pad);
        p.drawLine(start, end);
    }
}

void TPCPadDisplay::DrawColorBar(QPainter& p) const {
    float barX = width() - 50;
    float barY = 60;
    float barW = 20;
    float barH = height() - 120;

    // Draw gradient
    for (int i = 0; i < static_cast<int>(barH); i++) {
        float value = 1.0f - static_cast<float>(i) / barH;
        float charge = m_useLogScale ?
            std::pow(10.0f, std::log10(m_minCharge) + value * (std::log10(m_maxCharge) - std::log10(m_minCharge))) :
            m_minCharge + value * (m_maxCharge - m_minCharge);
        p.setPen(ChargeToColor(charge));
        p.drawLine(static_cast<int>(barX), static_cast<int>(barY + i),
                   static_cast<int>(barX + barW), static_cast<int>(barY + i));
    }

    // Draw border
    p.setPen(Qt::lightGray);
    p.drawRect(static_cast<int>(barX), static_cast<int>(barY),
               static_cast<int>(barW), static_cast<int>(barH));

    // Draw labels
    p.setFont(QFont("Arial", 8));
    QString maxLabel = QString::number(m_maxCharge / 1000, 'f', 0) + "k";
    QString minLabel = QString::number(m_minCharge, 'f', 0);
    p.drawText(static_cast<int>(barX - 5), static_cast<int>(barY - 5), maxLabel);
    p.drawText(static_cast<int>(barX - 5), static_cast<int>(barY + barH + 12), minLabel);

    // Label
    p.save();
    p.translate(static_cast<int>(barX + barW + 15), static_cast<int>(barY + barH / 2));
    p.rotate(-90);
    p.drawText(-40, 0, "Charge (e)");
    p.restore();
}

// ============================================================================
// Event Handlers
// ============================================================================

void TPCPadDisplay::paintEvent(QPaintEvent*) {
    if (!m_cacheValid) {
        RenderToCache();
    }

    QPainter p(this);
    p.drawPixmap(0, 0, m_renderCache);
}

void TPCPadDisplay::mousePressEvent(QMouseEvent* event) {
    m_lastMousePos = event->pos();

    if (event->button() == Qt::LeftButton) {
        int row, pad;
        ScreenToPad(event->pos(), row, pad);
        if (row >= 0 && row < m_nRows && pad >= 0 && pad < m_nPads) {
            emit PadClicked(row, pad, m_currentSide);
        }
    }
}

void TPCPadDisplay::mouseMoveEvent(QMouseEvent* event) {
    if (event->buttons() & Qt::RightButton) {
        // Pan
        QPoint delta = event->pos() - m_lastMousePos;
        m_panX += delta.x();
        m_panY += delta.y();
        m_lastMousePos = event->pos();
        m_cacheValid = false;
        update();
    }

    // Update tooltip with pad info
    int row, pad;
    ScreenToPad(event->pos(), row, pad);
    if (row >= 0 && row < m_nRows && pad >= 0 && pad < m_nPads) {
        float charge = m_padData[m_currentSide][row][pad];
        setToolTip(QString("Row %1, Pad %2\nCharge: %3")
                   .arg(row).arg(pad).arg(charge, 0, 'f', 0));
    }
}

void TPCPadDisplay::wheelEvent(QWheelEvent* event) {
    float zoomFactor = 1.1f;
    if (event->angleDelta().y() > 0) {
        m_zoom *= zoomFactor;
    } else {
        m_zoom /= zoomFactor;
    }
    m_zoom = std::max(0.5f, std::min(10.0f, m_zoom));

    m_cacheValid = false;
    update();
}

void TPCPadDisplay::resizeEvent(QResizeEvent*) {
    m_cacheValid = false;
}

} // namespace nnbar

#endif // WITH_DASHBOARD
