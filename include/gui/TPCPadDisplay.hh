// ============================================================================
// TPCPadDisplay.hh
// Qt widget for TPC pad plane readout visualization
// Displays 2D pad occupancy map with color-coded charge
// ============================================================================

#ifndef TPC_PAD_DISPLAY_HH
#define TPC_PAD_DISPLAY_HH

#include "config.h"

#ifdef WITH_DASHBOARD

#include <QWidget>
#include <QPixmap>
#include <QTimer>
#include <vector>
#include <tuple>

namespace nnbar {

/**
 * @class TPCPadDisplay
 * @brief Qt widget showing TPC pad plane readout map
 *
 * Displays the TPC pad occupancy as a 2D heatmap, where:
 * - X-axis: pad number (azimuthal direction)
 * - Y-axis: pad row (radial direction)
 * - Color: charge collected (log scale, blue to red)
 *
 * Features:
 * - Real-time update during simulation
 * - Side selection (+z / -z anode)
 * - Zoom and pan
 * - Cluster markers
 */
class TPCPadDisplay : public QWidget {
    Q_OBJECT

public:
    explicit TPCPadDisplay(QWidget* parent = nullptr);
    ~TPCPadDisplay() override = default;

    /**
     * Set pad plane geometry
     * @param nRows Number of pad rows
     * @param nPads Number of pads per row
     */
    void SetGeometry(int nRows, int nPads);

    /**
     * Set pad occupancy data
     * @param data Vector of (row, pad, charge) tuples
     * @param side 0 = +z, 1 = -z
     */
    void SetPadData(const std::vector<std::tuple<int, int, float>>& data, int side);

    /**
     * Add cluster marker
     * @param row Center row
     * @param pad Center pad
     * @param charge Total charge
     */
    void AddCluster(int row, int pad, float charge);

    /**
     * Clear all data
     */
    void Clear();

    /**
     * Get current side being displayed
     */
    int GetCurrentSide() const { return m_currentSide; }

    /**
     * Set color scale range
     * @param minCharge Minimum charge (maps to blue)
     * @param maxCharge Maximum charge (maps to red)
     */
    void SetColorRange(float minCharge, float maxCharge);

    /**
     * Enable/disable log scale for colors
     */
    void SetLogScale(bool enable) { m_useLogScale = enable; update(); }

public slots:
    /**
     * Switch displayed side (0 = +z, 1 = -z)
     */
    void SetSide(int side);

    /**
     * Toggle between +z and -z
     */
    void ToggleSide();

signals:
    /**
     * Emitted when user clicks on a pad
     */
    void PadClicked(int row, int pad, int side);

protected:
    void paintEvent(QPaintEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;

private:
    // Geometry
    int m_nRows = 85;
    int m_nPads = 250;

    // Pad data: [side][row][pad] -> charge
    std::vector<std::vector<std::vector<float>>> m_padData;

    // Cluster markers: (row, pad, charge)
    std::vector<std::tuple<int, int, float>> m_clusters;

    // Display state
    int m_currentSide = 0;
    float m_minCharge = 100.0f;
    float m_maxCharge = 100000.0f;
    bool m_useLogScale = true;

    // View transform
    float m_zoom = 1.0f;
    float m_panX = 0.0f;
    float m_panY = 0.0f;
    QPoint m_lastMousePos;

    // Rendering cache
    QPixmap m_renderCache;
    bool m_cacheValid = false;

    // Helper methods
    QColor ChargeToColor(float charge) const;
    QPointF PadToScreen(int row, int pad) const;
    void ScreenToPad(const QPoint& screen, int& row, int& pad) const;
    void RenderToCache();
    void DrawColorBar(QPainter& p) const;
    void DrawTitle(QPainter& p) const;
    void DrawGrid(QPainter& p) const;
};

} // namespace nnbar

#endif // WITH_DASHBOARD
#endif // TPC_PAD_DISPLAY_HH
