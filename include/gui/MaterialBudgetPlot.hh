// ============================================================================
// MaterialBudgetPlot.hh
// Custom Qt widget for plotting material budget (X0, λ) vs pseudorapidity (η)
// ============================================================================

#ifndef MATERIAL_BUDGET_PLOT_HH
#define MATERIAL_BUDGET_PLOT_HH

#include "config.h"

#ifdef WITH_DASHBOARD

#include <QWidget>
#include <QPainter>
#include <QMouseEvent>
#include <vector>
#include <string>
#include <functional>

namespace nnbar {

// Data point for material budget
struct MaterialBudgetPoint {
    double eta;           // Pseudorapidity
    double x0;            // Radiation lengths
    double lambda;        // Interaction lengths
    std::string region;   // Detector region name
};

// Eta bin for stacked plot
struct EtaBin {
    double etaMin, etaMax;
    std::vector<std::pair<std::string, double>> contributions;  // region -> value
};

class MaterialBudgetPlot : public QWidget {
    Q_OBJECT

public:
    enum class PlotType {
        RadiationLength,    // X0 vs η
        InteractionLength   // λ vs η
    };

    explicit MaterialBudgetPlot(QWidget* parent = nullptr);
    ~MaterialBudgetPlot() = default;

    // Set plot type
    void SetPlotType(PlotType type);
    PlotType GetPlotType() const { return m_plotType; }

    // Data management
    void SetData(const std::vector<EtaBin>& bins);
    void ClearData();

    // Display options
    void SetShowLegend(bool show) { m_showLegend = show; update(); }
    void SetShowGrid(bool show) { m_showGrid = show; update(); }
    void SetStacked(bool stacked) { m_stacked = stacked; update(); }
    void SetTitle(const QString& title) { m_title = title; update(); }

    // Export data
    QString ExportToCSV() const;
    QString ExportToJSON() const;

    // Get data at eta
    QString GetDataAtEta(double eta) const;

signals:
    void etaHovered(double eta, const QString& info);
    void plotClicked(double eta);

protected:
    void paintEvent(QPaintEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;

private:
    void DrawAxes(QPainter& p);
    void DrawGrid(QPainter& p);
    void DrawData(QPainter& p);
    void DrawLegend(QPainter& p);
    void DrawTitle(QPainter& p);
    void DrawHoverInfo(QPainter& p);

    // Coordinate transformations
    QPointF DataToScreen(double eta, double value);
    double ScreenToEta(int x);
    double ScreenToValue(int y);

    // Color for detector region
    QColor ColorForRegion(const std::string& region);

    // Initialize with NNBAR detector geometry data
    void InitializeNNBARData();
    void UpdateYScale();
    double GetLambdaFactor(const std::string& region) const;

    PlotType m_plotType = PlotType::RadiationLength;
    std::vector<EtaBin> m_data;

    // Display options
    bool m_showLegend = true;
    bool m_showGrid = true;
    bool m_stacked = true;
    QString m_title;

    // Plot margins
    int m_marginLeft = 70;
    int m_marginRight = 20;
    int m_marginTop = 40;
    int m_marginBottom = 50;

    // Axis ranges
    double m_etaMin = -3.0;
    double m_etaMax = 3.0;
    double m_valueMin = 0.0;
    double m_valueMax = 20.0;  // Will be auto-scaled

    // Hover info
    int m_hoverX = -1;
    int m_hoverY = -1;
    double m_hoverEta = 0.0;

    // Region colors (consistent with ATLAS style)
    static const std::vector<std::pair<std::string, QColor>> s_regionColors;
};

} // namespace nnbar

#endif // WITH_DASHBOARD
#endif // MATERIAL_BUDGET_PLOT_HH
