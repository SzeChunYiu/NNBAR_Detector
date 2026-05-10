// ============================================================================
// EventDisplay.hh
// Premium-grade ATLAS/CMS-style event display for NNBAR detector simulation
// Accurate geometry rendering from actual detector construction data
// ============================================================================

#ifndef EVENT_DISPLAY_HH
#define EVENT_DISPLAY_HH

#include "config.h"

#ifdef WITH_DASHBOARD

#include <QWidget>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QTimer>
#include <QMatrix4x4>
#include <QVector3D>
#include <vector>
#include <mutex>
#include <cstdint>
#include <map>
#include "util/GeometryParameters.hh"

namespace nnbar {

// ============================================================================
// Enhanced Track Data Structures
// ============================================================================

struct TrackPoint {
    float x, y, z;          // Position in cm
    float t;                // Time in ns
    float ke;               // Kinetic energy at this point
    std::string process;    // Process that created this step
    std::string volume;     // Volume name
};

struct Track {
    int32_t trackId;
    int32_t parentId;
    int32_t pdg;
    float initialEnergy;
    std::string particleName;
    std::string creatorProcess;
    std::vector<TrackPoint> points;

    // Helper to get track direction at a point
    QVector3D DirectionAt(size_t index) const {
        if (index + 1 >= points.size()) return QVector3D(0, 0, 1);
        const auto& p1 = points[index];
        const auto& p2 = points[index + 1];
        QVector3D dir(p2.x - p1.x, p2.y - p1.y, p2.z - p1.z);
        return dir.normalized();
    }
};

// Interaction vertex where particles are created
struct InteractionVertex {
    float x, y, z;
    float t;
    std::string process;        // e.g., "nCapture", "pi-Inelastic", "compt"
    int32_t parentTrackId;
    std::vector<int32_t> daughterTrackIds;
    std::vector<int32_t> daughterPDGs;
    float totalDaughterEnergy;

    int Multiplicity() const { return daughterTrackIds.size(); }
};

// Calorimeter hit for visualization
struct CaloHit {
    float x, y, z;
    float energy;           // MeV
    float time;             // ns
    int detectorType;       // 0=Scintillator, 1=LeadGlass, 2=TPC
    int moduleId;
    std::string volumeName;
};

// Optical photon density in detector regions
struct OpticalPhotonRegion {
    int regionId;
    int photonCount;
    float x, y, z;
    std::string volumeName;  // Volume where photons are detected
};

// Energy deposition in a specific volume/region
struct EnergyDeposit {
    float x, y, z;          // Position in cm
    float energy;           // Energy deposited in MeV
    int detectorType;       // 0=Beampipe, 1=TPC, 2=Scintillator, 3=LeadGlass, 4=Shield
    int moduleId;           // Module/block index
    std::string volumeName;
};

// Detector type enumeration for energy tracking
enum class DetectorType {
    Beampipe = 0,
    TPC = 1,
    Scintillator = 2,
    LeadGlass = 3,
    Shield = 4,
    Other = 5
};

// ============================================================================
// Geometry Data Structures (loaded from CSV/construction)
// ============================================================================

struct LeadGlassBlock {
    double x, y, z;         // Position in cm
    double rotX, rotY, rotZ;// Rotation angles in degrees
    int surface;            // 0-3 for sides, 4-5 for front/back
    int index;
};

struct ScintillatorModule {
    double x, y, z;
    int surface;            // Which surface (top/bottom/left/right/front/back)
    int index;
};

struct TPCModule {
    double cx, cy, cz;      // Center position
    double hx, hy, hz;      // Half-dimensions
    bool isTypeII;          // Type II (horizontal) or Type I (vertical)
    int zSide;              // -1 for front, +1 for back
};

// ============================================================================
// Event Data Container
// ============================================================================

struct EventDisplayData {
    int32_t eventId = 0;
    int32_t runId = 0;
    std::vector<Track> tracks;
    std::vector<CaloHit> caloHits;
    std::vector<InteractionVertex> vertices;
    std::vector<OpticalPhotonRegion> opticalPhotonRegions;
    std::vector<EnergyDeposit> energyDeposits;
    int totalOpticalPhotons = 0;
    double totalEnergy = 0.0;
    double primaryEnergy = 0.0;
    std::string primaryParticle;
    std::string timestamp;

    // Particle Statistics
    int nPrimaries = 0;
    int nSecondaries = 0;
    int nGammas = 0;
    int nElectrons = 0;
    int nNeutrons = 0;
    int nProtons = 0;
    int nPions = 0;
    int nMuons = 0;
    int nOther = 0;

    // Energy deposited per detector (MeV)
    double edepBeampipe = 0.0;
    double edepTPC = 0.0;
    double edepScintillator = 0.0;
    double edepLeadGlass = 0.0;
    double edepShield = 0.0;
    double edepOther = 0.0;

    // Optical photon counts per detector
    int photonsScintillator = 0;
    int photonsLeadGlass = 0;
    int photonsTPC = 0;

    // Optical photon counts by process
    int scintPhotonCount = 0;
    int cerenkovPhotonCount = 0;
};

// View mode enumeration
enum class ViewMode {
    Transverse,     // XY view (beam axis into screen)
    Longitudinal,   // RZ view (side view)
    Perspective3D   // Full 3D rotatable view
};

// ============================================================================
// Main EventDisplay Class
// ============================================================================

class EventDisplay : public QWidget {
    Q_OBJECT

public:
    explicit EventDisplay(QWidget* parent = nullptr);
    ~EventDisplay();

    // Geometry loading
    void LoadGeometryData();
    bool IsGeometryLoaded() const { return m_geometryLoaded; }

    // View mode control
    void SetViewMode(ViewMode mode);
    ViewMode GetViewMode() const { return m_viewMode; }

    // Data update methods (thread-safe)
    void AddTrack(const Track& track);
    void AddCaloHit(const CaloHit& hit);
    void AddInteractionVertex(const InteractionVertex& vertex);
    void SetEventData(const EventDisplayData& data);
    void ClearEvent();

    // Display options
    void SetShowTracks(bool show) { m_showTracks = show; update(); }
    void SetShowCaloHits(bool show) { m_showCaloHits = show; update(); }
    void SetShowDetector(bool show) { m_showDetector = show; update(); }
    void SetShowGrid(bool show) { m_showGrid = show; update(); }
    void SetShowAxes(bool show) { m_showAxes = show; update(); }
    void SetShowOpticalPhotons(bool show) { m_showOpticalPhotons = show; update(); }
    void SetShowVertices(bool show) { m_showVertices = show; update(); }
    void SetShowTrackArrows(bool show) { m_showTrackArrows = show; update(); }
    void SetShowLabels(bool show) { m_showLabels = show; update(); }
    void SetShowLegend(bool show) { m_showLegend = show; update(); }
    void SetShowScaleBar(bool show) { m_showScaleBar = show; update(); }
    void SetShowRadiationLength(bool show) { m_showRadiationLength = show; update(); }
    void SetShowEnergyHeatmap(bool show) { m_showEnergyHeatmap = show; update(); }

    // Track filtering
    void SetMinTrackEnergy(float mev) { m_minTrackEnergy = mev; update(); }
    void SetShowElectrons(bool show) { m_showElectrons = show; update(); }
    void SetShowGammas(bool show) { m_showGammas = show; update(); }
    void SetShowNeutrons(bool show) { m_showNeutrons = show; update(); }
    void SetShowProtons(bool show) { m_showProtons = show; update(); }
    void SetShowPions(bool show) { m_showPions = show; update(); }
    void SetShowMuons(bool show) { m_showMuons = show; update(); }
    void SetShowOther(bool show) { m_showOther = show; update(); }

    // Optical photon data
    void AddOpticalPhotonRegion(const OpticalPhotonRegion& region);
    void SetTotalOpticalPhotons(int count);

    // Energy deposition data
    void AddEnergyDeposit(const EnergyDeposit& deposit);

    // Camera controls
    void ResetView();
    void ZoomIn();
    void ZoomOut();
    void SetZoom(float zoom) { m_zoom = zoom; update(); }
    float GetZoom() const { return m_zoom; }

    // Selection
    Track* GetSelectedTrack() { return m_selectedTrack; }
    InteractionVertex* GetSelectedVertex() { return m_selectedVertex; }

signals:
    void viewModeChanged(ViewMode mode);
    void eventDataUpdated();
    void trackSelected(const Track* track);
    void vertexSelected(const InteractionVertex* vertex);
    void geometryLoaded();

protected:
    void paintEvent(QPaintEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;
    void wheelEvent(QWheelEvent* event) override;
    void resizeEvent(QResizeEvent* event) override;

private:
    // ========== Geometry Loading ==========
    void LoadLeadGlassPositions();
    void LoadScintillatorPositions();
    void CalculateTPCModules();
    void CalculateShieldingGeometry();
    bool LoadCSV(const std::string& filename, std::vector<std::vector<double>>& data);

    // ========== Drawing Functions ==========
    void DrawBackground(QPainter& p);
    void DrawGrid(QPainter& p);
    void DrawDetectorGeometry(QPainter& p);
    void DrawDetectorXY(QPainter& p);
    void DrawDetectorRZ(QPainter& p);
    void DrawDetector3D(QPainter& p);

    // Individual detector components
    void DrawBeampipe(QPainter& p);
    void DrawTPCModules(QPainter& p);
    void DrawScintillatorModules(QPainter& p);
    void DrawLeadGlassBlocks(QPainter& p);
    void DrawCosmicShielding(QPainter& p);

    // Event data visualization
    void DrawTracks(QPainter& p);
    void DrawTrackWithArrows(QPainter& p, const Track& track);
    void DrawCaloHits(QPainter& p);
    void DrawInteractionVertices(QPainter& p);
    void DrawOpticalPhotonDensity(QPainter& p);
    void DrawEnergyHeatmap(QPainter& p);

    // Overlays and information
    void DrawAxes(QPainter& p);
    void DrawEventInfo(QPainter& p);
    void DrawStatisticsPanel(QPainter& p);
    void DrawLegend(QPainter& p);
    void DrawScaleBar(QPainter& p);
    void DrawLogo(QPainter& p);
    void DrawTooltip(QPainter& p);
    void DrawSelectedInfo(QPainter& p);

    // ========== Coordinate Transformation ==========
    QPointF ToScreen(float x, float y);
    QPointF ToScreen3D(float x, float y, float z);
    float ToScreenSize(float size);
    QPointF WorldToScreen(const QVector3D& world);

    // ========== Hit Testing ==========
    Track* HitTestTracks(const QPoint& pos);
    InteractionVertex* HitTestVertices(const QPoint& pos);
    CaloHit* HitTestHits(const QPoint& pos);

    // ========== Color Utilities ==========
    QColor ColorForPDG(int pdg, float alpha = 1.0f);
    QColor ColorForEnergy(float energy, float maxEnergy);
    QColor ColorForProcess(const std::string& process);
    QString NameForPDG(int pdg);
    QString SymbolForPDG(int pdg);

    // ========== Arrow Drawing ==========
    void DrawArrow(QPainter& p, const QPointF& from, const QPointF& to,
                   float arrowSize = 8.0f);
    void DrawArrowHead(QPainter& p, const QPointF& tip, const QPointF& direction,
                       float size = 8.0f);
    void DrawStarMarker(QPainter& p, const QPointF& center, float size);

    // ========== View State ==========
    ViewMode m_viewMode = ViewMode::Transverse;
    float m_zoom = 1.0f;
    float m_rotationX = 30.0f;
    float m_rotationY = -45.0f;
    float m_panX = 0.0f;
    float m_panY = 0.0f;
    QPoint m_lastMousePos;
    bool m_isDragging = false;
    bool m_isRotating = false;

    // ========== Display Options ==========
    bool m_showTracks = true;
    bool m_showCaloHits = true;
    bool m_showDetector = true;
    bool m_showGrid = true;
    bool m_showAxes = true;
    bool m_showOpticalPhotons = false;
    bool m_showVertices = true;
    bool m_showTrackArrows = true;
    bool m_showLabels = true;
    bool m_showLegend = true;
    bool m_showScaleBar = true;
    bool m_showRadiationLength = false;
    bool m_showEnergyHeatmap = false;

    // Track filtering
    float m_minTrackEnergy = 0.1f;  // MeV
    bool m_showElectrons = true;
    bool m_showGammas = true;
    bool m_showNeutrons = true;
    bool m_showProtons = true;
    bool m_showPions = true;
    bool m_showMuons = true;
    bool m_showOther = true;

    // ========== Geometry Data ==========
    bool m_geometryLoaded = false;
    std::vector<LeadGlassBlock> m_leadGlassBlocks;
    std::vector<LeadGlassBlock> m_leadGlassFB;
    std::vector<ScintillatorModule> m_scintModules;
    std::vector<TPCModule> m_tpcModules;

    // Shielding geometry (calculated)
    float m_shieldTopX, m_shieldTopY, m_shieldTopZ;
    float m_shieldSideX, m_shieldSideY, m_shieldSideZ;
    float m_shieldFBX, m_shieldFBY, m_shieldFBZ;

    // ========== Event Data ==========
    std::mutex m_mutex;
    EventDisplayData m_eventData;
    int m_pendingUpdates = 0;
    QTimer* m_refreshTimer = nullptr;

    // ========== Selection State ==========
    Track* m_selectedTrack = nullptr;
    InteractionVertex* m_selectedVertex = nullptr;
    CaloHit* m_selectedHit = nullptr;
    QPoint m_tooltipPos;
    QString m_tooltipText;

    // ========== Constants ==========
    static constexpr float SCALE = 0.003f;  // cm to display units
    static constexpr int ARROW_SPACING = 50;  // Pixels between arrows on tracks
    static constexpr float VERTEX_MARKER_SIZE = 8.0f;

    // ========== Geometry Parameters (from GeometryParameters singleton) ==========
    float BEAMPIPE_INNER() const { return GeometryParameters::Instance().BeampipeInnerRadius(); }
    float BEAMPIPE_OUTER() const { return GeometryParameters::Instance().BeampipeOuterRadius(); }
    float BEAMPIPE_HALF_Z() const { return GeometryParameters::Instance().BeampipeHalfZ(); }
    float TPC_DRIFT() const { return GeometryParameters::Instance().TPCDriftLength(); }
    float TPC_HALF_Z() const { return GeometryParameters::Instance().TPCHalfZ(); }
    float TPC_TYPE_II_WIDTH() const { return GeometryParameters::Instance().TPCTypeIIWidth(); }
    float TPC_TYPE_II_HEIGHT() const { return GeometryParameters::Instance().TPCTypeIIHeight(); }
    float TPC_TYPE_II_Y() const { return GeometryParameters::Instance().TPCTypeIIY(); }
    float TPC_TYPE_I_WIDTH() const { return GeometryParameters::Instance().TPCTypeIWidth(); }
    float TPC_TYPE_I_HEIGHT() const { return GeometryParameters::Instance().TPCTypeIHeight(); }
    float TPC_TYPE_I_X() const { return GeometryParameters::Instance().TPCTypeIX(); }
    float TPC_TYPE_I_Y() const { return GeometryParameters::Instance().TPCTypeIY(); }
    float SCINT_DIST() const { return GeometryParameters::Instance().ScintillatorDistance(); }
    float SCINT_HALF_Z() const { return GeometryParameters::Instance().ScintillatorHalfZ(); }
    float CALO_SIDE_DIST() const { return GeometryParameters::Instance().CaloSideDistance(); }
    float CALO_FB_DIST() const { return GeometryParameters::Instance().CaloFBDistance(); }
    float CALO_HALF_Y() const { return GeometryParameters::Instance().CaloHalfY(); }
    float SHIELD_HALF_X() const { return GeometryParameters::Instance().ShieldHalfX(); }
    float SHIELD_HALF_Y() const { return GeometryParameters::Instance().ShieldHalfY(); }
    float SHIELD_HALF_Z() const { return GeometryParameters::Instance().ShieldHalfZ(); }
};

} // namespace nnbar

#endif // WITH_DASHBOARD
#endif // EVENT_DISPLAY_HH
