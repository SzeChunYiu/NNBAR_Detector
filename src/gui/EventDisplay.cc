// ============================================================================
// EventDisplay.cc
// Premium-grade ATLAS/CMS-style event display implementation
// Accurate geometry from CSV data, tracks with arrows, interaction vertices
// ============================================================================

#include "gui/EventDisplay.hh"
#include "config.h"

#if WITH_DASHBOARD

#include <QPainter>
#include <QPainterPath>
#include <QFontMetrics>
#include <QPaintEvent>
#include <QtMath>
#include <QToolTip>
#include <QApplication>
#include <QDir>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "util/GeometryManager.hh"

namespace nnbar {

// ============================================================================
// Constructor / Destructor
// ============================================================================

EventDisplay::EventDisplay(QWidget* parent)
    : QWidget(parent)
{
    setMinimumSize(500, 400);
    setMouseTracking(true);
    setAttribute(Qt::WA_OpaquePaintEvent);
    setFocusPolicy(Qt::StrongFocus);

    // Setup adaptive refresh timer
    m_refreshTimer = new QTimer(this);
    connect(m_refreshTimer, &QTimer::timeout, this, [this]() {
        if (m_pendingUpdates > 0) {
            update();
            m_pendingUpdates = 0;
        }
    });
    m_refreshTimer->start(33);  // ~30 FPS max

    // Load geometry data from CSV files
    LoadGeometryData();
}

EventDisplay::~EventDisplay() {
    if (m_refreshTimer) {
        m_refreshTimer->stop();
        delete m_refreshTimer;
    }
}

// ============================================================================
// Geometry Loading from CSV Files
// ============================================================================

void EventDisplay::LoadGeometryData() {
    LoadLeadGlassPositions();
    LoadScintillatorPositions();
    CalculateTPCModules();
    CalculateShieldingGeometry();

    m_geometryLoaded = true;
    emit geometryLoaded();
}

bool EventDisplay::LoadCSV(const std::string& filename, std::vector<std::vector<double>>& csvData) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        qWarning() << "EventDisplay: Cannot open CSV file:" << QString::fromStdString(filename);
        return false;
    }

    csvData.clear();
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::vector<double> row;
        std::istringstream iss(line);
        std::string token;

        while (std::getline(iss, token, ',')) {
            try {
                row.push_back(std::stod(token));
            } catch (...) {
                row.push_back(0.0);
            }
        }
        if (!row.empty()) {
            csvData.push_back(row);
        }
    }

    file.close();
    return !csvData.empty();
}

void EventDisplay::LoadLeadGlassPositions() {
    m_leadGlassBlocks.clear();
    m_leadGlassFB.clear();

    // ========================================================================
    // NEW ALGORITHM: Generate lead glass positions algorithmically
    // All blocks on 6 sides of detector, all pointing toward center
    // ========================================================================

    // Lead glass block dimensions (cm)
    const double blockW = 8.0;   // x and z dimensions
    const double blockD = 25.0;  // depth (y dimension - radial)

    // Detector geometry parameters
    const double sideDistance = 262.0;    // Radial distance for side blocks (cm)
    const double fbDistance = 322.0;      // Z distance for front/back blocks (cm)
    const double zHalfLength = 260.0;     // Half-length of side coverage in Z (cm)
    const double xyHalfWidth = 260.0;     // Half-width of FB coverage in X/Y (cm)

    // Spacing between blocks (center-to-center)
    const double blockSpacing = 8.5;      // cm (slightly more than block width for gaps)

    int globalIndex = 0;

    // ========================================================================
    // SIDE SURFACES (4 surfaces: +Y, -X, -Y, +X rotated by 90° each)
    // ========================================================================
    for (int surface = 0; surface < 4; surface++) {
        double surfaceAngle = surface * 90.0;  // degrees
        double surfaceAngleRad = surfaceAngle * M_PI / 180.0;

        // Number of blocks along X (tangential) and Z directions
        int nBlocksX = static_cast<int>(2.0 * xyHalfWidth / blockSpacing);
        int nBlocksZ = static_cast<int>(2.0 * zHalfLength / blockSpacing);

        for (int ix = 0; ix < nBlocksX; ix++) {
            // X position in local coordinates (tangent to surface)
            double localX = -xyHalfWidth + blockSpacing/2 + ix * blockSpacing;

            for (int iz = 0; iz < nBlocksZ; iz++) {
                // Z position
                double z = -zHalfLength + blockSpacing/2 + iz * blockSpacing;

                // Local position: x is along surface, y is radial distance
                double localY = sideDistance + blockD / 2.0;  // Center of block

                // Transform to global coordinates by rotating around Z axis
                double x = localX * std::cos(surfaceAngleRad) - localY * std::sin(surfaceAngleRad);
                double y = localX * std::sin(surfaceAngleRad) + localY * std::cos(surfaceAngleRad);

                // Calculate rotation to point toward center
                // The block's "front" (PMT side) should face outward, body points inward
                double rotZ = -surfaceAngle;  // Rotate to face the surface direction

                // Calculate tilt angle based on position (blocks tilt toward center)
                double tiltX = std::atan2(localX, localY) * 180.0 / M_PI;
                double tiltZ = std::atan2(z, localY) * 180.0 / M_PI;

                LeadGlassBlock block;
                block.x = x;
                block.y = y;
                block.z = z;
                block.rotX = -tiltZ;      // Tilt in XZ plane toward center
                block.rotY = 0;
                block.rotZ = rotZ + tiltX; // Rotation in XY plane
                block.surface = surface;
                block.index = globalIndex++;

                m_leadGlassBlocks.push_back(block);
            }
        }
    }

    // ========================================================================
    // FRONT/BACK SURFACES (2 surfaces: -Z and +Z)
    // ========================================================================
    // Number of blocks in X and Y for front/back
    int nBlocksXY = static_cast<int>(2.0 * xyHalfWidth / blockSpacing);

    for (int iz = 0; iz < 2; iz++) {  // 0 = front (-Z), 1 = back (+Z)
        double zSign = (iz == 0) ? -1.0 : 1.0;
        double zPos = zSign * (fbDistance + blockD / 2.0);

        for (int ix = 0; ix < nBlocksXY; ix++) {
            double x = -xyHalfWidth + blockSpacing/2 + ix * blockSpacing;

            for (int iy = 0; iy < nBlocksXY; iy++) {
                double y = -xyHalfWidth + blockSpacing/2 + iy * blockSpacing;

                // Skip blocks that are inside the beampipe region
                double r = std::sqrt(x * x + y * y);
                if (r < 120.0) continue;  // Skip inner region (beampipe area)

                // Calculate rotation to point toward center (along Z axis)
                double tiltX = std::atan2(y, fbDistance) * 180.0 / M_PI;
                double tiltY = std::atan2(x, fbDistance) * 180.0 / M_PI;

                LeadGlassBlock block;
                block.x = x;
                block.y = y;
                block.z = zPos;
                block.rotX = zSign * (90.0 - tiltX);  // Point toward center
                block.rotY = 0;
                block.rotZ = tiltY;
                block.surface = (iz == 0) ? 4 : 5;  // 4 = front, 5 = back
                block.index = globalIndex++;

                m_leadGlassFB.push_back(block);
            }
        }
    }

    qDebug() << "EventDisplay: Generated" << m_leadGlassBlocks.size() << "side lead glass blocks";
    qDebug() << "EventDisplay: Generated" << m_leadGlassFB.size() << "front/back lead glass blocks";
    qDebug() << "EventDisplay: Total lead glass blocks:" << (m_leadGlassBlocks.size() + m_leadGlassFB.size());
}

void EventDisplay::LoadScintillatorPositions() {
    // Scintillator positions are calculated from geometry parameters
    // Matching Scintillator_geometry.cc calculations
    m_scintModules.clear();

    // Parameters from geometry code
    const double beampipeRadius = 114.0;  // cm (Beampipe_5_radius_2)
    const double tpcDrift = 85.0;         // cm
    const double tpcWall = 0.2;           // cm
    const double scintModuleX = 40.0;     // cm (4 bars * 10cm)
    const double scintModuleY = 30.0;     // cm (10 layers * 3cm)
    const double scintModuleZ = 40.0;     // cm
    const double dy = 20.0;               // cm gap

    const double scintBaseDist = beampipeRadius + tpcDrift + 2.0 * tpcWall + dy + scintModuleY / 2.0;
    const int nModulesX = 10;
    const int nModulesZ = 11;

    const double totalWidthX = 2.0 * beampipeRadius + 2.0 * tpcDrift + 4.0 * tpcWall + dy + scintModuleY;
    const double totalWidthZ = 500.0 + 2.0 * 2.0;  // Beampipe_5_len + 2*Beampipe_6_len

    const double dz = (totalWidthZ - nModulesZ * scintModuleZ) / (nModulesZ - 1.0);
    const double dx = (totalWidthX - nModulesX * scintModuleX) / (nModulesX - 1.0);

    int index = 0;

    // Generate positions for 4 surfaces (top, right, bottom, left)
    for (int surface = 0; surface < 4; surface++) {
        for (int ix = 0; ix < nModulesX; ix++) {
            for (int iz = 0; iz < nModulesZ; iz++) {
                double localX = -(nModulesX * scintModuleX + (nModulesX - 1) * dx) / 2.0
                               + 0.5 * scintModuleX + ix * (scintModuleX + dx);
                double localZ = -(nModulesZ * scintModuleZ + (nModulesZ - 1) * dz) / 2.0
                               + 0.5 * scintModuleZ + iz * (scintModuleZ + dz);

                ScintillatorModule mod;
                mod.index = index++;
                mod.surface = surface;

                // Apply rotation for each surface
                switch (surface) {
                    case 0:  // Top
                        mod.x = localX;
                        mod.y = scintBaseDist;
                        mod.z = localZ;
                        break;
                    case 1:  // Right
                        mod.x = scintBaseDist;
                        mod.y = localX;
                        mod.z = localZ;
                        break;
                    case 2:  // Bottom
                        mod.x = localX;
                        mod.y = -scintBaseDist;
                        mod.z = localZ;
                        break;
                    case 3:  // Left
                        mod.x = -scintBaseDist;
                        mod.y = localX;
                        mod.z = localZ;
                        break;
                }

                m_scintModules.push_back(mod);
            }
        }
    }

    qDebug() << "EventDisplay: Generated" << m_scintModules.size() << "scintillator module positions";
}

void EventDisplay::CalculateTPCModules() {
    m_tpcModules.clear();

    // TPC parameters from geometry (TPC_geometry.cc)
    // Beampipe_5_radius_2 = 114 cm
    // TPC_drift_len = 85 cm
    // TPC_wall_thickness = 0.2 cm
    const double beampipeR = 114.0;
    const double tpcDrift = 85.0;
    const double tpcWall = 0.2;
    const double tpcHalfZ = 252.0;     // cm (Beampipe_5_len + 2*Beampipe_6_len)/2

    // Type II (horizontal) - TPC_x_2 = 2*beampipeR + 2*wall, TPC_y_2 = drift + 2*wall
    const double typeII_hx = beampipeR + tpcWall;         // 114.2 cm half-width
    const double typeII_hy = (tpcDrift + 2*tpcWall) / 2;  // 42.7 cm half-height
    const double typeII_hz = tpcHalfZ;
    const double typeII_Y = beampipeR + typeII_hy;        // 156.7 cm center position

    // Type I (vertical) - TPC_x_1 = drift + 2*wall, TPC_y_1 = (2*beampipeR + 2*typeII_h)/2
    const double typeI_hx = (tpcDrift + 2*tpcWall) / 2;   // 42.7 cm half-width
    const double typeI_hy = (2*beampipeR + 2*(tpcDrift + 2*tpcWall)) / 4;  // 99.7 cm half-height
    const double typeI_hz = tpcHalfZ;
    const double typeI_X = typeII_hx + typeI_hx;          // 156.9 cm X position
    const double typeI_Y = typeI_hy;                      // 99.7 cm Y position

    // Add Type II modules (horizontal)
    m_tpcModules.push_back({0, typeII_Y, 0, typeII_hx, typeII_hy, typeII_hz, true, 0});
    m_tpcModules.push_back({0, -typeII_Y, 0, typeII_hx, typeII_hy, typeII_hz, true, 0});

    // Add Type I modules (vertical corners)
    m_tpcModules.push_back({-typeI_X, typeI_Y, 0, typeI_hx, typeI_hy, typeI_hz, false, 0});
    m_tpcModules.push_back({typeI_X, typeI_Y, 0, typeI_hx, typeI_hy, typeI_hz, false, 0});
    m_tpcModules.push_back({-typeI_X, -typeI_Y, 0, typeI_hx, typeI_hy, typeI_hz, false, 0});
    m_tpcModules.push_back({typeI_X, -typeI_Y, 0, typeI_hx, typeI_hy, typeI_hz, false, 0});

    qDebug() << "EventDisplay: Configured" << m_tpcModules.size() << "TPC modules";
}

void EventDisplay::CalculateShieldingGeometry() {
    // From Cosmic_Shielding_geometry.cc
    const double shieldOffset = 50.0;      // cm
    const double shieldThickness = 200.0;  // cm (2m)
    const double distCenterLG = 275.0;     // cm (2.75m)
    const double distCenterLGFB = 335.0;   // cm (3.35m)

    m_shieldTopX = 2.0 * distCenterLG + 2.0 * shieldOffset;
    m_shieldTopY = shieldThickness;
    m_shieldTopZ = 2.0 * distCenterLGFB + 2.0 * shieldOffset;

    m_shieldSideX = shieldThickness;
    m_shieldSideY = 2.0 * (distCenterLG + shieldOffset + m_shieldTopY);
    m_shieldSideZ = m_shieldTopZ;

    m_shieldFBX = m_shieldTopX + 2.0 * m_shieldSideX;
    m_shieldFBY = m_shieldSideY;
    m_shieldFBZ = shieldThickness;
}

// ============================================================================
// Coordinate Transformation
// ============================================================================

QPointF EventDisplay::ToScreen(float x, float y) {
    float cx = width() / 2.0f + m_panX;
    float cy = height() / 2.0f + m_panY;
    float scale = std::min(width(), height()) * 0.4f * m_zoom * SCALE;

    return QPointF(cx + x * scale, cy - y * scale);
}

QPointF EventDisplay::ToScreen3D(float x, float y, float z) {
    float cx = width() / 2.0f + m_panX;
    float cy = height() / 2.0f + m_panY;
    float scale = std::min(width(), height()) * 0.3f * m_zoom * SCALE;

    float radX = m_rotationX * M_PI / 180.0f;
    float radY = m_rotationY * M_PI / 180.0f;

    // Rotate around Y axis
    float x1 = x * std::cos(radY) - z * std::sin(radY);
    float z1 = x * std::sin(radY) + z * std::cos(radY);

    // Rotate around X axis
    float y1 = y * std::cos(radX) - z1 * std::sin(radX);
    float z2 = y * std::sin(radX) + z1 * std::cos(radX);

    // Perspective projection
    float perspective = 1.0f / (1.0f + z2 * 0.0002f);
    return QPointF(cx + x1 * scale * perspective, cy - y1 * scale * perspective);
}

float EventDisplay::ToScreenSize(float size) {
    float scale = std::min(width(), height()) * 0.4f * m_zoom * SCALE;
    return size * scale;
}

// ============================================================================
// Main Paint Function
// ============================================================================

void EventDisplay::paintEvent(QPaintEvent* event) {
    Q_UNUSED(event);

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::TextAntialiasing);

    // Draw layers in order
    DrawBackground(painter);
    if (m_showGrid) DrawGrid(painter);
    if (m_showEnergyHeatmap) DrawEnergyHeatmap(painter);
    if (m_showOpticalPhotons) DrawOpticalPhotonDensity(painter);
    if (m_showDetector) DrawDetectorGeometry(painter);
    if (m_showCaloHits) DrawCaloHits(painter);
    if (m_showVertices) DrawInteractionVertices(painter);
    if (m_showTracks) DrawTracks(painter);
    if (m_showAxes) DrawAxes(painter);

    // Information overlays
    DrawEventInfo(painter);
    DrawStatisticsPanel(painter);
    if (m_showLegend) DrawLegend(painter);
    if (m_showScaleBar) DrawScaleBar(painter);
    DrawLogo(painter);
    DrawSelectedInfo(painter);

    if (!m_tooltipText.isEmpty()) {
        DrawTooltip(painter);
    }
}

void EventDisplay::DrawBackground(QPainter& p) {
    // Dark gradient background
    QLinearGradient gradient(0, 0, 0, height());
    gradient.setColorAt(0, QColor(8, 8, 18));
    gradient.setColorAt(1, QColor(2, 2, 8));
    p.fillRect(rect(), gradient);
}

void EventDisplay::DrawGrid(QPainter& p) {
    p.setPen(QPen(QColor(25, 25, 35), 1));

    float cx = width() / 2.0f + m_panX;
    float cy = height() / 2.0f + m_panY;

    // Adaptive grid size based on zoom
    float gridSize = 50.0f;
    if (m_zoom > 2.0f) gridSize = 25.0f;
    if (m_zoom < 0.5f) gridSize = 100.0f;

    for (float x = cx; x < width(); x += gridSize) {
        p.drawLine(QPointF(x, 0), QPointF(x, height()));
        p.drawLine(QPointF(2*cx - x, 0), QPointF(2*cx - x, height()));
    }
    for (float y = cy; y < height(); y += gridSize) {
        p.drawLine(QPointF(0, y), QPointF(width(), y));
        p.drawLine(QPointF(0, 2*cy - y), QPointF(width(), 2*cy - y));
    }
}

// ============================================================================
// Detector Geometry Drawing
// ============================================================================

void EventDisplay::DrawDetectorGeometry(QPainter& p) {
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_viewMode == ViewMode::Transverse) {
        DrawDetectorXY(p);
    } else if (m_viewMode == ViewMode::Longitudinal) {
        DrawDetectorRZ(p);
    } else {
        DrawDetector3D(p);
    }
}

void EventDisplay::DrawDetectorXY(QPainter& p) {
    // Draw from outside in
    DrawCosmicShielding(p);
    DrawLeadGlassBlocks(p);
    DrawScintillatorModules(p);
    DrawTPCModules(p);
    DrawBeampipe(p);
}

void EventDisplay::DrawDetectorRZ(QPainter& p) {
    DrawCosmicShielding(p);
    DrawLeadGlassBlocks(p);
    DrawScintillatorModules(p);
    DrawTPCModules(p);
    DrawBeampipe(p);
}

void EventDisplay::DrawBeampipe(QPainter& p) {
    QPointF center = ToScreen(0, 0);
    float outerR = ToScreenSize(BEAMPIPE_OUTER());
    float innerR = ToScreenSize(BEAMPIPE_INNER());

    if (m_viewMode == ViewMode::Transverse) {
        // Transparent outline only - no fill
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(QColor(150, 150, 160), 2));
        p.drawEllipse(center, outerR, outerR);

        // Inner boundary
        p.setPen(QPen(QColor(100, 100, 110), 1));
        p.drawEllipse(center, innerR, innerR);
    } else if (m_viewMode == ViewMode::Longitudinal) {
        float halfZ = BEAMPIPE_HALF_Z();
        float outer = BEAMPIPE_OUTER();
        float inner = BEAMPIPE_INNER();

        // Transparent outline only
        p.setBrush(Qt::NoBrush);
        p.setPen(QPen(QColor(140, 140, 150), 1));

        // Top wall
        QPointF tl = ToScreen(-halfZ, outer);
        QPointF br = ToScreen(halfZ, inner);
        p.drawRect(QRectF(tl, br));

        // Bottom wall
        tl = ToScreen(-halfZ, -inner);
        br = ToScreen(halfZ, -outer);
        p.drawRect(QRectF(tl, br));
    }
}

void EventDisplay::DrawTPCModules(QPainter& p) {
    // Calculate fill color based on optical photon count (log scale)
    int photonCount = m_showOpticalPhotons ? m_eventData.photonsTPC : 0;
    QColor fillColor = Qt::transparent;
    if (photonCount > 0) {
        // Log-scale brightness: intensity = log10(count + 1) / log10(max + 1)
        float logCount = std::log10(static_cast<float>(photonCount) + 1.0f);
        float logMax = std::log10(10000.0f);  // Assume max ~10000 photons
        float intensity = std::min(1.0f, logCount / logMax);
        int alpha = static_cast<int>(30 + 150 * intensity);
        int green = static_cast<int>(100 + 155 * intensity);
        fillColor = QColor(0, green, green, alpha);
    }

    p.setPen(QPen(QColor(0, 220, 220), 1));

    for (const auto& tpc : m_tpcModules) {
        if (m_viewMode == ViewMode::Transverse) {
            QPointF tl = ToScreen(tpc.cx - tpc.hx, tpc.cy + tpc.hy);
            QPointF br = ToScreen(tpc.cx + tpc.hx, tpc.cy - tpc.hy);
            QRectF rect(tl, br);

            // Fill with optical photon color if enabled
            if (photonCount > 0) {
                p.setBrush(fillColor);
            } else {
                p.setBrush(Qt::NoBrush);
            }
            p.drawRect(rect);

            // Draw drift layers as faint lines
            p.setBrush(Qt::NoBrush);
            p.setPen(QPen(QColor(0, 180, 180, 50), 1));
            int nLayers = 6;
            for (int i = 1; i < nLayers; i++) {
                float frac = i / float(nLayers);
                if (tpc.isTypeII) {
                    float y = tpc.cy + tpc.hy - frac * 2 * tpc.hy;
                    p.drawLine(ToScreen(tpc.cx - tpc.hx, y), ToScreen(tpc.cx + tpc.hx, y));
                } else {
                    float x = tpc.cx - tpc.hx + frac * 2 * tpc.hx;
                    p.drawLine(ToScreen(x, tpc.cy + tpc.hy), ToScreen(x, tpc.cy - tpc.hy));
                }
            }
            p.setPen(QPen(QColor(0, 220, 220), 1));

            // Label
            if (m_showLabels && m_zoom > 0.8f) {
                p.setPen(QColor(0, 200, 200));
                p.setFont(QFont("Sans", 7));
                QString label = tpc.isTypeII ? "TPC-II" : "TPC-I";
                p.drawText(rect.center() + QPointF(-15, 4), label);
            }
        } else if (m_viewMode == ViewMode::Longitudinal) {
            // RZ view - draw as rectangles at top and bottom
            float r = std::sqrt(tpc.cx * tpc.cx + tpc.cy * tpc.cy);
            if (photonCount > 0) {
                p.setBrush(fillColor);
            } else {
                p.setBrush(Qt::NoBrush);
            }
            if (tpc.cy > 0) {
                QPointF tl = ToScreen(-tpc.hz, r + tpc.hy);
                QPointF br = ToScreen(tpc.hz, r - tpc.hy);
                p.drawRect(QRectF(tl, br));
            } else if (tpc.cy < 0) {
                QPointF tl = ToScreen(-tpc.hz, -r + tpc.hy);
                QPointF br = ToScreen(tpc.hz, -r - tpc.hy);
                p.drawRect(QRectF(tl, br));
            }
        }
    }
}

void EventDisplay::DrawScintillatorModules(QPainter& p) {
    // Exact scintillator geometry from Scintillator_geometry.cc
    const double beampipeR = 114.0;       // Beampipe_5_radius_2
    const double tpcDrift = 85.0;         // TPC_drift_len
    const double tpcWall = 0.2;           // TPC_wall_thickness
    const double scintModX = 40.0;        // scint_module_x (4 bars * 10cm)
    const double scintModY = 30.0;        // scint_module_y (10 layers * 3cm)
    const double dy = 20.0;               // gap between TPC and scintillator
    const int nModX = 10;                 // n_scint_module_x
    const int nModZ = 11;                 // n_scint_module_z

    // Individual layer dimensions (stacked in radial direction)
    const double layerHeight = 3.0;       // scint_y = 3 cm
    const int nLayersPerModule = 10;      // 10 layers per module depth

    // Calculate derived values exactly as in G4 code
    const double scintBaseDist = beampipeR + tpcDrift + 2.0 * tpcWall;  // 199.4 cm
    const double totalWidthX = 2.0 * beampipeR + 2.0 * tpcDrift + 4.0 * tpcWall + dy + scintModY;  // 448.8 cm
    const double widthX = 2.0 * beampipeR + 2.0 * tpcDrift + 4.0 * tpcWall;  // 398.8 cm

    const double dx = (totalWidthX - nModX * scintModX) / (nModX - 1.0);  // gap in X

    // Group position offsets
    const double xOffset = (totalWidthX - widthX) / 2.0;  // 25 cm
    const double scintY = scintBaseDist + scintModY / 2.0 + dy;  // 234.4 cm

    // Calculate fill color based on optical photon count (log scale)
    int photonCount = m_showOpticalPhotons ? m_eventData.photonsScintillator : 0;
    QColor fillColor = Qt::transparent;
    if (photonCount > 0) {
        // Log-scale brightness: intensity = log10(count + 1) / log10(max + 1)
        float logCount = std::log10(static_cast<float>(photonCount) + 1.0f);
        float logMax = std::log10(100000.0f);  // Scintillators produce more photons
        float intensity = std::min(1.0f, logCount / logMax);
        int alpha = static_cast<int>(30 + 150 * intensity);
        int green = static_cast<int>(100 + 155 * intensity);
        fillColor = QColor(30, green, 50, alpha);
    }

    if (photonCount > 0) {
        p.setBrush(fillColor);
    } else {
        p.setBrush(Qt::NoBrush);
    }
    p.setPen(QPen(QColor(60, 200, 60, 150), 1));

    if (m_viewMode == ViewMode::Transverse) {
        // Draw individual layers (bars stacked in radial direction)
        for (int i = 0; i < nModX; i++) {
            for (int j = 0; j < nModZ; j++) {
                // Local position within the group
                double localX = -(nModX * scintModX + (nModX - 1) * dx) / 2.0
                              + 0.5 * scintModX + i * (scintModX + dx);

                // Draw layers for each surface
                for (int layer = 0; layer < nLayersPerModule; layer++) {
                    double layerOffset = -scintModY/2 + layerHeight/2 + layer * layerHeight;

                    // Surface 0: TOP - layers stack upward
                    double modX = localX + xOffset;
                    double modY = scintY + layerOffset;
                    QPointF tl = ToScreen(modX - scintModX/2, modY + layerHeight/2);
                    QPointF br = ToScreen(modX + scintModX/2, modY - layerHeight/2);
                    p.drawRect(QRectF(tl, br));

                    // Surface 1: LEFT - layers stack leftward
                    modX = -scintY - layerOffset;
                    modY = localX + xOffset;
                    tl = ToScreen(modX - layerHeight/2, modY + scintModX/2);
                    br = ToScreen(modX + layerHeight/2, modY - scintModX/2);
                    p.drawRect(QRectF(tl, br));

                    // Surface 2: BOTTOM - layers stack downward
                    modX = localX - xOffset;
                    modY = -scintY - layerOffset;
                    tl = ToScreen(modX - scintModX/2, modY + layerHeight/2);
                    br = ToScreen(modX + scintModX/2, modY - layerHeight/2);
                    p.drawRect(QRectF(tl, br));

                    // Surface 3: RIGHT - layers stack rightward
                    modX = scintY + layerOffset;
                    modY = localX - xOffset;
                    tl = ToScreen(modX - layerHeight/2, modY + scintModX/2);
                    br = ToScreen(modX + layerHeight/2, modY - scintModX/2);
                    p.drawRect(QRectF(tl, br));
                }
            }
        }
    } else if (m_viewMode == ViewMode::Longitudinal) {
        // Draw individual layers in RZ view
        const float halfZ = 252.0f;

        for (int layer = 0; layer < nLayersPerModule; layer++) {
            double layerOffset = -scintModY/2 + layerHeight/2 + layer * layerHeight;

            // Top layers
            QPointF tl = ToScreen(-halfZ, scintY + layerOffset + layerHeight/2);
            QPointF br = ToScreen(halfZ, scintY + layerOffset - layerHeight/2);
            p.drawRect(QRectF(tl, br));

            // Bottom layers
            tl = ToScreen(-halfZ, -scintY - layerOffset + layerHeight/2);
            br = ToScreen(halfZ, -scintY - layerOffset - layerHeight/2);
            p.drawRect(QRectF(tl, br));
        }
    }
}

void EventDisplay::DrawLeadGlassBlocks(QPainter& p) {
    const float blockW = 8.0f;   // cm (x and z dimensions)
    const float blockD = 25.0f;  // cm depth (y dimension - along radial direction)

    auto& geoMgr = GeometryManager::Instance();

    // Helper to get color for a specific block based on its energy deposit
    auto getBlockColor = [&](int blockIndex) -> QColor {
        // Default: no fill unless showing heatmap or optical photons
        if (!m_showEnergyHeatmap && !m_showOpticalPhotons) {
            return Qt::transparent;
        }

        const LeadGlassInfo* info = geoMgr.GetLeadGlassBlock(blockIndex);
        if (!info) return Qt::transparent;

        double energy = info->energyDeposit;
        int photons = info->nPhotons;

        if (m_showEnergyHeatmap && energy > 0.01) {
            // Blue-to-Red colormap for energy (log scale)
            float logE = std::log10(static_cast<float>(energy) + 1.0f);
            float logMax = std::log10(100.0f);  // Max ~100 MeV per block
            float intensity = std::min(1.0f, logE / logMax);

            int r, g, b;
            if (intensity < 0.25f) {
                float t = intensity / 0.25f;
                r = 0; g = static_cast<int>(255 * t); b = 255;
            } else if (intensity < 0.5f) {
                float t = (intensity - 0.25f) / 0.25f;
                r = 0; g = 255; b = static_cast<int>(255 * (1 - t));
            } else if (intensity < 0.75f) {
                float t = (intensity - 0.5f) / 0.25f;
                r = static_cast<int>(255 * t); g = 255; b = 0;
            } else {
                float t = (intensity - 0.75f) / 0.25f;
                r = 255; g = static_cast<int>(255 * (1 - t)); b = 0;
            }
            return QColor(r, g, b, 180);
        } else if (m_showOpticalPhotons && photons > 0) {
            // Yellow-orange for photon count
            float logCount = std::log10(static_cast<float>(photons) + 1.0f);
            float logMax = std::log10(5000.0f);
            float intensity = std::min(1.0f, logCount / logMax);
            int alpha = static_cast<int>(50 + 200 * intensity);
            int red = static_cast<int>(150 + 105 * intensity);
            int green = static_cast<int>(100 + 80 * intensity);
            return QColor(red, green, 0, alpha);
        }
        return Qt::transparent;
    };

    p.setPen(QPen(QColor(200, 160, 0, 150), 1));

    if (m_viewMode == ViewMode::Transverse) {
        // Draw side lead glass blocks from loaded CSV data
        for (const auto& block : m_leadGlassBlocks) {
            QColor fillColor = getBlockColor(block.index);
            p.setBrush(fillColor.alpha() > 0 ? fillColor : Qt::NoBrush);

            QPointF screenPos = ToScreen(block.x, block.y);
            float screenSize = ToScreenSize(blockW);

            p.save();
            p.translate(screenPos);
            p.rotate(-block.rotZ);
            p.drawRect(QRectF(-screenSize/2, -screenSize/2, screenSize, screenSize));
            p.restore();
        }
    } else if (m_viewMode == ViewMode::Longitudinal) {
        // Draw in RZ view
        for (const auto& block : m_leadGlassBlocks) {
            QColor fillColor = getBlockColor(block.index);
            p.setBrush(fillColor.alpha() > 0 ? fillColor : Qt::NoBrush);

            float r = std::sqrt(block.x * block.x + block.y * block.y);
            float sign = (block.y >= 0) ? 1.0f : -1.0f;

            QPointF screenPos = ToScreen(block.z, r * sign);
            float screenW = ToScreenSize(blockW);
            float screenD = ToScreenSize(blockD);

            // Draw block (depth along radial direction)
            if (sign > 0) {
                p.drawRect(QRectF(screenPos.x() - screenW/2, screenPos.y() - screenD, screenW, screenD));
            } else {
                p.drawRect(QRectF(screenPos.x() - screenW/2, screenPos.y(), screenW, screenD));
            }
        }

        // Draw front/back lead glass
        for (const auto& block : m_leadGlassFB) {
            QColor fillColor = getBlockColor(block.index);
            p.setBrush(fillColor.alpha() > 0 ? fillColor : Qt::NoBrush);

            QPointF screenPos = ToScreen(block.z, block.y);
            float screenW = ToScreenSize(blockW);
            float screenD = ToScreenSize(blockD);

            p.drawRect(QRectF(screenPos.x() - screenD/2, screenPos.y() - screenW/2, screenD, screenW));
        }
    }
}

void EventDisplay::DrawCosmicShielding(QPainter& p) {
    // Outer lead/concrete shield - transparent outline only
    const float distLG = 275.0f;
    const float offset = 50.0f;
    const float thickness = 200.0f;

    p.setBrush(Qt::NoBrush);

    if (m_viewMode == ViewMode::Transverse) {
        float outerX = distLG + offset + thickness;
        float outerY = distLG + offset + thickness;

        // Outer boundary
        QPointF tl = ToScreen(-outerX, outerY);
        QPointF br = ToScreen(outerX, -outerY);
        p.setPen(QPen(QColor(150, 60, 60, 150), 2));
        p.drawRect(QRectF(tl, br));

        // Inner steel layer (30cm inside)
        float innerX = outerX - 30.0f;
        float innerY = outerY - 30.0f;
        tl = ToScreen(-innerX, innerY);
        br = ToScreen(innerX, -innerY);
        p.setPen(QPen(QColor(100, 100, 110, 100), 1));
        p.drawRect(QRectF(tl, br));

        // Label
        if (m_showLabels) {
            p.setPen(QColor(130, 70, 70));
            p.setFont(QFont("Sans", 8));
            p.drawText(ToScreen(-outerX + 10, outerY - 10), "Cosmic Shield");
        }
    } else if (m_viewMode == ViewMode::Longitudinal) {
        float halfZ = 335.0f + offset + thickness;
        float halfY = distLG + offset + thickness;

        QPointF tl = ToScreen(-halfZ, halfY);
        QPointF br = ToScreen(halfZ, -halfY);
        p.setPen(QPen(QColor(150, 60, 60, 150), 2));
        p.drawRect(QRectF(tl, br));
    }
}

void EventDisplay::DrawDetector3D(QPainter& p) {
    // Simplified 3D wireframe view
    auto drawBox3D = [&](float cx, float cy, float cz, float hx, float hy, float hz, QColor color) {
        QPointF corners[8] = {
            ToScreen3D(cx - hx, cy - hy, cz - hz),
            ToScreen3D(cx + hx, cy - hy, cz - hz),
            ToScreen3D(cx + hx, cy + hy, cz - hz),
            ToScreen3D(cx - hx, cy + hy, cz - hz),
            ToScreen3D(cx - hx, cy - hy, cz + hz),
            ToScreen3D(cx + hx, cy - hy, cz + hz),
            ToScreen3D(cx + hx, cy + hy, cz + hz),
            ToScreen3D(cx - hx, cy + hy, cz + hz),
        };

        p.setPen(QPen(color, 1));
        p.setBrush(QColor(color.red(), color.green(), color.blue(), 20));

        // Draw faces
        QPolygonF front; front << corners[0] << corners[1] << corners[2] << corners[3];
        p.drawPolygon(front);

        // Connect front to back
        for (int i = 0; i < 4; i++) {
            p.drawLine(corners[i], corners[i + 4]);
        }
    };

    // Draw detector components
    drawBox3D(0, 0, 0, 450, 350, 400, QColor(100, 50, 50));  // Shield
    drawBox3D(0, 200, 0, 150, 40, 130, QColor(0, 200, 200));  // TPC top
    drawBox3D(0, -200, 0, 150, 40, 130, QColor(0, 200, 200)); // TPC bottom

    // 3D axes
    p.setPen(QPen(QColor(255, 80, 80), 2));
    p.drawLine(ToScreen3D(0, 0, 0), ToScreen3D(100, 0, 0));
    p.drawText(ToScreen3D(110, 0, 0), "X");

    p.setPen(QPen(QColor(80, 255, 80), 2));
    p.drawLine(ToScreen3D(0, 0, 0), ToScreen3D(0, 100, 0));
    p.drawText(ToScreen3D(0, 110, 0), "Y");

    p.setPen(QPen(QColor(80, 80, 255), 2));
    p.drawLine(ToScreen3D(0, 0, 0), ToScreen3D(0, 0, 100));
    p.drawText(ToScreen3D(0, 0, 110), "Z (beam)");
}

// ============================================================================
// Track Drawing with Directional Arrows
// ============================================================================

void EventDisplay::DrawTracks(QPainter& p) {
    std::lock_guard<std::mutex> lock(m_mutex);

    for (const auto& track : m_eventData.tracks) {
        // Apply filters
        if (track.initialEnergy < m_minTrackEnergy) continue;
        if (track.points.size() < 2) continue;

        int absPdg = std::abs(track.pdg);
        if (absPdg == 11 && !m_showElectrons) continue;
        if (absPdg == 22 && !m_showGammas) continue;
        if (absPdg == 2112 && !m_showNeutrons) continue;
        if (absPdg == 2212 && !m_showProtons) continue;
        if ((absPdg == 211 || absPdg == 111) && !m_showPions) continue;
        if (absPdg == 13 && !m_showMuons) continue;

        // Skip optical photons
        if (track.pdg == 0 || track.pdg == -22) continue;

        DrawTrackWithArrows(p, track);
    }
}

void EventDisplay::DrawTrackWithArrows(QPainter& p, const Track& track) {
    QColor color = ColorForPDG(track.pdg);
    float lineWidth = 1.5f;  // Thinner lines

    // Slightly thicker lines for primaries
    if (track.parentId == 0) {
        lineWidth = 2.0f;
        color = color.lighter(120);
    }

    // Ensure no fill - only stroke the path
    p.setBrush(Qt::NoBrush);

    // Dashed lines for neutrals (gammas, neutrons)
    bool isNeutral = (std::abs(track.pdg) == 22 || std::abs(track.pdg) == 2112);
    if (isNeutral) {
        p.setPen(QPen(color, lineWidth, Qt::DashLine));
    } else {
        p.setPen(QPen(color, lineWidth));
    }

    // Simplify trajectory by removing collinear points
    // This makes straight-line tracks render as clean straight lines
    auto simplifyTrajectory = [](const std::vector<TrackPoint>& points, float tolerance = 0.001f)
        -> std::vector<TrackPoint> {
        if (points.size() <= 2) return points;

        std::vector<TrackPoint> simplified;
        simplified.push_back(points[0]);

        for (size_t i = 1; i < points.size() - 1; i++) {
            const auto& prev = simplified.back();
            const auto& curr = points[i];
            const auto& next = points[i + 1];

            // Calculate direction vectors
            float dx1 = curr.x - prev.x;
            float dy1 = curr.y - prev.y;
            float dz1 = curr.z - prev.z;
            float dx2 = next.x - curr.x;
            float dy2 = next.y - curr.y;
            float dz2 = next.z - curr.z;

            // Normalize directions
            float len1 = std::sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
            float len2 = std::sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);

            if (len1 > 1e-6 && len2 > 1e-6) {
                dx1 /= len1; dy1 /= len1; dz1 /= len1;
                dx2 /= len2; dy2 /= len2; dz2 /= len2;

                // Calculate cross product magnitude (sin of angle)
                float crossX = dy1*dz2 - dz1*dy2;
                float crossY = dz1*dx2 - dx1*dz2;
                float crossZ = dx1*dy2 - dy1*dx2;
                float crossMag = std::sqrt(crossX*crossX + crossY*crossY + crossZ*crossZ);

                // If directions are not collinear (cross product > tolerance), keep the point
                if (crossMag > tolerance) {
                    simplified.push_back(curr);
                }
            }
        }

        simplified.push_back(points.back());
        return simplified;
    };

    // Simplify the track points to remove collinear intermediate points
    std::vector<TrackPoint> simplifiedPoints = simplifyTrajectory(track.points);

    // Convert to screen coordinates
    std::vector<QPointF> screenPoints;
    for (const auto& pt : simplifiedPoints) {
        QPointF screenPt;
        if (m_viewMode == ViewMode::Transverse) {
            screenPt = ToScreen(pt.x, pt.y);
        } else if (m_viewMode == ViewMode::Longitudinal) {
            float r = std::sqrt(pt.x * pt.x + pt.y * pt.y);
            screenPt = ToScreen(pt.z, r * (pt.y >= 0 ? 1 : -1));
        } else {
            screenPt = ToScreen3D(pt.x, pt.y, pt.z);
        }
        screenPoints.push_back(screenPt);
    }

    // Draw track path
    QPainterPath path;
    if (!screenPoints.empty()) {
        path.moveTo(screenPoints[0]);
        for (size_t i = 1; i < screenPoints.size(); i++) {
            path.lineTo(screenPoints[i]);
        }
    }

    p.drawPath(path);

    // Draw directional arrows along the track
    if (m_showTrackArrows && screenPoints.size() >= 2) {
        float accumulatedDist = 0;
        const float arrowInterval = ARROW_SPACING / m_zoom;  // Adapt to zoom

        for (size_t i = 1; i < screenPoints.size(); i++) {
            QPointF delta = screenPoints[i] - screenPoints[i-1];
            float segmentLen = std::sqrt(delta.x() * delta.x() + delta.y() * delta.y());

            accumulatedDist += segmentLen;

            if (accumulatedDist >= arrowInterval && segmentLen > 5) {
                // Draw arrow at midpoint of this segment
                QPointF midPt = (screenPoints[i-1] + screenPoints[i]) / 2.0;
                QPointF dir = delta / segmentLen;

                DrawArrowHead(p, midPt, dir, 6.0f);
                accumulatedDist = 0;
            }
        }

        // Always draw arrow at track end
        if (screenPoints.size() >= 2) {
            QPointF last = screenPoints.back();
            QPointF prev = screenPoints[screenPoints.size() - 2];
            QPointF delta = last - prev;
            float len = std::sqrt(delta.x() * delta.x() + delta.y() * delta.y());
            if (len > 3) {
                DrawArrowHead(p, last, delta / len, 8.0f);
            }
        }
    }

    // Draw start marker for primaries
    if (track.parentId == 0 && !screenPoints.empty()) {
        p.setBrush(color);
        p.setPen(QPen(Qt::white, 1));
        p.drawEllipse(screenPoints[0], 5, 5);
    }

    // Highlight selected track
    if (m_selectedTrack && m_selectedTrack->trackId == track.trackId) {
        p.setPen(QPen(Qt::yellow, lineWidth + 2));
        p.drawPath(path);
    }
}

void EventDisplay::DrawArrowHead(QPainter& p, const QPointF& tip, const QPointF& direction, float size) {
    // Calculate perpendicular direction
    QPointF perp(-direction.y(), direction.x());

    QPointF base = tip - direction * size;
    QPointF left = base + perp * size * 0.5f;
    QPointF right = base - perp * size * 0.5f;

    QPolygonF arrow;
    arrow << tip << left << right;

    p.save();
    p.setBrush(p.pen().color());
    p.setPen(Qt::NoPen);
    p.drawPolygon(arrow);
    p.restore();
}

// ============================================================================
// Interaction Vertices
// ============================================================================

void EventDisplay::DrawInteractionVertices(QPainter& p) {
    std::lock_guard<std::mutex> lock(m_mutex);

    for (const auto& vtx : m_eventData.vertices) {
        QPointF screenPos;
        if (m_viewMode == ViewMode::Transverse) {
            screenPos = ToScreen(vtx.x, vtx.y);
        } else if (m_viewMode == ViewMode::Longitudinal) {
            float r = std::sqrt(vtx.x * vtx.x + vtx.y * vtx.y);
            screenPos = ToScreen(vtx.z, r * (vtx.y >= 0 ? 1 : -1));
        } else {
            screenPos = ToScreen3D(vtx.x, vtx.y, vtx.z);
        }

        // Size based on multiplicity
        float markerSize = VERTEX_MARKER_SIZE + vtx.Multiplicity() * 1.5f;

        // Color based on process
        QColor color = ColorForProcess(vtx.process);

        // Draw starburst for high multiplicity
        if (vtx.Multiplicity() > 3) {
            int numRays = vtx.Multiplicity();
            p.setPen(QPen(color.lighter(130), 1));
            for (int i = 0; i < numRays; i++) {
                float angle = i * 2.0f * M_PI / numRays;
                QPointF rayEnd = screenPos + QPointF(std::cos(angle), std::sin(angle)) * markerSize * 1.5f;
                p.drawLine(screenPos, rayEnd);
            }
        }

        // Draw vertex marker (star shape)
        p.setPen(QPen(Qt::white, 1));
        p.setBrush(color);
        DrawStarMarker(p, screenPos, markerSize);

        // Label for significant vertices
        if (m_showLabels && vtx.Multiplicity() > 2 && m_zoom > 0.7f) {
            p.setPen(Qt::white);
            p.setFont(QFont("Sans", 7));
            QString label = QString::fromStdString(vtx.process);
            if (label.length() > 10) label = label.left(8) + "..";
            p.drawText(screenPos + QPointF(markerSize + 3, 4), label);
        }
    }
}

void EventDisplay::DrawStarMarker(QPainter& p, const QPointF& center, float size) {
    // 4-pointed star
    QPolygonF star;
    for (int i = 0; i < 8; i++) {
        float angle = i * M_PI / 4.0f - M_PI / 2.0f;
        float r = (i % 2 == 0) ? size : size * 0.4f;
        star << center + QPointF(std::cos(angle) * r, std::sin(angle) * r);
    }
    p.drawPolygon(star);
}

// ============================================================================
// Calorimeter Hits
// ============================================================================

void EventDisplay::DrawCaloHits(QPainter& p) {
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_eventData.caloHits.empty()) return;

    // Find max energy for normalization
    float maxE = 0.1f;
    for (const auto& hit : m_eventData.caloHits) {
        maxE = std::max(maxE, hit.energy);
    }

    for (const auto& hit : m_eventData.caloHits) {
        QColor color = ColorForEnergy(hit.energy, maxE);

        QPointF screenPt;
        if (m_viewMode == ViewMode::Transverse) {
            screenPt = ToScreen(hit.x, hit.y);
        } else if (m_viewMode == ViewMode::Longitudinal) {
            float r = std::sqrt(hit.x * hit.x + hit.y * hit.y);
            screenPt = ToScreen(hit.z, r * (hit.y >= 0 ? 1 : -1));
        } else {
            screenPt = ToScreen3D(hit.x, hit.y, hit.z);
        }

        // Size proportional to energy (log scale for visibility)
        float size = 3.0f + 8.0f * std::log10(1.0f + hit.energy / maxE * 10.0f) / std::log10(11.0f);

        // Different shape for different detectors
        p.setPen(Qt::NoPen);
        p.setBrush(color);

        switch (hit.detectorType) {
            case 0:  // Scintillator - square
                p.drawRect(QRectF(screenPt.x() - size/2, screenPt.y() - size/2, size, size));
                break;
            case 1:  // Lead glass - diamond
                {
                    QPolygonF diamond;
                    diamond << screenPt + QPointF(0, -size)
                           << screenPt + QPointF(size, 0)
                           << screenPt + QPointF(0, size)
                           << screenPt + QPointF(-size, 0);
                    p.drawPolygon(diamond);
                }
                break;
            case 2:  // TPC - circle
            default:
                p.drawEllipse(screenPt, size, size);
                break;
        }
    }
}

void EventDisplay::DrawOpticalPhotonDensity(QPainter& p) {
    std::lock_guard<std::mutex> lock(m_mutex);

    if (m_eventData.opticalPhotonRegions.empty()) return;

    // Find max count for log-scale normalization
    int maxCount = 1;
    for (const auto& region : m_eventData.opticalPhotonRegions) {
        maxCount = std::max(maxCount, region.photonCount);
    }

    // Use log10 scale for intensity
    float logMax = std::log10(static_cast<float>(maxCount) + 1.0f);

    for (const auto& region : m_eventData.opticalPhotonRegions) {
        // Log-scale intensity (brightness)
        float logCount = std::log10(static_cast<float>(region.photonCount) + 1.0f);
        float intensity = logMax > 0 ? logCount / logMax : 0.0f;

        // Green color with brightness based on log-scale photon count
        int green = static_cast<int>(80 + 175 * intensity);
        int alpha = static_cast<int>(40 + 160 * intensity);

        p.setPen(Qt::NoPen);
        p.setBrush(QColor(20, green, 50, alpha));

        QPointF screenPt;
        if (m_viewMode == ViewMode::Transverse) {
            screenPt = ToScreen(region.x, region.y);
        } else if (m_viewMode == ViewMode::Longitudinal) {
            float r = std::sqrt(region.x * region.x + region.y * region.y);
            screenPt = ToScreen(region.z, r * (region.y >= 0 ? 1 : -1));
        } else {
            screenPt = ToScreen3D(region.x, region.y, region.z);
        }

        // Size based on detector module size (approximate)
        float size = ToScreenSize(20.0f);
        p.drawRect(QRectF(screenPt.x() - size/2, screenPt.y() - size/2, size, size));
    }
}

void EventDisplay::DrawEnergyHeatmap(QPainter& p) {
    std::lock_guard<std::mutex> lock(m_mutex);

    // Helper function: Blue-to-Red colormap based on log-scale intensity
    auto blueToRed = [](float intensity, int alpha = 150) -> QColor {
        // intensity: 0.0 (cold/blue) to 1.0 (hot/red)
        // Uses smooth transition: blue -> cyan -> green -> yellow -> red
        intensity = std::max(0.0f, std::min(1.0f, intensity));
        int r, g, b;
        if (intensity < 0.25f) {
            // Blue to Cyan
            float t = intensity / 0.25f;
            r = 0;
            g = static_cast<int>(255 * t);
            b = 255;
        } else if (intensity < 0.5f) {
            // Cyan to Green
            float t = (intensity - 0.25f) / 0.25f;
            r = 0;
            g = 255;
            b = static_cast<int>(255 * (1 - t));
        } else if (intensity < 0.75f) {
            // Green to Yellow
            float t = (intensity - 0.5f) / 0.25f;
            r = static_cast<int>(255 * t);
            g = 255;
            b = 0;
        } else {
            // Yellow to Red
            float t = (intensity - 0.75f) / 0.25f;
            r = 255;
            g = static_cast<int>(255 * (1 - t));
            b = 0;
        }
        return QColor(r, g, b, alpha);
    };

    // Get total energy depositions per detector
    double edepTPC = m_eventData.edepTPC;
    double edepScint = m_eventData.edepScintillator;
    double edepLG = m_eventData.edepLeadGlass;
    double totalEdep = edepTPC + edepScint + edepLG + m_eventData.edepOther;

    if (totalEdep < 0.01) return;  // No significant energy deposition

    // Find max for normalization (log scale)
    double maxEdep = std::max({edepTPC, edepScint, edepLG, 1.0});
    float logMax = std::log10(static_cast<float>(maxEdep) + 1.0f);

    p.setPen(Qt::NoPen);

    // Draw TPC heatmap fill
    if (edepTPC > 0.01) {
        float logE = std::log10(static_cast<float>(edepTPC) + 1.0f);
        float intensity = logMax > 0 ? logE / logMax : 0.0f;
        QColor color = blueToRed(intensity, 100);

        for (const auto& tpc : m_tpcModules) {
            if (m_viewMode == ViewMode::Transverse) {
                QPointF tl = ToScreen(tpc.cx - tpc.hx, tpc.cy + tpc.hy);
                QPointF br = ToScreen(tpc.cx + tpc.hx, tpc.cy - tpc.hy);
                p.setBrush(color);
                p.drawRect(QRectF(tl, br));
            } else if (m_viewMode == ViewMode::Longitudinal) {
                float r = std::sqrt(tpc.cx * tpc.cx + tpc.cy * tpc.cy);
                if (tpc.cy > 0) {
                    QPointF tl = ToScreen(-tpc.hz, r + tpc.hy);
                    QPointF br = ToScreen(tpc.hz, r - tpc.hy);
                    p.setBrush(color);
                    p.drawRect(QRectF(tl, br));
                } else if (tpc.cy < 0) {
                    QPointF tl = ToScreen(-tpc.hz, -r + tpc.hy);
                    QPointF br = ToScreen(tpc.hz, -r - tpc.hy);
                    p.setBrush(color);
                    p.drawRect(QRectF(tl, br));
                }
            }
        }
    }

    // Draw Scintillator heatmap fill
    if (edepScint > 0.01) {
        float logE = std::log10(static_cast<float>(edepScint) + 1.0f);
        float intensity = logMax > 0 ? logE / logMax : 0.0f;
        QColor color = blueToRed(intensity, 120);
        p.setBrush(color);

        // Approximate scintillator regions
        const double scintDist = 234.4;  // cm from center
        const double scintThick = 30.0;
        if (m_viewMode == ViewMode::Transverse) {
            // Top/Bottom
            QPointF tl = ToScreen(-220, scintDist + scintThick/2);
            QPointF br = ToScreen(220, scintDist - scintThick/2);
            p.drawRect(QRectF(tl, br));
            tl = ToScreen(-220, -scintDist + scintThick/2);
            br = ToScreen(220, -scintDist - scintThick/2);
            p.drawRect(QRectF(tl, br));
        }
    }

    // Draw Lead Glass heatmap fill
    if (edepLG > 0.01) {
        float logE = std::log10(static_cast<float>(edepLG) + 1.0f);
        float intensity = logMax > 0 ? logE / logMax : 0.0f;
        QColor color = blueToRed(intensity, 130);
        p.setBrush(color);

        // Draw filled regions for lead glass blocks
        for (const auto& block : m_leadGlassBlocks) {
            if (m_viewMode == ViewMode::Transverse) {
                QPointF pt = ToScreen(block.x, block.y);
                float size = ToScreenSize(12.0f);
                p.drawRect(QRectF(pt.x() - size/2, pt.y() - size/2, size, size));
            }
        }
        for (const auto& block : m_leadGlassFB) {
            if (m_viewMode == ViewMode::Longitudinal) {
                float r = std::sqrt(block.x * block.x + block.y * block.y);
                QPointF pt = ToScreen(block.z, r * (block.y >= 0 ? 1 : -1));
                float size = ToScreenSize(12.0f);
                p.drawRect(QRectF(pt.x() - size/2, pt.y() - size/2, size, size));
            }
        }
    }

    // Draw energy deposition colorbar legend
    int barX = width() - 30;
    int barY = 100;
    int barH = 150;
    int barW = 15;

    // Background
    p.setBrush(QColor(20, 20, 30, 200));
    p.setPen(QPen(QColor(60, 60, 70), 1));
    p.drawRect(barX - 5, barY - 25, barW + 50, barH + 50);

    // Title
    p.setPen(QColor(200, 200, 200));
    p.setFont(QFont("Sans", 8, QFont::Bold));
    p.drawText(barX - 5, barY - 8, "Edep");

    // Draw colorbar
    for (int i = 0; i < barH; i++) {
        float intensity = 1.0f - static_cast<float>(i) / barH;
        p.setPen(blueToRed(intensity, 255));
        p.drawLine(barX, barY + i, barX + barW, barY + i);
    }

    // Labels
    p.setPen(QColor(180, 180, 180));
    p.setFont(QFont("Sans", 7));
    p.drawText(barX + barW + 3, barY + 5, QString("%1").arg(maxEdep, 0, 'f', 0));
    p.drawText(barX + barW + 3, barY + barH/2, QString("%1").arg(maxEdep/2, 0, 'f', 0));
    p.drawText(barX + barW + 3, barY + barH, "0 MeV");
}

// ============================================================================
// Information Overlays
// ============================================================================

void EventDisplay::DrawAxes(QPainter& p) {
    if (m_viewMode == ViewMode::Perspective3D) return;  // 3D has its own axes

    float cx = 60.0f;
    float cy = height() - 60.0f;
    float len = 40.0f;

    // Background circle
    p.setBrush(QColor(20, 20, 30, 180));
    p.setPen(Qt::NoPen);
    p.drawEllipse(QPointF(cx, cy), 50, 50);

    // X/Z axis (red)
    p.setPen(QPen(QColor(255, 80, 80), 2));
    p.drawLine(QPointF(cx, cy), QPointF(cx + len, cy));
    QString xLabel = (m_viewMode == ViewMode::Longitudinal) ? "Z" : "X";
    p.drawText(QPointF(cx + len + 5, cy + 5), xLabel);

    // Y/R axis (green)
    p.setPen(QPen(QColor(80, 255, 80), 2));
    p.drawLine(QPointF(cx, cy), QPointF(cx, cy - len));
    QString yLabel = (m_viewMode == ViewMode::Longitudinal) ? "R" : "Y";
    p.drawText(QPointF(cx + 5, cy - len - 5), yLabel);

    // View label
    p.setPen(QColor(150, 150, 150));
    p.setFont(QFont("Sans", 8));
    QString viewName = (m_viewMode == ViewMode::Transverse) ? "X-Y View" : "R-Z View";
    p.drawText(QPointF(cx - 25, cy + 55), viewName);
}

void EventDisplay::DrawEventInfo(QPainter& p) {
    std::lock_guard<std::mutex> lock(m_mutex);

    p.setPen(QColor(220, 220, 220));
    QFont font("Sans", 10);
    font.setBold(true);
    p.setFont(font);

    int x = 15, y = 25;
    int lineH = 18;

    // Event identification
    p.drawText(x, y, QString("Run: %1  Event: %2")
        .arg(m_eventData.runId)
        .arg(m_eventData.eventId));
    y += lineH;

    // Primary particle info
    font.setBold(false);
    font.setPointSize(9);
    p.setFont(font);
    p.setPen(QColor(180, 180, 180));

    if (!m_eventData.primaryParticle.empty()) {
        p.drawText(x, y, QString("Primary: %1 @ %2 MeV")
            .arg(QString::fromStdString(m_eventData.primaryParticle))
            .arg(m_eventData.primaryEnergy, 0, 'f', 1));
        y += lineH;
    }

    // Track and hit counts
    p.drawText(x, y, QString("Tracks: %1  Hits: %2  Vertices: %3")
        .arg(m_eventData.tracks.size())
        .arg(m_eventData.caloHits.size())
        .arg(m_eventData.vertices.size()));
    y += lineH;

    // Total energy
    p.drawText(x, y, QString("Total E_dep: %1 MeV")
        .arg(m_eventData.totalEnergy, 0, 'f', 2));
}

void EventDisplay::DrawStatisticsPanel(QPainter& p) {
    std::lock_guard<std::mutex> lock(m_mutex);

    // Draw panel in top-right corner - expanded for energy info
    int panelW = 180;
    int panelH = 320;
    int x = width() - panelW - 15;
    int y = 50;

    // Panel background
    p.setBrush(QColor(15, 15, 25, 220));
    p.setPen(QPen(QColor(60, 60, 80), 1));
    p.drawRoundedRect(x, y, panelW, panelH, 5, 5);

    // Title - Particle Statistics
    p.setPen(QColor(100, 200, 255));
    QFont titleFont("Sans", 9);
    titleFont.setBold(true);
    p.setFont(titleFont);
    p.drawText(x + 10, y + 18, "Particle Statistics");

    // Separator line
    p.setPen(QPen(QColor(60, 60, 80), 1));
    p.drawLine(x + 10, y + 25, x + panelW - 10, y + 25);

    // Particle Statistics
    QFont statFont("Sans", 8);
    p.setFont(statFont);

    int row = y + 42;
    int rowH = 15;

    auto drawStatRow = [&](const QString& label, int count, QColor color) {
        p.setPen(color);
        p.setBrush(color);
        p.drawEllipse(x + 12, row - 4, 8, 8);
        p.setPen(QColor(200, 200, 200));
        p.drawText(x + 25, row, label);
        p.drawText(x + panelW - 40, row, QString::number(count));
        row += rowH;
    };

    drawStatRow("Primaries", m_eventData.nPrimaries, QColor(255, 255, 255));
    drawStatRow("Gammas", m_eventData.nGammas, ColorForPDG(22));
    drawStatRow("Electrons", m_eventData.nElectrons, ColorForPDG(11));
    drawStatRow("Protons", m_eventData.nProtons, ColorForPDG(2212));
    drawStatRow("Neutrons", m_eventData.nNeutrons, ColorForPDG(2112));
    drawStatRow("Pions", m_eventData.nPions, ColorForPDG(211));
    drawStatRow("Muons", m_eventData.nMuons, ColorForPDG(13));
    drawStatRow("Other", m_eventData.nOther, ColorForPDG(0));

    // Section separator
    row += 5;
    p.setPen(QPen(QColor(60, 60, 80), 1));
    p.drawLine(x + 10, row, x + panelW - 10, row);
    row += 15;

    // Title - Energy Deposition
    p.setPen(QColor(255, 180, 100));
    p.setFont(titleFont);
    p.drawText(x + 10, row, "Energy Deposition");
    row += 15;

    p.setFont(statFont);

    auto drawEnergyRow = [&](const QString& label, double energy, QColor color) {
        p.setPen(color);
        p.setBrush(color);
        p.drawRect(x + 12, row - 6, 8, 8);
        p.setPen(QColor(200, 200, 200));
        p.drawText(x + 25, row, label);
        // Format energy with appropriate units
        QString energyStr;
        if (energy < 0.01) {
            energyStr = QString("%1 keV").arg(energy * 1000, 0, 'f', 1);
        } else if (energy < 1000) {
            energyStr = QString("%1 MeV").arg(energy, 0, 'f', 2);
        } else {
            energyStr = QString("%1 GeV").arg(energy / 1000, 0, 'f', 2);
        }
        p.drawText(x + panelW - 70, row, energyStr);
        row += rowH;
    };

    drawEnergyRow("TPC", m_eventData.edepTPC, QColor(0, 220, 220));
    drawEnergyRow("Scintillator", m_eventData.edepScintillator, QColor(60, 200, 60));
    drawEnergyRow("Lead Glass", m_eventData.edepLeadGlass, QColor(255, 180, 0));
    drawEnergyRow("Beampipe", m_eventData.edepBeampipe, QColor(150, 150, 160));
    drawEnergyRow("Shield", m_eventData.edepShield, QColor(150, 60, 60));
    drawEnergyRow("Other", m_eventData.edepOther, QColor(150, 100, 200));

    // Total
    row += 3;
    p.setPen(QPen(QColor(60, 60, 80), 1));
    p.drawLine(x + 10, row, x + panelW - 10, row);
    row += 12;
    p.setPen(QColor(255, 255, 200));
    p.drawText(x + 12, row, "Total:");
    QString totalStr;
    double total = m_eventData.totalEnergy;
    if (total < 1000) {
        totalStr = QString("%1 MeV").arg(total, 0, 'f', 2);
    } else {
        totalStr = QString("%1 GeV").arg(total / 1000, 0, 'f', 2);
    }
    p.drawText(x + panelW - 70, row, totalStr);
}

void EventDisplay::DrawLegend(QPainter& p) {
    // Draw in bottom-right
    int x = width() - 180;
    int y = height() - 130;

    // Background
    p.setBrush(QColor(15, 15, 25, 200));
    p.setPen(QPen(QColor(60, 60, 80), 1));
    p.drawRoundedRect(x, y, 165, 120, 5, 5);

    // Title
    p.setPen(QColor(180, 180, 180));
    p.setFont(QFont("Sans", 8, QFont::Bold));
    p.drawText(x + 10, y + 15, "Particle Legend");

    // Legend items
    struct LegendItem { const char* name; int pdg; };
    LegendItem items[] = {
        {"e-/e+", 11},
        {"gamma", 22},
        {"mu-/mu+", 13},
        {"proton", 2212},
        {"neutron", 2112},
        {"pi+/pi-", 211},
    };

    p.setFont(QFont("Sans", 8));
    int row = y + 32;

    for (const auto& item : items) {
        QColor color = ColorForPDG(item.pdg);
        p.setPen(color);
        p.drawLine(x + 10, row, x + 35, row);

        // Dashed for neutrals
        if (item.pdg == 22 || item.pdg == 2112) {
            p.setPen(QPen(color, 2, Qt::DashLine));
            p.drawLine(x + 10, row, x + 35, row);
        }

        p.setPen(QColor(200, 200, 200));
        p.drawText(x + 42, row + 4, item.name);
        row += 15;
    }
}

void EventDisplay::DrawScaleBar(QPainter& p) {
    // Draw scale bar at bottom
    int x = 15;
    int y = height() - 25;

    // Calculate scale in cm
    float pixelsPerCm = ToScreenSize(1.0f);
    float scaleLength = 100.0f;  // cm

    // Adjust scale to nice round number
    float targetPixels = 100.0f;  // Target bar length in pixels
    float cmPerBar = targetPixels / pixelsPerCm;

    // Round to nice number
    if (cmPerBar > 500) scaleLength = 500;
    else if (cmPerBar > 200) scaleLength = 200;
    else if (cmPerBar > 100) scaleLength = 100;
    else if (cmPerBar > 50) scaleLength = 50;
    else if (cmPerBar > 20) scaleLength = 20;
    else scaleLength = 10;

    float barPixels = scaleLength * pixelsPerCm;

    // Draw bar
    p.setPen(QPen(QColor(200, 200, 200), 2));
    p.drawLine(QPointF(x, y), QPointF(x + barPixels, y));

    // End caps
    p.drawLine(QPointF(x, y - 5), QPointF(x, y + 5));
    p.drawLine(QPointF(x + barPixels, y - 5), QPointF(x + barPixels, y + 5));

    // Label
    p.setFont(QFont("Sans", 8));
    QString label = QString("%1 cm").arg(scaleLength, 0, 'f', 0);
    if (scaleLength >= 100) label = QString("%1 m").arg(scaleLength / 100.0f, 0, 'f', 1);
    p.drawText(QPointF(x + barPixels/2 - 15, y - 8), label);
}

void EventDisplay::DrawLogo(QPainter& p) {
    QFont font("Sans", 18);
    font.setBold(true);
    p.setFont(font);

    p.setPen(QColor(255, 215, 0));
    p.drawText(width() - 110, 30, "NNBAR");

    font.setPointSize(9);
    font.setBold(false);
    p.setFont(font);
    p.setPen(QColor(140, 140, 140));
    p.drawText(width() - 110, 47, "Event Display");
}

void EventDisplay::DrawSelectedInfo(QPainter& p) {
    if (!m_selectedTrack && !m_selectedVertex) return;

    // Draw info panel for selected object
    int panelW = 200;
    int panelH = 140;
    int x = 15;
    int y = height() - panelH - 50;

    p.setBrush(QColor(20, 20, 35, 230));
    p.setPen(QPen(QColor(255, 215, 0), 2));
    p.drawRoundedRect(x, y, panelW, panelH, 5, 5);

    p.setFont(QFont("Sans", 9, QFont::Bold));
    p.setPen(QColor(255, 215, 0));

    if (m_selectedTrack) {
        p.drawText(x + 10, y + 18, "Selected Track");

        p.setFont(QFont("Sans", 8));
        p.setPen(QColor(200, 200, 200));

        int row = y + 35;
        p.drawText(x + 10, row, QString("Track ID: %1").arg(m_selectedTrack->trackId));
        row += 15;
        p.drawText(x + 10, row, QString("Particle: %1").arg(QString::fromStdString(m_selectedTrack->particleName)));
        row += 15;
        p.drawText(x + 10, row, QString("PDG: %1").arg(m_selectedTrack->pdg));
        row += 15;
        p.drawText(x + 10, row, QString("Energy: %1 MeV").arg(m_selectedTrack->initialEnergy, 0, 'f', 2));
        row += 15;
        p.drawText(x + 10, row, QString("Parent: %1").arg(m_selectedTrack->parentId));
        row += 15;
        p.drawText(x + 10, row, QString("Creator: %1").arg(QString::fromStdString(m_selectedTrack->creatorProcess)));
        row += 15;
        p.drawText(x + 10, row, QString("Steps: %1").arg(m_selectedTrack->points.size()));
    } else if (m_selectedVertex) {
        p.drawText(x + 10, y + 18, "Selected Vertex");

        p.setFont(QFont("Sans", 8));
        p.setPen(QColor(200, 200, 200));

        int row = y + 35;
        p.drawText(x + 10, row, QString("Process: %1").arg(QString::fromStdString(m_selectedVertex->process)));
        row += 15;
        p.drawText(x + 10, row, QString("Position: (%.1f, %.1f, %.1f) cm")
            .arg(m_selectedVertex->x).arg(m_selectedVertex->y).arg(m_selectedVertex->z));
        row += 15;
        p.drawText(x + 10, row, QString("Time: %1 ns").arg(m_selectedVertex->t, 0, 'f', 2));
        row += 15;
        p.drawText(x + 10, row, QString("Multiplicity: %1").arg(m_selectedVertex->Multiplicity()));
        row += 15;
        p.drawText(x + 10, row, QString("Total E: %1 MeV").arg(m_selectedVertex->totalDaughterEnergy, 0, 'f', 2));
        row += 15;
        p.drawText(x + 10, row, QString("Parent Track: %1").arg(m_selectedVertex->parentTrackId));
    }
}

void EventDisplay::DrawTooltip(QPainter& p) {
    if (m_tooltipText.isEmpty()) return;

    p.setFont(QFont("Sans", 8));
    QFontMetrics fm(p.font());
    QRect textRect = fm.boundingRect(m_tooltipText);

    int padding = 5;
    QRect bgRect(m_tooltipPos.x() + 15,
                 m_tooltipPos.y() - textRect.height() - padding * 2,
                 textRect.width() + padding * 2,
                 textRect.height() + padding * 2);

    p.setBrush(QColor(40, 40, 50, 230));
    p.setPen(QPen(QColor(100, 100, 120), 1));
    p.drawRoundedRect(bgRect, 3, 3);

    p.setPen(QColor(220, 220, 220));
    p.drawText(bgRect.adjusted(padding, padding, -padding, -padding), m_tooltipText);
}

// ============================================================================
// Color Utilities
// ============================================================================

QColor EventDisplay::ColorForPDG(int pdg, float alpha) {
    int absPdg = std::abs(pdg);
    QColor color;

    if (absPdg == 11) color = QColor(0, 220, 255);          // Electron - cyan
    else if (absPdg == 22) color = QColor(50, 255, 50);     // Photon - green
    else if (absPdg == 13) color = QColor(255, 100, 100);   // Muon - red
    else if (absPdg == 2212) color = QColor(255, 165, 0);   // Proton - orange
    else if (absPdg == 2112) color = QColor(180, 180, 200); // Neutron - gray
    else if (absPdg == 211 || absPdg == 111) color = QColor(255, 255, 0);  // Pion - yellow
    else if (absPdg == 321 || absPdg == 311) color = QColor(255, 0, 255);  // Kaon - magenta
    else if (absPdg == 3122) color = QColor(128, 180, 255); // Lambda - light blue
    else color = QColor(255, 200, 100);  // Default - gold

    color.setAlphaF(alpha);
    return color;
}

QColor EventDisplay::ColorForEnergy(float energy, float maxEnergy) {
    float t = std::min(1.0f, energy / maxEnergy);

    // Blue -> Cyan -> Green -> Yellow -> Red
    if (t < 0.25f) {
        float s = t / 0.25f;
        return QColor(0, static_cast<int>(255 * s), 255);
    } else if (t < 0.5f) {
        float s = (t - 0.25f) / 0.25f;
        return QColor(0, 255, static_cast<int>(255 * (1 - s)));
    } else if (t < 0.75f) {
        float s = (t - 0.5f) / 0.25f;
        return QColor(static_cast<int>(255 * s), 255, 0);
    } else {
        float s = (t - 0.75f) / 0.25f;
        return QColor(255, static_cast<int>(255 * (1 - s)), 0);
    }
}

QColor EventDisplay::ColorForProcess(const std::string& process) {
    if (process.find("Inelastic") != std::string::npos) return QColor(255, 100, 100);
    if (process.find("Elastic") != std::string::npos) return QColor(100, 255, 100);
    if (process.find("Capture") != std::string::npos) return QColor(255, 255, 100);
    if (process.find("Decay") != std::string::npos) return QColor(255, 100, 255);
    if (process.find("compt") != std::string::npos) return QColor(100, 200, 255);
    if (process.find("phot") != std::string::npos) return QColor(100, 255, 200);
    if (process.find("conv") != std::string::npos) return QColor(255, 200, 100);
    return QColor(200, 200, 200);
}

QString EventDisplay::NameForPDG(int pdg) {
    int absPdg = std::abs(pdg);
    bool antiParticle = (pdg < 0 && absPdg != 22 && absPdg != 111 && absPdg != 2112);

    QString name;
    switch (absPdg) {
        case 11: name = "electron"; break;
        case 12: name = "nu_e"; break;
        case 13: name = "muon"; break;
        case 14: name = "nu_mu"; break;
        case 22: name = "gamma"; break;
        case 111: name = "pi0"; break;
        case 211: name = "pi+"; break;
        case 321: name = "K+"; break;
        case 2112: name = "neutron"; break;
        case 2212: name = "proton"; break;
        case 3122: name = "Lambda"; break;
        default: name = QString("PDG:%1").arg(pdg);
    }

    if (antiParticle && absPdg == 11) name = "positron";
    else if (antiParticle && absPdg == 13) name = "anti-muon";
    else if (antiParticle && absPdg == 211) name = "pi-";
    else if (antiParticle && absPdg == 321) name = "K-";
    else if (antiParticle && absPdg == 2212) name = "anti-proton";

    return name;
}

// ============================================================================
// Mouse Events
// ============================================================================

void EventDisplay::mousePressEvent(QMouseEvent* event) {
    m_lastMousePos = event->pos();

    if (event->button() == Qt::LeftButton) {
        // Check for track/vertex selection
        m_selectedTrack = HitTestTracks(event->pos());
        m_selectedVertex = HitTestVertices(event->pos());

        if (m_selectedTrack) {
            emit trackSelected(m_selectedTrack);
        } else if (m_selectedVertex) {
            emit vertexSelected(m_selectedVertex);
        }

        m_isDragging = true;
    } else if (event->button() == Qt::RightButton) {
        m_isRotating = true;
    }

    update();
}

void EventDisplay::mouseMoveEvent(QMouseEvent* event) {
    QPoint delta = event->pos() - m_lastMousePos;

    if (m_isDragging) {
        m_panX += delta.x();
        m_panY += delta.y();
        update();
    } else if (m_isRotating && m_viewMode == ViewMode::Perspective3D) {
        m_rotationY += delta.x() * 0.5f;
        m_rotationX += delta.y() * 0.5f;
        update();
    } else {
        // Hover tooltip
        Track* hoverTrack = HitTestTracks(event->pos());
        InteractionVertex* hoverVtx = HitTestVertices(event->pos());

        if (hoverTrack) {
            m_tooltipPos = event->pos();
            m_tooltipText = QString("%1 (E=%2 MeV)")
                .arg(QString::fromStdString(hoverTrack->particleName))
                .arg(hoverTrack->initialEnergy, 0, 'f', 1);
            update();
        } else if (hoverVtx) {
            m_tooltipPos = event->pos();
            m_tooltipText = QString("%1 (n=%2)")
                .arg(QString::fromStdString(hoverVtx->process))
                .arg(hoverVtx->Multiplicity());
            update();
        } else {
            if (!m_tooltipText.isEmpty()) {
                m_tooltipText.clear();
                update();
            }
        }
    }

    m_lastMousePos = event->pos();
}

void EventDisplay::mouseReleaseEvent(QMouseEvent* event) {
    Q_UNUSED(event);
    m_isDragging = false;
    m_isRotating = false;
}

void EventDisplay::wheelEvent(QWheelEvent* event) {
    if (event->angleDelta().y() > 0) {
        m_zoom *= 1.15f;
    } else {
        m_zoom /= 1.15f;
        if (m_zoom < 0.1f) m_zoom = 0.1f;
    }
    update();
}

void EventDisplay::resizeEvent(QResizeEvent* event) {
    QWidget::resizeEvent(event);
    update();
}

// ============================================================================
// Hit Testing
// ============================================================================

Track* EventDisplay::HitTestTracks(const QPoint& pos) {
    std::lock_guard<std::mutex> lock(m_mutex);

    const float hitRadius = 10.0f;

    for (auto& track : m_eventData.tracks) {
        if (track.points.size() < 2) continue;

        for (const auto& pt : track.points) {
            QPointF screenPt;
            if (m_viewMode == ViewMode::Transverse) {
                screenPt = ToScreen(pt.x, pt.y);
            } else if (m_viewMode == ViewMode::Longitudinal) {
                float r = std::sqrt(pt.x * pt.x + pt.y * pt.y);
                screenPt = ToScreen(pt.z, r * (pt.y >= 0 ? 1 : -1));
            } else {
                screenPt = ToScreen3D(pt.x, pt.y, pt.z);
            }

            float dx = pos.x() - screenPt.x();
            float dy = pos.y() - screenPt.y();
            if (dx*dx + dy*dy < hitRadius*hitRadius) {
                return &track;
            }
        }
    }
    return nullptr;
}

InteractionVertex* EventDisplay::HitTestVertices(const QPoint& pos) {
    std::lock_guard<std::mutex> lock(m_mutex);

    const float hitRadius = 15.0f;

    for (auto& vtx : m_eventData.vertices) {
        QPointF screenPt;
        if (m_viewMode == ViewMode::Transverse) {
            screenPt = ToScreen(vtx.x, vtx.y);
        } else if (m_viewMode == ViewMode::Longitudinal) {
            float r = std::sqrt(vtx.x * vtx.x + vtx.y * vtx.y);
            screenPt = ToScreen(vtx.z, r * (vtx.y >= 0 ? 1 : -1));
        } else {
            screenPt = ToScreen3D(vtx.x, vtx.y, vtx.z);
        }

        float dx = pos.x() - screenPt.x();
        float dy = pos.y() - screenPt.y();
        if (dx*dx + dy*dy < hitRadius*hitRadius) {
            return &vtx;
        }
    }
    return nullptr;
}

// ============================================================================
// View Mode Control
// ============================================================================

void EventDisplay::SetViewMode(ViewMode mode) {
    m_viewMode = mode;
    update();
    emit viewModeChanged(mode);
}

void EventDisplay::ResetView() {
    m_zoom = 1.0f;
    m_panX = 0.0f;
    m_panY = 0.0f;
    m_rotationX = 30.0f;
    m_rotationY = -45.0f;
    m_selectedTrack = nullptr;
    m_selectedVertex = nullptr;
    update();
}

void EventDisplay::ZoomIn() {
    m_zoom *= 1.2f;
    update();
}

void EventDisplay::ZoomOut() {
    m_zoom /= 1.2f;
    if (m_zoom < 0.1f) m_zoom = 0.1f;
    update();
}

// ============================================================================
// Data Methods
// ============================================================================

void EventDisplay::AddTrack(const Track& track) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_eventData.tracks.push_back(track);
    m_pendingUpdates++;
}

void EventDisplay::AddCaloHit(const CaloHit& hit) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_eventData.caloHits.push_back(hit);
    m_pendingUpdates++;
}

void EventDisplay::AddInteractionVertex(const InteractionVertex& vertex) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_eventData.vertices.push_back(vertex);
    m_pendingUpdates++;
}

void EventDisplay::SetEventData(const EventDisplayData& evtData) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_eventData = evtData;
    m_pendingUpdates++;
    emit eventDataUpdated();
}

void EventDisplay::ClearEvent() {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_eventData.tracks.clear();
    m_eventData.caloHits.clear();
    m_eventData.vertices.clear();
    m_eventData.opticalPhotonRegions.clear();
    m_eventData.energyDeposits.clear();
    m_eventData.totalOpticalPhotons = 0;
    m_eventData.eventId = 0;
    m_eventData.totalEnergy = 0.0;

    // Particle counts
    m_eventData.nPrimaries = 0;
    m_eventData.nSecondaries = 0;
    m_eventData.nGammas = 0;
    m_eventData.nElectrons = 0;
    m_eventData.nNeutrons = 0;
    m_eventData.nProtons = 0;
    m_eventData.nPions = 0;
    m_eventData.nMuons = 0;
    m_eventData.nOther = 0;

    // Energy per detector
    m_eventData.edepBeampipe = 0.0;
    m_eventData.edepTPC = 0.0;
    m_eventData.edepScintillator = 0.0;
    m_eventData.edepLeadGlass = 0.0;
    m_eventData.edepShield = 0.0;
    m_eventData.edepOther = 0.0;

    // Optical photon counts
    m_eventData.photonsScintillator = 0;
    m_eventData.photonsLeadGlass = 0;
    m_eventData.photonsTPC = 0;

    m_selectedTrack = nullptr;
    m_selectedVertex = nullptr;
    update();
}

void EventDisplay::AddOpticalPhotonRegion(const OpticalPhotonRegion& region) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_eventData.opticalPhotonRegions.push_back(region);
    m_pendingUpdates++;
}

void EventDisplay::SetTotalOpticalPhotons(int count) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_eventData.totalOpticalPhotons = count;
    m_pendingUpdates++;
}

void EventDisplay::AddEnergyDeposit(const EnergyDeposit& deposit) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_eventData.energyDeposits.push_back(deposit);

    // Update energy totals per detector
    switch (deposit.detectorType) {
        case 0:  // Beampipe
            m_eventData.edepBeampipe += deposit.energy;
            break;
        case 1:  // TPC
            m_eventData.edepTPC += deposit.energy;
            break;
        case 2:  // Scintillator
            m_eventData.edepScintillator += deposit.energy;
            break;
        case 3:  // LeadGlass
            m_eventData.edepLeadGlass += deposit.energy;
            break;
        case 4:  // Shield
            m_eventData.edepShield += deposit.energy;
            break;
        default:
            m_eventData.edepOther += deposit.energy;
            break;
    }

    m_eventData.totalEnergy += deposit.energy;
    m_pendingUpdates++;
}

} // namespace nnbar

#endif // WITH_DASHBOARD
