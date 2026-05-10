// ============================================================================
// PhoenixDisplay.hh
// Phoenix event display integration for NNBAR detector simulation
// Embeds the CERN Phoenix framework via QWebEngineView
// Requires Qt5WebEngineWidgets - disabled if not available
// ============================================================================

#ifndef PHOENIX_DISPLAY_HH
#define PHOENIX_DISPLAY_HH

#include "config.h"

// Only compile Phoenix display if both dashboard and WebEngine are available
#if WITH_DASHBOARD && WITH_PHOENIX

#include <QWidget>
#include <QWebEngineView>
#include <QWebChannel>
#include <QVBoxLayout>
#include <QJsonObject>
#include <QJsonArray>
#include <QJsonDocument>
#include <vector>
#include <mutex>
#include <cstdint>

namespace nnbar {

// Track data for Phoenix visualization
struct PhoenixTrack {
    std::vector<std::array<double, 3>> positions;  // x, y, z points
    int32_t pdg;
    double energy;
    std::string color;
};

// Calorimeter cluster for Phoenix
struct PhoenixCluster {
    double x, y, z;
    double energy;
    double eta, phi;
    std::string type;  // "ecal", "hcal", etc.
};

// Bridge object for JavaScript communication
class PhoenixBridge : public QObject {
    Q_OBJECT

public:
    explicit PhoenixBridge(QObject* parent = nullptr) : QObject(parent) {}

public slots:
    void onPhoenixReady() { emit phoenixReady(); }
    void onEventLoaded(const QString& info) { emit eventLoaded(info); }

signals:
    void phoenixReady();
    void eventLoaded(const QString& info);
    void loadEventData(const QString& jsonData);
    void clearEvent();
    void setViewAngle(double theta, double phi);
};

class PhoenixDisplay : public QWidget {
    Q_OBJECT

public:
    explicit PhoenixDisplay(QWidget* parent = nullptr);
    ~PhoenixDisplay();

    // Data update methods (thread-safe)
    void AddTrack(const PhoenixTrack& track);
    void AddCluster(const PhoenixCluster& cluster);
    void SetEventInfo(int32_t runId, int32_t eventId, const QString& timestamp);
    void ClearEvent();
    void LoadEvent();

    // View controls
    void SetViewAngle(double theta, double phi);
    void ResetView();

signals:
    void eventDataReady(const QString& json);

private slots:
    void onPhoenixReady();
    void onLoadFinished(bool ok);

private:
    void InitializePhoenix();
    QString GeneratePhoenixHTML();
    QString GenerateEventJSON();

    QWebEngineView* m_webView = nullptr;
    QWebChannel* m_channel = nullptr;
    PhoenixBridge* m_bridge = nullptr;

    std::mutex m_mutex;
    std::vector<PhoenixTrack> m_tracks;
    std::vector<PhoenixCluster> m_clusters;
    int32_t m_runId = 0;
    int32_t m_eventId = 0;
    QString m_timestamp;

    bool m_phoenixReady = false;
};

} // namespace nnbar

#endif // WITH_DASHBOARD && WITH_PHOENIX
#endif // PHOENIX_DISPLAY_HH
