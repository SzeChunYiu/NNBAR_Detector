// ============================================================================
// PhoenixDisplay.cc
// Phoenix event display integration implementation
// ============================================================================

#include "config.h"
#include "gui/PhoenixDisplay.hh"

#if WITH_DASHBOARD && WITH_PHOENIX

#include <QWebEngineSettings>
#include <QUrl>
#include <QDebug>

namespace nnbar {

// ============================================================================
// Constructor
// ============================================================================

PhoenixDisplay::PhoenixDisplay(QWidget* parent) : QWidget(parent) {
    QVBoxLayout* layout = new QVBoxLayout(this);
    layout->setContentsMargins(0, 0, 0, 0);

    m_webView = new QWebEngineView(this);
    m_webView->setMinimumSize(600, 400);

    // Enable WebGL and local storage
    QWebEngineSettings* settings = m_webView->settings();
    settings->setAttribute(QWebEngineSettings::WebGLEnabled, true);
    settings->setAttribute(QWebEngineSettings::LocalStorageEnabled, true);
    settings->setAttribute(QWebEngineSettings::JavascriptEnabled, true);
    settings->setAttribute(QWebEngineSettings::LocalContentCanAccessRemoteUrls, true);

    layout->addWidget(m_webView);

    // Set up web channel for JavaScript communication
    m_channel = new QWebChannel(this);
    m_bridge = new PhoenixBridge(this);
    m_channel->registerObject("phoenixBridge", m_bridge);
    m_webView->page()->setWebChannel(m_channel);

    connect(m_bridge, &PhoenixBridge::phoenixReady, this, &PhoenixDisplay::onPhoenixReady);
    connect(m_webView, &QWebEngineView::loadFinished, this, &PhoenixDisplay::onLoadFinished);

    InitializePhoenix();
}

PhoenixDisplay::~PhoenixDisplay() = default;

// ============================================================================
// Initialize Phoenix
// ============================================================================

void PhoenixDisplay::InitializePhoenix() {
    QString html = GeneratePhoenixHTML();
    m_webView->setHtml(html, QUrl("qrc:/phoenix/"));
}

// ============================================================================
// Generate Phoenix HTML with embedded viewer
// ============================================================================

QString PhoenixDisplay::GeneratePhoenixHTML() {
    return R"HTML(
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>NNBAR Event Display</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/three@0.128.0/examples/js/controls/OrbitControls.js"></script>
    <script src="qrc:///qtwebchannel/qwebchannel.js"></script>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            background: #0a0a12;
            overflow: hidden;
            font-family: 'Segoe UI', Arial, sans-serif;
        }
        #canvas-container {
            width: 100vw;
            height: 100vh;
            position: relative;
        }
        #info-overlay {
            position: absolute;
            top: 10px;
            left: 10px;
            color: #00d4ff;
            font-size: 12px;
            background: rgba(0,0,0,0.7);
            padding: 10px 15px;
            border-radius: 5px;
            border: 1px solid #333;
            z-index: 100;
        }
        #info-overlay h2 {
            color: #ffd700;
            font-size: 16px;
            margin-bottom: 8px;
        }
        #info-overlay .label { color: #888; }
        #info-overlay .value { color: #00ff88; font-family: monospace; }
        #legend {
            position: absolute;
            bottom: 10px;
            left: 10px;
            color: white;
            font-size: 11px;
            background: rgba(0,0,0,0.7);
            padding: 10px;
            border-radius: 5px;
            border: 1px solid #333;
        }
        .legend-item { display: flex; align-items: center; margin: 3px 0; }
        .legend-color { width: 20px; height: 3px; margin-right: 8px; border-radius: 1px; }
        #controls {
            position: absolute;
            top: 10px;
            right: 10px;
            display: flex;
            gap: 5px;
        }
        .ctrl-btn {
            background: rgba(0,120,212,0.8);
            color: white;
            border: none;
            padding: 8px 12px;
            border-radius: 4px;
            cursor: pointer;
            font-size: 11px;
        }
        .ctrl-btn:hover { background: rgba(0,140,232,0.9); }
    </style>
</head>
<body>
    <div id="canvas-container"></div>

    <div id="info-overlay">
        <h2>NNBAR Event Display</h2>
        <div><span class="label">Run:</span> <span class="value" id="run-id">-</span></div>
        <div><span class="label">Event:</span> <span class="value" id="event-id">-</span></div>
        <div><span class="label">Tracks:</span> <span class="value" id="track-count">0</span></div>
        <div><span class="label">Clusters:</span> <span class="value" id="cluster-count">0</span></div>
    </div>

    <div id="legend">
        <div class="legend-item"><div class="legend-color" style="background:#00ffff"></div>Electrons</div>
        <div class="legend-item"><div class="legend-color" style="background:#00ff00"></div>Photons</div>
        <div class="legend-item"><div class="legend-color" style="background:#ff6666"></div>Muons</div>
        <div class="legend-item"><div class="legend-color" style="background:#ffcc00"></div>Hadrons</div>
        <div class="legend-item"><div class="legend-color" style="background:#ff00ff;height:8px;width:8px;border-radius:50%"></div>Clusters</div>
    </div>

    <div id="controls">
        <button class="ctrl-btn" onclick="resetView()">Reset View</button>
        <button class="ctrl-btn" onclick="viewXY()">X-Y</button>
        <button class="ctrl-btn" onclick="viewRZ()">R-Z</button>
        <button class="ctrl-btn" onclick="view3D()">3D</button>
    </div>

    <script>
        // Three.js setup
        let scene, camera, renderer, controls;
        let detectorGroup, tracksGroup, clustersGroup;

        const SCALE = 0.01;  // cm to display units

        function init() {
            // Scene
            scene = new THREE.Scene();
            scene.background = new THREE.Color(0x0a0a12);

            // Camera
            camera = new THREE.PerspectiveCamera(45, window.innerWidth / window.innerHeight, 0.1, 1000);
            camera.position.set(8, 6, 8);

            // Renderer
            renderer = new THREE.WebGLRenderer({ antialias: true });
            renderer.setSize(window.innerWidth, window.innerHeight);
            renderer.setPixelRatio(window.devicePixelRatio);
            document.getElementById('canvas-container').appendChild(renderer.domElement);

            // Controls
            controls = new THREE.OrbitControls(camera, renderer.domElement);
            controls.enableDamping = true;
            controls.dampingFactor = 0.05;

            // Groups
            detectorGroup = new THREE.Group();
            tracksGroup = new THREE.Group();
            clustersGroup = new THREE.Group();
            scene.add(detectorGroup);
            scene.add(tracksGroup);
            scene.add(clustersGroup);

            // Create detector geometry
            createDetector();

            // Lighting
            const ambientLight = new THREE.AmbientLight(0x404040, 0.5);
            scene.add(ambientLight);
            const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
            directionalLight.position.set(5, 10, 5);
            scene.add(directionalLight);

            // Grid
            const gridHelper = new THREE.GridHelper(10, 20, 0x333366, 0x222244);
            scene.add(gridHelper);

            // Axes
            const axesHelper = new THREE.AxesHelper(2);
            scene.add(axesHelper);

            window.addEventListener('resize', onWindowResize);
            animate();

            // Setup Qt WebChannel
            if (typeof QWebChannel !== 'undefined') {
                new QWebChannel(qt.webChannelTransport, function(channel) {
                    window.phoenixBridge = channel.objects.phoenixBridge;
                    phoenixBridge.phoenixReady();

                    phoenixBridge.loadEventData.connect(function(jsonData) {
                        loadEventFromJSON(jsonData);
                    });

                    phoenixBridge.clearEvent.connect(function() {
                        clearTracks();
                        clearClusters();
                    });
                });
            }
        }

        function createDetector() {
            // TPC - inner cylinder
            const tpcInnerGeom = new THREE.CylinderGeometry(0.5, 0.5, 6, 32, 1, true);
            const tpcMat = new THREE.MeshBasicMaterial({
                color: 0x00aaaa,
                transparent: true,
                opacity: 0.2,
                wireframe: true
            });
            const tpcInner = new THREE.Mesh(tpcInnerGeom, tpcMat);
            tpcInner.rotation.x = Math.PI / 2;
            detectorGroup.add(tpcInner);

            // TPC - outer cylinder
            const tpcOuterGeom = new THREE.CylinderGeometry(2, 2, 6, 32, 1, true);
            const tpcOuter = new THREE.Mesh(tpcOuterGeom, tpcMat);
            tpcOuter.rotation.x = Math.PI / 2;
            detectorGroup.add(tpcOuter);

            // Scintillator layer
            const scintGeom = new THREE.CylinderGeometry(2.5, 2.5, 6, 32, 1, true);
            const scintMat = new THREE.MeshBasicMaterial({
                color: 0x22aa22,
                transparent: true,
                opacity: 0.15,
                wireframe: true
            });
            const scint = new THREE.Mesh(scintGeom, scintMat);
            scint.rotation.x = Math.PI / 2;
            detectorGroup.add(scint);

            // Calorimeter (Lead Glass)
            const caloGeom = new THREE.CylinderGeometry(3.5, 3.5, 8, 32, 1, true);
            const caloMat = new THREE.MeshBasicMaterial({
                color: 0xaaaa00,
                transparent: true,
                opacity: 0.1,
                wireframe: true
            });
            const calo = new THREE.Mesh(caloGeom, caloMat);
            calo.rotation.x = Math.PI / 2;
            detectorGroup.add(calo);

            // Endcaps
            const endcapGeom = new THREE.RingGeometry(0.5, 3.5, 32);
            const endcapMat = new THREE.MeshBasicMaterial({
                color: 0xaa8800,
                transparent: true,
                opacity: 0.15,
                side: THREE.DoubleSide
            });
            const endcap1 = new THREE.Mesh(endcapGeom, endcapMat);
            endcap1.position.z = 4;
            detectorGroup.add(endcap1);
            const endcap2 = new THREE.Mesh(endcapGeom, endcapMat);
            endcap2.position.z = -4;
            detectorGroup.add(endcap2);
        }

        function addTrack(points, color) {
            if (points.length < 2) return;

            const geometry = new THREE.BufferGeometry();
            const positions = new Float32Array(points.length * 3);

            for (let i = 0; i < points.length; i++) {
                positions[i * 3] = points[i][0] * SCALE;
                positions[i * 3 + 1] = points[i][1] * SCALE;
                positions[i * 3 + 2] = points[i][2] * SCALE;
            }

            geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));

            const material = new THREE.LineBasicMaterial({
                color: color,
                linewidth: 2
            });

            const line = new THREE.Line(geometry, material);
            tracksGroup.add(line);
        }

        function addCluster(x, y, z, energy, color) {
            const size = 0.05 + Math.min(0.3, energy * 0.001);
            const geometry = new THREE.SphereGeometry(size, 16, 16);
            const material = new THREE.MeshPhongMaterial({
                color: color,
                transparent: true,
                opacity: 0.8,
                emissive: color,
                emissiveIntensity: 0.3
            });
            const sphere = new THREE.Mesh(geometry, material);
            sphere.position.set(x * SCALE, y * SCALE, z * SCALE);
            clustersGroup.add(sphere);
        }

        function clearTracks() {
            while(tracksGroup.children.length > 0) {
                tracksGroup.remove(tracksGroup.children[0]);
            }
        }

        function clearClusters() {
            while(clustersGroup.children.length > 0) {
                clustersGroup.remove(clustersGroup.children[0]);
            }
        }

        function loadEventFromJSON(jsonString) {
            try {
                const data = JSON.parse(jsonString);

                clearTracks();
                clearClusters();

                // Update info
                document.getElementById('run-id').textContent = data.runId || '-';
                document.getElementById('event-id').textContent = data.eventId || '-';

                // Add tracks
                if (data.tracks) {
                    data.tracks.forEach(track => {
                        const color = getColorForPDG(track.pdg);
                        addTrack(track.positions, color);
                    });
                    document.getElementById('track-count').textContent = data.tracks.length;
                }

                // Add clusters
                if (data.clusters) {
                    data.clusters.forEach(cluster => {
                        addCluster(cluster.x, cluster.y, cluster.z, cluster.energy, 0xff00ff);
                    });
                    document.getElementById('cluster-count').textContent = data.clusters.length;
                }

            } catch (e) {
                console.error('Error loading event:', e);
            }
        }

        function getColorForPDG(pdg) {
            const absPdg = Math.abs(pdg);
            if (absPdg === 11) return 0x00ffff;  // electron - cyan
            if (absPdg === 22) return 0x00ff00;  // photon - green
            if (absPdg === 13) return 0xff6666;  // muon - red
            if (absPdg === 2212) return 0xff9900; // proton - orange
            if (absPdg === 2112) return 0xcccccc; // neutron - gray
            if (absPdg === 211 || absPdg === 111) return 0xffff00; // pion - yellow
            return 0xffcc00; // default hadron - gold
        }

        function resetView() {
            camera.position.set(8, 6, 8);
            camera.lookAt(0, 0, 0);
            controls.reset();
        }

        function viewXY() {
            camera.position.set(0, 0, 12);
            camera.lookAt(0, 0, 0);
        }

        function viewRZ() {
            camera.position.set(12, 0, 0);
            camera.lookAt(0, 0, 0);
        }

        function view3D() {
            camera.position.set(8, 6, 8);
            camera.lookAt(0, 0, 0);
        }

        function onWindowResize() {
            camera.aspect = window.innerWidth / window.innerHeight;
            camera.updateProjectionMatrix();
            renderer.setSize(window.innerWidth, window.innerHeight);
        }

        function animate() {
            requestAnimationFrame(animate);
            controls.update();
            renderer.render(scene, camera);
        }

        // Start
        init();
    </script>
</body>
</html>
)HTML";
}

// ============================================================================
// Event Handlers
// ============================================================================

void PhoenixDisplay::onLoadFinished(bool ok) {
    if (ok) {
        qDebug() << "Phoenix HTML loaded successfully";
    } else {
        qDebug() << "Failed to load Phoenix HTML";
    }
}

void PhoenixDisplay::onPhoenixReady() {
    m_phoenixReady = true;
    qDebug() << "Phoenix is ready for events";
}

// ============================================================================
// Data Methods
// ============================================================================

void PhoenixDisplay::AddTrack(const PhoenixTrack& track) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_tracks.push_back(track);
}

void PhoenixDisplay::AddCluster(const PhoenixCluster& cluster) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_clusters.push_back(cluster);
}

void PhoenixDisplay::SetEventInfo(int32_t runId, int32_t eventId, const QString& timestamp) {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_runId = runId;
    m_eventId = eventId;
    m_timestamp = timestamp;
}

void PhoenixDisplay::ClearEvent() {
    std::lock_guard<std::mutex> lock(m_mutex);
    m_tracks.clear();
    m_clusters.clear();

    if (m_phoenixReady && m_bridge) {
        emit m_bridge->clearEvent();
    }
}

void PhoenixDisplay::LoadEvent() {
    if (!m_phoenixReady) return;

    QString json = GenerateEventJSON();
    if (m_bridge) {
        emit m_bridge->loadEventData(json);
    }
}

QString PhoenixDisplay::GenerateEventJSON() {
    std::lock_guard<std::mutex> lock(m_mutex);

    QJsonObject root;
    root["runId"] = m_runId;
    root["eventId"] = m_eventId;
    root["timestamp"] = m_timestamp;

    // Tracks
    QJsonArray tracksArray;
    for (const auto& track : m_tracks) {
        QJsonObject trackObj;
        trackObj["pdg"] = track.pdg;
        trackObj["energy"] = track.energy;

        QJsonArray posArray;
        for (const auto& pos : track.positions) {
            QJsonArray point;
            point.append(pos[0]);
            point.append(pos[1]);
            point.append(pos[2]);
            posArray.append(point);
        }
        trackObj["positions"] = posArray;
        tracksArray.append(trackObj);
    }
    root["tracks"] = tracksArray;

    // Clusters
    QJsonArray clustersArray;
    for (const auto& cluster : m_clusters) {
        QJsonObject clusterObj;
        clusterObj["x"] = cluster.x;
        clusterObj["y"] = cluster.y;
        clusterObj["z"] = cluster.z;
        clusterObj["energy"] = cluster.energy;
        clusterObj["type"] = QString::fromStdString(cluster.type);
        clustersArray.append(clusterObj);
    }
    root["clusters"] = clustersArray;

    return QString(QJsonDocument(root).toJson(QJsonDocument::Compact));
}

// ============================================================================
// View Controls
// ============================================================================

void PhoenixDisplay::SetViewAngle(double theta, double phi) {
    if (m_phoenixReady && m_bridge) {
        emit m_bridge->setViewAngle(theta, phi);
    }
}

void PhoenixDisplay::ResetView() {
    if (m_webView) {
        m_webView->page()->runJavaScript("resetView();");
    }
}

} // namespace nnbar

#endif // WITH_DASHBOARD && WITH_PHOENIX
