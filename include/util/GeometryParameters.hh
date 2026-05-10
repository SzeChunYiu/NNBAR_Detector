// ============================================================================
// GeometryParameters.hh
// Singleton to share geometry parameters between DetectorConstruction and GUI
// Enables dynamic geometry updates in event display
// ============================================================================

#ifndef GEOMETRY_PARAMETERS_HH
#define GEOMETRY_PARAMETERS_HH

#include <mutex>
#include <map>
#include <string>

namespace nnbar {

class GeometryParameters {
public:
    // Singleton access
    static GeometryParameters& Instance() {
        static GeometryParameters instance;
        return instance;
    }

    // Prevent copying
    GeometryParameters(const GeometryParameters&) = delete;
    GeometryParameters& operator=(const GeometryParameters&) = delete;

    // Thread-safe setters (called from DetectorConstruction)
    void Set(const std::string& name, double value) {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_params[name] = value;
    }

    // Thread-safe getters (called from EventDisplay)
    double Get(const std::string& name, double defaultValue = 0.0) const {
        std::lock_guard<std::mutex> lock(m_mutex);
        auto it = m_params.find(name);
        return (it != m_params.end()) ? it->second : defaultValue;
    }

    bool Has(const std::string& name) const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return m_params.find(name) != m_params.end();
    }

    // Convenience accessors for common geometry values (in cm for display)
    // Beampipe
    double BeampipeInnerRadius() const { return Get("beampipe.inner_radius", 112.0); }
    double BeampipeOuterRadius() const { return Get("beampipe.outer_radius", 114.0); }
    double BeampipeThickness() const { return Get("beampipe.thickness", 2.0); }
    double BeampipeHalfZ() const { return Get("beampipe.half_z", 250.0); }

    // TPC Type I (vertical)
    double TPCTypeIWidth() const { return Get("tpc.type1.width", 85.4); }
    double TPCTypeIHeight() const { return Get("tpc.type1.height", 199.4); }
    double TPCTypeIX() const { return Get("tpc.type1.x", 156.9); }
    double TPCTypeIY() const { return Get("tpc.type1.y", 99.7); }

    // TPC Type II (horizontal)
    double TPCTypeIIWidth() const { return Get("tpc.type2.width", 228.4); }
    double TPCTypeIIHeight() const { return Get("tpc.type2.height", 85.4); }
    double TPCTypeIIY() const { return Get("tpc.type2.y", 156.7); }

    // TPC common
    double TPCDriftLength() const { return Get("tpc.drift_length", 85.0); }
    double TPCHalfZ() const { return Get("tpc.half_z", 126.0); }

    // Scintillator
    double ScintillatorDistance() const { return Get("scint.distance", 275.0); }
    double ScintillatorHalfZ() const { return Get("scint.half_z", 275.0); }

    // Lead Glass Calorimeter
    double CaloSideDistance() const { return Get("calo.side_distance", 275.0); }
    double CaloFBDistance() const { return Get("calo.fb_distance", 335.0); }
    double CaloHalfY() const { return Get("calo.half_y", 300.0); }

    // Cosmic Shielding
    double ShieldHalfX() const { return Get("shield.half_x", 550.0); }
    double ShieldHalfY() const { return Get("shield.half_y", 275.0); }
    double ShieldHalfZ() const { return Get("shield.half_z", 335.0); }

    // Check if geometry has been registered
    bool IsInitialized() const {
        std::lock_guard<std::mutex> lock(m_mutex);
        return !m_params.empty();
    }

    // Clear all parameters
    void Clear() {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_params.clear();
    }

private:
    GeometryParameters() = default;
    ~GeometryParameters() = default;

    mutable std::mutex m_mutex;
    std::map<std::string, double> m_params;
};

} // namespace nnbar

#endif // GEOMETRY_PARAMETERS_HH
