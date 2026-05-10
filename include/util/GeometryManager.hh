// ============================================================================
// GeometryManager.hh
// Singleton manager for detector geometry information
// Maps replica/copy numbers to positions, orientations, and dimensions
// Enables lookup of specific volumes for visualization and analysis
// ============================================================================

#ifndef GEOMETRY_MANAGER_HH
#define GEOMETRY_MANAGER_HH

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include <vector>
#include <map>
#include <string>
#include <cstdint>

namespace nnbar {

// ============================================================================
// Data structures for volume information
// ============================================================================

/**
 * @struct VolumeInfo
 * @brief Complete information about a detector volume
 */
struct VolumeInfo {
    int32_t copyNumber;         // Replica/copy number from G4PVPlacement
    std::string name;           // Volume name (e.g., "LeadGlassPV", "ScintPV")
    std::string type;           // Volume type (e.g., "scint_side", "lg_front")

    // Position in global coordinates (mm)
    double x, y, z;

    // Orientation (rotation angles in degrees)
    double rotX, rotY, rotZ;

    // Dimensions (mm)
    double sizeX, sizeY, sizeZ;

    // Hierarchical info
    int32_t surfaceIndex;       // Which surface (0=top, 1=right, 2=bottom, 3=left, 4=front, 5=back)
    int32_t layerIndex;         // Layer within module (for scintillator)
    int32_t staveIndex;         // Stave/bar index within layer

    // For visualization
    double energyDeposit;       // Current energy deposit (MeV)
    int32_t nOpticalPhotons;    // Number of optical photons detected
};

/**
 * @struct ScintillatorInfo
 * @brief Specific information for scintillator bars
 */
struct ScintillatorInfo {
    int32_t copyNumber;         // Module copy number
    int32_t barCopyNumber;      // Bar copy number within module
    std::string barType;        // "H" (horizontal) or "V" (vertical)

    // Module position
    double moduleX, moduleY, moduleZ;

    // Bar position within module (local coordinates)
    double localX, localY, localZ;

    // Global position of bar center
    double globalX, globalY, globalZ;

    // Surface and orientation
    int32_t surfaceIndex;       // 0-5 for top/right/bottom/left/front/back
    double rotationAngle;       // Rotation around Z (for side surfaces)

    // Bar dimensions
    double barX, barY, barZ;

    // Physics data (updated during simulation)
    double energyDeposit;
    int32_t nPhotons;
};

/**
 * @struct LeadGlassInfo
 * @brief Specific information for lead glass blocks
 */
struct LeadGlassInfo {
    int32_t copyNumber;         // Copy number from placement

    // Position from CSV data
    double csvX, csvY, csvZ;    // Original CSV coordinates (cm)
    double csvRotX, csvRotZ;    // Original CSV rotation angles

    // Global position (mm)
    double globalX, globalY, globalZ;

    // Orientation
    int32_t surfaceIndex;       // 0-3 for sides, 4=front, 5=back
    double rotationAngle;       // Surface rotation (0, 90, 180, 270 deg)

    // Block dimensions
    static constexpr double blockX = 80.0;   // mm
    static constexpr double blockY = 250.0;  // mm (PMT direction)
    static constexpr double blockZ = 80.0;   // mm

    // Physics data
    double energyDeposit;
    int32_t nPhotons;
    int32_t nPMTHits;
};

/**
 * @struct TPCModuleInfo
 * @brief Information for TPC modules
 */
struct TPCModuleInfo {
    int32_t moduleIndex;        // 0-11
    int32_t moduleType;         // 1 = Type I (vertical), 2 = Type II (horizontal)

    // Position (mm)
    double centerX, centerY, centerZ;

    // Dimensions (mm)
    double sizeX, sizeY, sizeZ;

    // Drift direction
    int32_t driftAxis;          // 0=X, 1=Y
    int32_t driftSign;          // +1 or -1

    // Location
    std::string location;       // "front" or "back"
    std::string position;       // "top", "bottom", "top-left", etc.
};

/**
 * @struct BeampipeInfo
 * @brief Information for beampipe sections
 */
struct BeampipeInfo {
    int32_t sectionIndex;       // 1-8 for main sections
    std::string name;           // "Beampipe_1", "Beampipe_5", etc.
    std::string type;           // "wall", "coating", "cap", "beamstop"

    // Position (mm)
    double posZ;                // Z position of section center

    // Dimensions (mm)
    double innerRadius1;        // Inner radius at -z end
    double outerRadius1;        // Outer radius at -z end
    double innerRadius2;        // Inner radius at +z end
    double outerRadius2;        // Outer radius at +z end
    double length;              // Half-length

    // Material
    std::string material;       // "Aluminum", "B4C", "Copper"

    // Physics data
    double energyDeposit;
    int32_t nNeutronCaptures;
};

// ============================================================================
// GeometryManager Singleton
// ============================================================================

class GeometryManager {
public:
    static GeometryManager& Instance();

    // Initialization - call after geometry construction
    void Initialize();
    void Clear();

    // ========================================================================
    // Registration methods (called during geometry construction)
    // ========================================================================

    /**
     * Register a scintillator module
     */
    void RegisterScintillatorModule(int32_t copyNumber,
                                    double x, double y, double z,
                                    int32_t surfaceIndex,
                                    double rotationAngle);

    /**
     * Register a scintillator bar within a module
     */
    void RegisterScintillatorBar(int32_t moduleCopyNumber,
                                 int32_t barCopyNumber,
                                 const std::string& barType,
                                 double localX, double localY, double localZ,
                                 double barSizeX, double barSizeY, double barSizeZ);

    /**
     * Register a lead glass block
     */
    void RegisterLeadGlassBlock(int32_t copyNumber,
                                double x, double y, double z,
                                double rotX, double rotZ,
                                int32_t surfaceIndex);

    /**
     * Register a TPC module
     */
    void RegisterTPCModule(int32_t moduleIndex,
                           int32_t moduleType,
                           double x, double y, double z,
                           double sizeX, double sizeY, double sizeZ,
                           int32_t driftAxis, int32_t driftSign);

    /**
     * Register a beampipe section
     */
    void RegisterBeampipeSection(int32_t sectionIndex,
                                 const std::string& name,
                                 const std::string& type,
                                 double posZ,
                                 double innerR1, double outerR1,
                                 double innerR2, double outerR2,
                                 double length,
                                 const std::string& material);

    // ========================================================================
    // Lookup methods
    // ========================================================================

    /**
     * Get scintillator bar info by copy numbers
     */
    const ScintillatorInfo* GetScintillatorBar(int32_t moduleCopyNumber,
                                               int32_t barCopyNumber) const;

    /**
     * Get all scintillator bars on a surface
     */
    std::vector<const ScintillatorInfo*> GetScintillatorBarsOnSurface(int32_t surfaceIndex) const;

    /**
     * Get lead glass block by copy number
     */
    const LeadGlassInfo* GetLeadGlassBlock(int32_t copyNumber) const;

    /**
     * Get all lead glass blocks on a surface
     */
    std::vector<const LeadGlassInfo*> GetLeadGlassBlocksOnSurface(int32_t surfaceIndex) const;

    /**
     * Get TPC module by index
     */
    const TPCModuleInfo* GetTPCModule(int32_t moduleIndex) const;

    /**
     * Find which TPC module contains a point
     */
    int32_t FindTPCModule(double x, double y, double z) const;

    /**
     * Get beampipe section by index
     */
    const BeampipeInfo* GetBeampipeSection(int32_t sectionIndex) const;

    /**
     * Get all beampipe sections
     */
    const std::vector<BeampipeInfo>& GetAllBeampipeSections() const { return m_beampipeSections; }

    /**
     * Find which lead glass block is at a position (within tolerance)
     */
    int32_t FindLeadGlassBlock(double x, double y, double z, double tolerance = 50.0) const;

    /**
     * Find which scintillator bar is at a position
     */
    std::pair<int32_t, int32_t> FindScintillatorBar(double x, double y, double z,
                                                     double tolerance = 20.0) const;

    // ========================================================================
    // Energy deposit recording (for visualization)
    // ========================================================================

    /**
     * Add energy deposit to a lead glass block
     */
    void AddLeadGlassEnergyDeposit(int32_t copyNumber, double energy);

    /**
     * Add energy deposit to a scintillator bar
     */
    void AddScintillatorEnergyDeposit(int32_t moduleCopyNumber,
                                      int32_t barCopyNumber,
                                      double energy);

    /**
     * Add optical photon count
     */
    void AddLeadGlassPhotons(int32_t copyNumber, int32_t nPhotons);
    void AddScintillatorPhotons(int32_t moduleCopyNumber,
                                int32_t barCopyNumber,
                                int32_t nPhotons);

    /**
     * Clear all energy deposits (call at start of event)
     */
    void ClearEventData();

    // ========================================================================
    // Statistics and reporting
    // ========================================================================

    int32_t GetNumScintillatorModules() const { return m_scintModules.size(); }
    int32_t GetNumScintillatorBars() const { return m_scintBars.size(); }
    int32_t GetNumLeadGlassBlocks() const { return m_leadGlassBlocks.size(); }
    int32_t GetNumTPCModules() const { return m_tpcModules.size(); }
    int32_t GetNumBeampipeSections() const { return m_beampipeSections.size(); }

    void PrintStatistics() const;
    void PrintLeadGlassMap() const;
    void PrintScintillatorMap() const;

    // Get all data for visualization
    const std::map<int32_t, LeadGlassInfo>& GetAllLeadGlass() const { return m_leadGlassBlocks; }
    const std::vector<ScintillatorInfo>& GetAllScintillatorBars() const { return m_scintBars; }
    const std::vector<TPCModuleInfo>& GetAllTPCModules() const { return m_tpcModules; }

private:
    GeometryManager() = default;
    ~GeometryManager() = default;
    GeometryManager(const GeometryManager&) = delete;
    GeometryManager& operator=(const GeometryManager&) = delete;

    // Storage
    std::vector<VolumeInfo> m_scintModules;              // Scintillator module positions
    std::vector<ScintillatorInfo> m_scintBars;           // Individual bar info
    std::map<int32_t, LeadGlassInfo> m_leadGlassBlocks;  // By copy number
    std::vector<TPCModuleInfo> m_tpcModules;             // TPC modules
    std::vector<BeampipeInfo> m_beampipeSections;        // Beampipe sections

    // Lookup maps
    std::map<std::pair<int32_t, int32_t>, size_t> m_scintBarIndex;  // (module, bar) -> index

    bool m_initialized = false;
};

// ============================================================================
// Surface index constants
// ============================================================================

namespace Surface {
    constexpr int32_t TOP = 0;
    constexpr int32_t RIGHT = 1;
    constexpr int32_t BOTTOM = 2;
    constexpr int32_t LEFT = 3;
    constexpr int32_t FRONT = 4;
    constexpr int32_t BACK = 5;

    inline const char* GetName(int32_t index) {
        switch (index) {
            case TOP: return "Top";
            case RIGHT: return "Right";
            case BOTTOM: return "Bottom";
            case LEFT: return "Left";
            case FRONT: return "Front";
            case BACK: return "Back";
            default: return "Unknown";
        }
    }
}

} // namespace nnbar

#endif // GEOMETRY_MANAGER_HH
