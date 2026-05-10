// ============================================================================
// TPCFieldMap.hh
// Electric field map for MWPC-type TPC with Ar/CO2 gas mixture
// Provides field values and gas transport properties for electron drift
// ============================================================================

#ifndef TPC_FIELD_MAP_HH
#define TPC_FIELD_MAP_HH

#include "config.h"
#include <vector>
#include <cmath>

namespace nnbar {

/**
 * @struct MWPCParameters
 * @brief Parameters for Multi-Wire Proportional Chamber configuration
 */
struct MWPCParameters {
    // Wire geometry
    float wireSpacing = 2.0f;       // Wire pitch in mm (center-to-center)
    float wireRadius = 0.010f;      // Wire radius in mm (10 um typical)
    float cathodeWireSpacing = 4.0f; // Cathode wire pitch in mm
    float gapHeight = 4.0f;         // Gap between anode and cathode planes (mm)

    // Voltages (in V)
    float anodeVoltage = 2000.0f;   // Anode wire voltage
    float cathodeVoltage = 0.0f;    // Cathode wire voltage
    float driftVoltage = -5000.0f;  // Drift electrode voltage

    // Drift region
    float driftLength = 850.0f;     // Drift length in mm
    float driftField = 250.0f;      // Nominal drift field V/cm
};

/**
 * @struct GasProperties
 * @brief Transport properties for Ar/CO2 gas mixture
 * Pre-computed from Magboltz for various E-field values
 */
struct GasProperties {
    float arFraction = 0.80f;       // Argon fraction
    float co2Fraction = 0.20f;      // CO2 fraction
    float pressure = 760.0f;        // Pressure in Torr
    float temperature = 293.15f;    // Temperature in Kelvin

    // Transport parameters (at nominal drift field)
    float driftVelocity = 0.0055f;  // cm/ns
    float diffusionL = 0.023f;      // cm/sqrt(cm) longitudinal
    float diffusionT = 0.030f;      // cm/sqrt(cm) transverse
    float townsendCoeff = 0.0f;     // 1/cm (for avalanche)
    float attachmentCoeff = 0.0f;   // 1/cm
    float ionMobility = 1.5e-4f;    // cm^2/(V*s) for Ar+

    // Ionization parameters
    float wValue = 26.0f;           // eV per ion pair
    float fanoFactor = 0.17f;       // Fano factor for Ar
};

/**
 * @class TPCFieldMap
 * @brief Electric field calculator for MWPC TPC
 *
 * Calculates E-field at any point in the TPC volume including:
 * - Uniform drift field in the bulk
 * - Non-uniform field near MWPC wires
 * - Field distortions near boundaries
 *
 * Also provides gas transport properties as function of E-field.
 */
class TPCFieldMap {
public:
    TPCFieldMap();
    ~TPCFieldMap() = default;

    // ========== Configuration ==========

    /**
     * Set MWPC geometry parameters
     */
    void SetMWPCParameters(const MWPCParameters& params);

    /**
     * Set gas mixture properties
     */
    void SetGasProperties(const GasProperties& props);

    /**
     * Set TPC geometry (all in mm)
     * @param innerRadius Inner radius (beampipe)
     * @param outerRadius Outer radius
     * @param halfLength Half-length in z
     */
    void SetTPCGeometry(float innerRadius, float outerRadius, float halfLength);

    /**
     * Initialize field map (pre-compute lookup tables)
     */
    void Initialize();

    // ========== Field Calculation ==========

    /**
     * Get electric field at a point (x, y, z in mm)
     * Returns field vector (Ex, Ey, Ez) in V/cm
     */
    void GetField(float x, float y, float z,
                  float& Ex, float& Ey, float& Ez) const;

    /**
     * Get field magnitude at a point (mm -> V/cm)
     */
    float GetFieldMagnitude(float x, float y, float z) const;

    /**
     * Check if point is in active drift region
     */
    bool IsInDriftRegion(float x, float y, float z) const;

    /**
     * Check if point is in amplification region (near wires)
     */
    bool IsInAmplificationRegion(float x, float y, float z) const;

    // ========== Transport Properties ==========

    /**
     * Get drift velocity at given E-field (cm/ns)
     * Uses Magboltz-based parametrization for Ar/CO2
     */
    float GetDriftVelocity(float E) const;

    /**
     * Get longitudinal diffusion coefficient (cm/sqrt(cm))
     */
    float GetDiffusionL(float E) const;

    /**
     * Get transverse diffusion coefficient (cm/sqrt(cm))
     */
    float GetDiffusionT(float E) const;

    /**
     * Get Townsend coefficient for avalanche (1/cm)
     */
    float GetTownsendCoefficient(float E) const;

    /**
     * Get attachment coefficient (1/cm)
     */
    float GetAttachmentCoefficient(float E) const;

    // ========== Accessors ==========

    const MWPCParameters& GetMWPCParams() const { return m_mwpc; }
    const GasProperties& GetGasProps() const { return m_gas; }
    float GetNominalDriftField() const { return m_mwpc.driftField; }

private:
    MWPCParameters m_mwpc;
    GasProperties m_gas;

    // TPC geometry (mm)
    float m_innerRadius = 1140.0f;
    float m_outerRadius = 2000.0f;
    float m_halfLength = 3000.0f;

    bool m_initialized = false;

    // Lookup tables for gas transport (indexed by E-field in V/cm)
    static constexpr int N_FIELD_POINTS = 100;
    static constexpr float MIN_FIELD = 50.0f;   // V/cm
    static constexpr float MAX_FIELD = 10000.0f; // V/cm

    std::vector<float> m_driftVelocityTable;
    std::vector<float> m_diffusionLTable;
    std::vector<float> m_diffusionTTable;
    std::vector<float> m_townsendTable;
    std::vector<float> m_attachmentTable;

    // Helper methods
    void BuildTransportTables();
    float InterpolateTable(const std::vector<float>& table, float E) const;
    int FieldToIndex(float E) const;
};

// ============================================================================
// Magboltz-based parametrization functions for Ar/CO2 80/20
// These are empirical fits to Magboltz simulation data
// ============================================================================

/**
 * Calculate drift velocity for Ar/CO2 mixture
 * @param E Electric field in V/cm
 * @param co2Frac CO2 fraction (0-1)
 * @param pressure Pressure in Torr
 * @param temperature Temperature in K
 * @return Drift velocity in cm/us
 */
float MagboltzDriftVelocity(float E, float co2Frac, float pressure, float temperature);

/**
 * Calculate longitudinal diffusion coefficient
 * @return Diffusion in cm/sqrt(cm)
 */
float MagboltzDiffusionL(float E, float co2Frac, float pressure, float temperature);

/**
 * Calculate transverse diffusion coefficient
 * @return Diffusion in cm/sqrt(cm)
 */
float MagboltzDiffusionT(float E, float co2Frac, float pressure, float temperature);

/**
 * Calculate Townsend coefficient (for avalanche near wires)
 * @return Alpha in 1/cm
 */
float MagboltzTownsend(float E, float co2Frac, float pressure);

} // namespace nnbar

#endif // TPC_FIELD_MAP_HH
