// ============================================================================
// TPCPadReadout.hh
// TPC pad plane readout for MWPC-type rectangular TPC (ALICE TPC style)
// Simulates induced signals on readout pads from drifting electrons
// ============================================================================

#ifndef TPC_PAD_READOUT_HH
#define TPC_PAD_READOUT_HH

#include "config.h"
#include <vector>
#include <array>
#include <cstdint>
#include <memory>

namespace nnbar {

/**
 * @struct TPCPadParameters
 * @brief Configuration for TPC pad plane geometry
 * Based on ALICE TPC-style rectangular pad layout
 */
struct TPCPadParameters {
    // Pad dimensions (mm)
    float padWidth = 4.0f;          // Pad width in r-phi direction
    float padLength = 7.5f;         // Pad length in z direction (short pads)
    float padLengthLong = 15.0f;    // Longer pads in outer regions

    // Pad plane layout
    int nPadRows = 159;             // Number of pad rows (radially)
    int nPadsPerRow = 128;          // Number of pads per row (azimuthally)
    int nInnerRows = 63;            // Short pad region (inner rows)
    int nOuterRows = 96;            // Long pad region (outer rows)

    // Pad plane position (mm from center)
    float innerRadius = 850.0f;     // Inner radius of pad plane
    float outerRadius = 2500.0f;    // Outer radius of pad plane

    // For rectangular TPC: define x/y extent instead of cylindrical
    float xMin = -2000.0f;          // Min x coordinate
    float xMax = 2000.0f;           // Max x coordinate
    float yMin = -2000.0f;          // Min y coordinate
    float yMax = 2000.0f;           // Max y coordinate

    // Readout sides (z positions)
    float anodePosZ = 850.0f;       // +z anode position (mm)
    float anodeNegZ = -850.0f;      // -z anode position (mm)

    // Signal parameters
    float shaperPeakingTime = 160.0f;  // ns (ALICE PASA shaper)
    float samplingRate = 10.0f;        // MHz (time bins per us)
    int nTimeBins = 1000;              // Total time bins per drift
    float electronGain = 8000.0f;      // Gas gain at wires
    float padResponseWidth = 3.0f;     // PRF sigma in mm
};

/**
 * @struct PadSignal
 * @brief Signal on a single pad
 */
struct PadSignal {
    int padRow;                     // Pad row index (0 = innermost)
    int padNum;                     // Pad number in row
    int side;                       // 0 = +z side, 1 = -z side
    float x, y;                     // Pad center position (mm)
    std::vector<float> timeSamples; // ADC samples over time bins
    float totalCharge;              // Integrated charge
    float maxAmplitude;             // Peak amplitude
    int peakTimeBin;                // Time bin of peak
};

/**
 * @struct TPCCluster
 * @brief Reconstructed space point from pad signals
 */
struct TPCCluster {
    float x, y, z;                  // Reconstructed position (mm)
    float t;                        // Reconstructed time (ns)
    float charge;                   // Total charge (ADC counts)
    int padRow;                     // Pad row
    int nPads;                      // Number of pads in cluster
    int nTimeBins;                  // Number of time bins
    float sigmaX, sigmaY, sigmaZ;   // Cluster widths
};

/**
 * @class TPCPadReadout
 * @brief TPC pad plane readout simulation
 *
 * Simulates the signal induction on rectangular readout pads from
 * electrons arriving at the MWPC anode plane.
 *
 * Features:
 * - Pad response function (PRF) for charge sharing
 * - Shaper response for time development
 * - Digitization to ADC counts
 * - Cluster finding
 */
class TPCPadReadout {
public:
    TPCPadReadout();
    ~TPCPadReadout() = default;

    // ========== Configuration ==========

    /**
     * Set pad plane parameters
     */
    void SetParameters(const TPCPadParameters& params);

    /**
     * Initialize pad plane (create pad map)
     */
    void Initialize();

    // ========== Signal Simulation ==========

    /**
     * Add electron arrival at anode
     * @param x, y Position on anode plane (mm)
     * @param t Arrival time (ns)
     * @param nElectrons Number of electrons (after gas gain)
     * @param side 0 = +z anode, 1 = -z anode
     */
    void AddElectronSignal(float x, float y, float t, int nElectrons, int side);

    /**
     * Process all signals and digitize
     * Called at end of event
     */
    void Digitize();

    /**
     * Find clusters from pad signals
     */
    void FindClusters();

    /**
     * Clear all signals for new event
     */
    void ClearEvent();

    // ========== Data Access ==========

    /**
     * Get pad signals with non-zero charge
     */
    const std::vector<PadSignal>& GetPadSignals() const { return m_padSignals; }

    /**
     * Get reconstructed clusters
     */
    const std::vector<TPCCluster>& GetClusters() const { return m_clusters; }

    /**
     * Get 2D pad occupancy map (for visualization)
     * Returns: vector of (row, pad, charge) for all hit pads
     */
    std::vector<std::tuple<int, int, float>> GetPadOccupancy(int side) const;

    /**
     * Get pad position in local coordinates
     * @param row Pad row index
     * @param pad Pad number in row
     * @param x, y Output: pad center position (mm)
     */
    void GetPadPosition(int row, int pad, float& x, float& y) const;

    /**
     * Get pad indices from position
     * @return true if position is on a valid pad
     */
    bool GetPadFromPosition(float x, float y, int& row, int& pad) const;

    // ========== Statistics ==========

    int GetNPadRows() const { return m_params.nPadRows; }
    int GetNPadsPerRow() const { return m_params.nPadsPerRow; }
    int GetNTimeBins() const { return m_params.nTimeBins; }
    int GetNHitPads(int side) const;
    int GetNClusters() const { return static_cast<int>(m_clusters.size()); }

private:
    TPCPadParameters m_params;
    bool m_initialized = false;

    // Pad geometry lookup
    struct PadInfo {
        float x, y;         // Center position
        float width, length; // Dimensions
    };
    std::vector<std::vector<PadInfo>> m_padMap;  // [row][pad]

    // Signal buffers
    // [side][row][pad] -> time samples
    std::array<std::vector<std::vector<std::vector<float>>>, 2> m_signalBuffer;

    // Output data
    std::vector<PadSignal> m_padSignals;
    std::vector<TPCCluster> m_clusters;

    // Helper methods
    void CreatePadMap();
    float PadResponseFunction(float dx, float dy) const;
    float ShaperResponse(float dt) const;
    void AddChargeToSignal(int row, int pad, int side, float t, float charge);
};

} // namespace nnbar

#endif // TPC_PAD_READOUT_HH
