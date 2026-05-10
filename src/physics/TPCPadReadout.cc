// ============================================================================
// TPCPadReadout.cc
// TPC pad plane readout implementation for MWPC rectangular TPC
// ============================================================================

#include "physics/TPCPadReadout.hh"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <numeric>

namespace nnbar {

// ============================================================================
// Constructor
// ============================================================================

TPCPadReadout::TPCPadReadout() {
    // Default parameters based on NNBAR TPC geometry
    // TPC drift length: 850 mm, TPC z-extent: ~2500 mm
    m_params.padWidth = 4.0f;       // 4mm pads (similar to ALICE)
    m_params.padLength = 7.5f;      // 7.5mm for inner region
    m_params.padLengthLong = 15.0f; // 15mm for outer region
    m_params.nPadRows = 85;         // Match TPC layers (drift/10mm)
    m_params.nPadsPerRow = 250;     // Cover ~1000mm width with 4mm pads
    m_params.nInnerRows = 42;       // Inner half uses short pads
    m_params.nOuterRows = 43;       // Outer half uses long pads
    m_params.innerRadius = 1120.0f; // Beampipe radius
    m_params.outerRadius = 1970.0f; // After drift
    m_params.xMin = -1000.0f;
    m_params.xMax = 1000.0f;
    m_params.yMin = -1000.0f;
    m_params.yMax = 1000.0f;
    m_params.anodePosZ = 850.0f;
    m_params.anodeNegZ = -850.0f;
    m_params.shaperPeakingTime = 160.0f;
    m_params.samplingRate = 10.0f;
    m_params.nTimeBins = 1000;
    m_params.electronGain = 8000.0f;
    m_params.padResponseWidth = 3.0f;
}

// ============================================================================
// Configuration
// ============================================================================

void TPCPadReadout::SetParameters(const TPCPadParameters& params) {
    m_params = params;
    m_initialized = false;
}

void TPCPadReadout::Initialize() {
    if (m_initialized) return;

    std::cout << "TPCPadReadout: Initializing..." << std::endl;
    std::cout << "  Pad rows: " << m_params.nPadRows << std::endl;
    std::cout << "  Pads per row: " << m_params.nPadsPerRow << std::endl;
    std::cout << "  Pad size: " << m_params.padWidth << " x " << m_params.padLength << " mm" << std::endl;
    std::cout << "  Time bins: " << m_params.nTimeBins << std::endl;

    CreatePadMap();

    // Initialize signal buffers for both sides
    for (int side = 0; side < 2; side++) {
        m_signalBuffer[side].resize(m_params.nPadRows);
        for (int row = 0; row < m_params.nPadRows; row++) {
            m_signalBuffer[side][row].resize(m_params.nPadsPerRow);
            for (int pad = 0; pad < m_params.nPadsPerRow; pad++) {
                m_signalBuffer[side][row][pad].resize(m_params.nTimeBins, 0.0f);
            }
        }
    }

    m_initialized = true;
    std::cout << "TPCPadReadout: Initialized" << std::endl;
}

// ============================================================================
// Pad Map Creation
// ============================================================================

void TPCPadReadout::CreatePadMap() {
    m_padMap.resize(m_params.nPadRows);

    // Calculate pad positions for rectangular geometry
    float totalWidth = m_params.xMax - m_params.xMin;
    float padSpacing = totalWidth / m_params.nPadsPerRow;

    for (int row = 0; row < m_params.nPadRows; row++) {
        m_padMap[row].resize(m_params.nPadsPerRow);

        // Row position (radial distance from center)
        float rowDist = m_params.innerRadius +
                        (row + 0.5f) * (m_params.outerRadius - m_params.innerRadius) / m_params.nPadRows;

        // Pad length depends on region
        float padLen = (row < m_params.nInnerRows) ? m_params.padLength : m_params.padLengthLong;

        for (int pad = 0; pad < m_params.nPadsPerRow; pad++) {
            // Pad center position
            float x = m_params.xMin + (pad + 0.5f) * padSpacing;
            float y = rowDist;  // In local coordinates

            m_padMap[row][pad].x = x;
            m_padMap[row][pad].y = y;
            m_padMap[row][pad].width = m_params.padWidth;
            m_padMap[row][pad].length = padLen;
        }
    }
}

// ============================================================================
// Signal Simulation
// ============================================================================

void TPCPadReadout::AddElectronSignal(float x, float y, float t, int nElectrons, int side) {
    if (!m_initialized) {
        Initialize();
    }

    // Find affected pads using pad response function
    int centerRow, centerPad;
    if (!GetPadFromPosition(x, y, centerRow, centerPad)) {
        return;  // Outside pad plane
    }

    // Calculate charge including gas gain
    float charge = static_cast<float>(nElectrons) * m_params.electronGain;

    // Spread charge over nearby pads using PRF
    int rowRange = 2;  // Affect ±2 rows
    int padRange = 3;  // Affect ±3 pads

    for (int dr = -rowRange; dr <= rowRange; dr++) {
        int row = centerRow + dr;
        if (row < 0 || row >= m_params.nPadRows) continue;

        for (int dp = -padRange; dp <= padRange; dp++) {
            int pad = centerPad + dp;
            if (pad < 0 || pad >= m_params.nPadsPerRow) continue;

            // Calculate position difference
            float padX, padY;
            GetPadPosition(row, pad, padX, padY);
            float dx = x - padX;
            float dy = y - padY;

            // Calculate pad response
            float prfWeight = PadResponseFunction(dx, dy);
            if (prfWeight < 0.001f) continue;

            // Add charge with shaper response over time
            float weightedCharge = charge * prfWeight;
            AddChargeToSignal(row, pad, side, t, weightedCharge);
        }
    }
}

float TPCPadReadout::PadResponseFunction(float dx, float dy) const {
    // 2D Gaussian pad response function
    float sigma = m_params.padResponseWidth;
    float r2 = (dx * dx + dy * dy) / (2.0f * sigma * sigma);
    return std::exp(-r2);
}

float TPCPadReadout::ShaperResponse(float dt) const {
    // Semi-Gaussian shaper response (ALICE PASA style)
    // Peaking time = 160 ns
    if (dt < 0) return 0.0f;

    float tau = m_params.shaperPeakingTime / 2.5f;  // Time constant
    float t_tau = dt / tau;

    // Response: (t/tau)^4 * exp(-4*t/tau)
    float x = t_tau;
    float response = x * x * x * x * std::exp(-4.0f * x);

    return response * 4.0f;  // Normalize peak to ~1
}

void TPCPadReadout::AddChargeToSignal(int row, int pad, int side, float t, float charge) {
    // Convert time to time bin
    float timeBinSize = 1000.0f / m_params.samplingRate;  // ns per bin
    int centerBin = static_cast<int>(t / timeBinSize);

    // Add shaped signal over time bins
    int binRange = static_cast<int>(5.0f * m_params.shaperPeakingTime / timeBinSize);

    for (int db = -binRange; db <= binRange; db++) {
        int bin = centerBin + db;
        if (bin < 0 || bin >= m_params.nTimeBins) continue;

        float dt = db * timeBinSize;
        float shaperWeight = ShaperResponse(std::abs(dt));

        m_signalBuffer[side][row][pad][bin] += charge * shaperWeight;
    }
}

// ============================================================================
// Digitization and Cluster Finding
// ============================================================================

void TPCPadReadout::Digitize() {
    m_padSignals.clear();

    // ADC threshold (in electrons * gain)
    float threshold = 100.0f;  // Minimum signal to record

    for (int side = 0; side < 2; side++) {
        for (int row = 0; row < m_params.nPadRows; row++) {
            for (int pad = 0; pad < m_params.nPadsPerRow; pad++) {
                const auto& samples = m_signalBuffer[side][row][pad];

                // Find max amplitude
                float maxAmp = 0.0f;
                int peakBin = 0;
                for (int bin = 0; bin < m_params.nTimeBins; bin++) {
                    if (samples[bin] > maxAmp) {
                        maxAmp = samples[bin];
                        peakBin = bin;
                    }
                }

                if (maxAmp < threshold) continue;

                // Calculate total charge
                float totalCharge = 0.0f;
                for (int bin = 0; bin < m_params.nTimeBins; bin++) {
                    totalCharge += samples[bin];
                }

                // Create pad signal
                PadSignal sig;
                sig.padRow = row;
                sig.padNum = pad;
                sig.side = side;
                GetPadPosition(row, pad, sig.x, sig.y);
                sig.timeSamples = samples;
                sig.totalCharge = totalCharge;
                sig.maxAmplitude = maxAmp;
                sig.peakTimeBin = peakBin;

                m_padSignals.push_back(sig);
            }
        }
    }

    std::cout << "TPCPadReadout: Digitized " << m_padSignals.size() << " pad signals" << std::endl;
}

void TPCPadReadout::FindClusters() {
    m_clusters.clear();

    // Simple clustering: group adjacent pad signals
    // More sophisticated clustering would use connected-component analysis

    std::vector<bool> used(m_padSignals.size(), false);

    for (size_t i = 0; i < m_padSignals.size(); i++) {
        if (used[i]) continue;

        const auto& seed = m_padSignals[i];

        TPCCluster cluster;
        cluster.x = seed.x * seed.totalCharge;
        cluster.y = seed.y * seed.totalCharge;
        cluster.charge = seed.totalCharge;
        cluster.padRow = seed.padRow;
        cluster.nPads = 1;
        cluster.nTimeBins = 0;

        used[i] = true;

        // Find neighboring signals
        for (size_t j = i + 1; j < m_padSignals.size(); j++) {
            if (used[j]) continue;

            const auto& neighbor = m_padSignals[j];

            // Check if adjacent
            if (neighbor.side != seed.side) continue;
            if (std::abs(neighbor.padRow - seed.padRow) > 2) continue;
            if (std::abs(neighbor.padNum - seed.padNum) > 3) continue;
            if (std::abs(neighbor.peakTimeBin - seed.peakTimeBin) > 10) continue;

            // Add to cluster
            cluster.x += neighbor.x * neighbor.totalCharge;
            cluster.y += neighbor.y * neighbor.totalCharge;
            cluster.charge += neighbor.totalCharge;
            cluster.nPads++;

            used[j] = true;
        }

        // Calculate centroid
        cluster.x /= cluster.charge;
        cluster.y /= cluster.charge;

        // Calculate z from drift time
        float timeBinSize = 1000.0f / m_params.samplingRate;
        float driftVelocity = 0.0055f * 10.0f;  // cm/ns -> mm/ns
        cluster.z = seed.peakTimeBin * timeBinSize * driftVelocity;
        if (seed.side == 1) cluster.z = -cluster.z;

        cluster.t = seed.peakTimeBin * timeBinSize;

        m_clusters.push_back(cluster);
    }

    std::cout << "TPCPadReadout: Found " << m_clusters.size() << " clusters" << std::endl;
}

void TPCPadReadout::ClearEvent() {
    // Clear signal buffers
    for (int side = 0; side < 2; side++) {
        for (auto& row : m_signalBuffer[side]) {
            for (auto& pad : row) {
                std::fill(pad.begin(), pad.end(), 0.0f);
            }
        }
    }

    m_padSignals.clear();
    m_clusters.clear();
}

// ============================================================================
// Data Access
// ============================================================================

std::vector<std::tuple<int, int, float>> TPCPadReadout::GetPadOccupancy(int side) const {
    std::vector<std::tuple<int, int, float>> occupancy;

    for (const auto& sig : m_padSignals) {
        if (sig.side == side) {
            occupancy.emplace_back(sig.padRow, sig.padNum, sig.totalCharge);
        }
    }

    return occupancy;
}

void TPCPadReadout::GetPadPosition(int row, int pad, float& x, float& y) const {
    if (row >= 0 && row < m_params.nPadRows &&
        pad >= 0 && pad < m_params.nPadsPerRow) {
        x = m_padMap[row][pad].x;
        y = m_padMap[row][pad].y;
    } else {
        x = y = 0.0f;
    }
}

bool TPCPadReadout::GetPadFromPosition(float x, float y, int& row, int& pad) const {
    // Convert position to pad indices
    float totalWidth = m_params.xMax - m_params.xMin;
    float padSpacing = totalWidth / m_params.nPadsPerRow;

    pad = static_cast<int>((x - m_params.xMin) / padSpacing);
    if (pad < 0 || pad >= m_params.nPadsPerRow) return false;

    // Row from radial distance (simplified)
    row = static_cast<int>((y - m_params.innerRadius) /
                           (m_params.outerRadius - m_params.innerRadius) * m_params.nPadRows);
    if (row < 0 || row >= m_params.nPadRows) return false;

    return true;
}

int TPCPadReadout::GetNHitPads(int side) const {
    int count = 0;
    for (const auto& sig : m_padSignals) {
        if (sig.side == side) count++;
    }
    return count;
}

} // namespace nnbar
