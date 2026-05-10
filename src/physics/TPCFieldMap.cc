// ============================================================================
// TPCFieldMap.cc
// Electric field map implementation for MWPC TPC
// ============================================================================

#include "physics/TPCFieldMap.hh"
#include <cmath>
#include <iostream>
#include <algorithm>

namespace nnbar {

// ============================================================================
// Constructor
// ============================================================================

TPCFieldMap::TPCFieldMap() {
    m_driftVelocityTable.resize(N_FIELD_POINTS);
    m_diffusionLTable.resize(N_FIELD_POINTS);
    m_diffusionTTable.resize(N_FIELD_POINTS);
    m_townsendTable.resize(N_FIELD_POINTS);
    m_attachmentTable.resize(N_FIELD_POINTS);
}

// ============================================================================
// Configuration
// ============================================================================

void TPCFieldMap::SetMWPCParameters(const MWPCParameters& params) {
    m_mwpc = params;
    m_initialized = false;
}

void TPCFieldMap::SetGasProperties(const GasProperties& props) {
    m_gas = props;
    m_initialized = false;
}

void TPCFieldMap::SetTPCGeometry(float innerRadius, float outerRadius, float halfLength) {
    m_innerRadius = innerRadius;
    m_outerRadius = outerRadius;
    m_halfLength = halfLength;
    m_initialized = false;
}

void TPCFieldMap::Initialize() {
    if (m_initialized) return;

    std::cout << "TPCFieldMap: Initializing..." << std::endl;
    std::cout << "  Gas: Ar/CO2 " << (m_gas.arFraction * 100) << "/"
              << (m_gas.co2Fraction * 100) << std::endl;
    std::cout << "  Pressure: " << m_gas.pressure << " Torr" << std::endl;
    std::cout << "  Temperature: " << m_gas.temperature << " K" << std::endl;
    std::cout << "  Drift field: " << m_mwpc.driftField << " V/cm" << std::endl;
    std::cout << "  Drift length: " << m_mwpc.driftLength << " mm" << std::endl;

    BuildTransportTables();

    // Print nominal transport properties
    float nominalE = m_mwpc.driftField;
    std::cout << "  Drift velocity @ " << nominalE << " V/cm: "
              << GetDriftVelocity(nominalE) * 1000.0f << " um/ns" << std::endl;
    std::cout << "  Diffusion L: " << GetDiffusionL(nominalE) * 10000.0f << " um/sqrt(cm)" << std::endl;
    std::cout << "  Diffusion T: " << GetDiffusionT(nominalE) * 10000.0f << " um/sqrt(cm)" << std::endl;

    m_initialized = true;
    std::cout << "TPCFieldMap: Initialized" << std::endl;
}

// ============================================================================
// Field Calculation
// ============================================================================

void TPCFieldMap::GetField(float x, float y, float z,
                           float& Ex, float& Ey, float& Ez) const {
    // Convert from mm to cm for internal calculations
    float x_cm = x / 10.0f;
    float y_cm = y / 10.0f;
    float z_cm = z / 10.0f;

    // Check if in active region
    float r = std::sqrt(x_cm * x_cm + y_cm * y_cm);

    // Default: uniform drift field in z direction
    Ex = 0.0f;
    Ey = 0.0f;
    Ez = m_mwpc.driftField;  // V/cm, pointing toward anode

    // Near wire region - more complex field
    if (IsInAmplificationRegion(x, y, z)) {
        // MWPC field calculation near anode wires
        // Simplified model: radial field increase near wires
        float distFromAnode = std::abs(z_cm - m_mwpc.driftLength / 10.0f);

        if (distFromAnode < m_mwpc.gapHeight / 10.0f) {
            // Field enhancement factor near wires
            float enhancement = 1.0f + 10.0f * (1.0f - distFromAnode / (m_mwpc.gapHeight / 10.0f));
            Ez *= enhancement;
        }
    }

    // Boundary effects near inner/outer radius
    float r_inner = m_innerRadius / 10.0f;
    float r_outer = m_outerRadius / 10.0f;

    if (r < r_inner + 5.0f) {  // Within 5cm of inner boundary
        // Field lines curve near beampipe
        float radialFactor = 0.1f * (1.0f - (r - r_inner) / 5.0f);
        Ex = m_mwpc.driftField * radialFactor * x_cm / r;
        Ey = m_mwpc.driftField * radialFactor * y_cm / r;
    }
}

float TPCFieldMap::GetFieldMagnitude(float x, float y, float z) const {
    float Ex, Ey, Ez;
    GetField(x, y, z, Ex, Ey, Ez);
    return std::sqrt(Ex * Ex + Ey * Ey + Ez * Ez);
}

bool TPCFieldMap::IsInDriftRegion(float x, float y, float z) const {
    float r = std::sqrt(x * x + y * y);
    return (r > m_innerRadius && r < m_outerRadius &&
            std::abs(z) < m_halfLength);
}

bool TPCFieldMap::IsInAmplificationRegion(float x, float y, float z) const {
    // Amplification region is within gapHeight of anode plane
    float anodeZ = m_mwpc.driftLength;  // mm
    float distFromAnode = std::abs(std::abs(z) - anodeZ);
    return (distFromAnode < m_mwpc.gapHeight && IsInDriftRegion(x, y, z));
}

// ============================================================================
// Transport Properties
// ============================================================================

float TPCFieldMap::GetDriftVelocity(float E) const {
    return InterpolateTable(m_driftVelocityTable, E);
}

float TPCFieldMap::GetDiffusionL(float E) const {
    return InterpolateTable(m_diffusionLTable, E);
}

float TPCFieldMap::GetDiffusionT(float E) const {
    return InterpolateTable(m_diffusionTTable, E);
}

float TPCFieldMap::GetTownsendCoefficient(float E) const {
    return InterpolateTable(m_townsendTable, E);
}

float TPCFieldMap::GetAttachmentCoefficient(float E) const {
    return InterpolateTable(m_attachmentTable, E);
}

// ============================================================================
// Lookup Table Helpers
// ============================================================================

void TPCFieldMap::BuildTransportTables() {
    float co2Frac = m_gas.co2Fraction;
    float pressure = m_gas.pressure;
    float temperature = m_gas.temperature;

    for (int i = 0; i < N_FIELD_POINTS; i++) {
        float E = MIN_FIELD + (MAX_FIELD - MIN_FIELD) * i / (N_FIELD_POINTS - 1);

        // Convert from cm/us to cm/ns (divide by 1000)
        m_driftVelocityTable[i] = MagboltzDriftVelocity(E, co2Frac, pressure, temperature) / 1000.0f;
        m_diffusionLTable[i] = MagboltzDiffusionL(E, co2Frac, pressure, temperature);
        m_diffusionTTable[i] = MagboltzDiffusionT(E, co2Frac, pressure, temperature);
        m_townsendTable[i] = MagboltzTownsend(E, co2Frac, pressure);
        m_attachmentTable[i] = 0.0f;  // Negligible for Ar/CO2
    }
}

int TPCFieldMap::FieldToIndex(float E) const {
    E = std::max(MIN_FIELD, std::min(MAX_FIELD, E));
    return static_cast<int>((E - MIN_FIELD) / (MAX_FIELD - MIN_FIELD) * (N_FIELD_POINTS - 1));
}

float TPCFieldMap::InterpolateTable(const std::vector<float>& table, float E) const {
    E = std::max(MIN_FIELD, std::min(MAX_FIELD, E));

    float frac = (E - MIN_FIELD) / (MAX_FIELD - MIN_FIELD) * (N_FIELD_POINTS - 1);
    int idx = static_cast<int>(frac);
    float t = frac - idx;

    if (idx >= N_FIELD_POINTS - 1) {
        return table[N_FIELD_POINTS - 1];
    }

    return table[idx] * (1.0f - t) + table[idx + 1] * t;
}

// ============================================================================
// Magboltz Parametrizations for Ar/CO2
// Based on fits to Magboltz 11.6 simulation data
// ============================================================================

float MagboltzDriftVelocity(float E, float co2Frac, float pressure, float temperature) {
    // Reduced field: E/N (N is number density, proportional to pressure)
    float reducedE = E * (760.0f / pressure) * (temperature / 293.15f);

    // Empirical fit for Ar/CO2 mixtures
    // v = a * E_red / (1 + b * E_red + c * E_red^2)
    // Coefficients depend on CO2 fraction

    float a = 4.5f + 1.5f * co2Frac;      // Low field mobility factor
    float b = 0.008f + 0.004f * co2Frac;  // Medium field saturation
    float c = 1.0e-6f;                     // High field term

    float vDrift = a * reducedE / (1.0f + b * reducedE + c * reducedE * reducedE);

    // Typical range: 1-10 cm/us for 100-1000 V/cm
    return vDrift;  // cm/us
}

float MagboltzDiffusionL(float E, float co2Frac, float pressure, float temperature) {
    // Longitudinal diffusion coefficient
    // D_L = sigma_L^2 / (2 * drift_distance)
    // sigma_L ~ sqrt(2 * D * t) where D is diffusion coefficient

    // Empirical: D_L decreases with increasing E-field
    float reducedE = E * (760.0f / pressure) * (temperature / 293.15f);

    // Base diffusion (thermal limit)
    float D0 = 0.04f;  // cm/sqrt(cm) at low field

    // Field dependence: D_L ~ D0 / sqrt(1 + E/E0)
    float E0 = 500.0f;  // Reference field
    float diffL = D0 / std::sqrt(1.0f + reducedE / E0);

    // CO2 quenching effect (CO2 reduces diffusion)
    diffL *= (1.0f - 0.3f * co2Frac);

    return diffL;  // cm/sqrt(cm)
}

float MagboltzDiffusionT(float E, float co2Frac, float pressure, float temperature) {
    // Transverse diffusion is typically 1.2-1.5x longitudinal
    float diffL = MagboltzDiffusionL(E, co2Frac, pressure, temperature);
    return 1.3f * diffL;  // cm/sqrt(cm)
}

float MagboltzTownsend(float E, float co2Frac, float pressure) {
    // Townsend first ionization coefficient (alpha)
    // Only significant at high fields (near wires in MWPC)

    float reducedE = E * (760.0f / pressure);

    // Threshold field for ionization in Ar/CO2
    float Ethresh = 10000.0f;  // V/cm (effective threshold)

    if (reducedE < Ethresh) {
        return 0.0f;
    }

    // Empirical: alpha/N ~ A * exp(-B * N / E)
    // For Ar/CO2, approximate:
    float A = 30.0f * (1.0f - 0.5f * co2Frac);  // 1/cm at high field
    float B = 20000.0f;  // V/cm characteristic field

    float alpha = A * std::exp(-B / reducedE);

    return alpha;  // 1/cm
}

} // namespace nnbar
