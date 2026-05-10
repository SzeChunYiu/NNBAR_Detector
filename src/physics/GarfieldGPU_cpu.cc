// ============================================================================
// GarfieldGPU_cpu.cc
// CPU fallback implementation when CUDA is not available
// Uses OpenMP for parallelization on multi-core CPUs
// ============================================================================

#include "physics/GarfieldGPU.hh"
#include <iostream>
#include <random>
#include <cmath>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace nnbar {

// ============================================================================
// GarfieldGPU Implementation (CPU Fallback)
// ============================================================================

GarfieldGPU::GarfieldGPU() {
    // Default TPC geometry
    m_geometry.innerRadius = 114.0f;
    m_geometry.outerRadius = 200.0f;
    m_geometry.halfLength = 300.0f;
    m_geometry.driftLength = 85.0f;
    m_geometry.cathodeZ = 0.0f;
    m_geometry.anodeZ = 85.0f;

    // Default gas properties (Ar/CO2 80/20)
    m_gasData.electricField = 250.0f;
    m_gasData.driftVelocity = 0.0055f;
    m_gasData.diffusionL = 0.023f;
    m_gasData.diffusionT = 0.030f;
    m_gasData.townsendCoeff = 0.0f;
    m_gasData.attachmentCoeff = 0.0f;
}

GarfieldGPU::~GarfieldGPU() {
    // No GPU resources to free in CPU mode
}

bool GarfieldGPU::Initialize() {
    std::cout << "=============================================" << std::endl;
    std::cout << "GarfieldGPU: CPU Mode (CUDA not available)" << std::endl;
#ifdef _OPENMP
    int nThreads = omp_get_max_threads();
    std::cout << "  OpenMP threads: " << nThreads << std::endl;
    m_deviceName = "CPU (" + std::to_string(nThreads) + " threads)";
#else
    std::cout << "  Single-threaded mode" << std::endl;
    m_deviceName = "CPU (single-threaded)";
#endif
    std::cout << "=============================================" << std::endl;

    m_gpuAvailable = false;  // No GPU
    m_initialized = true;
    return true;
}

void GarfieldGPU::SetGeometry(const TPCGeometry& geom) {
    m_geometry = geom;
}

void GarfieldGPU::SetGasProperties(float arFrac, float co2Frac,
                                    float pressure, float temperature) {
    float E = m_gasData.electricField;
    float reducedField = E / (pressure / 760.0f);

    float a = 0.00015f + 0.00003f * co2Frac;
    float b = 0.00001f;
    float v0 = a * reducedField / (1.0f + b * reducedField);

    float tempCorrection = std::sqrt(temperature / 293.15f);
    m_gasData.driftVelocity = v0 * tempCorrection;

    float diffBase = 0.02f + 0.01f * co2Frac;
    m_gasData.diffusionL = diffBase * std::sqrt(293.15f / temperature);
    m_gasData.diffusionT = 1.3f * m_gasData.diffusionL;

    std::cout << "GarfieldGPU: Gas properties updated" << std::endl;
    std::cout << "  Drift velocity: " << m_gasData.driftVelocity << " cm/ns" << std::endl;
}

void GarfieldGPU::SetElectricField(float field) {
    m_gasData.electricField = field;
}

void GarfieldGPU::AddClusters(const std::vector<IonizationCluster>& clusters) {
    for (const auto& cluster : clusters) {
        for (int i = 0; i < cluster.nElectrons; i++) {
            AddElectron(cluster.x, cluster.y, cluster.z, cluster.t);
        }
    }
}

void GarfieldGPU::AddElectron(float x, float y, float z, float t) {
    ElectronState e;
    e.x = x;
    e.y = y;
    e.z = z;
    e.t = t;
    e.vx = e.vy = e.vz = 0.0f;
    e.status = 0;
    m_electrons.push_back(e);
}

void GarfieldGPU::ClearElectrons() {
    m_electrons.clear();
    m_results.clear();
    m_nCollected = 0;
    m_totalCharge = 0.0f;
}

// Single electron drift simulation
static DriftResult DriftSingleElectron(
    const ElectronState& e,
    const GasTransportData& gasData,
    const TPCGeometry& geometry,
    float maxTime,
    bool avalanche,
    std::mt19937& rng
) {
    DriftResult result;
    result.driftTime = 0.0f;
    result.pathLength = 0.0f;
    result.nAvalanche = 1;
    result.status = 0;

    std::normal_distribution<float> normal(0.0f, 1.0f);
    std::uniform_real_distribution<float> uniform(0.0f, 1.0f);

    float vDrift = gasData.driftVelocity;
    float diffL = gasData.diffusionL;
    float diffT = gasData.diffusionT;
    float alpha = gasData.townsendCoeff;
    float eta = gasData.attachmentCoeff;

    float rInner = geometry.innerRadius;
    float rOuter = geometry.outerRadius;
    float halfZ = geometry.halfLength;
    float anodeZ = geometry.anodeZ;
    float cathodeZ = geometry.cathodeZ;

    float driftDir = (anodeZ > cathodeZ) ? 1.0f : -1.0f;
    float dt = 1.0f;  // ns

    float t = e.t;
    float x = e.x;
    float y = e.y;
    float z = e.z;

    int nSteps = 0;
    const int maxSteps = 100000;

    while (t < maxTime && nSteps < maxSteps) {
        float dz = vDrift * dt * driftDir;
        float stepLength = std::abs(dz);

        float sigmaL = diffL * std::sqrt(stepLength);
        dz += sigmaL * normal(rng);

        float sigmaT = diffT * std::sqrt(stepLength);
        float dx = sigmaT * normal(rng);
        float dy = sigmaT * normal(rng);

        x += dx;
        y += dy;
        z += dz;
        t += dt;
        result.pathLength += std::sqrt(dx*dx + dy*dy + dz*dz);

        float r = std::sqrt(x*x + y*y);

        if (r < rInner) {
            result.status = 2;
            break;
        }
        if (r > rOuter || std::abs(z) > halfZ) {
            result.status = 3;
            break;
        }

        if ((driftDir > 0 && z >= anodeZ) || (driftDir < 0 && z <= anodeZ)) {
            result.status = 1;
            break;
        }

        if (avalanche && alpha > 0) {
            float distToAnode = std::abs(z - anodeZ);
            if (distToAnode < 1.0f) {
                float pAvalanche = 1.0f - std::exp(-alpha * stepLength);
                if (uniform(rng) < pAvalanche) {
                    result.nAvalanche++;
                }
            }
        }

        if (eta > 0) {
            float pAttach = 1.0f - std::exp(-eta * stepLength);
            if (uniform(rng) < pAttach) {
                result.status = 2;
                break;
            }
        }

        nSteps++;
    }

    result.finalX = x;
    result.finalY = y;
    result.finalZ = z;
    result.driftTime = t - e.t;

    return result;
}

int GarfieldGPU::DriftElectrons() {
    if (m_electrons.empty()) {
        return 0;
    }

    int nElectrons = static_cast<int>(m_electrons.size());
    m_results.resize(nElectrons);

    auto startTime = std::chrono::high_resolution_clock::now();

#ifdef _OPENMP
    #pragma omp parallel
    {
        // Each thread gets its own RNG
        std::mt19937 rng(42 + omp_get_thread_num());

        #pragma omp for
        for (int i = 0; i < nElectrons; i++) {
            m_results[i] = DriftSingleElectron(
                m_electrons[i],
                m_gasData,
                m_geometry,
                m_maxDriftTime,
                m_avalancheEnabled,
                rng
            );
        }
    }
#else
    std::mt19937 rng(42);
    for (int i = 0; i < nElectrons; i++) {
        m_results[i] = DriftSingleElectron(
            m_electrons[i],
            m_gasData,
            m_geometry,
            m_maxDriftTime,
            m_avalancheEnabled,
            rng
        );
    }
#endif

    auto endTime = std::chrono::high_resolution_clock::now();
    m_kernelTimeMs = std::chrono::duration<float, std::milli>(endTime - startTime).count();

    // Count results
    m_nCollected = 0;
    m_totalCharge = 0.0f;
    for (const auto& r : m_results) {
        if (r.status == 1) {
            m_nCollected++;
            m_totalCharge += static_cast<float>(r.nAvalanche);
        }
    }

    std::cout << "GarfieldGPU (CPU): Drifted " << nElectrons << " electrons in "
              << m_kernelTimeMs << " ms" << std::endl;
    std::cout << "  Collected: " << m_nCollected
              << " (" << (100.0f * m_nCollected / nElectrons) << "%)" << std::endl;

    return m_nCollected;
}

// ============================================================================
// Utility Functions
// ============================================================================

float CalculateDriftVelocity(float electricField, float pressure,
                              float temperature, float arFrac) {
    float reducedField = electricField / (pressure / 760.0f);
    float co2Frac = 1.0f - arFrac;

    float a = 0.00015f + 0.00003f * co2Frac;
    float b = 0.00001f;
    float v0 = a * reducedField / (1.0f + b * reducedField);

    return v0 * std::sqrt(temperature / 293.15f);
}

void CalculateDiffusion(float electricField, float pressure,
                        float& diffL, float& diffT, float arFrac) {
    float co2Frac = 1.0f - arFrac;
    float diffBase = 0.02f + 0.01f * co2Frac;
    float fieldFactor = 1.0f + 100.0f / electricField;

    diffL = diffBase * std::sqrt(fieldFactor);
    diffT = 1.3f * diffL;
}

std::vector<IonizationCluster> GenerateIonization(
    float x0, float y0, float z0,
    float dx, float dy, float dz,
    float energy,
    float pathLength,
    float dEdx
) {
    std::vector<IonizationCluster> clusters;

    const float ionizationEnergy = 26.0e-6f;  // MeV (for Ar)

    float totalEdep = dEdx * pathLength;
    if (totalEdep > energy) totalEdep = energy;

    int totalElectrons = static_cast<int>(totalEdep / ionizationEnergy);
    if (totalElectrons < 1) totalElectrons = 1;

    float meanClusterSize = 3.0f;
    int nClusters = std::max(1, static_cast<int>(totalElectrons / meanClusterSize));

    for (int i = 0; i < nClusters; i++) {
        float frac = (i + 0.5f) / nClusters;

        IonizationCluster cluster;
        cluster.x = x0 + dx * pathLength * frac;
        cluster.y = y0 + dy * pathLength * frac;
        cluster.z = z0 + dz * pathLength * frac;
        cluster.t = 0.0f;
        cluster.nElectrons = totalElectrons / nClusters;
        cluster.energy = totalEdep / nClusters * 1e6f;

        clusters.push_back(cluster);
    }

    return clusters;
}

} // namespace nnbar
