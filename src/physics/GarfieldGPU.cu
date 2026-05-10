// ============================================================================
// GarfieldGPU.cu
// CUDA implementation of GPU-accelerated electron drift simulation
// ============================================================================

#include "physics/GarfieldGPU.hh"
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cmath>
#include <iostream>

namespace nnbar {

// ============================================================================
// CUDA Error Checking Macro
// ============================================================================
#define CUDA_CHECK(call) do { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        std::cerr << "CUDA Error: " << cudaGetErrorString(err) \
                  << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
        return; \
    } \
} while(0)

// ============================================================================
// Device Constants
// ============================================================================
__constant__ GasTransportData d_gasData;
__constant__ TPCGeometry d_geometry;

// ============================================================================
// CUDA Kernel: Initialize Random States
// ============================================================================
__global__ void InitRandomStates(curandState* states, unsigned long seed, int n) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        curand_init(seed, idx, 0, &states[idx]);
    }
}

// ============================================================================
// CUDA Kernel: Electron Drift with Diffusion
// ============================================================================
/**
 * Each thread handles one electron's drift trajectory
 * Uses Langevin equation with thermal diffusion:
 *   dx = v_drift * dt + sqrt(2*D*dt) * N(0,1)
 */
__global__ void DriftKernel(
    ElectronState* electrons,
    DriftResult* results,
    curandState* randStates,
    int nElectrons,
    float dt,           // Time step (ns)
    float maxTime,      // Maximum simulation time (ns)
    bool avalanche      // Enable avalanche multiplication
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nElectrons) return;

    // Get electron state
    ElectronState e = electrons[idx];
    if (e.status != 0) return;  // Skip inactive electrons

    // Get local random state
    curandState localRand = randStates[idx];

    // Initialize result
    DriftResult result;
    result.driftTime = 0.0f;
    result.pathLength = 0.0f;
    result.nAvalanche = 1;  // Start with 1 electron
    result.status = 0;

    // Drift parameters from gas data
    float vDrift = d_gasData.driftVelocity;   // cm/ns
    float diffL = d_gasData.diffusionL;       // cm/sqrt(cm)
    float diffT = d_gasData.diffusionT;       // cm/sqrt(cm)
    float alpha = d_gasData.townsendCoeff;    // 1/cm (avalanche)
    float eta = d_gasData.attachmentCoeff;    // 1/cm (attachment)

    // Geometry
    float rInner = d_geometry.innerRadius;
    float rOuter = d_geometry.outerRadius;
    float halfZ = d_geometry.halfLength;
    float anodeZ = d_geometry.anodeZ;
    float cathodeZ = d_geometry.cathodeZ;

    // Determine drift direction (toward anode)
    float driftDir = (anodeZ > cathodeZ) ? 1.0f : -1.0f;

    // Drift loop
    float t = e.t;
    float x = e.x;
    float y = e.y;
    float z = e.z;

    int nSteps = 0;
    const int maxSteps = 100000;

    while (t < maxTime && nSteps < maxSteps) {
        // Calculate step length
        float dz = vDrift * dt * driftDir;
        float stepLength = fabsf(dz);

        // Add diffusion (Gaussian random walk)
        // Longitudinal diffusion (along drift)
        float sigmaL = diffL * sqrtf(stepLength);
        float randL = curand_normal(&localRand);
        dz += sigmaL * randL;

        // Transverse diffusion (perpendicular to drift)
        float sigmaT = diffT * sqrtf(stepLength);
        float randTx = curand_normal(&localRand);
        float randTy = curand_normal(&localRand);
        float dx = sigmaT * randTx;
        float dy = sigmaT * randTy;

        // Update position
        x += dx;
        y += dy;
        z += dz;
        t += dt;
        result.pathLength += sqrtf(dx*dx + dy*dy + dz*dz);

        // Check boundaries
        float r = sqrtf(x*x + y*y);

        // Absorbed at inner boundary (beampipe)
        if (r < rInner) {
            result.status = 2;  // Absorbed
            break;
        }

        // Out of bounds radially
        if (r > rOuter) {
            result.status = 3;  // Out of bounds
            break;
        }

        // Out of bounds in z
        if (fabsf(z) > halfZ) {
            result.status = 3;  // Out of bounds
            break;
        }

        // Check if reached anode
        if ((driftDir > 0 && z >= anodeZ) || (driftDir < 0 && z <= anodeZ)) {
            result.status = 1;  // Collected!
            break;
        }

        // Avalanche multiplication (if enabled and near anode)
        if (avalanche && alpha > 0) {
            float distToAnode = fabsf(z - anodeZ);
            if (distToAnode < 1.0f) {  // Within 1cm of anode
                // Probability of creating secondary electron
                float pAvalanche = 1.0f - expf(-alpha * stepLength);
                if (curand_uniform(&localRand) < pAvalanche) {
                    result.nAvalanche++;
                }
            }
        }

        // Attachment (electron loss)
        if (eta > 0) {
            float pAttach = 1.0f - expf(-eta * stepLength);
            if (curand_uniform(&localRand) < pAttach) {
                result.status = 2;  // Absorbed by attachment
                break;
            }
        }

        nSteps++;
    }

    // Store final position
    result.finalX = x;
    result.finalY = y;
    result.finalZ = z;
    result.driftTime = t - e.t;

    // Write result
    results[idx] = result;

    // Save random state
    randStates[idx] = localRand;
}

// ============================================================================
// GarfieldGPU Implementation
// ============================================================================

GarfieldGPU::GarfieldGPU() {
    // Default TPC geometry
    m_geometry.innerRadius = 114.0f;   // cm (beampipe)
    m_geometry.outerRadius = 200.0f;   // cm
    m_geometry.halfLength = 300.0f;    // cm
    m_geometry.driftLength = 85.0f;    // cm
    m_geometry.cathodeZ = 0.0f;
    m_geometry.anodeZ = 85.0f;

    // Default gas properties (Ar/CO2 80/20)
    m_gasData.electricField = 250.0f;  // V/cm
    m_gasData.driftVelocity = 0.0055f; // cm/ns (typical for Ar/CO2)
    m_gasData.diffusionL = 0.023f;     // cm/sqrt(cm)
    m_gasData.diffusionT = 0.030f;     // cm/sqrt(cm)
    m_gasData.townsendCoeff = 0.0f;    // No avalanche by default
    m_gasData.attachmentCoeff = 0.0f;  // No attachment
}

GarfieldGPU::~GarfieldGPU() {
    FreeGPUMemory();
}

bool GarfieldGPU::Initialize() {
    // Check for CUDA devices
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);

    if (err != cudaSuccess || deviceCount == 0) {
        std::cerr << "GarfieldGPU: No CUDA devices available" << std::endl;
        m_gpuAvailable = false;
        return false;
    }

    // Get device properties
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, 0);
    m_deviceName = prop.name;

    std::cout << "=============================================" << std::endl;
    std::cout << "GarfieldGPU: GPU Acceleration Initialized" << std::endl;
    std::cout << "  Device: " << m_deviceName << std::endl;
    std::cout << "  Compute Capability: " << prop.major << "." << prop.minor << std::endl;
    std::cout << "  Global Memory: " << (prop.totalGlobalMem / 1024 / 1024) << " MB" << std::endl;
    std::cout << "  SM Count: " << prop.multiProcessorCount << std::endl;
    std::cout << "=============================================" << std::endl;

    m_gpuAvailable = true;
    m_initialized = true;
    return true;
}

void GarfieldGPU::SetGeometry(const TPCGeometry& geom) {
    m_geometry = geom;
}

void GarfieldGPU::SetGasProperties(float arFrac, float co2Frac,
                                    float pressure, float temperature) {
    // Calculate transport properties from gas composition
    float E = m_gasData.electricField;

    // Empirical drift velocity formula for Ar/CO2
    // Based on Magboltz calculations
    float reducedField = E / (pressure / 760.0f);  // V/cm at 1 atm

    // Drift velocity (simplified empirical model)
    // v = a * E / (1 + b * E) with pressure/temperature correction
    float a = 0.00015f + 0.00003f * co2Frac;  // Mobility coefficient
    float b = 0.00001f;
    float v0 = a * reducedField / (1.0f + b * reducedField);

    // Temperature correction
    float tempCorrection = sqrtf(temperature / 293.15f);
    m_gasData.driftVelocity = v0 * tempCorrection;

    // Diffusion coefficients (empirical)
    // Increase with lower E/p
    float diffBase = 0.02f + 0.01f * co2Frac;
    m_gasData.diffusionL = diffBase * sqrtf(293.15f / temperature);
    m_gasData.diffusionT = 1.3f * m_gasData.diffusionL;  // Typically DT > DL

    std::cout << "GarfieldGPU: Gas properties updated" << std::endl;
    std::cout << "  Ar/CO2: " << (arFrac*100) << "/" << (co2Frac*100) << std::endl;
    std::cout << "  Drift velocity: " << m_gasData.driftVelocity << " cm/ns" << std::endl;
    std::cout << "  Diffusion L/T: " << m_gasData.diffusionL << "/" << m_gasData.diffusionT << " cm/sqrt(cm)" << std::endl;
}

void GarfieldGPU::SetElectricField(float field) {
    m_gasData.electricField = field;
    // Recalculate transport with new field
    // (Would need to store gas composition to do this properly)
}

void GarfieldGPU::AddClusters(const std::vector<IonizationCluster>& clusters) {
    for (const auto& cluster : clusters) {
        for (int i = 0; i < cluster.nElectrons; i++) {
            // Add electron with small random offset within cluster
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
    e.status = 0;  // Active
    m_electrons.push_back(e);
}

void GarfieldGPU::ClearElectrons() {
    m_electrons.clear();
    m_results.clear();
    m_nCollected = 0;
    m_totalCharge = 0.0f;
}

void GarfieldGPU::AllocateGPUMemory(int nElectrons) {
    // Free any existing memory
    FreeGPUMemory();

    // Allocate device memory
    cudaMalloc(&m_d_electrons, nElectrons * sizeof(ElectronState));
    cudaMalloc(&m_d_results, nElectrons * sizeof(DriftResult));
    cudaMalloc(&m_d_randStates, nElectrons * sizeof(curandState));

    // Copy constants to device
    cudaMemcpyToSymbol(d_gasData, &m_gasData, sizeof(GasTransportData));
    cudaMemcpyToSymbol(d_geometry, &m_geometry, sizeof(TPCGeometry));
}

void GarfieldGPU::FreeGPUMemory() {
    if (m_d_electrons) { cudaFree(m_d_electrons); m_d_electrons = nullptr; }
    if (m_d_results) { cudaFree(m_d_results); m_d_results = nullptr; }
    if (m_d_randStates) { cudaFree(m_d_randStates); m_d_randStates = nullptr; }
}

void GarfieldGPU::LaunchDriftKernel(int nElectrons) {
    // Configure kernel launch
    int threadsPerBlock = 256;
    int numBlocks = (nElectrons + threadsPerBlock - 1) / threadsPerBlock;

    // Time step for drift (ns) - smaller = more accurate, slower
    float dt = 1.0f;  // 1 ns time step

    // Create CUDA events for timing
    cudaEvent_t start, stop;
    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    // Initialize random states
    cudaEventRecord(start);
    InitRandomStates<<<numBlocks, threadsPerBlock>>>(
        (curandState*)m_d_randStates,
        12345,  // Seed
        nElectrons
    );

    // Launch drift kernel
    DriftKernel<<<numBlocks, threadsPerBlock>>>(
        (ElectronState*)m_d_electrons,
        (DriftResult*)m_d_results,
        (curandState*)m_d_randStates,
        nElectrons,
        dt,
        m_maxDriftTime,
        m_avalancheEnabled
    );

    cudaEventRecord(stop);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&m_kernelTimeMs, start, stop);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}

int GarfieldGPU::DriftElectrons() {
    if (!m_gpuAvailable || m_electrons.empty()) {
        return 0;
    }

    int nElectrons = static_cast<int>(m_electrons.size());

    // Allocate GPU memory
    AllocateGPUMemory(nElectrons);

    // Copy electrons to device
    cudaMemcpy(m_d_electrons, m_electrons.data(),
               nElectrons * sizeof(ElectronState), cudaMemcpyHostToDevice);

    // Launch kernel
    LaunchDriftKernel(nElectrons);

    // Copy results back
    m_results.resize(nElectrons);
    cudaMemcpy(m_results.data(), m_d_results,
               nElectrons * sizeof(DriftResult), cudaMemcpyDeviceToHost);

    // Count collected electrons and total charge
    m_nCollected = 0;
    m_totalCharge = 0.0f;
    for (const auto& r : m_results) {
        if (r.status == 1) {  // Collected
            m_nCollected++;
            m_totalCharge += static_cast<float>(r.nAvalanche);
        }
    }

    std::cout << "GarfieldGPU: Drifted " << nElectrons << " electrons in "
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
    // Empirical formula for Ar/CO2 mixtures
    float reducedField = electricField / (pressure / 760.0f);
    float co2Frac = 1.0f - arFrac;

    float a = 0.00015f + 0.00003f * co2Frac;
    float b = 0.00001f;
    float v0 = a * reducedField / (1.0f + b * reducedField);

    return v0 * sqrtf(temperature / 293.15f);
}

void CalculateDiffusion(float electricField, float pressure,
                        float& diffL, float& diffT, float arFrac) {
    float co2Frac = 1.0f - arFrac;
    float diffBase = 0.02f + 0.01f * co2Frac;

    // Diffusion increases at lower fields
    float fieldFactor = 1.0f + 100.0f / electricField;

    diffL = diffBase * sqrtf(fieldFactor);
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

    // Average ionization energy for Ar (26 eV)
    const float ionizationEnergy = 26.0e-6f;  // MeV

    // Total energy deposited
    float totalEdep = dEdx * pathLength;
    if (totalEdep > energy) totalEdep = energy;

    // Total number of electrons
    int totalElectrons = static_cast<int>(totalEdep / ionizationEnergy);
    if (totalElectrons < 1) totalElectrons = 1;

    // Average cluster size (~3 electrons for minimum ionizing particle)
    float meanClusterSize = 3.0f;
    int nClusters = std::max(1, static_cast<int>(totalElectrons / meanClusterSize));

    // Generate clusters along path
    for (int i = 0; i < nClusters; i++) {
        float frac = (i + 0.5f) / nClusters;

        IonizationCluster cluster;
        cluster.x = x0 + dx * pathLength * frac;
        cluster.y = y0 + dy * pathLength * frac;
        cluster.z = z0 + dz * pathLength * frac;
        cluster.t = 0.0f;
        cluster.nElectrons = totalElectrons / nClusters;
        cluster.energy = totalEdep / nClusters * 1e6f;  // Convert to eV

        clusters.push_back(cluster);
    }

    return clusters;
}

} // namespace nnbar
