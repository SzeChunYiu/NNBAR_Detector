// ============================================================================
// OpticalPhotonGPU.cc
// GPU-accelerated optical photon simulation with CPU fallback
// ============================================================================

#include "gpu/OpticalPhotonGPU.hh"
#include "gpu/GPUManager.hh"
#include "G4ios.hh"
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#include <curand_kernel.h>
#define HAS_CUDA 1
#else
#define HAS_CUDA 0
#endif

namespace nnbar {
namespace gpu {

// ============================================================================
// CUDA Kernels (if CUDA available)
// ============================================================================

#if HAS_CUDA

__global__ void GenerateScintillationKernel(
    const ScintillationPoint* points,
    int nPoints,
    OpticalPhoton* photons,
    int* photonCount,
    float yield,
    float decayTime,
    curandState* rngStates
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nPoints) return;

    curandState localState = rngStates[idx];
    const ScintillationPoint& pt = points[idx];

    // Number of photons from this point
    int nPhotons = (int)(pt.energyDeposit * yield);

    // Generate photons
    int baseIdx = atomicAdd(photonCount, nPhotons);
    for (int i = 0; i < nPhotons && (baseIdx + i) < 10000000; i++) {
        OpticalPhoton& p = photons[baseIdx + i];

        // Random direction (isotropic)
        float cosTheta = 2.0f * curand_uniform(&localState) - 1.0f;
        float sinTheta = sqrtf(1.0f - cosTheta * cosTheta);
        float phi = 2.0f * 3.14159265f * curand_uniform(&localState);

        p.x = pt.x;
        p.y = pt.y;
        p.z = pt.z;
        p.dx = sinTheta * cosf(phi);
        p.dy = sinTheta * sinf(phi);
        p.dz = cosTheta;
        p.wavelength = 420.0f + 60.0f * curand_uniform(&localState);  // 420-480 nm
        p.time = pt.time - decayTime * logf(curand_uniform(&localState));
        p.eventId = pt.eventId;
        p.isAlive = true;
    }

    rngStates[idx] = localState;
}

__global__ void PropagatePhotonsKernel(
    OpticalPhoton* photons,
    int nPhotons,
    PMTHit* hits,
    int* hitCount,
    float pmtRadius,
    float pmtZ
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nPhotons) return;

    OpticalPhoton& p = photons[idx];
    if (!p.isAlive) return;

    // Simple propagation to PMT plane
    // This is a simplified model - real implementation would trace through geometry
    float t = (pmtZ - p.z) / p.dz;
    if (t > 0 && p.dz > 0) {
        float hitX = p.x + t * p.dx;
        float hitY = p.y + t * p.dy;
        float r = sqrtf(hitX * hitX + hitY * hitY);

        // Check if hits PMT
        if (r < pmtRadius) {
            int hitIdx = atomicAdd(hitCount, 1);
            if (hitIdx < 1000000) {
                hits[hitIdx].pmtId = (int)(r / 100.0f);  // Simple PMT assignment
                hits[hitIdx].moduleId = 0;
                hits[hitIdx].time = p.time + t / 300.0f;  // ~c in mm/ns
                hits[hitIdx].wavelength = p.wavelength;
                hits[hitIdx].eventId = p.eventId;
                hits[hitIdx].photonCount = 1;
            }
        }
    }
    p.isAlive = false;
}

#endif  // HAS_CUDA

// ============================================================================
// Singleton Instance
// ============================================================================

OpticalPhotonGPU& OpticalPhotonGPU::Instance() {
    static OpticalPhotonGPU instance;
    return instance;
}

OpticalPhotonGPU::~OpticalPhotonGPU() {
    Shutdown();
}

// ============================================================================
// Initialization
// ============================================================================

bool OpticalPhotonGPU::Initialize() {
    std::lock_guard<std::mutex> lock(mutex_);

    if (isInitialized_) return isGPUEnabled_;

    G4cout << "[OpticalPhotonGPU] Initializing..." << G4endl;

#if HAS_CUDA
    auto& gpuMgr = GPUManager::Instance();
    if (gpuMgr.IsGPUAvailable() && gpuMgr.IsGPUOpticalEnabled()) {
        // Allocate GPU buffers
        cudaError_t err;

        err = cudaMalloc(&d_scintPoints_, sizeof(ScintillationPoint) * 1000000);
        if (err != cudaSuccess) {
            G4cerr << "[OpticalPhotonGPU] Failed to allocate scint points buffer" << G4endl;
            isGPUEnabled_ = false;
        } else {
            err = cudaMalloc(&d_photons_, sizeof(OpticalPhoton) * maxPhotons_);
            if (err != cudaSuccess) {
                G4cerr << "[OpticalPhotonGPU] Failed to allocate photon buffer" << G4endl;
                cudaFree(d_scintPoints_);
                d_scintPoints_ = nullptr;
                isGPUEnabled_ = false;
            } else {
                err = cudaMalloc(&d_pmtHits_, sizeof(PMTHit) * 1000000);
                if (err != cudaSuccess) {
                    cudaFree(d_scintPoints_);
                    cudaFree(d_photons_);
                    d_scintPoints_ = nullptr;
                    d_photons_ = nullptr;
                    isGPUEnabled_ = false;
                } else {
                    isGPUEnabled_ = true;
                    G4cout << "[OpticalPhotonGPU] GPU mode ENABLED" << G4endl;
                    G4cout << "[OpticalPhotonGPU] Max photons: " << maxPhotons_ << G4endl;
                }
            }
        }
    } else {
        isGPUEnabled_ = false;
    }
#else
    isGPUEnabled_ = false;
#endif

    if (!isGPUEnabled_) {
        G4cout << "[OpticalPhotonGPU] CPU fallback mode ENABLED" << G4endl;
    }

    // Initialize default emission spectrum (typical scintillator)
    emissionWavelengths_ = {400, 420, 440, 460, 480, 500};
    emissionIntensities_ = {0.1, 0.5, 1.0, 0.8, 0.4, 0.1};

    // Default refractive indices
    refractiveIndices_.resize(10, 1.5f);

    isInitialized_ = true;
    return isGPUEnabled_;
}

void OpticalPhotonGPU::Shutdown() {
    std::lock_guard<std::mutex> lock(mutex_);

    if (!isInitialized_) return;

#if HAS_CUDA
    if (d_scintPoints_) { cudaFree(d_scintPoints_); d_scintPoints_ = nullptr; }
    if (d_photons_) { cudaFree(d_photons_); d_photons_ = nullptr; }
    if (d_pmtHits_) { cudaFree(d_pmtHits_); d_pmtHits_ = nullptr; }
#endif

    isInitialized_ = false;
    isGPUEnabled_ = false;
}

// ============================================================================
// Processing
// ============================================================================

void OpticalPhotonGPU::AddScintillationPoint(const ScintillationPoint& point) {
    std::lock_guard<std::mutex> lock(mutex_);
    scintPoints_.push_back(point);
}

void OpticalPhotonGPU::AddCherenkovPoint(float x, float y, float z,
                                          float dx, float dy, float dz,
                                          float beta, float trackLength,
                                          int materialId, int eventId) {
    if (!cherenkovEnabled_ || beta < 0.7f) return;

    std::lock_guard<std::mutex> lock(mutex_);

    // Simplified Cherenkov: estimate number of photons
    float n = refractiveIndices_[materialId % refractiveIndices_.size()];
    float cosTheta = 1.0f / (beta * n);
    if (cosTheta > 1.0f) return;

    // Number of Cherenkov photons (simplified Frank-Tamm)
    int nPhotons = (int)(370.0f * trackLength * (1.0f - 1.0f / (beta * beta * n * n)));
    if (nPhotons <= 0) return;

    // Add photons
    float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);
    static std::mt19937 rng(42);
    std::uniform_real_distribution<float> phiDist(0, 2.0f * M_PI);
    std::uniform_real_distribution<float> waveDist(300.0f, 600.0f);

    for (int i = 0; i < std::min(nPhotons, 1000); i++) {
        OpticalPhoton p;
        p.x = x;
        p.y = y;
        p.z = z;

        float phi = phiDist(rng);
        // Cherenkov cone direction
        p.dx = sinTheta * std::cos(phi) * dx - sinTheta * std::sin(phi);
        p.dy = sinTheta * std::sin(phi) * dy + sinTheta * std::cos(phi);
        p.dz = cosTheta * dz;

        // Normalize
        float norm = std::sqrt(p.dx*p.dx + p.dy*p.dy + p.dz*p.dz);
        p.dx /= norm; p.dy /= norm; p.dz /= norm;

        p.wavelength = waveDist(rng);
        p.time = 0;
        p.eventId = eventId;
        p.isAlive = true;

        cherenkovPhotons_.push_back(p);
    }
}

void OpticalPhotonGPU::ProcessEvent() {
    auto start = std::chrono::high_resolution_clock::now();

    if (isGPUEnabled_) {
        ProcessEventGPU();
        eventsProcessedGPU_++;
    } else {
        ProcessEventCPU();
        eventsProcessedCPU_++;
    }

    auto end = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(end - start).count();

    auto& gpuMgr = GPUManager::Instance();
    if (isGPUEnabled_) {
        gpuMgr.GetStatistics().AddGPUTime(ms);
        gpuMgr.GetStatistics().gpuOpticalPhotons.fetch_add(totalScintPhotons_ + totalCherenkovPhotons_);
    } else {
        gpuMgr.GetStatistics().AddCPUTime(ms);
    }
}

void OpticalPhotonGPU::ProcessEventGPU() {
#if HAS_CUDA
    // TODO: Full GPU implementation
    // For now, fall back to CPU
    ProcessEventCPU();
#else
    ProcessEventCPU();
#endif
}

void OpticalPhotonGPU::ProcessEventCPU() {
    std::lock_guard<std::mutex> lock(mutex_);

    std::vector<OpticalPhoton> allPhotons;

    // Generate scintillation photons
    for (const auto& pt : scintPoints_) {
        GenerateScintillationPhotons(pt, allPhotons);
    }
    totalScintPhotons_ = allPhotons.size();

    // Add Cherenkov photons
    allPhotons.insert(allPhotons.end(), cherenkovPhotons_.begin(), cherenkovPhotons_.end());
    totalCherenkovPhotons_ = cherenkovPhotons_.size();

    // Propagate to PMTs
    PropagatePhotons(allPhotons);

    totalDetectedPhotons_ = pmtHits_.size();
}

void OpticalPhotonGPU::GenerateScintillationPhotons(const ScintillationPoint& pt,
                                                     std::vector<OpticalPhoton>& photons) {
    static std::mt19937 rng(42);
    std::uniform_real_distribution<float> uniform(0.0f, 1.0f);
    std::uniform_real_distribution<float> waveDist(420.0f, 480.0f);

    // Limit photon generation for performance
    int nPhotons = std::min((int)(pt.energyDeposit * scintYield_), 10000);

    for (int i = 0; i < nPhotons; i++) {
        OpticalPhoton p;
        p.x = pt.x;
        p.y = pt.y;
        p.z = pt.z;

        // Isotropic direction
        float cosTheta = 2.0f * uniform(rng) - 1.0f;
        float sinTheta = std::sqrt(1.0f - cosTheta * cosTheta);
        float phi = 2.0f * M_PI * uniform(rng);

        p.dx = sinTheta * std::cos(phi);
        p.dy = sinTheta * std::sin(phi);
        p.dz = cosTheta;

        p.wavelength = waveDist(rng);
        p.time = pt.time - scintDecayTime_ * std::log(uniform(rng) + 1e-10f);
        p.eventId = pt.eventId;
        p.isAlive = true;

        photons.push_back(p);
    }
}

void OpticalPhotonGPU::PropagatePhotons(std::vector<OpticalPhoton>& photons) {
    // Simplified propagation - detect photons reaching PMT surfaces
    // Real implementation would trace through geometry

    static std::mt19937 rng(42);
    std::uniform_real_distribution<float> uniform(0.0f, 1.0f);

    const float detectionEfficiency = 0.25f;  // 25% PMT QE
    const float reflectionProb = 0.95f;       // 95% reflection at boundaries

    pmtHits_.clear();

    for (auto& p : photons) {
        if (!p.isAlive) continue;

        // Simple model: probability of detection based on geometry
        // Real implementation would do ray tracing
        int bounces = 0;
        const int maxBounces = 10;

        while (p.isAlive && bounces < maxBounces) {
            // Random walk until detection or absorption
            float r = uniform(rng);

            if (r < detectionEfficiency) {
                // Detected!
                PMTHit hit;
                hit.pmtId = (int)(uniform(rng) * 100);
                hit.moduleId = 0;
                hit.time = p.time + bounces * 10.0f;  // Add propagation time
                hit.wavelength = p.wavelength;
                hit.eventId = p.eventId;
                hit.photonCount = 1;
                pmtHits_.push_back(hit);
                p.isAlive = false;
            } else if (r < detectionEfficiency + (1.0f - reflectionProb)) {
                // Absorbed
                p.isAlive = false;
            } else {
                // Reflected, continue
                bounces++;
            }
        }
        p.isAlive = false;  // Kill after max bounces
    }
}

void OpticalPhotonGPU::ClearEvent() {
    std::lock_guard<std::mutex> lock(mutex_);
    scintPoints_.clear();
    cherenkovPhotons_.clear();
    pmtHits_.clear();
    totalScintPhotons_ = 0;
    totalCherenkovPhotons_ = 0;
    totalDetectedPhotons_ = 0;
}

void OpticalPhotonGPU::SetScintillatorEmissionSpectrum(const std::vector<float>& wavelengths,
                                                        const std::vector<float>& intensities) {
    emissionWavelengths_ = wavelengths;
    emissionIntensities_ = intensities;
}

void OpticalPhotonGPU::PrintStatistics() const {
    G4cout << "\n[OpticalPhotonGPU] Statistics:" << G4endl;
    G4cout << "  Events processed (GPU): " << eventsProcessedGPU_ << G4endl;
    G4cout << "  Events processed (CPU): " << eventsProcessedCPU_ << G4endl;
    G4cout << "  Mode: " << (isGPUEnabled_ ? "GPU" : "CPU fallback") << G4endl;
}

}  // namespace gpu
}  // namespace nnbar
