// ============================================================================
// TPCDriftGPU.cc
// GPU-accelerated TPC electron drift simulation with CPU fallback
// ============================================================================

#include "gpu/TPCDriftGPU.hh"
#include "gpu/GPUManager.hh"
#include "G4ios.hh"
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>
#include <fstream>

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

__global__ void InitRNGKernel(curandState* states, int nStates, unsigned long seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < nStates) {
        curand_init(seed, idx, 0, &states[idx]);
    }
}

__global__ void DriftElectronsKernel(
    const IonizationCluster* clusters,
    int nClusters,
    DriftedElectron* electrons,
    int* electronCount,
    float driftVelocity,
    float diffusionL,
    float diffusionT,
    float driftLength,
    curandState* rngStates
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nClusters) return;

    curandState localState = rngStates[idx];
    const IonizationCluster& cluster = clusters[idx];

    // Drift distance
    float driftDist = driftLength - cluster.z;  // Assuming readout at z=driftLength
    if (driftDist < 0) return;

    // Drift time
    float driftTime = driftDist / driftVelocity;

    // Diffusion sigmas
    float sigmaL = diffusionL * sqrtf(driftDist / 10.0f);  // Convert to sqrt(cm)
    float sigmaT = diffusionT * sqrtf(driftDist / 10.0f);

    // Process each electron in cluster
    int baseIdx = atomicAdd(electronCount, cluster.electrons);
    for (int i = 0; i < cluster.electrons && (baseIdx + i) < 10000000; i++) {
        DriftedElectron& e = electrons[baseIdx + i];

        // Add diffusion
        float dx = sigmaT * curand_normal(&localState);
        float dy = sigmaT * curand_normal(&localState);
        float dz = sigmaL * curand_normal(&localState);

        e.x = cluster.x + dx;
        e.y = cluster.y + dy;
        e.z = driftLength + dz;  // At readout plane
        e.time = cluster.time + driftTime + dz / driftVelocity;
        e.clusterId = idx;
        e.eventId = cluster.eventId;
    }

    rngStates[idx] = localState;
}

__global__ void CollectOnPadsKernel(
    const DriftedElectron* electrons,
    int nElectrons,
    TPCPadHit* hits,
    int* hitCount,
    float padSizeX,
    float padSizeY,
    int nPadsX,
    int nPadsY
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= nElectrons) return;

    const DriftedElectron& e = electrons[idx];

    // Find pad
    int padX = (int)((e.x + nPadsX * padSizeX / 2) / padSizeX);
    int padY = (int)((e.y + nPadsY * padSizeY / 2) / padSizeY);

    if (padX >= 0 && padX < nPadsX && padY >= 0 && padY < nPadsY) {
        int padId = padY * nPadsX + padX;
        int hitIdx = atomicAdd(hitCount, 1);
        if (hitIdx < 1000000) {
            hits[hitIdx].padId = padId;
            hits[hitIdx].moduleId = 0;  // Would need geometry info
            hits[hitIdx].layerId = 0;
            hits[hitIdx].charge = 1.0f;
            hits[hitIdx].time = e.time;
            hits[hitIdx].eventId = e.eventId;
        }
    }
}

#endif  // HAS_CUDA

// ============================================================================
// Singleton Instance
// ============================================================================

TPCDriftGPU& TPCDriftGPU::Instance() {
    static TPCDriftGPU instance;
    return instance;
}

TPCDriftGPU::~TPCDriftGPU() {
    Shutdown();
}

// ============================================================================
// Initialization
// ============================================================================

bool TPCDriftGPU::Initialize() {
    std::lock_guard<std::mutex> lock(mutex_);

    if (isInitialized_) return isGPUEnabled_;

    G4cout << "[TPCDriftGPU] Initializing..." << G4endl;

#if HAS_CUDA
    auto& gpuMgr = GPUManager::Instance();
    if (gpuMgr.IsGPUAvailable() && gpuMgr.IsGPUTPCDriftEnabled()) {
        cudaError_t err;
        bool allOk = true;

        err = cudaMalloc(&d_clusters_, sizeof(IonizationCluster) * 1000000);
        if (err != cudaSuccess) { allOk = false; }

        if (allOk) {
            err = cudaMalloc(&d_electrons_, sizeof(DriftedElectron) * maxElectrons_);
            if (err != cudaSuccess) { allOk = false; }
        }

        if (allOk) {
            err = cudaMalloc(&d_padHits_, sizeof(TPCPadHit) * 1000000);
            if (err != cudaSuccess) { allOk = false; }
        }

        if (allOk) {
            err = cudaMalloc(&d_rngStates_, sizeof(curandState) * 1000000);
            if (err != cudaSuccess) { allOk = false; }
        }

        if (allOk) {
            // Initialize RNG states
            int blockSize = 256;
            int numBlocks = (1000000 + blockSize - 1) / blockSize;
            InitRNGKernel<<<numBlocks, blockSize>>>((curandState*)d_rngStates_, 1000000, 42);
            cudaDeviceSynchronize();

            isGPUEnabled_ = true;
            G4cout << "[TPCDriftGPU] GPU mode ENABLED" << G4endl;
            G4cout << "[TPCDriftGPU] Max electrons: " << maxElectrons_ << G4endl;
        } else {
            // Cleanup on failure
            if (d_clusters_) cudaFree(d_clusters_);
            if (d_electrons_) cudaFree(d_electrons_);
            if (d_padHits_) cudaFree(d_padHits_);
            if (d_rngStates_) cudaFree(d_rngStates_);
            d_clusters_ = d_electrons_ = d_padHits_ = d_rngStates_ = nullptr;
            isGPUEnabled_ = false;
        }
    } else {
        isGPUEnabled_ = false;
    }
#else
    isGPUEnabled_ = false;
#endif

    if (!isGPUEnabled_) {
        G4cout << "[TPCDriftGPU] CPU fallback mode ENABLED" << G4endl;
    }

    isInitialized_ = true;
    return isGPUEnabled_;
}

void TPCDriftGPU::Shutdown() {
    std::lock_guard<std::mutex> lock(mutex_);

    if (!isInitialized_) return;

#if HAS_CUDA
    if (d_clusters_) { cudaFree(d_clusters_); d_clusters_ = nullptr; }
    if (d_electrons_) { cudaFree(d_electrons_); d_electrons_ = nullptr; }
    if (d_padHits_) { cudaFree(d_padHits_); d_padHits_ = nullptr; }
    if (d_rngStates_) { cudaFree(d_rngStates_); d_rngStates_ = nullptr; }
    if (d_fieldMap_) { cudaFree(d_fieldMap_); d_fieldMap_ = nullptr; }
#endif

    isInitialized_ = false;
    isGPUEnabled_ = false;
}

// ============================================================================
// Field Map
// ============================================================================

bool TPCDriftGPU::LoadFieldMap(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        G4cerr << "[TPCDriftGPU] Failed to open field map: " << filename << G4endl;
        return false;
    }

    // Read header
    file >> fieldMapNx_ >> fieldMapNy_ >> fieldMapNz_;
    file >> fieldMapDx_ >> fieldMapDy_ >> fieldMapDz_;

    int totalPoints = fieldMapNx_ * fieldMapNy_ * fieldMapNz_;
    fieldMap_.resize(totalPoints);

    for (int i = 0; i < totalPoints; i++) {
        file >> fieldMap_[i].Ex >> fieldMap_[i].Ey >> fieldMap_[i].Ez;
    }

    hasFieldMap_ = true;

#if HAS_CUDA
    if (isGPUEnabled_) {
        if (d_fieldMap_) cudaFree(d_fieldMap_);
        cudaMalloc(&d_fieldMap_, sizeof(TPCFieldPoint) * totalPoints);
        cudaMemcpy(d_fieldMap_, fieldMap_.data(), sizeof(TPCFieldPoint) * totalPoints,
                   cudaMemcpyHostToDevice);
    }
#endif

    G4cout << "[TPCDriftGPU] Loaded field map: " << fieldMapNx_ << "x"
           << fieldMapNy_ << "x" << fieldMapNz_ << " points" << G4endl;

    return true;
}

void TPCDriftGPU::SetUniformField(float Ex, float Ey, float Ez) {
    hasFieldMap_ = false;
    electricField_ = std::sqrt(Ex*Ex + Ey*Ey + Ez*Ez);
}

// ============================================================================
// Processing
// ============================================================================

void TPCDriftGPU::AddIonizationCluster(const IonizationCluster& cluster) {
    std::lock_guard<std::mutex> lock(mutex_);
    clusters_.push_back(cluster);
    totalInputElectrons_ += cluster.electrons;
}

void TPCDriftGPU::ProcessEvent() {
    auto start = std::chrono::high_resolution_clock::now();

    nClusters_ = clusters_.size();

    if (nClusters_ == 0) {
        processingTimeMs_ = 0;
        return;
    }

    if (isGPUEnabled_) {
        ProcessEventGPU();
        eventsProcessedGPU_++;
    } else {
        ProcessEventCPU();
        eventsProcessedCPU_++;
    }

    auto end = std::chrono::high_resolution_clock::now();
    processingTimeMs_ = std::chrono::duration<float, std::milli>(end - start).count();

    auto& gpuMgr = GPUManager::Instance();
    if (isGPUEnabled_) {
        gpuMgr.GetStatistics().AddGPUTime(processingTimeMs_);
        gpuMgr.GetStatistics().gpuTPCElectrons.fetch_add(totalCollectedElectrons_);
    } else {
        gpuMgr.GetStatistics().AddCPUTime(processingTimeMs_);
    }
}

void TPCDriftGPU::ProcessEventGPU() {
#if HAS_CUDA
    int nClusters = clusters_.size();
    if (nClusters == 0) return;

    // Copy clusters to GPU
    cudaMemcpy(d_clusters_, clusters_.data(), sizeof(IonizationCluster) * nClusters,
               cudaMemcpyHostToDevice);

    // Allocate electron count on GPU
    int* d_electronCount;
    cudaMalloc(&d_electronCount, sizeof(int));
    cudaMemset(d_electronCount, 0, sizeof(int));

    // Launch drift kernel
    int blockSize = 256;
    int numBlocks = (nClusters + blockSize - 1) / blockSize;
    DriftElectronsKernel<<<numBlocks, blockSize>>>(
        (IonizationCluster*)d_clusters_,
        nClusters,
        (DriftedElectron*)d_electrons_,
        d_electronCount,
        driftVelocity_,
        diffusionL_,
        diffusionT_,
        driftLength_,
        (curandState*)d_rngStates_
    );

    // Get electron count
    int nElectrons;
    cudaMemcpy(&nElectrons, d_electronCount, sizeof(int), cudaMemcpyDeviceToHost);
    totalCollectedElectrons_ = nElectrons;

    // Allocate hit count on GPU
    int* d_hitCount;
    cudaMalloc(&d_hitCount, sizeof(int));
    cudaMemset(d_hitCount, 0, sizeof(int));

    // Launch pad collection kernel
    numBlocks = (nElectrons + blockSize - 1) / blockSize;
    CollectOnPadsKernel<<<numBlocks, blockSize>>>(
        (DriftedElectron*)d_electrons_,
        nElectrons,
        (TPCPadHit*)d_padHits_,
        d_hitCount,
        5.0f, 5.0f,  // Pad size (mm)
        100, 100     // Number of pads
    );

    // Get hit count and copy hits back
    int nHits;
    cudaMemcpy(&nHits, d_hitCount, sizeof(int), cudaMemcpyDeviceToHost);
    padHits_.resize(nHits);
    cudaMemcpy(padHits_.data(), d_padHits_, sizeof(TPCPadHit) * nHits, cudaMemcpyDeviceToHost);

    cudaFree(d_electronCount);
    cudaFree(d_hitCount);
#else
    ProcessEventCPU();
#endif
}

void TPCDriftGPU::ProcessEventCPU() {
    std::lock_guard<std::mutex> lock(mutex_);

    driftedElectrons_.clear();

    // Drift electrons from each cluster
    int clusterId = 0;
    for (const auto& cluster : clusters_) {
        for (int i = 0; i < cluster.electrons; i++) {
            DriftedElectron e = DriftElectron(
                cluster.x, cluster.y, cluster.z, cluster.time,
                clusterId, cluster.eventId
            );
            driftedElectrons_.push_back(e);
        }
        clusterId++;
    }

    totalCollectedElectrons_ = driftedElectrons_.size();

    // Collect on pads
    CollectOnPads(driftedElectrons_);
}

DriftedElectron TPCDriftGPU::DriftElectron(float x, float y, float z, float t,
                                            int clusterId, int eventId) {
    static std::mt19937 rng(42);
    std::normal_distribution<float> gauss(0.0f, 1.0f);

    DriftedElectron e;

    // Drift distance (assuming readout at z = driftLength_)
    float driftDist = driftLength_ - z;
    if (driftDist < 0) driftDist = 0;

    // Drift time
    float driftTime = driftDist / driftVelocity_;

    // Diffusion
    float sigmaL = diffusionL_ * std::sqrt(driftDist / 10.0f);  // sqrt(cm) to sqrt(mm)
    float sigmaT = diffusionT_ * std::sqrt(driftDist / 10.0f);

    e.x = x + sigmaT * gauss(rng);
    e.y = y + sigmaT * gauss(rng);
    e.z = driftLength_ + sigmaL * gauss(rng);
    e.time = t + driftTime;
    e.clusterId = clusterId;
    e.eventId = eventId;

    return e;
}

void TPCDriftGPU::CollectOnPads(const std::vector<DriftedElectron>& electrons) {
    padHits_.clear();

    const float padSizeX = 5.0f;  // mm
    const float padSizeY = 5.0f;  // mm
    const int nPadsX = 100;
    const int nPadsY = 100;

    for (const auto& e : electrons) {
        int padX = (int)((e.x + nPadsX * padSizeX / 2) / padSizeX);
        int padY = (int)((e.y + nPadsY * padSizeY / 2) / padSizeY);

        if (padX >= 0 && padX < nPadsX && padY >= 0 && padY < nPadsY) {
            TPCPadHit hit;
            hit.padId = padY * nPadsX + padX;
            hit.moduleId = 0;
            hit.layerId = 0;
            hit.charge = gasGain_;  // Apply gas gain
            hit.time = e.time;
            hit.eventId = e.eventId;
            padHits_.push_back(hit);
        }
    }
}

void TPCDriftGPU::ClearEvent() {
    std::lock_guard<std::mutex> lock(mutex_);
    clusters_.clear();
    driftedElectrons_.clear();
    padHits_.clear();
    nClusters_ = 0;
    totalInputElectrons_ = 0;
    totalCollectedElectrons_ = 0;
}

void TPCDriftGPU::PrintStatistics() const {
    G4cout << "\n[TPCDriftGPU] Statistics:" << G4endl;
    G4cout << "  Events processed (GPU): " << eventsProcessedGPU_ << G4endl;
    G4cout << "  Events processed (CPU): " << eventsProcessedCPU_ << G4endl;
    G4cout << "  Mode: " << (isGPUEnabled_ ? "GPU" : "CPU fallback") << G4endl;
}

}  // namespace gpu
}  // namespace nnbar
