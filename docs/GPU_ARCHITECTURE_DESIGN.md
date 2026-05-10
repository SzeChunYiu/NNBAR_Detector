# NNBAR Full GPU Simulation Architecture

## Executive Summary

This document outlines the architecture for a fully GPU-accelerated NNBAR detector simulation. The goal is to move all computationally intensive physics from CPU to GPU while maintaining physics accuracy.

## Current Architecture (Hybrid CPU/GPU)

```
┌─────────────────────────────────────────────────────────────────────┐
│                         GEANT4 (CPU)                                │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐                 │
│  │  Primary    │  │  Hadronic   │  │  Optical    │                 │
│  │  Generator  │  │  Physics    │  │  Photons    │                 │
│  └─────────────┘  └─────────────┘  └─────────────┘                 │
│         │               │               │                           │
│         ▼               ▼               ▼                           │
│  ┌─────────────────────────────────────────────────────────────┐   │
│  │              Geometry Navigation (CPU)                       │   │
│  └─────────────────────────────────────────────────────────────┘   │
│         │                                                           │
│         ▼                                                           │
│  ┌─────────────┐                                                   │
│  │ EM Offload  │ ──────────────────────────────────────────────┐   │
│  └─────────────┘                                               │   │
└─────────────────────────────────────────────────────────────────│───┘
                                                                  │
┌─────────────────────────────────────────────────────────────────│───┐
│                      CELERITAS (GPU)                            │   │
│  ┌─────────────────────────────────────────────────────────────┐│   │
│  │  EM Physics: e-, e+, gamma                                  ││   │
│  │  - Bremsstrahlung, Pair Production                          ││   │
│  │  - Compton Scattering, Photoelectric                        ││   │
│  │  - Multiple Scattering, Ionization                          │◄───┘
│  └─────────────────────────────────────────────────────────────┘│
│                              │                                  │
│                              ▼                                  │
│  ┌─────────────────────────────────────────────────────────────┐│
│  │  GeantSimpleCalo: Energy Deposits → Lead Glass, Scint       ││
│  └─────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────────┘
```

## Proposed Full GPU Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                    CONTROL LAYER (CPU - Minimal)                    │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐                 │
│  │   Event     │  │   I/O &     │  │  Analysis   │                 │
│  │  Dispatch   │  │  Parquet    │  │  Interface  │                 │
│  └──────┬──────┘  └──────┬──────┘  └──────┬──────┘                 │
└─────────┼────────────────┼────────────────┼─────────────────────────┘
          │                │                │
          ▼                ▼                ▼
┌─────────────────────────────────────────────────────────────────────┐
│                         GPU MEMORY (Unified)                        │
│  ┌─────────────────────────────────────────────────────────────────┐│
│  │  Primary Particles Buffer    │  Hit Collection Buffers          ││
│  │  Secondary Stack             │  Energy Deposit Maps             ││
│  │  Track State Arrays          │  Optical Photon Buffer           ││
│  └─────────────────────────────────────────────────────────────────┘│
└─────────────────────────────────────────────────────────────────────┘
          │
          ▼
┌─────────────────────────────────────────────────────────────────────┐
│                    GPU COMPUTE KERNELS                              │
│                                                                     │
│  ┌───────────────────────────────────────────────────────────────┐ │
│  │ LAYER 1: GEOMETRY ENGINE (VecGeom GPU)                        │ │
│  │  - BVH/Octree navigation                                      │ │
│  │  - Distance-to-boundary calculation                           │ │
│  │  - Material lookup                                            │ │
│  └───────────────────────────────────────────────────────────────┘ │
│                              │                                      │
│                              ▼                                      │
│  ┌───────────────────────────────────────────────────────────────┐ │
│  │ LAYER 2: PHYSICS ENGINES                                      │ │
│  │                                                               │ │
│  │  ┌─────────────────┐  ┌─────────────────┐  ┌───────────────┐ │ │
│  │  │ EM PHYSICS      │  │ HADRONIC        │  │ OPTICAL       │ │ │
│  │  │ (Celeritas)     │  │ (AdePT/Custom)  │  │ (Opticks)     │ │ │
│  │  │                 │  │                 │  │               │ │ │
│  │  │ • e-/e+/gamma   │  │ • Proton        │  │ • Scintillation│ │ │
│  │  │ • Brems/Pair    │  │ • Neutron       │  │ • Cherenkov   │ │ │
│  │  │ • Compton/PE    │  │ • Pion/Kaon     │  │ • Reflection  │ │ │
│  │  │ • MSC/Ioniz     │  │ • Nuclear rxns  │  │ • Absorption  │ │ │
│  │  └─────────────────┘  └─────────────────┘  └───────────────┘ │ │
│  └───────────────────────────────────────────────────────────────┘ │
│                              │                                      │
│                              ▼                                      │
│  ┌───────────────────────────────────────────────────────────────┐ │
│  │ LAYER 3: DETECTOR RESPONSE                                    │ │
│  │                                                               │ │
│  │  ┌─────────────────┐  ┌─────────────────┐  ┌───────────────┐ │ │
│  │  │ TPC DRIFT       │  │ SCINTILLATOR    │  │ LEAD GLASS    │ │ │
│  │  │ (CUDA Kernels)  │  │ (GPU Response)  │  │ (GPU Calo)    │ │ │
│  │  │                 │  │                 │  │               │ │ │
│  │  │ • E-field map   │  │ • Light yield   │  │ • Cherenkov   │ │ │
│  │  │ • Drift physics │  │ • Timing        │  │ • PMT response│ │ │
│  │  │ • Diffusion     │  │ • SiPM model    │  │ • Digitization│ │ │
│  │  │ • Amplification │  │ • Digitization  │  │               │ │ │
│  │  └─────────────────┘  └─────────────────┘  └───────────────┘ │ │
│  └───────────────────────────────────────────────────────────────┘ │
│                              │                                      │
│                              ▼                                      │
│  ┌───────────────────────────────────────────────────────────────┐ │
│  │ LAYER 4: HIT AGGREGATION & OUTPUT                             │ │
│  │  - Parallel reduction for energy sums                         │ │
│  │  - GPU-direct Parquet writing (RAPIDS cuDF)                   │ │
│  │  - Event reconstruction primitives                            │ │
│  └───────────────────────────────────────────────────────────────┘ │
└─────────────────────────────────────────────────────────────────────┘
```

## Implementation Phases

### Phase 1: Current (Completed)
- [x] Celeritas EM physics integration
- [x] GeantSimpleCalo for GPU energy recording
- [x] Per-event GPU energy tracking
- [x] Parquet output with GPU energy

### Phase 2: Optical Photon GPU (3-6 months)
**Goal:** Move optical photon simulation to GPU

**Technology:** Opticks or custom CUDA implementation

```cpp
// Proposed interface
namespace nnbar {
namespace gpu {

class OpticalPhotonEngine {
public:
    // Initialize with geometry and material properties
    void Initialize(const G4VPhysicalVolume* world);

    // Generate scintillation photons on GPU
    void GenerateScintillation(
        const thrust::device_vector<float4>& positions,  // x,y,z,edep
        const thrust::device_vector<int>& materials,
        thrust::device_vector<OpticalPhoton>& photons
    );

    // Generate Cherenkov photons on GPU
    void GenerateCherenkov(
        const thrust::device_vector<TrackState>& tracks,
        thrust::device_vector<OpticalPhoton>& photons
    );

    // Propagate photons to PMTs/SiPMs
    void PropagateToDetectors(
        thrust::device_vector<OpticalPhoton>& photons,
        thrust::device_vector<PMTHit>& hits
    );
};

}  // namespace gpu
}  // namespace nnbar
```

### Phase 3: TPC Drift GPU (3-6 months)
**Goal:** Full GPU TPC simulation

**Current:** Garfield++ with optional OpenMP
**Target:** Native CUDA kernels

```cpp
namespace nnbar {
namespace gpu {

class TPCDriftEngine {
public:
    // Load electric field map to GPU
    void LoadFieldMap(const std::string& fieldFile);

    // Process ionization clusters
    struct IonizationCluster {
        float3 position;
        float time;
        int electrons;
        int eventId;
    };

    // Drift electrons from ionization to readout
    __global__ void DriftElectrons(
        const IonizationCluster* clusters,
        int nClusters,
        const float* fieldMapX,
        const float* fieldMapY,
        const float* fieldMapZ,
        float3 fieldDim,
        DriftedElectron* output
    );

    // Apply diffusion
    __global__ void ApplyDiffusion(
        DriftedElectron* electrons,
        int nElectrons,
        float diffusionL,
        float diffusionT,
        curandState* rngStates
    );

    // Collect at readout pads
    __global__ void CollectAtPads(
        const DriftedElectron* electrons,
        int nElectrons,
        const PadGeometry* pads,
        PadSignal* signals
    );
};

}  // namespace gpu
}  // namespace nnbar
```

### Phase 4: Hadronic Physics GPU (12-24 months)
**Goal:** GPU-accelerated hadronic physics

**Challenge:** This is the hardest part - hadronic physics is complex

**Options:**
1. **AdePT** (CERN): Geant4-compatible GPU tracking, experimental hadronic support
2. **Extended Celeritas**: Future versions may add hadronic physics
3. **Custom Implementation**: CUDA kernels for specific processes

```cpp
// Conceptual hadronic GPU interface
namespace nnbar {
namespace gpu {

class HadronicPhysicsEngine {
public:
    // Simplified hadronic models for GPU
    enum class Model {
        BERTINI_SIMPLIFIED,  // Simplified Bertini cascade
        PARAMETRIC,          // Parametric model (fast but less accurate)
        NEURAL_NETWORK       // ML-based surrogate model
    };

    void SetModel(Model m);

    // Process hadronic interactions
    void ProcessInteractions(
        thrust::device_vector<HadronTrack>& tracks,
        thrust::device_vector<SecondaryParticle>& secondaries
    );

    // Neural network surrogate for complex physics
    class NNSurrogate {
        // Trained on Geant4 FTFP_BERT output
        void LoadModel(const std::string& onnxPath);
        void Infer(const HadronState& input, FinalState& output);
    };
};

}  // namespace gpu
}  // namespace nnbar
```

### Phase 5: Full Integration (6-12 months)
**Goal:** Unified GPU event processing

```cpp
namespace nnbar {
namespace gpu {

class UnifiedGPUSimulation {
public:
    // Single entry point for full GPU simulation
    void SimulateEvent(
        const PrimaryParticle* primaries,
        int nPrimaries,
        EventOutput& output
    );

private:
    // All engines on GPU
    std::unique_ptr<GeometryEngine> geometry_;      // VecGeom
    std::unique_ptr<EMPhysicsEngine> emPhysics_;    // Celeritas
    std::unique_ptr<HadronicEngine> hadPhysics_;    // AdePT/Custom
    std::unique_ptr<OpticalEngine> optical_;         // Opticks
    std::unique_ptr<TPCDriftEngine> tpcDrift_;      // Custom CUDA
    std::unique_ptr<DetectorResponse> response_;    // Custom CUDA

    // GPU memory management
    MemoryPool deviceMemory_;

    // Event processing pipeline
    void TransportPrimaries();
    void ProcessEMShowers();
    void ProcessHadronicCascades();
    void GenerateOpticalPhotons();
    void DriftTPCElectrons();
    void CollectDetectorHits();
    void WriteOutput();
};

}  // namespace gpu
}  // namespace nnbar
```

## Performance Estimates

| Component | Current (CPU) | GPU Target | Speedup |
|-----------|--------------|------------|---------|
| EM Physics | ~100 evt/s | ~10,000 evt/s | 100x |
| Hadronic | ~50 evt/s | ~1,000 evt/s | 20x |
| Optical Photons | ~10 evt/s | ~5,000 evt/s | 500x |
| TPC Drift | ~20 evt/s | ~2,000 evt/s | 100x |
| **Total** | ~5 evt/s | ~500 evt/s | **100x** |

## Hardware Requirements

### Minimum (Development)
- NVIDIA GPU: RTX 3080 or better (10+ GB VRAM)
- CUDA 11.0+
- 32 GB RAM

### Recommended (Production)
- NVIDIA GPU: A100 40GB or H100 80GB
- CUDA 12.0+
- 128 GB RAM
- NVMe storage for fast I/O

### HPC Cluster
- Multiple A100/H100 GPUs
- GPU-direct storage (GPUDirect Storage)
- High-bandwidth interconnect (NVLink, InfiniBand)

## Dependencies

### Required
- CUDA Toolkit 11.0+
- Celeritas (current)
- VecGeom (GPU geometry)
- Thrust (GPU algorithms)

### Optional
- Opticks (optical photons)
- AdePT (hadronic physics)
- RAPIDS cuDF (GPU Parquet)
- TensorRT (ML inference)

## Risks and Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| Hadronic GPU not mature | High | Use hybrid CPU/GPU, prioritize EM |
| Memory limitations | Medium | Implement memory pooling, event batching |
| Physics validation | High | Extensive comparison with CPU Geant4 |
| Complex geometry | Medium | Simplify geometry for GPU, use BVH |

## Conclusion

Full GPU simulation is achievable in phases:
1. **Now:** EM physics on GPU (done)
2. **Short-term:** Optical photons, TPC drift
3. **Medium-term:** Simplified hadronic physics
4. **Long-term:** Full hadronic on GPU

The architecture supports incremental migration while maintaining physics accuracy.
