# NNBAR Detector Calibration and Particle Identification Plan

## 1. Physics Context

### 1.1 Neutron-Antineutron Annihilation Products
- **Primary products**: ~5 pions per annihilation (π⁺, π⁻, π⁰)
- **π⁰ decay**: π⁰ → γγ (immediate, τ ~ 10⁻¹⁶ s)
- **Secondary products**: Protons, neutrons from nuclear breakup

### 1.2 Energy Spectra (from simulation/thesis)
| Particle | Peak Energy | Range |
|----------|-------------|-------|
| Charged pions (π±) | ~100 MeV | 50-600 MeV |
| Photons (γ from π⁰) | ~300 MeV | 50-800 MeV |
| Protons | ~200 MeV | 50-500 MeV |

---

## 2. Detector Subsystems and Calibration

### 2.1 Scintillator Calorimeter
**Material**: BC-408 plastic scintillator
**Light yield**: 10,000 photons/MeV
**Attenuation length**: 210 cm

#### Calibration Approach
1. **MIP Reference**: 300 MeV π⁺ (minimum ionizing)
   - Establishes baseline: ~2 MeV/cm × 10,000 photons/MeV = 20,000 photons/cm

2. **Energy Scan**: 50-600 MeV pions
   - Maps non-linear response due to Birks' law quenching
   - Birks' formula: dL/dx = S × (dE/dx) / (1 + kB × dE/dx)

3. **π⁻ vs π⁺ Comparison**
   - π⁻ nuclear capture at rest produces different signature
   - Important for charge-symmetric reconstruction

#### Key Calibration Files
- `calib_pion_mip.mac` - MIP reference
- `calib_pion_energy_scan.mac` - Energy response
- `calib_pion_minus.mac` - π⁻ comparison

### 2.2 Lead Glass Calorimeter
**Material**: Schott SF5 lead glass
**Cerenkov yield**: ~200 photons/MeV
**Radiation length**: X₀ = 1.76 cm
**Molière radius**: 3.5 cm

#### Calibration Approach
1. **Energy Scan**: 50-800 MeV photons
   - Maps Cerenkov yield vs true photon energy
   - Expected resolution: σ/E = 5%/√E ⊕ 2%

2. **Shower Containment**:
   - Lead glass is 25 cm deep = 14.2 X₀
   - 99% EM shower containment for E < 1 GeV

3. **Electron Validation**
   - Cross-check using electrons (same EM shower physics)
   - Validates pair-production vs direct Cerenkov

#### Key Calibration Files
- `calib_gamma_energy_scan.mac` - Photon response
- `calib_gamma_all_surfaces.mac` - Uniformity check
- `calib_electron_validation.mac` - Cross-validation

---

## 3. Particle Identification Strategy

### 3.1 Available Discriminating Variables

#### From TPC (Time Projection Chamber)
| Variable | Description | PID Power |
|----------|-------------|-----------|
| dE/dx | Energy loss per unit length | High |
| Track length | Total path before stopping/interaction | Medium |
| Track curvature | Momentum from magnetic field | High |
| Topology | Kinks, secondary vertices | High |

#### From Scintillator
| Variable | Description | PID Power |
|----------|-------------|-----------|
| Total light | Integrated scintillation | Medium |
| Hit pattern | Spatial distribution | Medium |
| Timing | Time-of-flight | High |

#### From Lead Glass
| Variable | Description | PID Power |
|----------|-------------|-----------|
| Shower shape | Longitudinal/lateral profile | High |
| E/p ratio | Energy vs momentum | High |

### 3.2 dE/dx Particle Identification

The Bethe-Bloch formula predicts:
```
dE/dx ∝ (Z²/β²) × [ln(2meγ²β²/I) - β²]
```

**Key insight**: At same momentum p, different masses give different β:
- β = p / √(p² + m²)

**Expected dE/dx bands** (at p = 200 MeV/c):
| Particle | Mass (MeV) | β | dE/dx (MeV/cm) |
|----------|-----------|-----|----------------|
| π± | 140 | 0.82 | ~2.5 |
| μ± | 106 | 0.88 | ~2.3 |
| p | 938 | 0.21 | ~15 |

**Separation power**:
- Protons are clearly separated from pions (factor ~6 in dE/dx)
- Pions vs muons require additional information (similar mass)

### 3.3 Pion vs Muon Discrimination

Since π± and μ± have similar masses, dE/dx alone cannot distinguish them.

**Additional discriminators:**

1. **Decay topology**: π⁺ → μ⁺ + νμ (τ = 26 ns)
   - Look for kinks in TPC tracks
   - Muon appears with different momentum direction

2. **Range**: At same energy, μ travels further (no hadronic interactions)
   - Pions can undergo nuclear interactions (π⁺ + n → p + π⁰)

3. **Scintillator light yield**:
   - Different quenching factors for π vs μ
   - Muons are cleaner (no nuclear breakup)

4. **Time structure**:
   - Pions from annihilation are prompt (t ~ 0)
   - Muons from pion decay appear delayed

### 3.4 Proton Identification

Protons are easy to identify from:
1. **High dE/dx**: Factor 5-10 higher than pions at same momentum
2. **Short range**: Protons stop quickly due to high dE/dx
3. **No decay**: Stable particle, straight tracks
4. **Annihilation signature**: Often accompanied by nuclear fragments

### 3.5 Photon Identification (in Lead Glass)

1. **No TPC track**: Photons don't ionize until pair production
2. **EM shower shape**: Narrow, well-contained
3. **Timing**: Speed of light (c), no mass delay
4. **π⁰ mass reconstruction**: Two photons with m_γγ ~ 135 MeV

---

## 4. Recommended Calibration Procedure

### Phase 1: Scintillator Calibration
```bash
./nnbar-calo-sim -m macro/calibration/scintillator/calib_pion_mip.mac
./nnbar-calo-sim -m macro/calibration/scintillator/calib_pion_energy_scan.mac
./nnbar-calo-sim -m macro/calibration/scintillator/calib_pion_minus.mac
```

### Phase 2: Lead Glass Calibration
```bash
./nnbar-calo-sim -m macro/calibration/leadglass/calib_gamma_energy_scan.mac
./nnbar-calo-sim -m macro/calibration/leadglass/calib_gamma_all_surfaces.mac
./nnbar-calo-sim -m macro/calibration/leadglass/calib_electron_validation.mac
```

### Phase 3: Run All
```bash
./nnbar-calo-sim -m macro/calibration/run_all_calibrations.mac
```

---

## 5. Expected Outputs

### 5.1 Calibration Constants
| Subsystem | Parameter | Expected Value |
|-----------|-----------|----------------|
| Scintillator | Light yield | 10,000 ± 500 ph/MeV |
| Scintillator | Birks' constant | 0.126 mm/MeV |
| Lead Glass | Cerenkov yield | 200 ± 20 ph/MeV |
| Lead Glass | Resolution | 5%/√E ⊕ 2% |

### 5.2 PID Performance Targets
| Separation | Target Efficiency | Target Purity |
|------------|-------------------|---------------|
| π± vs p | > 99% | > 98% |
| π± vs μ± | > 90% | > 85% |
| γ vs charged | > 99.9% | > 99% |

---

## 6. Additional Ideas for Improvement

### 6.1 Machine Learning PID
- Train neural network on combined TPC + calorimeter features
- Use simulation for training, validate with cosmic data
- Potential variables: dE/dx profile, shower shape, timing

### 6.2 Template Fitting
- Build probability density functions for each particle type
- Use likelihood ratio for classification
- Can combine multiple observables optimally

### 6.3 Topological Selection
- Vertex reconstruction quality
- Track multiplicity constraints
- Angular correlations (back-to-back for annihilation)

### 6.4 Energy Reconstruction Cross-Check
- Compare TPC momentum with calorimeter energy
- Helps identify misidentified particles
- E/p ratio is powerful discriminator

---

## 7. Implementation Notes

### 7.1 Code Changes Made
- Enhanced `PrimaryGeneratorAction` with calibration mode
- New Geant4 commands under `/calibration/` namespace
- Support for different particle types and energy ranges

### 7.2 Running Calibrations
```bash
cd /home/billy/nnbar/simulation/NNBAR_Detector/build
./nnbar-calo-sim -m macro/calibration/run_all_calibrations.mac
```

### 7.3 Analysis
After running calibrations, use the reconstruction framework:
```python
from nnbar_reconstruction.calibration import scintillator_calibration, leadglass_calibration
# Fit calibration curves from parquet output files
```
