# Offline Reconstruction

This package reconstructs analysis-level objects from the simulation Parquet
outputs without coupling the analysis to the Geant4 executable.

## Inputs

The reader expects files named like:

- `TPC_output_<run>.parquet`
- `Scintillator_output_<run>.parquet`
- `LeadGlass_output_<run>.parquet`
- `PMT_output_<run>.parquet`
- optional supporting tables such as `Particle`, `Silicon`, and `Interaction`

macOS `._*.parquet` sidecar files are ignored.

## Implemented Objects

- charged-track candidates from TPC hits
- electron/positron-pair candidates from the Chapter 8.2 TPC entry-point
  separation rule (`<=5 cm`)
- event vertices from the Chapter 7 TPC track-projection method onto the
  source-foil plane (`z=0`)
- preliminary charged pion/proton PID from TPC dE/dx plus scintillator range
- photon-like neutral objects from lead-glass energy deposits, with shower
  centroids and vertex-relative directions when a reconstructed event vertex is
  available
- pi0 candidates from photon pairs using those reconstructed photon directions
- event variables: calorimeter energy, object multiplicities, pion multiplicity,
  PMT photon counts, visible invariant mass, upper/lower scintillator and
  lead-glass energy, signed longitudinal calorimeter energy, transverse
  calorimeter energy, Chapter 7 in-time/out-of-time calorimeter energy, and
  sphericity
- thesis-threshold preliminary event selection and cumulative cut-flow columns

The pi0 defaults follow the thesis selection values available in Chapter 8:
mass in `100-180 MeV`, total energy `<=720 MeV`, scintillator energy
`<=250 MeV`, lead-glass energy `<=980 MeV`, lead-glass energy fraction `>=0.55`,
and opening angle `>=30 deg`.

The preliminary event selection defaults follow the visible thresholds in
Chapter 9: scintillator energy in `20-2000 MeV`, a TPC track with foil-like
origin provenance, pion multiplicity `>=1`, visible invariant mass `>=500 MeV`,
sphericity `>=0.2`, upper scintillator energy `<=320 MeV`, and lower
scintillator energy `<=930 MeV`. The visible invariant mass is reconstructed
from object directions and deposited/visible energies, so it is a first analysis
surface rather than a calibrated final mass estimator.

Charged-object directions are reconstructed from the ordered TPC hit positions
when those coordinates are available. The simulation momentum columns are used
only as a fallback for sparse tables without enough positional information.

The longitudinal and transverse energy variables follow the Chapter 9 event-shape
definitions using the detector z-axis as the beam axis: `EL = sum(E_i cos(alpha_i))`
and `ET = sum(E_i sin(alpha_i))`, evaluated separately for scintillator and
lead-glass hits and then combined for the full calorimeter.

The vertex table reports the mean of valid per-track projections to `z=0`, the
RMS radial spread of those projections, and the number of skipped tracks whose
TPC entry/exit points were insufficient or parallel to the foil plane.

Photon-like objects report the energy-weighted lead-glass shower centroid, the
vertex used for the neutral flight direction, and the corresponding path length.
When a reconstructed event vertex exists, the photon unit vector is computed from
`vertex -> shower centroid`; otherwise the package falls back to the historical
origin-to-centroid direction. This keeps pi0 masses and event-shape variables
consistent with the reconstructed interaction point instead of implicitly
assuming a source at `(0,0,0)`.

Charged-neutral separation for lead-glass clusters is geometric. TPC hits are
grouped into reconstructed track candidates, each candidate direction is taken
from the event vertex toward the farthest TPC hit, and a lead-glass cluster is
tagged as charged only when its vertex-to-centroid direction falls inside the
configurable matching cone. The simulation `Track_ID` remains available as
provenance in `source_track_id`, but it is not the charged-cluster decision rule.

The timing-window columns use the Chapter 7 formulas after a vertex and event
time are available. Scintillator hits are accepted between the flight time of a
`1000 MeV` charged pion minus `2 sigma_scint` and a `100 MeV` charged pion plus
`2 sigma_scint`; lead-glass hits are accepted within the photon flight time plus
or minus `2 sigma_lg`. The resolution values are configurable because detector
calibration constants are not fixed in this offline package.

## CLI

```bash
python3 -m nnbar_reconstruction.cli summarize build-codex-setup2/output --run 0
python3 -m nnbar_reconstruction.cli summarize build-codex-setup2/output --run 0 --tables-dir reconstruction_out
```

The CSV tables are intended as the handoff surface for plotting, cut-flow
studies, and future calibration-derived PID thresholds.
