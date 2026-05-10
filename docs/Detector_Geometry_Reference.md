# NNBAR Detector Geometry Reference
## Complete Geometry Documentation

**Version 1.0**

---

## Machine-Checkable Audit

Run the static geometry audit after geometry-source or reference-document
changes:

```bash
python3 -m nnbar_reconstruction.cli geometry-audit . --json output/studies/geometry_audit.json --fail-on-mismatch
```

The audit uses the active CMake source tree (`src/detector/*.cc`) as the
simulation geometry source, verifies that legacy `src/Detector_Module` sources
are not compiled, compares beampipe-5, TPC, scintillator, and lead-glass block
geometry against this reference, and scans compiled detector sources for
uninitialized loop counters.

---

# Table of Contents

1. [Coordinate System](#1-coordinate-system)
2. [Beampipe Geometry](#2-beampipe-geometry)
3. [TPC Geometry](#3-tpc-geometry)
4. [Scintillator Geometry](#4-scintillator-geometry)
5. [Lead Glass Calorimeter](#5-lead-glass-calorimeter)
6. [Complete Cross-Section Views](#6-complete-cross-section-views)

---

# 1. Coordinate System

```
    NNBAR Global Coordinate System

                           +Y (vertical up)
                            |
                            |
                            |
                            |
              +X ___________+___________ -X
             (right)       /|          (left)
                          / |
                         /  |
                        /   |
                      +Z    -Y (down)
                   (beam direction, downstream)


    Origin: Center of detector (center of Beampipe_5)

    Conventions:
    - +Z: Beam direction (toward beam stop)
    - +X: Horizontal right (when looking downstream)
    - +Y: Vertical upward
    - Azimuthal angle: measured from +X toward +Y

    Downstream observer looking along +Z:
    - +X is to the RIGHT
    - +Y is UP
    - -Z is toward observer (upstream)
```

---

# 2. Beampipe Geometry

## 2.1 Beampipe Section 5 (Detector Region)

The detector systems surround Beampipe Section 5.

```
    Beampipe_5 (Side View, XZ plane at y=0)

    z = -2500            z = 0             z = +2500
         |                 |                   |
         v                 v                   v
    +----+=======================================+----+
    |cap |             VACUUM                    |cap |
    |    |                                       |    |
    | 6  |          Beampipe_5                   | 7  |
    |    |         (r = 1120-1140 mm)            |    |
    +----+=======================================+----+

    Beampipe_5 dimensions:
    - Inner radius: 1120 mm
    - Outer radius: 1140 mm (with 20 mm Al wall)
    - Length: 5000 mm (z = -2500 to +2500)
    - Material: Aluminum with B4C coating inside
```

## 2.2 Cross-Section View

```
    Beampipe Cross-Section (XY plane at any z within detector)

                            +Y
                             |
                             |
                  ,----------+----------.
                ,'     B4C coating (1cm)  `.
              ,'                            `.
             /      ,------------------.      \
            /      /                    \      \
           |      |                      |      |  Al wall
           |      |       VACUUM         |      |  (20 mm)
           |      |                      |      |
            \      \                    /      /
             \      `------------------'      /
              `.                            ,'
                `.                        ,'
                  `----------+----------'
                             |
           -X ---------------+---------------- +X
                             |
                            -Y

    Radii:
    - Vacuum: r < 1100 mm
    - B4C coating: 1100-1110 mm
    - Al wall: 1120-1140 mm
```

---

# 3. TPC Geometry

## 3.1 Module Layout Overview

```
    TPC Module Arrangement (Cross-Section at z = 0)
    Looking DOWNSTREAM along +Z

                                    +Y
                                     |
                                     |
                        Module 0 (Type II, top)
                    +--------------------+  y = +1567
                    |     DRIFT -Y       |  854 mm tall
                    |    2284 x 854 mm   |  2284 mm wide
                    +--------------------+
                             |
    Module 1                 |                 Module 5
    (Type I)                 |                 (Type I)
    +-------+                |                 +-------+
    |       |                |                 |       |
    |DRIFT  |         ___________              |  DRIFT|
    | +X    |        /           \             |   -X  |
    |854x   |       |   BEAMPIPE  |            |   x854|
    |1994   |       |   r=1140mm  |            |  1994 |
    |       |        \___________/             |       |
    +-------+                |                 +-------+
    x=-1569                  |                  x=+1569
    y=+997                   |                  y=+997
                             |
    Module 2                 |                 Module 4
    (Type I)                 |                 (Type I)
    +-------+                |                 +-------+
    |       |                |                 |       |
    |DRIFT  |                |                 |  DRIFT|
    | +X    |                |                 |   -X  |
    |854x   |                |                 |   x854|
    |1994   |                |                 |  1994 |
    |       |                |                 |       |
    +-------+                |                 +-------+
    x=-1569                  |                  x=+1569
    y=-997                   |                  y=-997
                    +--------------------+
                    |    2284 x 854 mm   |
                    |     DRIFT +Y       |
                    +--------------------+  y = -1567
                        Module 3 (Type II, bottom)
                             |
           -X ---------------+---------------- +X


    KEY:
    - Type I (vertical): 854 mm wide x 1994 mm tall, drift along X
    - Type II (horizontal): 2284 mm wide x 854 mm tall, drift along Y
    - All modules extend 2520 mm in Z direction
    - Front modules (0-5) at z = -1260 mm
    - Back modules (6-11) at z = +1260 mm
```

## 3.2 Module Dimensions Table

| Module Type | Width (X) | Height (Y) | Length (Z) | Drift Dir |
|-------------|-----------|------------|------------|-----------|
| Type I | 854 mm | 1994 mm | 2520 mm | X |
| Type II | 2284 mm | 854 mm | 2520 mm | Y |

## 3.3 Module Positions

| Module | Type | X (mm) | Y (mm) | Z (mm) | Location |
|--------|------|--------|--------|--------|----------|
| 0 | II | 0 | +1567 | -1260 | Front Top |
| 1 | I | -1569 | +997 | -1260 | Front Top-Left |
| 2 | I | -1569 | -997 | -1260 | Front Bot-Left |
| 3 | II | 0 | -1567 | -1260 | Front Bottom |
| 4 | I | +1569 | -997 | -1260 | Front Bot-Right |
| 5 | I | +1569 | +997 | -1260 | Front Top-Right |
| 6 | II | 0 | +1567 | +1260 | Back Top |
| 7 | I | +1569 | +997 | +1260 | Back Top-Right |
| 8 | I | +1569 | -997 | +1260 | Back Bot-Right |
| 9 | II | 0 | -1567 | +1260 | Back Bottom |
| 10 | I | -1569 | -997 | +1260 | Back Bot-Left |
| 11 | I | -1569 | +997 | +1260 | Back Top-Left |

---

# 4. Scintillator Geometry

## 4.1 Overview

The scintillator system surrounds the TPC and consists of:
- **4 side surfaces** (top, bottom, left, right)
- **2 end surfaces** (front, back)

## 4.2 Key Dimensions

```
    Scintillator Reference Distances:

    scint_base_dist = Beampipe_5_radius_2 + TPC_drift_len + 2*TPC_wall_thickness
                    = 1140 mm + 850 mm + 4 mm
                    = 1994 mm (inner surface of scintillator layer)

    dy = 200 mm (gap between TPC and scintillator)

    Scintillator module surface position:
    = scint_base_dist + dy + scint_module_y/2
    = 1994 + 200 + 150 = 2344 mm (center of scintillator)
```

## 4.3 Side Scintillator Modules

```
    Side Scintillator Module Structure

    Module dimensions: 400 mm x 300 mm x 400 mm (X x Y x Z)

    Internal structure (10 alternating layers of 30 mm each):

    +--------------------------------------------+
    | Layer 9: Vertical bars (along X)           |
    +--------------------------------------------+
    | Layer 8: Horizontal bars (along Z)         |
    +--------------------------------------------+
    | Layer 7: Vertical bars                     |
    +--------------------------------------------+
    | Layer 6: Horizontal bars                   |
    +--------------------------------------------+
    | Layer 5: Vertical bars                     |
    +--------------------------------------------+
    | Layer 4: Horizontal bars                   |
    +--------------------------------------------+
    | Layer 3: Vertical bars                     |
    +--------------------------------------------+
    | Layer 2: Horizontal bars                   |
    +--------------------------------------------+
    | Layer 1: Vertical bars                     |
    +--------------------------------------------+
    | Layer 0: Horizontal bars                   |
    +--------------------------------------------+

    Each layer contains 4 scintillator bars
    - Horizontal bar: 100 mm x 30 mm x 400 mm
    - Vertical bar: 400 mm x 30 mm x 100 mm
```

## 4.4 Side Surface Arrangement

```
    Scintillator Side Modules (Cross-section at z=0)
    Looking DOWNSTREAM along +Z

                                 +Y
                                  |
                      TOP SURFACE (10 x 11 modules)
        +---+---+---+---+---+---+---+---+---+---+
        |   |   |   |   |   |   |   |   |   |   | y ~ +2344
        +---+---+---+---+---+---+---+---+---+---+

    +---+                                     +---+
    |   |                                     |   |
    +---+                                     +---+
    |   |      +---------------------+        |   |
    +---+      |                     |        +---+
    |   |      |        TPC          |        |   |
    +---+      |                     |        +---+
LEFT|   |      |      BEAMPIPE       |        |   | RIGHT
SIDE+---+      |                     |        +---+ SIDE
    |   |      |                     |        |   |
    +---+      +---------------------+        +---+
    |   |                                     |   |
    +---+                                     +---+
    |   |                                     |   |
    +---+                                     +---+
x~-2344                                       x~+2344

        +---+---+---+---+---+---+---+---+---+---+
        |   |   |   |   |   |   |   |   |   |   | y ~ -2344
        +---+---+---+---+---+---+---+---+---+---+
                      BOTTOM SURFACE
                                  |
          -X ---------------------+-------------------- +X


    Module count per surface:
    - Top: 10 x 11 = 110 modules
    - Bottom: 10 x 11 = 110 modules
    - Left: 10 x 11 = 110 modules (rotated 90 deg)
    - Right: 10 x 11 = 110 modules (rotated 90 deg)

    Total side modules: 440
```

## 4.5 Front/Back Scintillator Modules

```
    Front/Back Scintillator Module

    Module dimensions: 300 mm x 300 mm x 300 mm (X x Y x Z)

    Internal structure (10 layers, alternating orientation):
    - Horizontal bars: 300 mm x 50 mm x 30 mm (along X)
    - Vertical bars: 50 mm x 300 mm x 30 mm (along Y)

    z-position:
    - Front: z = -(Beampipe_5_len/2 + dz_fb + scint_module_fb_z/2)
           = -(2500 + 200 + 150) = -2850 mm
    - Back: z = +2850 mm
```

## 4.6 Front/Back Surface Arrangement

```
    Front Scintillator Surface (looking downstream, z = -2850 mm)

                                 +Y
                                  |
        Group 1 (top): 15 x 4 = 60 modules
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+

    +---+                                             +---+
    +---+                                             +---+
    +---+                  BEAMPIPE                   +---+
    +---+               (opening area)                +---+  Group 2 (sides)
    +---+                    ___                      +---+  4 x 7 = 28 each
    +---+                   /   \                     +---+
    +---+                  |     |                    +---+
    +---+                   \___/                     +---+

    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
    |   |   |   |   |   |   |   |   |   |   |   |   |   |   |   |
    +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+
        Group 1 (bottom): 15 x 4 = 60 modules
                                  |
          -X ---------------------+-------------------- +X

    Per front/back surface:
    - Group 1 (top + bottom): 2 x 60 = 120 modules
    - Group 2 (left + right): 2 x 28 = 56 modules
    - Total: 176 modules per surface

    Front + Back total: 352 modules
```

## 4.7 Scintillator Module Summary

| Surface | Rotation | Modules | Position |
|---------|----------|---------|----------|
| Top | 0 deg | 110 | y ~ +2344 mm |
| Left | 270 deg | 110 | x ~ -2344 mm |
| Bottom | 180 deg | 110 | y ~ -2344 mm |
| Right | 90 deg | 110 | x ~ +2344 mm |
| Front | rotX(0) | 176 | z ~ -2850 mm |
| Back | rotX(180) | 176 | z ~ +2850 mm |
| **Total** | | **792** | |

---

# 5. Lead Glass Calorimeter

## 5.1 Lead Glass Block Dimensions

```
    Single Lead Glass Block

    +--------+
    |        |
    |  80mm  |   PMT at top (+Y end)
    |  x     |<-- 250 mm length (Y direction)
    |  80mm  |
    |        |
    +--------+

    Dimensions: 80 mm x 250 mm x 80 mm (X x Y x Z)
    Material: Schott SF5 lead glass (density 6.22 g/cm3)

    Block structure:
    - Lead glass body: 80 x 250 x 80 mm
    - AlMgF2 coating on sides (0.01 mm)
    - PMT window at top (0.01 mm)
```

## 5.2 Lead Glass Placement

Lead glass positions are read from CSV files:
- `lead_glass_position.csv` - Side surface positions
- `lead_glass_position_fb.csv` - Front/back positions

The active Geant4 source reads these files from `data/lead_glass_position/`.
The current placement inputs are:

| Placement input | Active path | Rows | Surfaces using each row | Placed blocks |
|-----------------|-------------|------|--------------------------|---------------|
| Side surfaces | `data/lead_glass_position/lead_glass_position.csv` | 3233 | 4 | 12932 |
| Front/back surfaces | `data/lead_glass_position/lead_glass_position_fb.csv` | 2520 | 2 | 5040 |
| **Total placed blocks** | | | | **17972** |

The active importer clears prior rows before each load and raises a Geant4
`FatalException` if either placement file is missing or empty. This prevents a
run from silently constructing a detector with no lead-glass blocks.

## 5.3 Side Surface Arrangement

```
    Lead Glass Side Surfaces (Cross-section at z=0)
    Looking DOWNSTREAM along +Z

    The lead glass blocks are arranged on 4 surfaces around the
    scintillator layer, with the long axis (Y = 250 mm) pointing
    RADIALLY OUTWARD.

                                 +Y
                                  |
                      TOP SURFACE
         [][][][][][][][][][][][][][][][][][][][]
         [][][][][][][][][][][][][][][][][][][][]
         [][][][][][][][][][][][][][][][][][][][]
                    (radial out = +Y)
                             |
                             |
    [][]                                     [][]
    [][]                                     [][]
    [][]     +---------------------+         [][]
    [][]     |                     |         [][]
    [][]     |        TPC +        |         [][]
LEFT[][]     |     SCINTILLATOR    |         [][]RIGHT
SIDE[][]     |                     |         [][]SIDE
    [][]     |      BEAMPIPE       |         [][]
    [][]     |                     |         [][]
    [][]     +---------------------+         [][]
    [][]                                     [][]
    [][]                                     [][]
    (radial out = -X)               (radial out = +X)

         [][][][][][][][][][][][][][][][][][][][]
         [][][][][][][][][][][][][][][][][][][][]
         [][][][][][][][][][][][][][][][][][][][]
                      BOTTOM SURFACE
                    (radial out = -Y)
                                  |
          -X ---------------------+-------------------- +X

    Block orientation on each surface:
    - Top: PMT facing +Y (radially outward)
    - Bottom: PMT facing -Y (radially outward)
    - Left: PMT facing -X (radially outward)
    - Right: PMT facing +X (radially outward)

    Note: Blocks are rotated 90 deg between adjacent surfaces
```

## 5.4 Front/Back Surface Arrangement

```
    Front Lead Glass Surface (looking downstream, z < 0)

    Blocks on front/back surfaces have their long axis (250 mm)
    pointing along Z (beam direction), with PMT facing away from detector.

                                 +Y
                                  |
         [][][][][][][][][][][][][][][][][][][][]
         [][][][][][][][][][][][][][][][][][][][]
                    (PMT facing -Z)

    [][]                                     [][]
    [][]                                     [][]
    [][]                                     [][]
    [][]             BEAMPIPE                [][]
    [][]            (opening)                [][]
    [][]               ___                   [][]
    [][]              /   \                  [][]
    [][]             |     |                 [][]
    [][]              \___/                  [][]
    [][]                                     [][]
    [][]                                     [][]

         [][][][][][][][][][][][][][][][][][][][]
         [][][][][][][][][][][][][][][][][][][][]
                    (PMT facing -Z)
                                  |
          -X ---------------------+-------------------- +X

    Front blocks: PMT faces -Z (upstream, away from detector)
    Back blocks: PMT faces +Z (downstream, away from detector)
```

## 5.5 Lead Glass Index Assignment

The lead glass blocks are indexed as follows:

| Index Range | Surface | Notes |
|-------------|---------|-------|
| 0 - N1 | Top (i=0) | First surface |
| N1 - 2*N1 | Right (i=1) | 90 deg rotation |
| 2*N1 - 3*N1 | Bottom (i=2) | 180 deg rotation |
| 3*N1 - 4*N1 | Left (i=3) | 270 deg rotation |
| 4*N1 - 4*N1+N2 | Front | rotX(-90) |
| 4*N1+N2 - 4*N1+2*N2 | Back | rotX(+90) |

Where:
- N1 = number of blocks per side surface (from CSV)
- N2 = number of blocks per front/back surface (from CSV)

---

# 6. Complete Cross-Section Views

## 6.1 XY Cross-Section (at z = 0)

```
    Complete Detector Cross-Section (XY plane at z=0)
    Looking DOWNSTREAM along +Z axis
    (Not to scale - schematic representation)

                                          +Y
                                           |
    OUTER RADIUS                           |
    (~2600 mm)                             |
                    LEAD GLASS (outer)     |
            [][][][][][][][][][][][][][][][]
            [][][][][][][][][][][][][][][][]
                                           |
            +---+---+---+---+---+---+---+---+  SCINTILLATOR
            |   |   |   |   |   |   |   |   |  (top surface)
            +---+---+---+---+---+---+---+---+
                                           |
        []  +---------------------------+  []
        []  |       TPC Type II         |  []
        []  |      (Module 0)           |  []
        []  +---------------------------+  []
    []  []  +----+             +----+  []  []
    []  []  |    |   BEAMPIPE  |    |  []  []
LEAD[]  []  |TPC |   (vacuum)  |TPC |  []  []LEAD
GLASS   []  |    |     ___     |    |  []  GLASS
    []SCINT |Type|    /   \    |Type|SCINT[]
    []  []  | I  |   |     |   | I  |  []  []
    []  []  |(1,2|   |  O  |   |4,5)|  []  []
    []  []  |    |    \___/    |    |  []  []
    []  []  +----+             +----+  []  []
        []                             []
        []  +---------------------------+  []
        []  |       TPC Type II         |  []
        []  |      (Module 3)           |  []
        []  +---------------------------+  []
                                           |
            +---+---+---+---+---+---+---+---+
            |   |   |   |   |   |   |   |   |  SCINTILLATOR
            +---+---+---+---+---+---+---+---+  (bottom surface)
                                           |
            [][][][][][][][][][][][][][][][]
            [][][][][][][][][][][][][][][][]
                    LEAD GLASS (outer)     |
                                           |
    INNER RADIUS                           |
    (~1140 mm at beampipe)                 |
                                           |
          -X --------------------------+------------------------ +X

    Approximate radii from center:
    - Beampipe: 1120-1140 mm
    - TPC inner: 1140 mm
    - TPC outer: ~2000 mm
    - Scintillator: ~2000-2300 mm
    - Lead Glass: ~2300-2550 mm
```

## 6.2 XZ Cross-Section (at y = 0)

```
    Detector Side View (XZ plane at y=0)
    Looking from +Y (bird's eye view)

    z = -2850         z = -1260         z = +1260         z = +2850
        |                 |                 |                 |
        v                 v                 v                 v

    +-------+=========================================+-------+
    | FRONT |                                         | BACK  |
    | SCINT |          SIDE SCINTILLATOR              | SCINT |
    +-------+=========================================+-------+
    +-------+=========================================+-------+
    | FRONT |                                         | BACK  |
    | LG    |          SIDE LEAD GLASS                | LG    |
    +-------+=========================================+-------+

        +------+                             +------+
        | TPC  |                             | TPC  |
        |Front |                             |Back  |
        |0-5   |       BEAMPIPE              |6-11  |
        |      |      (z = -2500 to +2500)   |      |
        +------+                             +------+

    +-------+=========================================+-------+
    | FRONT |                                         | BACK  |
    | LG    |          SIDE LEAD GLASS                | LG    |
    +-------+=========================================+-------+
    +-------+=========================================+-------+
    | FRONT |                                         | BACK  |
    | SCINT |          SIDE SCINTILLATOR              | SCINT |
    +-------+=========================================+-------+

    -X -------------------- O -------------------- +X

    Z-positions:
    - Beampipe_5: z = -2500 to +2500 mm
    - TPC front modules: z = -1260 mm (center)
    - TPC back modules: z = +1260 mm (center)
    - Front scintillator: z ~ -2850 mm
    - Back scintillator: z ~ +2850 mm
```

## 6.3 YZ Cross-Section (at x = 0)

```
    Detector Front View (YZ plane at x=0)
    Looking from -X toward +X

    z = -2850         z = -1260         z = +1260         z = +2850
        |                 |                 |                 |
        v                 v                 v                 v
                                                          +Y
                                                           |
    +-------+         +------+         +------+         +-------+
    | FRONT |         | TPC  |         | TPC  |         | BACK  |
    | END   |         |Mod 0 |         |Mod 6 |         | END   |
    | CAP   |         +------+         +------+         | CAP   |
    |       |                                           |       |
    |       |              BEAMPIPE                     |       |
    |       |             (r=1140mm)                    |       |
    |       |               ___                         |       |
    |       |              /   \                        |       |
    |       |             |     |                       |       |
    |       |              \___/                        |       |
    |       |                                           |       |
    |       |         +------+         +------+         |       |
    |       |         | TPC  |         | TPC  |         |       |
    |       |         |Mod 3 |         |Mod 9 |         |       |
    +-------+         +------+         +------+         +-------+
                                                           |
          z <--------------------------------------------- +

    Note: Scintillator and Lead Glass layers wrap around the
    entire structure (not shown for clarity)
```

---

**Document Version History**

| Version | Date | Author | Changes |
|---------|------|--------|---------|
| 1.0 | 2026-01-09 | Claude | Initial release with accurate geometry |

---

*Generated for NNBAR Detector Simulation Project*
