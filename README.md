# vargrid

Variational gridding of radar polar volumes to Cartesian grids.

vargrid reads ODIM HDF5 radar volumes and produces CF-compliant NetCDF files on a projected Cartesian grid at multiple altitude levels. It supports two gridding methods: traditional CAPPI with inverse distance weighting, and a variational approach inspired by Brook et al. (2022) that minimises a regularised cost function balancing data fit against spatial smoothness.

## Features

- **Variational gridding** — second-order smoothness penalty with azimuthally-varying weights (Brook et al. 2022 Eq. 2-3), plus first-order Laplacian damping to suppress biharmonic ringing. Polak-Ribière conjugate gradient with Armijo backtracking line search. Cost-only evaluation in line search for performance.
- **Multi-field gridding** — automatically discovers all radar moments in the input volume (DBZH, VRADH, ZDR, RHOHV, etc.) and grids them all. Configurable via `include_fields` / `exclude_fields`.
- **Perona-Malik edge preservation** — optional first-order edge-preserving diffusion (activated via `vargrid_kappa > 0`) as an alternative to second-order smoothness.
- **O(1) gate-to-grid mapping** — radar gates are mapped to grid cells using precomputed bearing/range lookup tables (720 bearing bins × range bins at 250m resolution). No PROJ transforms needed for the observation operator.
- **Int16 bit-packing** — output data packed as CF-compliant `int16` with `scale_factor`/`add_offset` for ~50% file size reduction. Fixed physical ranges per known ODIM quantity.
- **CAPPI baseline** — traditional constant-altitude plan position indicator with IDW interpolation between sweeps.
- **CF-compliant output** — NetCDF output with CF-1.10 attributes, grid mapping, and per-field metadata mapped from ODIM quantity names.
- **Multi-threaded** — altitude layers are processed in parallel.

## Dependencies

- **bom-util** and **bom-core** — Bureau of Meteorology C++ utility and radar libraries
- **HDF5** — for reading ODIM polar volume files
- **NetCDF** — for writing gridded output
- **PROJ** — for map projection (grid geometry only; not used in gate-to-grid mapping)
- **C++20** compiler (GCC ≥ 11, Clang ≥ 14)

## Building

```bash
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/opt/bom
make
```

## Usage

```bash
vargrid [options] config.conf input.pvol.h5 output.nc
```

Generate a sample configuration file:

```bash
vargrid --generate > my_config.conf
```

### Options

| Flag | Description |
|---|---|
| `-h, --help` | Show usage and exit |
| `-g, --generate` | Print sample configuration to stdout |
| `-t, --trace=LEVEL` | Set log level: `none`, `status`, `error`, `warning`, `log`, `debug` |

## Configuration

### Grid geometry

```
proj4 "+proj=aea +lat_1=-32.2 +lat_2=-35.2 +lon_0=151.209 +lat_0=-33.7008 +a=6378137 +b=6356752.31414 +units=m"
size "301 301"
left_top "-150500 150500"
cell_delta "1000 -1000"
units m
origin xy
```

### Altitude layers

```
altitude_base 0.0
altitude_step 500.0
layer_count 13
```

### Field selection

By default, all fields found in the input volume are gridded. To restrict:

```
include_fields "DBZH VRADH ZDR"    # only grid these
exclude_fields "SQI CCOR"          # grid everything except these
```

`include_fields` takes precedence if both are set.

### Gridding method

```
gridding_method variational    # or "cappi"
```

### Variational parameters

| Parameter | Default | Description |
|---|---|---|
| `vargrid_lambda_h` | 0.01 | Horizontal smoothing weight. Higher = smoother. |
| `vargrid_max_alt_diff` | 2000 | Max altitude difference (m) for gate selection |
| `vargrid_max_iterations` | 200 | Maximum conjugate gradient iterations |
| `vargrid_tolerance` | 1e-5 | Convergence tolerance (relative gradient norm) |
| `vargrid_beam_power` | 2.0 | Exponent for range-dependent observation weight: w = (r₀/r)^p |
| `vargrid_ref_range` | 10000 | Reference range (m) where beam weight = 1.0 |
| `vargrid_min_weight` | 0.01 | Minimum observation weight floor |
| `vargrid_use_nearest_init` | true | Initialise solver with weighted nearest-gate values |
| `vargrid_kappa` | 0 | Perona-Malik edge threshold (data units). 0 = second-order smoothness |
| `vargrid_range_spacing` | 250 | Radar range gate spacing (m) for azimuthal weights. 0 = isotropic |
| `vargrid_mask_distance` | 3 | Max extrapolation distance from observations (grid cells) |

### Output options

| Parameter | Default | Description |
|---|---|---|
| `pack_output` | true | Pack data as int16 with CF scale_factor/add_offset (~50% smaller) |
| `output_obs_count` | false | Write `nobs_<FIELD>` diagnostic variables |
| `velocity` | VRADH | Velocity moment name (controls undetect/nodata handling) |

### CAPPI parameters

| Parameter | Default | Description |
|---|---|---|
| `max_alt_dist` | 20000 | Max altitude distance (m) for sweep selection |
| `idw_pwr` | 2.0 | Inverse distance weighting exponent |

## How it works

### Variational gridding

The variational method finds the gridded field **φ** that minimises (Brook et al. 2022):

**J(φ) = Jₒ + λH · (Js + Jd)**

where:
- **Jₒ = Σ wᵢ (φ[gᵢ] - dᵢ)²** — data fidelity. Observation weights model beam volume broadening: w = (r₀/r)^p × (1 - |Δz|/Δz_max).
- **Js = Σ [Wx·φ_xx² + Wy·φ_yy²]** — second-order smoothness (Eq. 2). Uses centered second differences with Neumann boundary conditions.
- **Jd = Σ_edges (φ_j - φ_k)²** — first-order Laplacian damping. Suppresses oscillatory ringing inherent to the biharmonic (fourth-order) operator produced by second-order smoothness.

**Azimuthal weights** (Eq. 3): Wx = C - A·cos(2φ), Wy = C + A·cos(2φ), where φ is the bearing from the radar, f = Δr/(beamwidth × range), A = |f-1|/2, C = (f+1)/2. This applies stronger smoothing along the azimuthal direction (where data spacing is coarser at long range). Weights are isotropic within 20km of the radar and precomputed once per field-layer.

Alternatively, **Perona-Malik edge-preserving diffusion** is available via `vargrid_kappa > 0`: Js = Σ κ² log(1 + d²/κ²) where d is the difference between adjacent cells. This preserves strong gradients while smoothing weak noise.

The solver uses **Polak-Ribière nonlinear conjugate gradient** with Armijo backtracking line search. The line search uses a cost-only evaluation (~3× cheaper than full gradient computation) for efficiency. Cells beyond `vargrid_mask_distance` grid steps from any observation are masked as NaN using BFS distance propagation.

### Gate-to-grid mapping

Radar gates are mapped to the Cartesian grid using a precomputed bearing/range lookup table. For each grid cell, the bearing and ground range from the radar are computed using `wgs84.latlon_to_bearing_range` (same geodetic function as the CAPPI method). These are binned into a 2D table (720 bearing bins × range bins at 250m resolution) for O(1) nearest-cell lookup per gate. No PROJ transforms are needed — the approach is coordinate-system agnostic.

### CAPPI

The CAPPI method finds the two nearest sweeps above and below the target altitude for each grid cell and interpolates using inverse distance weighting. Sweeps are pre-sorted by elevation for early termination.

## References

Brook, J. P., Protat, A., Soderholm, J., Dee, D., & Gregory, J. (2022). A variational interpolation method for gridding weather radar data. *Journal of Atmospheric and Oceanic Technology*, 39(12), 1853–1871.

## Repository structure

```
vargrid/
├── CMakeLists.txt
├── README.md
├── src/
│   ├── main.cc                        — CLI, config parsing, orchestration
│   ├── types.h                        — Data structures (sweep, volume) and BOM includes
│   ├── util.h                         — phase_timer, split_fields (header-only)
│   ├── reader.h / reader.cc           — ODIM reading, field discovery, metadata extraction
│   ├── writer.h / writer.cc           — CF-compliant NetCDF output with int16 packing
│   ├── grid.h / grid.cc              — Altitude init, flip/flipud utilities
│   └── gridding/
│       ├── config.h                   — vargrid_config parameter struct
│       ├── observation_operator.h/cc  — Bearing/range lookup, gate-to-grid mapping, Wx/Wy precomputation
│       ├── cost_function.h            — J(x), ∇J(x), evaluate_cost (header-only)
│       ├── solver.h                   — Polak-Ribière CG with cost-only line search (header-only)
│       ├── variational.h/cc           — Variational gridding driver, BFS masking
│       └── cappi.h/cc                 — CAPPI IDW gridding
└── test/
    ├── test_harness.h                 — Minimal test framework
    ├── test_main.cc                   — Test runner
    ├── test_laplacian.cc              — Laplacian operator verification
    ├── test_cost_function.cc          — Gradient vs finite differences
    ├── test_solver.cc                 — CG convergence tests
    └── test_field_selection.cc        — Include/exclude logic
```

## License

Copyright Commonwealth of Australia, Bureau of Meteorology. See LICENSE for details.