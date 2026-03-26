# vargrid

Variational gridding of radar polar volumes to Cartesian grids.

vargrid reads ODIM HDF5 radar volumes and produces CF-compliant NetCDF files on a projected Cartesian grid at multiple altitude levels. It supports two gridding methods: traditional CAPPI with inverse distance weighting, and a variational approach that minimises a regularised cost function balancing data fit against spatial smoothness.

## Features

- **Multi-field gridding** — automatically discovers all radar moments in the input volume (DBZH, VRADH, ZDR, RHOHV, etc.) and grids them all. Field selection is configurable via `include_fields` / `exclude_fields`.
- **Variational gridding** — minimises J(x) = Jₒ(x) + α·Jₛ(x) using Polak-Ribière conjugate gradient with Armijo backtracking line search. Observation weights model beam volume broadening with range.
- **CAPPI baseline** — traditional constant-altitude plan position indicator with IDW interpolation between sweeps, retained as a fast alternative and comparison baseline.
- **CF-compliant output** — NetCDF output with standard CF-1.10 attributes, grid mapping, and per-field metadata (units, standard_name, long_name) mapped from ODIM quantity names.
- **Multi-threaded** — altitude layers are processed in parallel via a bounded thread pool.

## Dependencies

- **[bom-util]** and **[bom-core]** — Bureau of Meteorology C++ utility and radar libraries
- **HDF5** — for reading ODIM polar volume files
- **NetCDF** — for writing gridded output
- **PROJ** — for map projection (geographic ↔ projected coordinate transforms)
- **C++20** compiler (GCC ≥ 11, Clang ≥ 14)

## Building

```bash
mkdir build && cd build
cmake ..
make
```

The BOM libraries are expected at the install prefix. System dependencies (HDF5, NetCDF, PROJ) are found via CMake config files or pkg-config.

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

The configuration file controls grid geometry, field selection, and gridding parameters. Key sections:

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
# Only grid these fields
include_fields "DBZH VRADH ZDR"

# Or grid everything except these
exclude_fields "SQI CCOR"
```

`include_fields` takes precedence if both are set. A warning is logged for any requested field not found in the input.

### Gridding method

```
gridding_method variational    # or "cappi"
```

### Variational parameters

| Parameter | Default | Description |
|---|---|---|
| `vargrid_alpha` | 1.0 | Regularisation weight (larger = smoother) |
| `vargrid_max_alt_diff` | 2000 | Max altitude difference (m) for gate selection |
| `vargrid_max_iterations` | 50 | Maximum conjugate gradient iterations |
| `vargrid_tolerance` | 1e-5 | Convergence tolerance (relative residual) |
| `vargrid_beam_power` | 2.0 | Exponent for range-dependent observation weight |
| `vargrid_ref_range` | 10000 | Reference range (m) where beam weight = 1.0 |
| `vargrid_min_weight` | 0.01 | Minimum observation weight floor |
| `vargrid_use_nearest_init` | true | Initialise solver with weighted nearest-gate values |

### CAPPI parameters

| Parameter | Default | Description |
|---|---|---|
| `max_alt_dist` | 20000 | Max altitude distance (m) for sweep selection |
| `idw_pwr` | 2.0 | Inverse distance weighting exponent |

### Other

| Parameter | Default | Description |
|---|---|---|
| `velocity` | VRADH | Velocity moment name (controls undetect handling) |
| `output_obs_count` | false | Write `nobs_<FIELD>` diagnostic variables |

## Repository structure

```
vargrid/
├── CMakeLists.txt
├── README.md
├── src/
│   ├── main.cc                        — CLI and orchestration
│   ├── types.h                        — Data structures (sweep, volume) and BOM includes
│   ├── util.h                         — phase_timer, split_fields (header-only)
│   ├── reader.h / reader.cc           — ODIM reading, field discovery, metadata extraction
│   ├── writer.h / writer.cc           — CF-compliant NetCDF output
│   ├── grid.h / grid.cc              — Altitude init, flip/flipud utilities
│   └── gridding/
│       ├── config.h                   — vargrid_config parameter struct
│       ├── observation_operator.h/cc  — Gate-to-grid mapping and projection
│       ├── cost_function.h            — J(x), ∇J(x), Laplacian (header-only)
│       ├── solver.h                   — Polak-Ribière CG solver (header-only)
│       ├── variational.h/cc           — Variational gridding driver
│       └── cappi.h/cc                 — CAPPI IDW gridding
```

## How it works

### Variational gridding

The variational method finds the gridded field **x** that minimises:

**J(x) = Jₒ(x) + α · Jₛ(x)**

where Jₒ is the weighted sum of squared differences between the gridded field and radar observations, and Jₛ is the Laplacian smoothness penalty. The observation weights model beam volume broadening: w = (r₀/r)^p, so nearby gates are trusted more than far-range gates whose pulse volumes average over larger physical regions.

The solver uses Polak-Ribière nonlinear conjugate gradient with Armijo backtracking line search. Convergence is typically achieved in 10–50 iterations.

### Gate-to-grid mapping

Radar gates are mapped to the Cartesian grid by projecting their geographic coordinates (lat/lon from geodetic beam propagation) into the grid's projected coordinate system using PROJ. This projection is precomputed once on the main thread for all gates, then shared read-only across worker threads to avoid PROJ thread-safety issues.

### CAPPI

The CAPPI method finds the two nearest sweeps above and below the target altitude for each grid cell and interpolates using inverse distance weighting. Sweeps are pre-sorted by elevation for early termination.

[bom-util]: https://gitlab.com/bom-radar/bom-util
[bom-core]: https://gitlab.com/bom-radar/bom-core