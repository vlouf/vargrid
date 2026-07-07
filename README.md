# vargrid

`vargrid` grids weather radar polar volumes onto a projected Cartesian grid and writes CF-compliant NetCDF output.

This README focuses on the production-ready methods:

- CAPPI (inverse-distance weighted vertical interpolation)
- LeROI (sweep-surface interpolation based on Dahl et al. 2019)

## What It Does

- Reads radar moments from ODIM HDF5 volumes
- Also supports CF/Radial-style NetCDF input
- Produces 3D Cartesian fields (`z`, `y`, `x`) in a user-defined map projection
- Supports selecting or excluding moments by name
- Optionally runs post-processing (currently Steiner convective/stratiform classification)
- Writes CF metadata and optional packed int16 output

## Build

### Dependencies

At configure/build time, `vargrid` uses:

- CMake >= 3.14
- C++20 compiler
- `bom-util`
- `bom-core`
- HDF5 (C + HL)
- netCDF
- PROJ
- Threads
- Optional: UDUNITS2

### Compile

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

On single-config generators (Ninja/Unix Makefiles), the executable is typically:

```text
build/vargrid
```

On multi-config generators (Visual Studio), use the selected configuration output.

## Command Line

```bash
vargrid [options] config.conf input_volume output.nc
```

Options:

| Flag | Description |
|---|---|
| `-h`, `--help` | Show usage and exit |
| `-v`, `--version` | Print version and exit |
| `-g`, `--generate [input]` | Print sample config and exit. If `input` is provided, projection center is initialized from radar location. |
| `-t`, `--trace=LEVEL` | Set log level: `none`, `status`, `error`, `warning`, `log`, `debug` |

Examples:

```bash
# Generate a template config
vargrid -g > flow.conf

# Generate template config centered on a specific ODIM file
vargrid -g 52_20240806_033500.pvol.h5 > flow.conf

# Run gridding
vargrid -t log flow.conf 52_20240806_033500.pvol.h5 out.nc
```

## Configuration Overview

Required geometry keys:

- `proj4`
- `size`
- `left_top`
- `cell_delta`
- `units`
- `origin` (`xy` or `ij`)
- `altitude_base`
- `altitude_step`
- `layer_count`

Common processing keys:

- `gridding_method` (`cappi` or `leroi`)
- `altitude_reference` (`sea` or `radar`)
- `include_fields` / `exclude_fields` (space-separated lists)
- `velocity` (default `VRADH`, used for velocity-specific undetect handling)
- `pack_output` (`true`/`false`)
- `output_obs_count` (only meaningful for methods not covered here; usually leave `false`)

### Minimal CAPPI Config

```conf
proj4 "+proj=aea +lat_1=-28.7 +lat_2=-38.7 +lon_0=151.20900 +lat_0=-33.70080 +units=m +ellps=GRS80"
size "301 301"
left_top "-150500 150500"
cell_delta "1000 -1000"
units m
origin xy

altitude_base 0.0
altitude_step 500.0
layer_count 13
altitude_reference sea

gridding_method cappi

# CAPPI parameters
max_alt_dist 20000
idw_pwr 2.0

pack_output true
output_obs_count false
```

### Minimal LeROI Config

```conf
proj4 "+proj=aea +lat_1=-28.7 +lat_2=-38.7 +lon_0=151.20900 +lat_0=-33.70080 +units=m +ellps=GRS80"
size "301 301"
left_top "-150500 150500"
cell_delta "1000 -1000"
units m
origin xy

altitude_base 0.0
altitude_step 500.0
layer_count 13
altitude_reference sea

gridding_method leroi

# LeROI parameters
leroi_weight barnes
leroi_roi 0
leroi_idw_pwr 2.0
leroi_ground_elevation -999

pack_output true
output_obs_count false
```

## CAPPI Notes

- Uses nearest lower and upper sweeps around each requested altitude, then blends with IDW.
- `idw_pwr` controls vertical weighting sharpness.
- `max_alt_dist` limits how far a sweep gate can be from requested altitude.
- If only one side (lower or upper) is available, CAPPI can still use it if within a fixed radius threshold.

## LeROI Notes

- Precomputes horizontally interpolated sweep surfaces for each field, then slices vertically to target altitudes.
- Excludes near-vertical birdbath scans.
- Requires at least two non-vertical sweeps.
- `leroi_weight` options: `barnes`, `cressman`, `idw`, `bilinear`.
- `leroi_roi=0` auto-estimates radius of influence (not used by `bilinear`).
- `leroi_ground_elevation` can anchor the lowest sweep at radar altitude for low-elevation scans.

## Field Selection

By default, all discovered moments are gridded.

Use either:

- `include_fields "DBZH VRADH ZDR"`
- `exclude_fields "SQI CCOR"`

If both are absent, all moments are processed.

## Post-Processing

Optional post-processors can be applied per altitude layer.

Current processor:

- `steiner` (outputs `STEINER_CLASS` from gridded `DBZH`)

Enable with:

```conf
post_processors "steiner"
```

## Testing

Unit tests are optional and disabled by default.

```bash
cmake -S . -B build -DVARGRID_BUILD_TESTS=ON
cmake --build build
ctest --test-dir build --output-on-failure
```

## License

See `LICENSE`.