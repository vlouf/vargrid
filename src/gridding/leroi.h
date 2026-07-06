#ifndef VARGRID_GRIDDING_LEROI_H
#define VARGRID_GRIDDING_LEROI_H

#include "leroi_core.h"
#include "observation_operator.h"  // grid_bearings
#include "../types.h"

#include <map>
#include <string>
#include <vector>

// Configuration for the leroi gridding method (Dahl et al. 2019).
struct leroi_config {
  leroi::weight_type weight = leroi::weight_type::barnes;
  float roi = 0.0f;                 // radius of influence (m); 0 = auto
  float idw_pwr = 2.0f;             // power for idw weighting
  float ground_elevation = -999.0f; // snap lowest sweep to h=0 if elev <= this (deg)
};

// Parse leroi_* options from the configuration file.
auto parse_leroi_config(io::configuration const& config) -> leroi_config;

// Volume-level precomputed sweep surfaces: for every used sweep (sorted by
// elevation, birdbaths excluded), the horizontally interpolated field values
// and the sweep surface heights on the 2D grid.
struct leroi_grids {
  size_t nx = 0, ny = 0;
  std::vector<std::vector<float>> heights;  // [nsweeps][ny*nx]
  std::map<std::string, std::vector<std::vector<float>>> surfaces;
};

// Horizontally interpolate all sweeps of all fields onto the grid.
// Called once per volume (main thread; parallelised internally over sweeps).
// max_altitude is the top output level, used for the auto ROI estimate.
auto leroi_precompute(
      std::map<std::string, volume> const& volumes
    , grid_bearings const& gb
    , leroi_config const& cfg
    , float max_altitude
    ) -> leroi_grids;

// Extract one altitude layer for one field by vertical interpolation
// between sweep surfaces. Thread-safe (read-only on pre).
auto leroi_slice(
      leroi_grids const& pre
    , std::string const& field
    , float altitude
    ) -> array2f;

#endif // VARGRID_GRIDDING_LEROI_H
