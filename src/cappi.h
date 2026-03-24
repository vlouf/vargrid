#ifndef CAPPI_H
#define CAPPI_H

#include "pch.h"
using namespace bom;

// Precomputed per-grid-point, per-sweep lookup indices.
// Computed once and reused across all altitude layers and CAPPI calls
// for the same volume, eliminating redundant geodetic and search computations.
struct sweep_lookup {
  size_t ibin;  // bin index for this grid point in this sweep
  size_t iray;  // ray index for this grid point in this sweep
};

struct grid_lookup {
  // For each grid point [y][x], stores per-sweep indices.
  // Layout: grid_ny * grid_nx * nsweeps
  std::vector<sweep_lookup> data;
  size_t grid_nx;
  size_t grid_ny;
  size_t nsweeps;

  sweep_lookup const& at(size_t y, size_t x, size_t isweep) const {
    return data[(y * grid_nx + x) * nsweeps + isweep];
  }
  sweep_lookup& at(size_t y, size_t x, size_t isweep) {
    return data[(y * grid_nx + x) * nsweeps + isweep];
  }
};

auto precompute_grid_lookup(
      volume const& vol
    , array2<latlon> const& latlons
    ) -> grid_lookup;

auto compute_roi(const float range, const float beamwidth, const float max_roi) -> float;
auto find_ground_range_bin(array1<bin_info> const& bins, float target) -> size_t;
auto find_ray(array1<angle> const& rays, angle target) -> size_t;

auto generate_cappi(
      volume const& vol
    , array2<latlon> const& latlons
    , float max_alt_diff
    , float idw_pwr
    , float altitude
    , float roi = 2500.f
    , const float beamwidth = 1.f
    ) -> array2f;

// Optimized version using precomputed lookup table
auto generate_cappi(
      volume const& vol
    , grid_lookup const& lookup
    , size_t grid_ny
    , size_t grid_nx
    , float max_alt_diff
    , float idw_pwr
    , float altitude
    , float roi = 2500.f
    ) -> array2f;

#endif