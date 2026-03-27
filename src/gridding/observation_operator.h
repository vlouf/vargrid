#ifndef VARGRID_GRIDDING_OBSERVATION_OPERATOR_H
#define VARGRID_GRIDDING_OBSERVATION_OPERATOR_H

#include <vector>
#include <string>
#include <cstddef>
#include <cmath>

namespace bom { class angle; struct latlon; struct latlonalt; template<typename T> class array2; }

// A single observation contribution to one grid cell.
struct observation {
  size_t grid_index;
  float value;
  float weight;
  float alt_dist;
};

// Sparse mapping from radar gates to grid cells.
struct observation_operator {
  std::vector<observation> obs;
  size_t grid_nx;
  size_t grid_ny;
  std::vector<int> obs_count;

  float kappa = 0.0f;

  // Per-cell azimuth from radar (degrees, [0,360)) — used for Wx/Wy weights.
  std::vector<float> cell_azimuth_deg;

  // Per-cell distance to nearest observation (in grid cells) — used for JB.
  // 0 = cell has observations, large = far from data.
  std::vector<float> dist_to_obs;

  auto grid_size() const -> size_t { return grid_nx * grid_ny; }
};

// Precomputed bearing/range lookup for mapping gates to grid cells.
struct grid_bearings {
  std::vector<float> cell_bearing_deg;  // [ny * nx]
  std::vector<float> cell_range;        // [ny * nx]

  std::vector<size_t> bin_to_grid;

  size_t n_bearing_bins;
  size_t n_range_bins;
  float bearing_bin_size;
  float range_bin_size;
  float max_range;

  size_t nx;
  size_t ny;
};

#endif // VARGRID_GRIDDING_OBSERVATION_OPERATOR_H
