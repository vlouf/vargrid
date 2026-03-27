#ifndef VARGRID_GRIDDING_OBSERVATION_OPERATOR_H
#define VARGRID_GRIDDING_OBSERVATION_OPERATOR_H

#include <vector>
#include <string>
#include <cstddef>
#include <cmath>

namespace bom { class angle; struct latlon; struct latlonalt; template<typename T> class array2; }

// A single observation contribution to one grid cell.
// With bilinear interpolation, a single radar gate produces up to 4 observations
// (one per surrounding grid cell), each with a bilinear weight multiplied into
// the observation weight.
struct observation {
  size_t grid_index;
  float value;
  float weight;      // observation quality weight × bilinear interpolation weight
  float alt_dist;
};

// Sparse mapping from radar gates to grid cells.
struct observation_operator {
  std::vector<observation> obs;
  size_t grid_nx;
  size_t grid_ny;
  std::vector<int> obs_count;

  float kappa = 0.0f;

  auto grid_size() const -> size_t { return grid_nx * grid_ny; }
};

// Precomputed bearing/range lookup for mapping gates to grid cells.
// Stores both a binned lookup table (for O(1) nearest-cell finding)
// and the per-cell bearing/range values (for bilinear weight computation).
struct grid_bearings {
  // Per-cell bearing (degrees, [0,360)) and ground range (m)
  std::vector<float> cell_bearing_deg;  // [ny * nx]
  std::vector<float> cell_range;        // [ny * nx]

  // Binned lookup: bin_to_grid[bearing_bin * n_range_bins + range_bin] = grid index
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
