#ifndef VARGRID_OBSERVATION_OPERATOR_H
#define VARGRID_OBSERVATION_OPERATOR_H

#include <vector>
#include <cstddef>
#include <cmath>

// A single observation: one radar gate contributing to the grid.
struct observation {
  size_t grid_index;  // Flattened index into the 2D grid (y * nx + x)
  float value;        // Observed value at this gate
  float weight;       // Observation quality weight (R⁻¹ component)
  float alt_dist;     // Altitude distance from target (for diagnostics)
};

// The observation operator H maps from radar polar coordinates to the
// Cartesian grid. For each altitude layer, it collects all radar gates
// that are close enough to the target altitude, computes which grid cell
// each gate falls into, and assigns a quality weight.
//
// This is a sparse operator: most grid cells have 0-3 contributing
// observations from different sweeps.
//
// The operator is built once per altitude layer and reused across
// solver iterations.
struct observation_operator {
  std::vector<observation> obs;  // All observations for this layer

  size_t grid_nx;
  size_t grid_ny;

  // Number of observations contributing to each grid cell.
  // Used to identify data-void regions.
  std::vector<int> obs_count;

  auto grid_size() const -> size_t { return grid_nx * grid_ny; }
};

#endif // VARGRID_OBSERVATION_OPERATOR_H
