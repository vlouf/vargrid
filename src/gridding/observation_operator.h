#ifndef VARGRID_GRIDDING_OBSERVATION_OPERATOR_H
#define VARGRID_GRIDDING_OBSERVATION_OPERATOR_H

#include <vector>
#include <string>
#include <cstddef>
#include <cmath>

// A single observation: one radar gate contributing to the grid.
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

  auto grid_size() const -> size_t { return grid_nx * grid_ny; }
};

// Precomputed projected coordinates for all gates in one sweep.
struct sweep_projections {
  std::vector<double> px;  // projected x (easting)
  std::vector<double> py;  // projected y (northing)
  size_t nrays;
  size_t nbins;
};

// Precomputed projections for all sweeps in a volume.
struct gate_projections {
  std::vector<sweep_projections> sweeps;
};

#endif // VARGRID_GRIDDING_OBSERVATION_OPERATOR_H
