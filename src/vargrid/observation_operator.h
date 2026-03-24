#ifndef VARGRID_OBSERVATION_OPERATOR_H
#define VARGRID_OBSERVATION_OPERATOR_H

#include <vector>
#include <string>
#include <cstddef>
#include <cmath>

// A single observation: one radar gate contributing to the grid.
struct observation {
  size_t grid_index;  // Flattened index into the 2D grid (y * nx + x)
  float value;        // Observed value at this gate
  float weight;       // Observation quality weight (R⁻¹ component)
  float alt_dist;     // Altitude distance from target (for diagnostics)
};

// The observation operator H: sparse mapping from radar gates to grid cells.
struct observation_operator {
  std::vector<observation> obs;
  size_t grid_nx;
  size_t grid_ny;
  std::vector<int> obs_count;

  auto grid_size() const -> size_t { return grid_nx * grid_ny; }
};

// Precomputed projected coordinates for all radar gates in all sweeps.
// Computed once on the main thread (PROJ is not thread-safe),
// then shared read-only across worker threads.
struct sweep_projections {
  std::vector<double> px;  // projected x (easting) per gate [nrays * nbins]
  std::vector<double> py;  // projected y (northing) per gate [nrays * nbins]
  size_t nrays;
  size_t nbins;
};

struct gate_projections {
  std::vector<sweep_projections> sweeps;
};

#endif // VARGRID_OBSERVATION_OPERATOR_H