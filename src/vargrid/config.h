#ifndef VARGRID_CONFIG_H
#define VARGRID_CONFIG_H

#include <cmath>
#include <limits>

// Configuration for the variational gridding solver.
// All parameters have sensible defaults for radar reflectivity gridding.
struct vargrid_config {
  // Regularisation weight (smoothness vs data fit).
  // Larger alpha = smoother output, less faithful to observations.
  // Typical range: 0.01 - 10.0
  float alpha = 1.0f;

  // Maximum altitude difference (m) between gate and target altitude
  // for a gate to contribute to the cost function.
  float max_alt_diff = 2000.0f;

  // Solver parameters
  int max_iterations = 50;     // Maximum CG iterations
  float tolerance = 1e-5f;     // Relative residual tolerance for convergence

  // Observation quality model
  float range_scale = 150000.f; // Range (m) at which observation weight halves
  float min_weight = 0.01f;     // Minimum observation weight (floor)

  // Initial guess strategy: "zero" or "nearest"
  // "nearest" uses nearest-gate values as starting point (faster convergence)
  bool use_nearest_init = true;

  // Background field value for areas with no observations.
  // NaN means no background constraint.
  float background = std::numeric_limits<float>::quiet_NaN();
};

#endif // VARGRID_CONFIG_H