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

  // Observation error model
  // The weight (inverse of observation error variance) for each gate is:
  //   w = w_beam * w_alt
  // where:
  //   w_beam = (ref_range / slant_range)^beam_power
  //     Models beam broadening: the pulse volume grows with range,
  //     so far gates average over larger volumes and are less
  //     representative of a single grid cell.
  //   w_alt = 1 - |alt_dist| / max_alt_diff
  //     Linear decay for gates further from the target altitude.
  float beam_power = 2.0f;       // Exponent for range-dependent weight (2 = inverse area)
  float ref_range = 10000.f;     // Reference range (m) where w_beam = 1.0
  float min_weight = 0.01f;      // Floor to prevent zero weights

  // Radar beamwidth in degrees (used for beam volume computation).
  // Overridden by the value read from the radar file if available.
  float beamwidth = 1.0f;

  // Initial guess strategy
  bool use_nearest_init = true;

  // Background field value for areas with no observations.
  float background = std::numeric_limits<float>::quiet_NaN();
};

#endif // VARGRID_CONFIG_H