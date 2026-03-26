#ifndef VARGRID_GRIDDING_CONFIG_H
#define VARGRID_GRIDDING_CONFIG_H

#include <cmath>
#include <limits>

// Configuration for the variational gridding solver.
struct vargrid_config {
  // Regularisation weight (smoothness vs data fit).
  float alpha = 1.0f;

  // Maximum altitude difference (m) for gate selection.
  float max_alt_diff = 2000.0f;

  // Solver parameters
  int max_iterations = 50;
  float tolerance = 1e-5f;

  // Observation error model
  float beam_power = 2.0f;       // Exponent for range-dependent weight
  float ref_range = 10000.f;     // Reference range (m) where w_beam = 1.0
  float min_weight = 0.01f;      // Floor to prevent zero weights
  float beamwidth = 1.0f;        // Radar beamwidth in degrees

  // Initial guess strategy
  bool use_nearest_init = true;

  // Background field value (NaN = no background constraint)
  float background = std::numeric_limits<float>::quiet_NaN();
};

#endif // VARGRID_GRIDDING_CONFIG_H
