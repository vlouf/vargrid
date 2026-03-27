#ifndef VARGRID_GRIDDING_CONFIG_H
#define VARGRID_GRIDDING_CONFIG_H

#include <cmath>
#include <limits>

// Configuration for the variational gridding solver.
struct vargrid_config {
  // Smoothing weights (Brook et al. 2022 formulation).
  // lambda_h: horizontal smoothing weight (second-order derivative penalty)
  // lambda_v: vertical smoothing weight (not used in 2D single-layer mode)
  // When lambda_h > 0, uses second-order smoothing ||phi_xx||^2 + ||phi_yy||^2
  // The old "alpha" parameter is mapped to lambda_h for backward compatibility.
  float lambda_h = 0.01f;

  // Maximum altitude difference (m) for gate selection.
  float max_alt_diff = 2000.0f;

  // Solver parameters
  int max_iterations = 200;
  float tolerance = 1e-5f;

  // Observation error model
  float beam_power = 2.0f;
  float ref_range = 10000.f;
  float min_weight = 0.01f;
  float beamwidth = 1.0f;

  // Perona-Malik edge threshold (data units). 0 = disabled.
  float kappa = 0.0f;

  // Initial guess strategy
  bool use_nearest_init = true;

  // Background field value (NaN = no background constraint)
  float background = std::numeric_limits<float>::quiet_NaN();
};

#endif // VARGRID_GRIDDING_CONFIG_H
