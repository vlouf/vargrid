#ifndef VARGRID_GRIDDING_CONFIG_H
#define VARGRID_GRIDDING_CONFIG_H

#include <cmath>
#include <limits>

struct vargrid_config {
  // Horizontal smoothing weight (Brook et al. 2022 Eq. 2).
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

  // Azimuthal weights (Brook et al. 2022 Eq. 3).
  // range_spacing = radar range gate spacing (m). 0 = disable azimuthal weights.
  float range_spacing = 250.0f;

  // Masking: cells further than this many grid cells from any observation
  // are set to NaN in the output. Controls how far smoothness can extrapolate.
  float mask_distance_cells = 3.0f;

  // Initial guess strategy
  bool use_nearest_init = true;

  // Background field value
  float background = std::numeric_limits<float>::quiet_NaN();
};

#endif // VARGRID_GRIDDING_CONFIG_H
