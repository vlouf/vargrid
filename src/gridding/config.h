#ifndef VARGRID_GRIDDING_CONFIG_H
#define VARGRID_GRIDDING_CONFIG_H

#include <cmath>
#include <limits>

// Configuration for the variational gridding solver.
struct vargrid_config {
  // Horizontal smoothing weight (Brook et al. 2022 Eq. 2).
  // Controls second-order derivative penalty: λH (||Wx φ_xx||² + ||Wy φ_yy||²)
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
  // When enabled, horizontal smoothing is weighted by azimuth to account
  // for the radar's sampling geometry: more smoothing along azimuths
  // (where data spacing is larger) than along range.
  // f = delta_r / delta_phi — ratio of range to azimuthal spacing.
  // Computed at runtime from range_spacing / (beamwidth * range).
  // Set to 0 to disable azimuthal weights (isotropic smoothing).
  float range_spacing = 250.0f;  // radar range gate spacing (m)

  // Background constraint (Brook et al. 2022 Eq. 4).
  // Penalises deviations from background in data voids:
  //   JB = ||exp(-rc²/r²) * φ||²
  // where r = distance to nearest observation, rc = cutoff radius.
  // Set lambda_b = 0 to disable.
  float lambda_b = 1.0f;          // background constraint weight
  float bg_cutoff_cells = 5.0f;   // rc in grid cells (auto-scaled to grid spacing)

  // Initial guess strategy
  bool use_nearest_init = true;

  // Background field value (NaN = no background constraint)
  float background = std::numeric_limits<float>::quiet_NaN();
};

#endif // VARGRID_GRIDDING_CONFIG_H
