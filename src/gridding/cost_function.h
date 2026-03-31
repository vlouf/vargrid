#ifndef VARGRID_COST_FUNCTION_H
#define VARGRID_COST_FUNCTION_H

#include "observation_operator.h"
#include "config.h"
#include <vector>
#include <cmath>

// Cost function:
//   J(φ) = Jo + λH · Js
//
//   Jo = Σ w_i (φ[g_i] - d_i)²                       (data fidelity)
//   Js = Σ [ Wx·φ_xx² + Wy·φ_yy² ]                   (smoothness, Eq. 2)
//
// Azimuthal weights (Eq. 3):
//   Wx = C - A·cos(2φ),  Wy = C + A·cos(2φ)

inline void compute_dxx(const float* x, float* dxx, size_t nx, size_t ny)
{
  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float left  = (i > 0)      ? x[idx - 1] : x[idx];
      float right = (i < nx - 1) ? x[idx + 1] : x[idx];
      dxx[idx] = left - 2.0f * x[idx] + right;
    }
}

inline void compute_dyy(const float* x, float* dyy, size_t nx, size_t ny)
{
  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float up   = (j > 0)      ? x[idx - nx] : x[idx];
      float down = (j < ny - 1) ? x[idx + nx] : x[idx];
      dyy[idx] = up - 2.0f * x[idx] + down;
    }
}

inline float pm_diffusivity(float diff_sq, float kappa_sq) {
  return 1.0f / (1.0f + diff_sq / kappa_sq);
}

inline auto evaluate_pm_smoothness(
    const float* x, float* grad_s,
    size_t nx, size_t ny, float kappa) -> float
{
  float kappa_sq = kappa * kappa;
  float Js = 0.0f;

  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float center = x[idx];
      float g_sum = 0.0f;

      if (i < nx - 1) { float d = center - x[idx+1]; g_sum += pm_diffusivity(d*d, kappa_sq) * d; Js += kappa_sq * std::log(1.0f + d*d / kappa_sq); }
      if (j < ny - 1) { float d = center - x[idx+nx]; g_sum += pm_diffusivity(d*d, kappa_sq) * d; Js += kappa_sq * std::log(1.0f + d*d / kappa_sq); }
      if (i > 0)      { float d = center - x[idx-1]; g_sum += pm_diffusivity(d*d, kappa_sq) * d; }
      if (j > 0)      { float d = center - x[idx-nx]; g_sum += pm_diffusivity(d*d, kappa_sq) * d; }

      grad_s[idx] = 2.0f * g_sum;
    }
  return Js;
}

// Main cost function and gradient evaluation.
inline auto evaluate_gradient(
    const float* x,
    float* grad,
    const observation_operator& H,
    const vargrid_config& cfg,
    float* work   // size >= 2 * nx * ny
    ) -> float
{
  size_t n = H.grid_size();
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;

  for (size_t i = 0; i < n; ++i)
    grad[i] = 0.0f;

  // Jo: data fidelity
  float Jo = 0.0f;
  for (auto& obs : H.obs) {
    float r = x[obs.grid_index] - obs.value;
    Jo += obs.weight * r * r;
    grad[obs.grid_index] += 2.0f * obs.weight * r;
  }

  // Js: smoothness
  float Js = 0.0f;

  if (H.kappa > 0.0f) {
    Js = evaluate_pm_smoothness(x, work, nx, ny, H.kappa);
    for (size_t i = 0; i < n; ++i)
      grad[i] += cfg.lambda_h * work[i];
  } else {
    // Second-order smoothness with azimuthal weights (Brook et al. 2022 Eq. 2-3)
    // plus a first-order Laplacian damping term to suppress biharmonic ringing.
    //
    //   Js = Σ [ Wx·φ_xx² + Wy·φ_yy² ]   (second-order, spreads information)
    //   Jd = Σ_edges (φ[j] - φ[k])²       (first-order, damps oscillations)
    //
    //   Total smoothness = λH · Js + (λH/10) · Jd
    //
    // The first-order term is weighted at 1/10 of λH — just enough to
    // suppress ringing without dominating the second-order behaviour.

    compute_dxx(x, work, nx, ny);
    compute_dyy(x, work + n, nx, ny);

    bool use_az = (cfg.range_spacing > 0.0f)
              && !H.cell_azimuth_deg.empty()
              && !H.cell_range.empty();

    // --- Second-order cost and weighted derivatives ---
    Js = 0.0f;
    for (size_t i = 0; i < n; ++i) {
      float Wx = 1.0f, Wy = 1.0f;

      if (use_az) {
        // Only apply anisotropic weights where azimuthal spacing >> range spacing.
        // Near the radar (< 20km), the grid is well-sampled in all directions
        // so isotropic smoothing is appropriate. This avoids rapidly varying
        // weights that cause solver oscillations at low altitudes.
        float cell_r = H.cell_range[i];
        if (cell_r < 20000.0f) {
          // Isotropic — leave Wx = Wy = 1
        } else {
          float delta_phi = cfg.beamwidth * (M_PI / 180.0f) * cell_r;
          float f = cfg.range_spacing / delta_phi;
          if (f > 1.0f) f = 1.0f;

          float Af = std::fabs(f - 1.0f) / 2.0f;
          float Cf = (f + 1.0f) / 2.0f;
          float cos2phi = std::cos(2.0f * H.cell_azimuth_deg[i] * (M_PI / 180.0f));

          Wx = Cf - Af * cos2phi;
          Wy = Cf + Af * cos2phi;
        }
      }

      Js += Wx * work[i] * work[i] + Wy * work[n + i] * work[n + i];
      work[i] *= Wx;
      work[n + i] *= Wy;
    }

    // --- Second-order gradient: 2 · (Dxx^T(Wx·φ_xx) + Dyy^T(Wy·φ_yy)) ---
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;

        float xx_l = (i > 0)      ? work[idx-1] : work[idx];
        float xx_r = (i < nx - 1) ? work[idx+1] : work[idx];
        float g_xx = xx_l - 2.0f * work[idx] + xx_r;

        float yy_u = (j > 0)      ? work[n+idx-nx] : work[n+idx];
        float yy_d = (j < ny - 1) ? work[n+idx+nx] : work[n+idx];
        float g_yy = yy_u - 2.0f * work[n+idx] + yy_d;

        grad[idx] += 2.0f * cfg.lambda_h * (g_xx + g_yy);
      }

    // --- First-order damping term: suppresses biharmonic ringing ---
    // Equal weight with the second-order term provides effective
    // oscillation damping while preserving the information-spreading
    // property of the second-order penalty.
    float lambda_damp = cfg.lambda_h;
    float Jd = 0.0f;
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        float center = x[idx];
        float g = 0.0f;

        if (i > 0)      { float d = center - x[idx-1];  g += d; }
        if (i < nx - 1) { float d = center - x[idx+1];  g += d; }
        if (j > 0)      { float d = center - x[idx-nx]; g += d; }
        if (j < ny - 1) { float d = center - x[idx+nx]; g += d; }

        // Cost: count each edge once (right and down only)
        if (i < nx - 1) { float d = center - x[idx+1];  Jd += d * d; }
        if (j < ny - 1) { float d = center - x[idx+nx]; Jd += d * d; }

        grad[idx] += 2.0f * lambda_damp * g;
      }

    Js += Jd;  // first-order damping included at same weight
  }

  return Jo + cfg.lambda_h * Js;
}

// Legacy — kept for unit tests
inline void apply_laplacian(const float* x, float* Lx, size_t nx, size_t ny)
{
  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float center = x[idx], sum = 0.0f;
      if (i > 0)      sum += center - x[idx - 1];
      if (i < nx - 1) sum += center - x[idx + 1];
      if (j > 0)      sum += center - x[idx - nx];
      if (j < ny - 1) sum += center - x[idx + nx];
      Lx[idx] = sum;
    }
}

#endif // VARGRID_COST_FUNCTION_H