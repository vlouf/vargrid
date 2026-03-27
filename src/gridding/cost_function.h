#ifndef VARGRID_COST_FUNCTION_H
#define VARGRID_COST_FUNCTION_H

#include "observation_operator.h"
#include "config.h"
#include <vector>
#include <cmath>

// Cost function following Brook et al. (2022):
//
//   J(φ) = Jo + λH·Js + λB·JB
//
//   Jo = sum_i w_i (Rφ - d_i)²                            (data fidelity)
//   Js = ||Wy·φ_yy||² + ||Wx·φ_xx||²                     (smoothness, Eq. 2)
//   JB = ||exp(-rc²/r²) · φ||²                            (background, Eq. 4)
//
// Azimuthal weights (Eq. 3):
//   Wx(φ) = C - A·cos(2φ),   Wy(φ) = C + A·cos(2φ)
//   where A = |f-1|/2, C = (f+1)/2, f = Δr/Δϕ

// Compute second derivative in x: φ_xx (centered, Neumann BC)
inline void compute_dxx(const float* x, float* dxx, size_t nx, size_t ny)
{
  for (size_t j = 0; j < ny; ++j) {
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float left  = (i > 0)      ? x[idx - 1] : x[idx];
      float right = (i < nx - 1) ? x[idx + 1] : x[idx];
      dxx[idx] = left - 2.0f * x[idx] + right;
    }
  }
}

// Compute second derivative in y: φ_yy (centered, Neumann BC)
inline void compute_dyy(const float* x, float* dyy, size_t nx, size_t ny)
{
  for (size_t j = 0; j < ny; ++j) {
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float up   = (j > 0)      ? x[idx - nx] : x[idx];
      float down = (j < ny - 1) ? x[idx + nx] : x[idx];
      dyy[idx] = up - 2.0f * x[idx] + down;
    }
  }
}

// Perona-Malik diffusivity
inline float pm_diffusivity(float diff_sq, float kappa_sq) {
  return 1.0f / (1.0f + diff_sq / kappa_sq);
}

// Perona-Malik first-order edge-preserving smoothness
inline auto evaluate_pm_smoothness(
    const float* x,
    float* grad_s,
    size_t nx,
    size_t ny,
    float kappa
    ) -> float
{
  float kappa_sq = kappa * kappa;
  float Js = 0.0f;

  for (size_t j = 0; j < ny; ++j) {
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float center = x[idx];
      float g_sum = 0.0f;

      if (i < nx - 1) {
        float d = center - x[idx + 1];
        g_sum += pm_diffusivity(d * d, kappa_sq) * d;
        Js += kappa_sq * std::log(1.0f + d * d / kappa_sq);
      }
      if (j < ny - 1) {
        float d = center - x[idx + nx];
        g_sum += pm_diffusivity(d * d, kappa_sq) * d;
        Js += kappa_sq * std::log(1.0f + d * d / kappa_sq);
      }
      if (i > 0) {
        float d = center - x[idx - 1];
        g_sum += pm_diffusivity(d * d, kappa_sq) * d;
      }
      if (j > 0) {
        float d = center - x[idx - nx];
        g_sum += pm_diffusivity(d * d, kappa_sq) * d;
      }

      grad_s[idx] = 2.0f * g_sum;
    }
  }

  return Js;
}

// Evaluate the full cost function and gradient.
//
//   J(φ) = Jo + λH·Js + λB·JB
//
// work buffer must be at least 2 * nx * ny floats.
inline auto evaluate_gradient(
    const float* x,
    float* grad,
    const observation_operator& H,
    const vargrid_config& cfg,
    float* work
    ) -> float
{
  size_t n = H.grid_size();
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;

  for (size_t i = 0; i < n; ++i)
    grad[i] = 0.0f;

  // === Jo: observation fidelity ===
  float Jo = 0.0f;
  for (auto& obs : H.obs) {
    float residual = x[obs.grid_index] - obs.value;
    Jo += obs.weight * residual * residual;
    grad[obs.grid_index] += 2.0f * obs.weight * residual;
  }

  // === Js: smoothness ===
  float Js = 0.0f;

  if (H.kappa > 0.0f) {
    // Perona-Malik (no azimuthal weights — they only apply to second-order)
    Js = evaluate_pm_smoothness(x, work, nx, ny, H.kappa);
    for (size_t i = 0; i < n; ++i)
      grad[i] += cfg.lambda_h * work[i];
  } else {
    // Second-order smoothness with azimuthal weights (Brook et al. 2022 Eq. 2-3).
    //
    // Js = sum_k [ Wx[k] * φ_xx[k]² + Wy[k] * φ_yy[k]² ]
    // ∇Js[k] = 2 * Dxx^T(Wx · φ_xx) + 2 * Dyy^T(Wy · φ_yy)

    compute_dxx(x, work, nx, ny);          // work[0..n) = φ_xx
    compute_dyy(x, work + n, nx, ny);      // work[n..2n) = φ_yy

    // Compute azimuthal weights per cell and apply to second derivatives.
    // f = range_spacing / (beamwidth_rad * cell_range)
    // If azimuthal weighting is disabled (range_spacing = 0), Wx = Wy = 1.
    bool use_azimuthal = (cfg.range_spacing > 0.0f)
                      && !H.cell_azimuth_deg.empty();

    Js = 0.0f;
    for (size_t i = 0; i < n; ++i) {
      float Wx = 1.0f, Wy = 1.0f;

      if (use_azimuthal) {
        float azimuth_rad = H.cell_azimuth_deg[i] * (M_PI / 180.0f);
        // f = Δr / Δϕ where Δϕ = beamwidth_rad * range
        // For cells near the radar, clamp range to avoid f->infinity
        float cell_r = (H.cell_azimuth_deg.size() > i) ? 0.0f : 0.0f;
        // We don't have cell_range in H, use a representative value.
        // Actually, we can get it from gb via the grid_bearings, but it's
        // not stored in H. Use the range_spacing and beamwidth to compute
        // a global f, or per-cell from stored bearing info.
        // For simplicity, use a global f based on mid-range of the grid.
        // TODO: store cell_range in H for per-cell f.
        float mid_range = 75000.0f;  // approximate mid-range
        float delta_phi = cfg.beamwidth * (M_PI / 180.0f) * mid_range;
        float f = cfg.range_spacing / delta_phi;
        if (f > 1.0f) f = 1.0f;  // clamp: range spacing < azimuthal spacing

        float A = std::fabs(f - 1.0f) / 2.0f;
        float C = (f + 1.0f) / 2.0f;
        float cos2phi = std::cos(2.0f * azimuth_rad);

        Wx = C - A * cos2phi;
        Wy = C + A * cos2phi;
      }

      // Apply weights to second derivatives
      float wx_dxx = Wx * work[i];
      float wy_dyy = Wy * work[n + i];

      Js += wx_dxx * work[i] + wy_dyy * work[n + i];

      // Store weighted derivatives for adjoint pass
      work[i] = wx_dxx;
      work[n + i] = wy_dyy;
    }

    // Gradient: 2 * (Dxx^T(Wx·φ_xx) + Dyy^T(Wy·φ_yy))
    for (size_t j = 0; j < ny; ++j) {
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;

        float xx_l = (i > 0)      ? work[idx - 1] : work[idx];
        float xx_r = (i < nx - 1) ? work[idx + 1] : work[idx];
        float g_xx = xx_l - 2.0f * work[idx] + xx_r;

        float yy_u = (j > 0)      ? work[n + idx - nx] : work[n + idx];
        float yy_d = (j < ny - 1) ? work[n + idx + nx] : work[n + idx];
        float g_yy = yy_u - 2.0f * work[n + idx] + yy_d;

        grad[idx] += 2.0f * cfg.lambda_h * (g_xx + g_yy);
      }
    }
  }

  // === JB: background constraint (Brook et al. 2022 Eq. 4) ===
  // JB = sum_k [ exp(-2 rc²/r²) * φ[k]² ]
  // ∇JB[k] = 2 * exp(-2 rc²/r²) * φ[k]
  // The exponential weighting means: cells near observations (r small) are
  // barely penalised, while cells far from data (r >> rc) are penalised
  // strongly toward zero (the background).
  float JB = 0.0f;

  if (cfg.lambda_b > 0.0f && !H.dist_to_obs.empty()) {
    float rc = cfg.bg_cutoff_cells;
    float rc_sq = rc * rc;
    float bg = std::isnan(cfg.background) ? 0.0f : cfg.background;

    for (size_t i = 0; i < n; ++i) {
      float r = H.dist_to_obs[i];
      if (r < 0.5f) continue;  // cell has observations — skip

      float r_sq = r * r;
      // Brook: exp(-rc²/r²) — large when r >> rc, small when r << rc
      float w_bg = std::exp(-rc_sq / r_sq);

      float dev = x[i] - bg;
      JB += w_bg * dev * dev;
      grad[i] += 2.0f * cfg.lambda_b * w_bg * dev;
    }
  }

  return Jo + cfg.lambda_h * Js + cfg.lambda_b * JB;
}

// Legacy Laplacian — kept for unit tests
inline void apply_laplacian(
    const float* x,
    float* Lx,
    size_t nx,
    size_t ny)
{
  for (size_t j = 0; j < ny; ++j) {
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float center = x[idx];
      float sum = 0.0f;
      if (i > 0)      sum += center - x[idx - 1];
      if (i < nx - 1) sum += center - x[idx + 1];
      if (j > 0)      sum += center - x[idx - nx];
      if (j < ny - 1) sum += center - x[idx + nx];
      Lx[idx] = sum;
    }
  }
}

#endif // VARGRID_COST_FUNCTION_H
