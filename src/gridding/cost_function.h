#ifndef VARGRID_COST_FUNCTION_H
#define VARGRID_COST_FUNCTION_H

#include "observation_operator.h"
#include "config.h"
#include <vector>
#include <cmath>

// Cost function:
//   J(φ) = Jo + λH · (Js2 + Js1)
//
//   Jo  = Σ w_i (φ[g_i] - d_i)²
//   Js2 = Σ [ Wx[k]·φ_xx[k]² + Wy[k]·φ_yy[k]² ]
//   Js1 = Σ_edges (φ[j] - φ[k])²

inline float pm_diffusivity(float diff_sq, float kappa_sq) {
  return 1.0f / (1.0f + diff_sq / kappa_sq);
}

// === Cost-only evaluation (no gradient) ===
// Used by the line search — ~3× cheaper than full gradient evaluation.
inline auto evaluate_cost(
    const float* __restrict__ x,
    const observation_operator& H,
    const vargrid_config& cfg
    ) -> float
{
  const size_t n = H.grid_size();
  const size_t nx = H.grid_nx;
  const size_t ny = H.grid_ny;

  // Jo
  float Jo = 0.0f;
  for (const auto& obs : H.obs) {
    float r = x[obs.grid_index] - obs.value;
    Jo += obs.weight * r * r;
  }

  float Js = 0.0f;

  if (H.kappa > 0.0f) {
    float kappa_sq = H.kappa * H.kappa;
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        float c = x[idx];
        if (i < nx - 1) { float d = c - x[idx+1]; Js += kappa_sq * std::log(1.0f + d*d / kappa_sq); }
        if (j < ny - 1) { float d = c - x[idx+nx]; Js += kappa_sq * std::log(1.0f + d*d / kappa_sq); }
      }
  } else {
    const bool has_az = !H.Wx.empty();
    const float* pWx = has_az ? H.Wx.data() : nullptr;
    const float* pWy = has_az ? H.Wy.data() : nullptr;

    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        float c = x[idx];
        float left  = (i > 0)      ? x[idx-1] : c;
        float right = (i < nx - 1) ? x[idx+1] : c;
        float up    = (j > 0)      ? x[idx-nx] : c;
        float down  = (j < ny - 1) ? x[idx+nx] : c;

        float dxx = left - 2.0f * c + right;
        float dyy = up - 2.0f * c + down;
        float wx = has_az ? pWx[idx] : 1.0f;
        float wy = has_az ? pWy[idx] : 1.0f;

        // Second-order cost
        Js += wx * dxx * dxx + wy * dyy * dyy;

        // First-order damping cost (right and down edges only)
        if (i < nx - 1) { float d = c - right; Js += d * d; }
        if (j < ny - 1) { float d = c - down;  Js += d * d; }
      }
  }

  return Jo + cfg.lambda_h * Js;
}

// === Full gradient + cost evaluation ===
inline auto evaluate_gradient(
    const float* __restrict__ x,
    float* __restrict__ grad,
    const observation_operator& H,
    const vargrid_config& cfg,
    float* __restrict__ work   // size >= 2 * nx * ny
    ) -> float
{
  const size_t n = H.grid_size();
  const size_t nx = H.grid_nx;
  const size_t ny = H.grid_ny;
  const float lh = cfg.lambda_h;

  for (size_t i = 0; i < n; ++i)
    grad[i] = 0.0f;

  // Jo
  float Jo = 0.0f;
  for (const auto& obs : H.obs) {
    float r = x[obs.grid_index] - obs.value;
    Jo += obs.weight * r * r;
    grad[obs.grid_index] += 2.0f * obs.weight * r;
  }

  float Js = 0.0f;

  if (H.kappa > 0.0f) {
    float kappa_sq = H.kappa * H.kappa;
    // Perona-Malik: combined cost + gradient in single pass
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        float c = x[idx];
        float g_sum = 0.0f;

        if (i < nx-1) { float d = c - x[idx+1]; g_sum += pm_diffusivity(d*d, kappa_sq) * d; Js += kappa_sq * std::log(1.0f + d*d/kappa_sq); }
        if (j < ny-1) { float d = c - x[idx+nx]; g_sum += pm_diffusivity(d*d, kappa_sq) * d; Js += kappa_sq * std::log(1.0f + d*d/kappa_sq); }
        if (i > 0)    { float d = c - x[idx-1]; g_sum += pm_diffusivity(d*d, kappa_sq) * d; }
        if (j > 0)    { float d = c - x[idx-nx]; g_sum += pm_diffusivity(d*d, kappa_sq) * d; }

        grad[idx] += lh * 2.0f * g_sum;
      }
  } else {
    float* __restrict__ wxx = work;
    float* __restrict__ wyy = work + n;

    const bool has_az = !H.Wx.empty();
    const float* __restrict__ pWx = has_az ? H.Wx.data() : nullptr;
    const float* __restrict__ pWy = has_az ? H.Wy.data() : nullptr;

    // --- Pass 1: dxx, dyy, weights, Js2 cost, first-order damping ---
    float Js2 = 0.0f, Jd = 0.0f;

    for (size_t j = 0; j < ny; ++j) {
      for (size_t i = 0; i < nx; ++i) {
        const size_t idx = j * nx + i;
        const float c = x[idx];

        float left  = (i > 0)      ? x[idx-1] : c;
        float right = (i < nx - 1) ? x[idx+1] : c;
        float up    = (j > 0)      ? x[idx-nx] : c;
        float down  = (j < ny - 1) ? x[idx+nx] : c;

        float dxx = left - 2.0f * c + right;
        float dyy = up - 2.0f * c + down;

        float wx = has_az ? pWx[idx] : 1.0f;
        float wy = has_az ? pWy[idx] : 1.0f;

        Js2 += wx * dxx * dxx + wy * dyy * dyy;
        wxx[idx] = wx * dxx;
        wyy[idx] = wy * dyy;

        // First-order damping: gradient + cost
        float g1 = (c - left) + (c - right) + (c - up) + (c - down);
        if (i < nx - 1) { float d = c - right; Jd += d * d; }
        if (j < ny - 1) { float d = c - down;  Jd += d * d; }

        grad[idx] += 2.0f * lh * g1;
      }
    }

    Js = Js2 + Jd;

    // --- Pass 2: Adjoint gradient for second-order term ---
    for (size_t j = 0; j < ny; ++j) {
      for (size_t i = 0; i < nx; ++i) {
        const size_t idx = j * nx + i;

        float xx_l = (i > 0)      ? wxx[idx-1] : wxx[idx];
        float xx_r = (i < nx - 1) ? wxx[idx+1] : wxx[idx];
        float g_xx = xx_l - 2.0f * wxx[idx] + xx_r;

        float yy_u = (j > 0)      ? wyy[idx-nx] : wyy[idx];
        float yy_d = (j < ny - 1) ? wyy[idx+nx] : wyy[idx];
        float g_yy = yy_u - 2.0f * wyy[idx] + yy_d;

        grad[idx] += 2.0f * lh * (g_xx + g_yy);
      }
    }
  }

  return Jo + lh * Js;
}

// Legacy — kept for unit tests
inline void apply_laplacian(const float* x, float* Lx, size_t nx, size_t ny)
{
  for (size_t j = 0; j < ny; ++j)
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float c = x[idx], sum = 0.0f;
      if (i > 0)      sum += c - x[idx-1];
      if (i < nx - 1) sum += c - x[idx+1];
      if (j > 0)      sum += c - x[idx-nx];
      if (j < ny - 1) sum += c - x[idx+nx];
      Lx[idx] = sum;
    }
}

#endif // VARGRID_COST_FUNCTION_H