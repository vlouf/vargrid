#ifndef VARGRID_COST_FUNCTION_H
#define VARGRID_COST_FUNCTION_H

#include "observation_operator.h"
#include "config.h"
#include <vector>
#include <cmath>

// Cost function:
//   J(φ) = Jo + λH · Js
//
//   Jo = Σ w_i (φ[g_i] - d_i)²           (data fidelity)
//   Js = Σ_edges w_e (φ_j - φ_k)²         (8-connected smoothness)
//
// The 8-connected stencil uses all 8 neighbours:
//   Cardinal (left/right/up/down): weight 1.0
//   Diagonal (corners):            weight 1/√2  (distance is √2 × cell size)
//
// This makes the smoothing nearly isotropic, eliminating the grid-axis-aligned
// star artefacts that a 4-connected stencil produces on radial data.

// Diagonal edge weight: 1/√2
constexpr float DIAG_W = 0.70710678f;

inline float pm_diffusivity(float diff_sq, float kappa_sq) {
  return 1.0f / (1.0f + diff_sq / kappa_sq);
}

// === Cost-only evaluation (no gradient) ===
inline auto evaluate_cost(
    const float* __restrict__ x,
    const observation_operator& H,
    const vargrid_config& cfg
    ) -> float
{
  const size_t nx = H.grid_nx;
  const size_t ny = H.grid_ny;

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
        // Cardinal edges (right and down)
        if (i < nx-1) { float d = c - x[idx+1];  Js += kappa_sq * std::log(1.0f + d*d/kappa_sq); }
        if (j < ny-1) { float d = c - x[idx+nx]; Js += kappa_sq * std::log(1.0f + d*d/kappa_sq); }
        // Diagonal edges (down-right and down-left)
        if (i < nx-1 && j < ny-1) { float d = c - x[idx+nx+1]; Js += DIAG_W * kappa_sq * std::log(1.0f + d*d/kappa_sq); }
        if (i > 0    && j < ny-1) { float d = c - x[idx+nx-1]; Js += DIAG_W * kappa_sq * std::log(1.0f + d*d/kappa_sq); }
      }
  } else {
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        float c = x[idx];
        // Cardinal edges (right and down)
        if (i < nx-1) { float d = c - x[idx+1];  Js += d * d; }
        if (j < ny-1) { float d = c - x[idx+nx]; Js += d * d; }
        // Diagonal edges (down-right and down-left)
        if (i < nx-1 && j < ny-1) { float d = c - x[idx+nx+1]; Js += DIAG_W * d * d; }
        if (i > 0    && j < ny-1) { float d = c - x[idx+nx-1]; Js += DIAG_W * d * d; }
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
    float* __restrict__ work   // unused but kept for API compatibility
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
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        float c = x[idx];
        float g = 0.0f;

        // Cardinal neighbours
        if (i < nx-1) { float d = c-x[idx+1];  g += pm_diffusivity(d*d,kappa_sq)*d; Js += kappa_sq*std::log(1.0f+d*d/kappa_sq); }
        if (i > 0)    { float d = c-x[idx-1];  g += pm_diffusivity(d*d,kappa_sq)*d; }
        if (j < ny-1) { float d = c-x[idx+nx]; g += pm_diffusivity(d*d,kappa_sq)*d; Js += kappa_sq*std::log(1.0f+d*d/kappa_sq); }
        if (j > 0)    { float d = c-x[idx-nx]; g += pm_diffusivity(d*d,kappa_sq)*d; }

        // Diagonal neighbours (weight DIAG_W)
        if (i<nx-1 && j<ny-1) { float d=c-x[idx+nx+1]; g+=DIAG_W*pm_diffusivity(d*d,kappa_sq)*d; Js+=DIAG_W*kappa_sq*std::log(1.0f+d*d/kappa_sq); }
        if (i>0    && j<ny-1) { float d=c-x[idx+nx-1]; g+=DIAG_W*pm_diffusivity(d*d,kappa_sq)*d; Js+=DIAG_W*kappa_sq*std::log(1.0f+d*d/kappa_sq); }
        if (i<nx-1 && j>0)    { float d=c-x[idx-nx+1]; g+=DIAG_W*pm_diffusivity(d*d,kappa_sq)*d; }
        if (i>0    && j>0)    { float d=c-x[idx-nx-1]; g+=DIAG_W*pm_diffusivity(d*d,kappa_sq)*d; }

        grad[idx] += 2.0f * lh * g;
      }
  } else {
    // 8-connected isotropic Laplacian smoothness.
    // Cardinal: weight 1.0,  Diagonal: weight 1/√2
    // Gradient: ∇Js[k] = 2 Σ_neighbours w_e (φ[k] - φ[neighbour])
    // Cost:     Js = Σ_edges w_e (φ_j - φ_k)² (each edge counted once)

    for (size_t j = 0; j < ny; ++j) {
      for (size_t i = 0; i < nx; ++i) {
        const size_t idx = j * nx + i;
        const float c = x[idx];
        float g = 0.0f;

        // Cardinal neighbours (gradient from all 4, cost from right+down)
        if (i > 0)    { g += c - x[idx-1]; }
        if (i < nx-1) { float d = c - x[idx+1];  g += d; Js += d * d; }
        if (j > 0)    { g += c - x[idx-nx]; }
        if (j < ny-1) { float d = c - x[idx+nx]; g += d; Js += d * d; }

        // Diagonal neighbours (gradient from all 4, cost from down-right+down-left)
        if (i > 0    && j > 0)    { g += DIAG_W * (c - x[idx-nx-1]); }
        if (i < nx-1 && j > 0)    { g += DIAG_W * (c - x[idx-nx+1]); }
        if (i > 0    && j < ny-1) { float d = c - x[idx+nx-1]; g += DIAG_W * d; Js += DIAG_W * d * d; }
        if (i < nx-1 && j < ny-1) { float d = c - x[idx+nx+1]; g += DIAG_W * d; Js += DIAG_W * d * d; }

        grad[idx] += 2.0f * lh * g;
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