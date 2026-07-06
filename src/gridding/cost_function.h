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

// Azimuthally-varying smoothness weights (Brook et al. 2022 Eq. 3).
// Wx/Wy are per-cell; an edge uses the mean of its two endpoints' weights so
// the weighting is symmetric (required for the gradient to be the exact
// derivative of the cost). Diagonal edges mix Wx and Wy evenly. Empty (or
// wrongly sized) Wx/Wy means isotropic weighting.
struct edge_weights {
  const float* Wx = nullptr;
  const float* Wy = nullptr;

  explicit edge_weights(const observation_operator& H) {
    if (H.Wx.size() == H.grid_size() && H.Wy.size() == H.grid_size()) {
      Wx = H.Wx.data();
      Wy = H.Wy.data();
    }
  }

  float x(size_t a, size_t b) const { return Wx ? 0.5f * (Wx[a] + Wx[b]) : 1.0f; }
  float y(size_t a, size_t b) const { return Wy ? 0.5f * (Wy[a] + Wy[b]) : 1.0f; }
  float d(size_t a, size_t b) const {
    return Wx ? 0.25f * (Wx[a] + Wy[a] + Wx[b] + Wy[b]) : 1.0f;
  }
};

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
  const edge_weights W{H};

  if (H.kappa > 0.0f) {
    float kappa_sq = H.kappa * H.kappa;
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        float c = x[idx];
        // Cardinal edges (right and down)
        if (i < nx-1) { float d = c - x[idx+1];  Js += W.x(idx, idx+1)  * kappa_sq * std::log(1.0f + d*d/kappa_sq); }
        if (j < ny-1) { float d = c - x[idx+nx]; Js += W.y(idx, idx+nx) * kappa_sq * std::log(1.0f + d*d/kappa_sq); }
        // Diagonal edges (down-right and down-left)
        if (i < nx-1 && j < ny-1) { float d = c - x[idx+nx+1]; Js += DIAG_W * W.d(idx, idx+nx+1) * kappa_sq * std::log(1.0f + d*d/kappa_sq); }
        if (i > 0    && j < ny-1) { float d = c - x[idx+nx-1]; Js += DIAG_W * W.d(idx, idx+nx-1) * kappa_sq * std::log(1.0f + d*d/kappa_sq); }
      }
  } else {
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        float c = x[idx];
        // Cardinal edges (right and down)
        if (i < nx-1) { float d = c - x[idx+1];  Js += W.x(idx, idx+1)  * d * d; }
        if (j < ny-1) { float d = c - x[idx+nx]; Js += W.y(idx, idx+nx) * d * d; }
        // Diagonal edges (down-right and down-left)
        if (i < nx-1 && j < ny-1) { float d = c - x[idx+nx+1]; Js += DIAG_W * W.d(idx, idx+nx+1) * d * d; }
        if (i > 0    && j < ny-1) { float d = c - x[idx+nx-1]; Js += DIAG_W * W.d(idx, idx+nx-1) * d * d; }
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
  const edge_weights W{H};

  if (H.kappa > 0.0f) {
    float kappa_sq = H.kappa * H.kappa;
    for (size_t j = 0; j < ny; ++j)
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        float c = x[idx];
        float g = 0.0f;

        // Cardinal neighbours
        if (i < nx-1) { float d = c-x[idx+1];  float w = W.x(idx, idx+1);  g += w*pm_diffusivity(d*d,kappa_sq)*d; Js += w*kappa_sq*std::log(1.0f+d*d/kappa_sq); }
        if (i > 0)    { float d = c-x[idx-1];  g += W.x(idx, idx-1)*pm_diffusivity(d*d,kappa_sq)*d; }
        if (j < ny-1) { float d = c-x[idx+nx]; float w = W.y(idx, idx+nx); g += w*pm_diffusivity(d*d,kappa_sq)*d; Js += w*kappa_sq*std::log(1.0f+d*d/kappa_sq); }
        if (j > 0)    { float d = c-x[idx-nx]; g += W.y(idx, idx-nx)*pm_diffusivity(d*d,kappa_sq)*d; }

        // Diagonal neighbours (weight DIAG_W)
        if (i<nx-1 && j<ny-1) { float d=c-x[idx+nx+1]; float w=W.d(idx,idx+nx+1); g+=DIAG_W*w*pm_diffusivity(d*d,kappa_sq)*d; Js+=DIAG_W*w*kappa_sq*std::log(1.0f+d*d/kappa_sq); }
        if (i>0    && j<ny-1) { float d=c-x[idx+nx-1]; float w=W.d(idx,idx+nx-1); g+=DIAG_W*w*pm_diffusivity(d*d,kappa_sq)*d; Js+=DIAG_W*w*kappa_sq*std::log(1.0f+d*d/kappa_sq); }
        if (i<nx-1 && j>0)    { float d=c-x[idx-nx+1]; g+=DIAG_W*W.d(idx,idx-nx+1)*pm_diffusivity(d*d,kappa_sq)*d; }
        if (i>0    && j>0)    { float d=c-x[idx-nx-1]; g+=DIAG_W*W.d(idx,idx-nx-1)*pm_diffusivity(d*d,kappa_sq)*d; }

        grad[idx] += 2.0f * lh * g;
      }
  } else {
    // 8-connected Laplacian smoothness with azimuthal weights (Brook Eq. 2-3).
    // Cardinal: weight Wx/Wy,  Diagonal: 1/√2 × mixed weight
    // Gradient: ∇Js[k] = 2 Σ_neighbours w_e (φ[k] - φ[neighbour])
    // Cost:     Js = Σ_edges w_e (φ_j - φ_k)² (each edge counted once)

    for (size_t j = 0; j < ny; ++j) {
      for (size_t i = 0; i < nx; ++i) {
        const size_t idx = j * nx + i;
        const float c = x[idx];
        float g = 0.0f;

        // Cardinal neighbours (gradient from all 4, cost from right+down)
        if (i > 0)    { g += W.x(idx, idx-1) * (c - x[idx-1]); }
        if (i < nx-1) { float d = c - x[idx+1];  float w = W.x(idx, idx+1);  g += w * d; Js += w * d * d; }
        if (j > 0)    { g += W.y(idx, idx-nx) * (c - x[idx-nx]); }
        if (j < ny-1) { float d = c - x[idx+nx]; float w = W.y(idx, idx+nx); g += w * d; Js += w * d * d; }

        // Diagonal neighbours (gradient from all 4, cost from down-right+down-left)
        if (i > 0    && j > 0)    { g += DIAG_W * W.d(idx, idx-nx-1) * (c - x[idx-nx-1]); }
        if (i < nx-1 && j > 0)    { g += DIAG_W * W.d(idx, idx-nx+1) * (c - x[idx-nx+1]); }
        if (i > 0    && j < ny-1) { float d = c - x[idx+nx-1]; float w = W.d(idx, idx+nx-1); g += DIAG_W * w * d; Js += DIAG_W * w * d * d; }
        if (i < nx-1 && j < ny-1) { float d = c - x[idx+nx+1]; float w = W.d(idx, idx+nx+1); g += DIAG_W * w * d; Js += DIAG_W * w * d * d; }

        grad[idx] += 2.0f * lh * g;
      }
    }
  }

  return Jo + lh * Js;
}

#endif // VARGRID_COST_FUNCTION_H