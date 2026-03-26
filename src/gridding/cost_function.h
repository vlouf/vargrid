#ifndef VARGRID_COST_FUNCTION_H
#define VARGRID_COST_FUNCTION_H

#include "observation_operator.h"
#include "config.h"
#include <vector>
#include <cmath>

// Cost function J(x) = Jo(x) + alpha * Js(x)
//
//   Jo(x) = sum_i w_i * (x[g_i] - y_i)^2           (observation fit)
//   Js(x) = sum_edges kappa^2 * log(1 + d^2/kappa^2) (Perona-Malik smoothness)
//
// where d = x[j] - x[k] is the difference across a grid edge.
//
// The Perona-Malik penalty behaves like:
//   - |d| << kappa: Js ~ d^2              (quadratic, strong smoothing)
//   - |d| >> kappa: Js ~ kappa^2 * log()  (sublinear, edges preserved)
//
// The edge threshold kappa controls the transition between smoothing
// (weak gradients) and preservation (strong gradients like storm edges).
// When kappa <= 0, the penalty falls back to the isotropic Laplacian.

// Apply the isotropic Laplacian operator L*x on a 2D grid.
// Uses 5-point stencil with Neumann (zero-gradient) boundary conditions.
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

// Perona-Malik diffusivity: g(s) = 1 / (1 + s / kappa^2)
// Controls edge sensitivity: strong gradients (s >> kappa^2) get
// low diffusivity (preserved), weak gradients get high diffusivity (smoothed).
inline float pm_diffusivity(float diff_sq, float kappa_sq) {
  return 1.0f / (1.0f + diff_sq / kappa_sq);
}

// Evaluate Perona-Malik smoothness cost and gradient.
//
// Cost:     Js = sum_edges kappa^2 * log(1 + d^2 / kappa^2)
// Gradient: dJs/dx[idx] = 2 * sum_neighbours g(d^2) * d
//
// Each edge is counted once for the cost (right and down directions),
// but both endpoints receive gradient contributions from all four neighbours.
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

      // Right neighbour
      if (i < nx - 1) {
        float d = center - x[idx + 1];
        float d_sq = d * d;
        g_sum += pm_diffusivity(d_sq, kappa_sq) * d;
        Js += kappa_sq * std::log(1.0f + d_sq / kappa_sq);
      }

      // Down neighbour
      if (j < ny - 1) {
        float d = center - x[idx + nx];
        float d_sq = d * d;
        g_sum += pm_diffusivity(d_sq, kappa_sq) * d;
        Js += kappa_sq * std::log(1.0f + d_sq / kappa_sq);
      }

      // Left neighbour (gradient only — cost counted from left cell)
      if (i > 0) {
        float d = center - x[idx - 1];
        g_sum += pm_diffusivity(d * d, kappa_sq) * d;
      }

      // Up neighbour (gradient only — cost counted from upper cell)
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
//   J(x) = Jo(x) + alpha * Js(x)
//
// When kappa > 0: Perona-Malik edge-preserving smoothness.
// When kappa <= 0: Isotropic Laplacian smoothness (original behaviour).
inline auto evaluate_gradient(
    const float* x,
    float* grad,
    const observation_operator& H,
    float alpha,
    float* work  // temporary buffer, size grid_nx * grid_ny
    ) -> float
{
  size_t n = H.grid_size();
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;

  for (size_t i = 0; i < n; ++i)
    grad[i] = 0.0f;

  // Observation term
  float Jo = 0.0f;
  for (auto& obs : H.obs) {
    float residual = x[obs.grid_index] - obs.value;
    Jo += obs.weight * residual * residual;
    grad[obs.grid_index] += 2.0f * obs.weight * residual;
  }

  // Smoothness term
  float Js = 0.0f;

  if (H.kappa > 0.0f) {
    // Perona-Malik edge-preserving
    Js = evaluate_pm_smoothness(x, work, nx, ny, H.kappa);
    for (size_t i = 0; i < n; ++i)
      grad[i] += alpha * work[i];
  } else {
    // Isotropic Laplacian fallback
    apply_laplacian(x, work, nx, ny);
    for (size_t i = 0; i < n; ++i) {
      Js += x[i] * work[i];
      grad[i] += 2.0f * alpha * work[i];
    }
  }

  return Jo + alpha * Js;
}

#endif // VARGRID_COST_FUNCTION_H