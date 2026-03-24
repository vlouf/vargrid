#ifndef VARGRID_COST_FUNCTION_H
#define VARGRID_COST_FUNCTION_H

#include "observation_operator.h"
#include "config.h"
#include <vector>
#include <cmath>

// Evaluate the cost function J(x) = Jo(x) + alpha * Js(x)
//
//   Jo(x) = sum_i w_i * (x[g_i] - y_i)^2    (observation fit)
//   Js(x) = sum_{j,k} (x[j] - x[k])^2        (Laplacian smoothness)
//
// Also computes the gradient ∇J(x) needed by the CG solver.

// Apply the Laplacian operator L*x on a 2D grid.
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
      int count = 0;

      // Neumann BC: at boundaries, the "missing" neighbour equals the center
      // so (center - neighbour) = 0, contributing nothing.
      if (i > 0)      { sum += center - x[idx - 1];  count++; }
      if (i < nx - 1) { sum += center - x[idx + 1];  count++; }
      if (j > 0)      { sum += center - x[idx - nx]; count++; }
      if (j < ny - 1) { sum += center - x[idx + nx]; count++; }

      Lx[idx] = sum;
    }
  }
}

// Compute gradient of J(x):
//   ∇J = ∇Jo + alpha * ∇Js
//
// ∇Jo[k] = 2 * sum_{obs i at grid k} w_i * (x[k] - y_i)
// ∇Js = 2 * L * x  (where L is the graph Laplacian)
//
// Returns the cost function value J(x) as well.
inline auto evaluate_gradient(
    const float* x,
    float* grad,
    const observation_operator& H,
    float alpha,
    float* work_laplacian  // temporary buffer, size grid_nx * grid_ny
    ) -> float
{
  size_t n = H.grid_size();
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;

  // Zero the gradient
  for (size_t i = 0; i < n; ++i)
    grad[i] = 0.0f;

  // Observation term: Jo and ∇Jo
  float Jo = 0.0f;
  for (auto& obs : H.obs) {
    float residual = x[obs.grid_index] - obs.value;
    Jo += obs.weight * residual * residual;
    grad[obs.grid_index] += 2.0f * obs.weight * residual;
  }

  // Smoothness term: Js and ∇Js
  apply_laplacian(x, work_laplacian, nx, ny);
  float Js = 0.0f;
  for (size_t i = 0; i < n; ++i) {
    Js += x[i] * work_laplacian[i];
    grad[i] += 2.0f * alpha * work_laplacian[i];
  }

  return Jo + alpha * Js;
}

#endif // VARGRID_COST_FUNCTION_H
