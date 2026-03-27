#ifndef VARGRID_COST_FUNCTION_H
#define VARGRID_COST_FUNCTION_H

#include "observation_operator.h"
#include "config.h"
#include <vector>
#include <cmath>

// Cost function following Brook et al. (2022):
//
//   J(φ) = ||d - Rφ||² + λH (||φ_xx||² + ||φ_yy||²)
//
// where R is the bilinear interpolation operator (encoded in H.obs),
// φ_xx and φ_yy are second-order centered finite differences with
// Neumann boundary conditions, and λH is the horizontal smoothing weight.
//
// Optionally, Perona-Malik edge-preserving diffusion can be applied
// to the first-order terms instead (when kappa > 0).

// Apply the isotropic Laplacian (first-order graph Laplacian).
// Kept for backward compatibility and unit tests.
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

// Second-order smoothness: ||φ_xx||² + ||φ_yy||²
// Uses centered second derivatives with Neumann BCs.
//
// φ_xx[j][i] = φ[j][i-1] - 2φ[j][i] + φ[j][i+1]
// At boundaries: Neumann => φ[-1] = φ[0], so φ_xx = -φ[0] + φ[1] at i=0
//
// Cost:  Js = sum (φ_xx)² + sum (φ_yy)²
// Gradient: ∇Js = 2 * (Dxx^T Dxx + Dyy^T Dyy) φ
//         = 2 * (L4x + L4y) φ  where L4 is the biharmonic-like operator
//
// For efficiency, we compute the second derivatives, then the cost and
// gradient in a single pass using the adjoint relationship.

// Compute second derivative in x: φ_xx
inline void compute_dxx(const float* x, float* dxx, size_t nx, size_t ny)
{
  for (size_t j = 0; j < ny; ++j) {
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;
      float left   = (i > 0)      ? x[idx - 1] : x[idx];  // Neumann BC
      float right  = (i < nx - 1) ? x[idx + 1] : x[idx];  // Neumann BC
      dxx[idx] = left - 2.0f * x[idx] + right;
    }
  }
}

// Compute second derivative in y: φ_yy
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

// Evaluate second-order smoothness cost and gradient.
// Js = ||φ_xx||² + ||φ_yy||²
// ∇Js = 2 * (Dxx^T * φ_xx + Dyy^T * φ_yy)
// The adjoint of Dxx is Dxx itself (symmetric stencil), so:
// ∇Js[idx] = 2 * (dxx_of_dxx[idx] + dyy_of_dyy[idx])
// But more efficiently: ∇Js = 2 * Dxx^T(φ_xx) + 2 * Dyy^T(φ_yy)
// Since Dxx is self-adjoint with Neumann BC, Dxx^T = Dxx.
inline auto evaluate_second_order_smoothness(
    const float* x,
    float* grad_s,     // output gradient of Js
    float* work1,      // scratch buffer (size nx*ny)
    float* work2,      // scratch buffer (size nx*ny)
    size_t nx,
    size_t ny
    ) -> float
{
  size_t n = nx * ny;

  // Compute φ_xx and φ_yy
  compute_dxx(x, work1, nx, ny);  // work1 = φ_xx
  compute_dyy(x, work2, nx, ny);  // work2 = φ_yy

  // Cost = ||φ_xx||² + ||φ_yy||²
  float Js = 0.0f;
  for (size_t i = 0; i < n; ++i)
    Js += work1[i] * work1[i] + work2[i] * work2[i];

  // Gradient = 2 * (Dxx^T φ_xx + Dyy^T φ_yy)
  // Since Dxx is self-adjoint: Dxx^T(φ_xx) = Dxx(φ_xx) = φ_xxxx ... no.
  // Actually: ∇_φ ||Dxx φ||² = 2 Dxx^T (Dxx φ)
  // The adjoint of the centered second difference is itself (it's symmetric).
  // So: grad_s = 2 * (Dxx(work1) + Dyy(work2))
  // But that gives fourth-order derivatives. Let's be more careful.
  //
  // Actually for a linear operator A:
  //   ∇_x ||Ax||² = 2 A^T A x
  // We already have Ax in work1/work2. We need A^T * (Ax).
  // For the centered second difference with Neumann BC, A^T = A.
  // So we need Dxx(work1) + Dyy(work2).

  // Compute Dxx^T(φ_xx) = Dxx(φ_xx) and Dyy^T(φ_yy) = Dyy(φ_yy)
  // We need two more buffers... but we can reuse: compute in-place.
  // Actually, we can compute the gradient directly by applying Dxx to work1
  // and Dyy to work2, summing into grad_s.

  // grad_s = 2 * (Dxx(work1) + Dyy(work2))
  for (size_t j = 0; j < ny; ++j) {
    for (size_t i = 0; i < nx; ++i) {
      size_t idx = j * nx + i;

      // Dxx applied to work1 (= φ_xx)
      float dxx_left  = (i > 0)      ? work1[idx - 1] : work1[idx];
      float dxx_right = (i < nx - 1) ? work1[idx + 1] : work1[idx];
      float grad_xx = dxx_left - 2.0f * work1[idx] + dxx_right;

      // Dyy applied to work2 (= φ_yy)
      float dyy_up   = (j > 0)      ? work2[idx - nx] : work2[idx];
      float dyy_down = (j < ny - 1) ? work2[idx + nx] : work2[idx];
      float grad_yy = dyy_up - 2.0f * work2[idx] + dyy_down;

      grad_s[idx] = 2.0f * (grad_xx + grad_yy);
    }
  }

  return Js;
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
// J(φ) = Jo(φ) + λH * Js(φ)
//
// Jo = sum_i w_i * (R*φ[i] - d_i)²   (bilinear R encoded in H.obs)
// Js = ||φ_xx||² + ||φ_yy||²         (second-order, Brook et al.)
//   or Perona-Malik (when kappa > 0)
inline auto evaluate_gradient(
    const float* x,
    float* grad,
    const observation_operator& H,
    float lambda_h,
    float* work  // temporary buffer, size >= 2 * grid_nx * grid_ny
    ) -> float
{
  size_t n = H.grid_size();
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;

  for (size_t i = 0; i < n; ++i)
    grad[i] = 0.0f;

  // Observation term: Jo = sum w_i (φ[g_i] - y_i)², ∇Jo[k] = 2 sum w_i (φ[k] - y_i)
  // With bilinear R, each observation already has the bilinear weight folded
  // into obs.weight, and obs.grid_index points to one of the 4 cells.
  // The cost and gradient accumulation is the same as before.
  float Jo = 0.0f;
  for (auto& obs : H.obs) {
    float residual = x[obs.grid_index] - obs.value;
    Jo += obs.weight * residual * residual;
    grad[obs.grid_index] += 2.0f * obs.weight * residual;
  }

  // Smoothness term
  float Js = 0.0f;

  if (H.kappa > 0.0f) {
    // Perona-Malik first-order edge-preserving
    Js = evaluate_pm_smoothness(x, work, nx, ny, H.kappa);
    for (size_t i = 0; i < n; ++i)
      grad[i] += lambda_h * work[i];
  } else {
    // Second-order smoothness (Brook et al. 2022)
    // Js = ||φ_xx||² + ||φ_yy||²
    // ∇Js = 2 * (Dxx^T Dxx φ + Dyy^T Dyy φ)
    // Uses work[0..n) for φ_xx and work[n..2n) for φ_yy

    compute_dxx(x, work, nx, ny);          // work[0..n) = φ_xx
    compute_dyy(x, work + n, nx, ny);      // work[n..2n) = φ_yy

    Js = 0.0f;
    for (size_t i = 0; i < n; ++i)
      Js += work[i] * work[i] + work[n + i] * work[n + i];

    // Gradient: 2 * (Dxx(φ_xx) + Dyy(φ_yy))
    // Dxx is self-adjoint with Neumann BC
    for (size_t j = 0; j < ny; ++j) {
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;

        float xx_l = (i > 0)      ? work[idx - 1] : work[idx];
        float xx_r = (i < nx - 1) ? work[idx + 1] : work[idx];
        float g_xx = xx_l - 2.0f * work[idx] + xx_r;

        float yy_u = (j > 0)      ? work[n + idx - nx] : work[n + idx];
        float yy_d = (j < ny - 1) ? work[n + idx + nx] : work[n + idx];
        float g_yy = yy_u - 2.0f * work[n + idx] + yy_d;

        grad[idx] += 2.0f * lambda_h * (g_xx + g_yy);
      }
    }
  }

  return Jo + lambda_h * Js;
}

#endif // VARGRID_COST_FUNCTION_H