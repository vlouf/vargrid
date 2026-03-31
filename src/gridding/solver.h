#ifndef VARGRID_SOLVER_H
#define VARGRID_SOLVER_H

#include "observation_operator.h"
#include "cost_function.h"
#include "config.h"

#include <vector>
#include <cmath>
#include <algorithm>

struct solver_result {
  int iterations;
  float final_cost;
  float final_residual;
  bool converged;
};

// Nonlinear conjugate gradient solver (Polak-Ribière).
//
// Minimises J(x) = Jo(x) + λH · Js(x).
//
// Key optimisation: the Armijo line search uses evaluate_cost (no gradient),
// which is ~3× cheaper than evaluate_gradient. The full gradient is computed
// only once per CG iteration after the step is accepted.
inline auto solve_cg(
    float* x,
    const observation_operator& H,
    const vargrid_config& cfg
    ) -> solver_result
{
  size_t n = H.grid_size();

  std::vector<float> grad(n);
  std::vector<float> grad_prev(n);
  std::vector<float> direction(n);
  std::vector<float> x_trial(n);
  std::vector<float> work(2 * n);

  // Initial gradient
  float cost = evaluate_gradient(x, grad.data(), H, cfg, work.data());

  float grad_norm_sq = 0.0f;
  for (size_t i = 0; i < n; ++i) {
    direction[i] = -grad[i];
    grad_norm_sq += grad[i] * grad[i];
  }

  float initial_grad_norm = std::sqrt(grad_norm_sq);
  if (initial_grad_norm < 1e-12f) {
    return {0, cost, 0.0f, true};
  }

  solver_result result;
  result.converged = false;

  for (int iter = 0; iter < cfg.max_iterations; ++iter) {

    // Directional derivative
    float dir_dot_grad = 0.0f;
    for (size_t i = 0; i < n; ++i)
      dir_dot_grad += direction[i] * grad[i];

    // Reset to steepest descent if not a descent direction
    if (dir_dot_grad >= 0.0f) {
      for (size_t i = 0; i < n; ++i)
        direction[i] = -grad[i];
      dir_dot_grad = -grad_norm_sq;
    }

    // Armijo backtracking line search — COST ONLY (no gradient)
    float step = 1.0f;
    constexpr float c1 = 1e-4f;
    constexpr float shrink = 0.5f;
    constexpr int max_ls = 20;

    float new_cost = cost;
    for (int ls = 0; ls < max_ls; ++ls) {
      for (size_t i = 0; i < n; ++i)
        x_trial[i] = x[i] + step * direction[i];

      // Cost-only evaluation — no gradient computation
      new_cost = evaluate_cost(x_trial.data(), H, cfg);

      if (new_cost <= cost + c1 * step * dir_dot_grad)
        break;

      step *= shrink;
    }

    // Accept the step
    for (size_t i = 0; i < n; ++i)
      x[i] = x_trial[i];
    cost = new_cost;

    // Compute gradient at the accepted point (once per CG iteration)
    evaluate_gradient(x, grad.data(), H, cfg, work.data());

    // Gradient norm
    float new_grad_norm_sq = 0.0f;
    for (size_t i = 0; i < n; ++i)
      new_grad_norm_sq += grad[i] * grad[i];

    float rel_residual = std::sqrt(new_grad_norm_sq) / initial_grad_norm;

    if (rel_residual < cfg.tolerance) {
      result.iterations = iter + 1;
      result.final_cost = cost;
      result.final_residual = rel_residual;
      result.converged = true;
      return result;
    }

    // Polak-Ribière beta
    float beta_num = 0.0f;
    if (iter > 0) {
      for (size_t i = 0; i < n; ++i)
        beta_num += grad[i] * (grad[i] - grad_prev[i]);
    }
    float beta = (iter == 0) ? 0.0f : std::max(0.0f, beta_num / grad_norm_sq);

    // Update search direction
    for (size_t i = 0; i < n; ++i)
      direction[i] = -grad[i] + beta * direction[i];

    std::copy(grad.begin(), grad.end(), grad_prev.begin());
    grad_norm_sq = new_grad_norm_sq;
  }

  float final_norm = std::sqrt(grad_norm_sq);
  result.iterations = cfg.max_iterations;
  result.final_cost = cost;
  result.final_residual = final_norm / initial_grad_norm;
  result.converged = false;
  return result;
}

#endif // VARGRID_SOLVER_H