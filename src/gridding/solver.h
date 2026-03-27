#ifndef VARGRID_SOLVER_H
#define VARGRID_SOLVER_H

#include "observation_operator.h"
#include "cost_function.h"
#include "config.h"

#include <vector>
#include <cmath>
#include <algorithm>

// Solver result with convergence diagnostics.
struct solver_result {
  int iterations;       // Number of CG iterations performed
  float final_cost;     // Final cost function value
  float final_residual; // Final gradient norm (relative)
  bool converged;       // Whether tolerance was reached
};

// Nonlinear conjugate gradient solver (Polak-Ribière variant).
//
// Minimises J(x) = Jo(x) + alpha * Js(x) where:
//   Jo = weighted sum of squared residuals to observations
//   Js = Laplacian smoothness penalty
//
// The problem is quadratic in x (both terms are quadratic), so CG
// converges in at most N iterations. In practice, 20-50 iterations
// suffice for a 301x301 grid.
//
// Uses backtracking line search along the CG direction to ensure
// sufficient decrease (Armijo condition).
inline auto solve_cg(
    float* x,                       // Initial guess / solution (modified in-place)
    const observation_operator& H,  // Observation operator
    const vargrid_config& cfg       // Solver configuration
    ) -> solver_result
{
  size_t n = H.grid_size();

  // Allocate working memory
  std::vector<float> grad(n);
  std::vector<float> grad_prev(n);
  std::vector<float> direction(n);
  std::vector<float> x_trial(n);
  std::vector<float> work(2 * n);  // scratch buffer for cost function (needs 2*n for second-order)

  // Evaluate initial gradient
  float cost = evaluate_gradient(x, grad.data(), H, cfg.lambda_h, work.data());

  // Initial search direction = -gradient (steepest descent)
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

    // Backtracking line search (Armijo condition)
    // Find step size t such that J(x + t*d) <= J(x) + c1*t*∇J·d
    float dir_dot_grad = 0.0f;
    for (size_t i = 0; i < n; ++i)
      dir_dot_grad += direction[i] * grad[i];

    // If search direction is not a descent direction, reset to steepest descent
    if (dir_dot_grad >= 0.0f) {
      for (size_t i = 0; i < n; ++i)
        direction[i] = -grad[i];
      dir_dot_grad = -grad_norm_sq;
    }

    float step = 1.0f;
    constexpr float c1 = 1e-4f;      // Armijo constant
    constexpr float shrink = 0.5f;    // Step reduction factor
    constexpr int max_ls = 20;        // Max line search iterations

    float new_cost = cost;
    for (int ls = 0; ls < max_ls; ++ls) {
      for (size_t i = 0; i < n; ++i)
        x_trial[i] = x[i] + step * direction[i];

      new_cost = evaluate_gradient(x_trial.data(), grad.data(), H, cfg.lambda_h, work.data());

      if (new_cost <= cost + c1 * step * dir_dot_grad)
        break;

      step *= shrink;
    }

    // Accept the step
    for (size_t i = 0; i < n; ++i)
      x[i] = x_trial[i];
    cost = new_cost;

    // Compute new gradient norm
    float new_grad_norm_sq = 0.0f;
    for (size_t i = 0; i < n; ++i)
      new_grad_norm_sq += grad[i] * grad[i];

    float rel_residual = std::sqrt(new_grad_norm_sq) / initial_grad_norm;

    // Check convergence
    if (rel_residual < cfg.tolerance) {
      result.iterations = iter + 1;
      result.final_cost = cost;
      result.final_residual = rel_residual;
      result.converged = true;
      return result;
    }

    // Polak-Ribière beta: use max(beta, 0) to ensure descent
    // beta = (∇J_new · (∇J_new - ∇J_old)) / (∇J_old · ∇J_old)
    float beta_num = 0.0f;
    if (iter > 0) {
      for (size_t i = 0; i < n; ++i)
        beta_num += grad[i] * (grad[i] - grad_prev[i]);
    }
    float beta = (iter == 0) ? 0.0f : std::max(0.0f, beta_num / grad_norm_sq);

    // Update search direction
    for (size_t i = 0; i < n; ++i)
      direction[i] = -grad[i] + beta * direction[i];

    // Store gradient for next iteration's beta computation
    std::copy(grad.begin(), grad.end(), grad_prev.begin());
    grad_norm_sq = new_grad_norm_sq;
  }

  // Did not converge within max_iterations
  float final_norm = std::sqrt(grad_norm_sq);
  result.iterations = cfg.max_iterations;
  result.final_cost = cost;
  result.final_residual = final_norm / initial_grad_norm;
  result.converged = false;
  return result;
}

#endif // VARGRID_SOLVER_H
