#include "vargrid.h"
#include "cost_function.h"

#include <cmath>
#include <algorithm>
#include <numeric>

using namespace bom;

// Compute a nearest-gate initial guess for faster CG convergence.
static void compute_initial_guess(
    float* x,
    const observation_operator& H,
    float background)
{
  size_t n = H.grid_size();

  std::vector<float> wsum(n, 0.0f);
  std::vector<float> wcount(n, 0.0f);

  for (auto& obs : H.obs) {
    wsum[obs.grid_index] += obs.weight * obs.value;
    wcount[obs.grid_index] += obs.weight;
  }

  for (size_t i = 0; i < n; ++i) {
    if (wcount[i] > 0.0f)
      x[i] = wsum[i] / wcount[i];
    else
      x[i] = background;
  }
}

// Mark unconstrained grid cells as nodata.
static void mask_unconstrained(
    float* x,
    const observation_operator& H)
{
  size_t n = H.grid_size();
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;

  constexpr int propagation_steps = 3;
  std::vector<bool> has_data(n, false);

  for (size_t i = 0; i < n; ++i)
    has_data[i] = (H.obs_count[i] > 0);

  std::vector<bool> next(n, false);
  for (int step = 0; step < propagation_steps; ++step) {
    for (size_t j = 0; j < ny; ++j) {
      for (size_t i = 0; i < nx; ++i) {
        size_t idx = j * nx + i;
        if (has_data[idx]) {
          next[idx] = true;
          continue;
        }
        bool neighbour_has_data = false;
        if (i > 0      && has_data[idx - 1])  neighbour_has_data = true;
        if (i < nx - 1 && has_data[idx + 1])  neighbour_has_data = true;
        if (j > 0      && has_data[idx - nx]) neighbour_has_data = true;
        if (j < ny - 1 && has_data[idx + nx]) neighbour_has_data = true;
        next[idx] = neighbour_has_data;
      }
    }
    std::swap(has_data, next);
  }

  for (size_t i = 0; i < n; ++i) {
    if (!has_data[i])
      x[i] = nodata;
  }
}

auto variational_grid(
      const observation_operator& H
    , const vargrid_config& cfg
    ) -> array2f
{
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;
  size_t n = H.grid_size();

  auto result = array2f{vec2z{nx, ny}};

  float bg = std::isnan(cfg.background) ? 0.0f : cfg.background;

  if (cfg.use_nearest_init) {
    compute_initial_guess(result.data(), H, bg);
  } else {
    for (size_t i = 0; i < n; ++i)
      result.data()[i] = bg;
  }

  if (!H.obs.empty()) {
    auto sr = solve_cg(result.data(), H, cfg);
    trace::log("  Variational grid: {} iterations, cost={:.4f}, residual={:.2e}, converged={}",
      sr.iterations, sr.final_cost, sr.final_residual, sr.converged ? "yes" : "no");
  } else {
    trace::warning("  Variational grid: no observations for this layer");
  }

  mask_unconstrained(result.data(), H);

  return result;
}

auto observation_density(const observation_operator& H) -> array2f
{
  auto result = array2f{vec2z{H.grid_nx, H.grid_ny}};
  for (size_t i = 0; i < H.grid_size(); ++i)
    result.data()[i] = static_cast<float>(H.obs_count[i]);
  return result;
}