#include "variational.h"
#include "cost_function.h"

static void compute_initial_guess(float* x, const observation_operator& H, float background)
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

// Mask cells that are very far from any observation as nodata.
// The background constraint JB handles smoothly decaying toward background
// in data voids, but cells beyond a large distance should be NaN in output.
// We use the precomputed dist_to_obs from the observation operator.
static void mask_distant_cells(float* x, const observation_operator& H, float max_dist_cells)
{
  size_t n = H.grid_size();

  if (H.dist_to_obs.empty()) return;

  for (size_t i = 0; i < n; ++i) {
    if (H.dist_to_obs[i] > max_dist_cells)
      x[i] = nodata;
  }
}

auto variational_grid(
      const observation_operator& H
    , const vargrid_config& cfg
    , size_t total_tasks
    , std::atomic<size_t>& completed_tasks
    ) -> array2f
{
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;
  size_t n = H.grid_size();

  auto result = array2f{vec2z{nx, ny}};
  float bg = std::isnan(cfg.background) ? 0.0f : cfg.background;

  if (cfg.use_nearest_init)
    compute_initial_guess(result.data(), H, bg);
  else
    for (size_t i = 0; i < n; ++i) result.data()[i] = bg;

  if (!H.obs.empty()) {
    auto sr = solve_cg(result.data(), H, cfg);

    auto done = completed_tasks.fetch_add(1) + 1;
    int pct = static_cast<int>(100.0 * done / total_tasks);
    trace::log("  [{:3d}%] Variational: {} iter, cost={:.1f}, residual={:.2e}, {}",
      pct, sr.iterations, sr.final_cost, sr.final_residual,
      sr.converged ? "converged" : "not converged");
  } else {
    completed_tasks.fetch_add(1);
    trace::warning("  Variational grid: no observations for this layer");
  }

  // Mask cells far from observations as nodata.
  // Use 2× the background cutoff as the mask threshold — cells beyond this
  // are fully dominated by the background constraint anyway.
  mask_distant_cells(result.data(), H, 2.0f * cfg.bg_cutoff_cells);
  return result;
}

auto observation_density(const observation_operator& H) -> array2f
{
  auto result = array2f{vec2z{H.grid_nx, H.grid_ny}};
  for (size_t i = 0; i < H.grid_size(); ++i)
    result.data()[i] = static_cast<float>(H.obs_count[i]);
  return result;
}
