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

// Mask cells far from any observation using BFS distance.
// Cells with no observation within mask_distance grid steps are set to NaN.
static void mask_by_distance(float* x, const observation_operator& H, float max_dist)
{
  size_t n = H.grid_size();
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;

  // BFS from cells that have observations
  std::vector<float> dist(n, std::numeric_limits<float>::max());
  std::vector<size_t> frontier;

  for (size_t i = 0; i < n; ++i) {
    if (H.obs_count[i] > 0) {
      dist[i] = 0.0f;
      frontier.push_back(i);
    }
  }

  float current_dist = 0.0f;
  while (!frontier.empty() && current_dist < max_dist) {
    current_dist += 1.0f;
    std::vector<size_t> next;
    for (auto idx : frontier) {
      size_t iy = idx / nx;
      size_t ix = idx % nx;
      auto try_cell = [&](size_t cx, size_t cy) {
        size_t cidx = cy * nx + cx;
        if (dist[cidx] > current_dist) {
          dist[cidx] = current_dist;
          next.push_back(cidx);
        }
      };
      if (ix > 0)       try_cell(ix - 1, iy);
      if (ix < nx - 1)  try_cell(ix + 1, iy);
      if (iy > 0)       try_cell(ix, iy - 1);
      if (iy < ny - 1)  try_cell(ix, iy + 1);
    }
    frontier = std::move(next);
  }

  for (size_t i = 0; i < n; ++i)
    if (dist[i] > max_dist)
      x[i] = nodata;
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
  }

  mask_by_distance(result.data(), H, cfg.mask_distance_cells);
  return result;
}

auto observation_density(const observation_operator& H) -> array2f
{
  auto result = array2f{vec2z{H.grid_nx, H.grid_ny}};
  for (size_t i = 0; i < H.grid_size(); ++i)
    result.data()[i] = static_cast<float>(H.obs_count[i]);
  return result;
}
