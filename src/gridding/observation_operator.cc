#include "observation_operator.h"
#include "config.h"
#include "../types.h"

static auto compute_weight(
      float slant_range
    , float alt_dist
    , float max_alt_diff
    , const vargrid_config& cfg
    ) -> float
{
  float effective_range = std::max(slant_range, cfg.ref_range);
  float w_beam = std::pow(cfg.ref_range / effective_range, cfg.beam_power);

  float w_alt = 1.0f - std::fabs(alt_dist) / max_alt_diff;
  if (w_alt < 0.0f) w_alt = 0.0f;

  float w = w_beam * w_alt;
  return std::max(w, cfg.min_weight);
}

auto precompute_grid_bearings(
      latlonalt const& radar_location
    , array2<latlon> const& latlons
    ) -> grid_bearings
{
  grid_bearings gb;
  auto ny = latlons.extents().y;
  auto nx = latlons.extents().x;
  gb.nx = nx;
  gb.ny = ny;

  gb.cell_bearing_deg.resize(ny * nx);
  gb.cell_range.resize(ny * nx);
  float max_range = 0;

  for (size_t y = 0; y < ny; ++y) {
    for (size_t x = 0; x < nx; ++x) {
      auto br = wgs84.latlon_to_bearing_range(radar_location, latlons[y][x]);
      size_t idx = y * nx + x;
      float b = br.first.normalize().degrees();
      if (b < 0) b += 360.0f;
      gb.cell_bearing_deg[idx] = b;
      gb.cell_range[idx] = br.second;
      if (br.second > max_range) max_range = br.second;
    }
  }

  // Build binned lookup table
  gb.n_bearing_bins = 720;
  gb.bearing_bin_size = 360.0f / gb.n_bearing_bins;
  gb.range_bin_size = 500.0f;
  gb.n_range_bins = static_cast<size_t>(max_range / gb.range_bin_size) + 2;
  gb.max_range = max_range;

  size_t total_bins = gb.n_bearing_bins * gb.n_range_bins;
  gb.bin_to_grid.assign(total_bins, static_cast<size_t>(-1));
  std::vector<float> bin_dist(total_bins, std::numeric_limits<float>::max());

  for (size_t idx = 0; idx < ny * nx; ++idx) {
    int ib = static_cast<int>(gb.cell_bearing_deg[idx] / gb.bearing_bin_size);
    if (ib < 0) ib = 0;
    if (ib >= static_cast<int>(gb.n_bearing_bins)) ib = gb.n_bearing_bins - 1;

    int ir = static_cast<int>(gb.cell_range[idx] / gb.range_bin_size);
    if (ir < 0) ir = 0;
    if (ir >= static_cast<int>(gb.n_range_bins)) ir = gb.n_range_bins - 1;

    size_t bin_idx = static_cast<size_t>(ib) * gb.n_range_bins + static_cast<size_t>(ir);

    float bearing_center = (ib + 0.5f) * gb.bearing_bin_size;
    float range_center = (ir + 0.5f) * gb.range_bin_size;
    float db = gb.cell_bearing_deg[idx] - bearing_center;
    float dr = gb.cell_range[idx] - range_center;
    float dist = db * db + (dr / gb.range_bin_size) * (dr / gb.range_bin_size);

    if (dist < bin_dist[bin_idx]) {
      bin_dist[bin_idx] = dist;
      gb.bin_to_grid[bin_idx] = idx;
    }
  }

  trace::debug("Grid bearing lookup: {} bearing bins x {} range bins, max range={:.0f}m",
    gb.n_bearing_bins, gb.n_range_bins, max_range);

  return gb;
}

static auto lookup_nearest(
      grid_bearings const& gb
    , float bearing_deg
    , float ground_range
    ) -> size_t
{
  if (ground_range > gb.max_range || ground_range < 0)
    return static_cast<size_t>(-1);

  if (bearing_deg < 0) bearing_deg += 360.0f;
  if (bearing_deg >= 360.0f) bearing_deg -= 360.0f;

  int ib = static_cast<int>(bearing_deg / gb.bearing_bin_size);
  if (ib < 0) ib = 0;
  if (ib >= static_cast<int>(gb.n_bearing_bins)) ib = gb.n_bearing_bins - 1;

  int ir = static_cast<int>(ground_range / gb.range_bin_size);
  if (ir < 0) ir = 0;
  if (ir >= static_cast<int>(gb.n_range_bins)) ir = gb.n_range_bins - 1;

  size_t bin_idx = static_cast<size_t>(ib) * gb.n_range_bins + static_cast<size_t>(ir);
  return gb.bin_to_grid[bin_idx];
}

// Compute squared physical distance between a gate (bearing_deg, ground_range)
// and a grid cell, in meters. Uses bearing difference × range as cross-track.
static inline float gate_cell_dist_sq(
      float gate_bearing_deg
    , float gate_range
    , float cell_bearing_deg
    , float cell_range
    )
{
  float db = gate_bearing_deg - cell_bearing_deg;
  if (db > 180.0f) db -= 360.0f;
  if (db < -180.0f) db += 360.0f;

  // Cross-track distance: bearing difference in radians × average range
  float avg_range = 0.5f * (gate_range + cell_range);
  float cross = db * (M_PI / 180.0f) * avg_range;

  // Along-track distance: range difference
  float along = gate_range - cell_range;

  return cross * cross + along * along;
}

// Add observation contributions using inverse-distance weighting among
// the nearest cell and its grid-axis neighbours. This properly spreads
// each gate's contribution to surrounding cells without making any
// assumptions about how bearing/range maps to grid x/y axes.
//
// We consider the nearest cell and its 4 direct neighbours (up/down/left/right).
// The gate's contribution is distributed proportionally to 1/d² where d is the
// physical distance from the gate to each cell.
static void add_interpolated_observation(
      observation_operator& H
    , grid_bearings const& gb
    , size_t nearest_idx
    , float gate_bearing_deg
    , float gate_range
    , float value
    , float obs_weight
    , float alt_dist
    )
{
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;
  size_t iy = nearest_idx / nx;
  size_t ix = nearest_idx % nx;

  // Collect the nearest cell and its 4 direct grid neighbours
  struct candidate { size_t idx; float inv_dist; };
  candidate candidates[5];
  int ncand = 0;

  auto add_candidate = [&](size_t cx, size_t cy) {
    size_t cidx = cy * nx + cx;
    float d2 = gate_cell_dist_sq(gate_bearing_deg, gate_range,
                                  gb.cell_bearing_deg[cidx], gb.cell_range[cidx]);
    // Avoid division by zero: if gate lands exactly on a cell, give it all the weight
    if (d2 < 1.0f) d2 = 1.0f;
    candidates[ncand++] = {cidx, 1.0f / d2};
  };

  // Center
  add_candidate(ix, iy);

  // 4 direct neighbours (if they exist)
  if (ix > 0)      add_candidate(ix - 1, iy);
  if (ix < nx - 1) add_candidate(ix + 1, iy);
  if (iy > 0)      add_candidate(ix, iy - 1);
  if (iy < ny - 1) add_candidate(ix, iy + 1);

  // Normalise weights
  float sum_w = 0.0f;
  for (int i = 0; i < ncand; ++i)
    sum_w += candidates[i].inv_dist;

  if (sum_w <= 0.0f) return;

  for (int i = 0; i < ncand; ++i) {
    float bw = candidates[i].inv_dist / sum_w;
    if (bw < 1e-6f) continue;
    H.obs.push_back({candidates[i].idx, value, obs_weight * bw, alt_dist});
    H.obs_count[candidates[i].idx]++;
  }
}

auto build_observation_operator(
      volume const& vol
    , grid_bearings const& gb
    , size_t grid_nx
    , size_t grid_ny
    , float altitude
    , float grid_spacing
    , const vargrid_config& cfg
    ) -> observation_operator
{
  observation_operator H;
  H.grid_nx = grid_nx;
  H.grid_ny = grid_ny;
  H.obs_count.resize(grid_nx * grid_ny, 0);
  H.kappa = cfg.kappa;

  // Copy per-cell azimuth for Wx/Wy weights
  H.cell_azimuth_deg = gb.cell_bearing_deg;

  size_t topscan = 0;
  for (size_t iscan = 1; iscan < vol.sweeps.size(); iscan++) {
    if (vol.sweeps[iscan].beam.elevation() > vol.sweeps[topscan].beam.elevation())
      topscan = iscan;
  }

  size_t total_obs = 0;
  size_t out_of_bounds = 0;

  for (size_t iscan = 0; iscan < vol.sweeps.size(); ++iscan) {
    if (iscan == topscan) continue;
    auto& scan = vol.sweeps[iscan];

    for (size_t ibin = 0; ibin < scan.bins.size(); ++ibin) {
      float alt_dist = scan.bins[ibin].altitude - altitude;
      if (std::fabs(alt_dist) > cfg.max_alt_diff) continue;

      for (size_t iray = 0; iray < scan.rays.size(); ++iray) {
        float val = scan.data[iray][ibin];
        if (std::isnan(val)) continue;
        if (std::fabs(val - undetect) < 0.1f) continue;

        float bearing_deg = scan.rays[iray].degrees();
        if (bearing_deg < 0) bearing_deg += 360.0f;
        float ground_range = scan.bins[ibin].ground_range;

        auto nearest = lookup_nearest(gb, bearing_deg, ground_range);
        if (nearest == static_cast<size_t>(-1)) {
          out_of_bounds++;
          continue;
        }

        float obs_w = compute_weight(scan.bins[ibin].slant_range, alt_dist, cfg.max_alt_diff, cfg);

        add_interpolated_observation(H, gb, nearest,
          bearing_deg, ground_range, val, obs_w, alt_dist);
        total_obs++;
      }
    }
  }

  trace::debug("  Observation operator: {} gates -> {} obs entries, {} out of bounds, alt={:.0f}m",
    total_obs, H.obs.size(), out_of_bounds, altitude);

  // Compute distance to nearest observation for each cell (BFS in grid cells).
  // Used by the background constraint JB (Brook et al. 2022 Eq. 4).
  {
    size_t n = grid_nx * grid_ny;
    H.dist_to_obs.assign(n, std::numeric_limits<float>::max());
    std::vector<size_t> frontier;

    for (size_t i = 0; i < n; ++i) {
      if (H.obs_count[i] > 0) {
        H.dist_to_obs[i] = 0.0f;
        frontier.push_back(i);
      }
    }

    // BFS: propagate distance in grid-cell units
    float current_dist = 0.0f;
    while (!frontier.empty()) {
      current_dist += 1.0f;
      std::vector<size_t> next_frontier;
      for (auto idx : frontier) {
        size_t iy = idx / grid_nx;
        size_t ix = idx % grid_nx;
        auto try_cell = [&](size_t cx, size_t cy) {
          size_t cidx = cy * grid_nx + cx;
          if (H.dist_to_obs[cidx] > current_dist) {
            H.dist_to_obs[cidx] = current_dist;
            next_frontier.push_back(cidx);
          }
        };
        if (ix > 0)            try_cell(ix - 1, iy);
        if (ix < grid_nx - 1)  try_cell(ix + 1, iy);
        if (iy > 0)            try_cell(ix, iy - 1);
        if (iy < grid_ny - 1)  try_cell(ix, iy + 1);
      }
      frontier = std::move(next_frontier);
    }
  }

  return H;
}
