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

  gb.n_bearing_bins = 720;
  gb.bearing_bin_size = 360.0f / gb.n_bearing_bins;
  gb.range_bin_size = 250.0f;  // match typical range gate spacing
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

// Build observation operator with nearest-neighbour mapping.
// Each gate maps to exactly 1 grid cell — fast and clean.
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
  H.cell_azimuth_deg = gb.cell_bearing_deg;
  H.cell_range = gb.cell_range;

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

        auto grid_idx = lookup_nearest(gb, bearing_deg, ground_range);
        if (grid_idx == static_cast<size_t>(-1)) {
          out_of_bounds++;
          continue;
        }

        float w = compute_weight(scan.bins[ibin].slant_range, alt_dist, cfg.max_alt_diff, cfg);

        H.obs.push_back({grid_idx, val, w, alt_dist});
        H.obs_count[grid_idx]++;
        total_obs++;
      }
    }
  }

  trace::debug("  Observation operator: {} obs, {} out of bounds, alt={:.0f}m",
    total_obs, out_of_bounds, altitude);

  return H;
}