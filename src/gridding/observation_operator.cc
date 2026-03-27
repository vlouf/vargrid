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

  // Compute bearing and range for every grid cell
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

// Find the nearest grid cell via binned lookup. Returns size_t(-1) if out of range.
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

// Add bilinear observation contributions for a single gate.
// The gate maps to a fractional position (fy, fx) in the grid.
// The 4 surrounding cells receive weighted contributions.
static void add_bilinear_observation(
      observation_operator& H
    , size_t ix0, size_t iy0   // lower-left grid cell
    , float wx, float wy       // fractional position within cell [0,1]
    , float value
    , float obs_weight
    , float alt_dist
    )
{
  size_t nx = H.grid_nx;
  size_t ny = H.grid_ny;

  // Bilinear weights for the 4 corners
  float w00 = (1.0f - wx) * (1.0f - wy);
  float w10 = wx * (1.0f - wy);
  float w01 = (1.0f - wx) * wy;
  float w11 = wx * wy;

  size_t ix1 = std::min(ix0 + 1, nx - 1);
  size_t iy1 = std::min(iy0 + 1, ny - 1);

  // Each corner gets the observation weighted by the bilinear coefficient
  auto add = [&](size_t ix, size_t iy, float bw) {
    if (bw < 1e-6f) return;
    size_t idx = iy * nx + ix;
    H.obs.push_back({idx, value, obs_weight * bw, alt_dist});
    H.obs_count[idx]++;
  };

  add(ix0, iy0, w00);
  add(ix1, iy0, w10);
  add(ix0, iy1, w01);
  add(ix1, iy1, w11);
}

// Build observation operator with bilinear interpolation.
// Each gate contributes to up to 4 grid cells proportionally.
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
        float ground_range = scan.bins[ibin].ground_range;

        // Find nearest grid cell
        auto nearest = lookup_nearest(gb, bearing_deg, ground_range);
        if (nearest == static_cast<size_t>(-1)) {
          out_of_bounds++;
          continue;
        }

        // Convert flat index to (ix, iy)
        size_t iy_near = nearest / grid_nx;
        size_t ix_near = nearest % grid_nx;

        // Compute fractional position relative to the nearest cell by
        // comparing the gate's bearing/range with neighbouring cells.
        // Find which quadrant the gate is in relative to the nearest cell.
        float cell_b = gb.cell_bearing_deg[nearest];
        float cell_r = gb.cell_range[nearest];

        // Bearing difference (handle 0/360 wrap)
        float db = bearing_deg - cell_b;
        if (db > 180.0f) db -= 360.0f;
        if (db < -180.0f) db += 360.0f;

        float dr = ground_range - cell_r;

        // Determine the lower-left cell of the bilinear quad.
        // If the gate is to the "right" (positive db in grid terms) of the cell,
        // the cell is the left edge. Otherwise, shift left by 1.
        // Same logic for range (up/down in the grid).
        //
        // The grid y-axis generally increases with range from radar,
        // and x-axis with bearing. But the exact mapping depends on
        // the grid orientation. For a regular grid, we approximate:
        // the bearing direction maps primarily to x, range to y.
        // Use the grid spacing to convert bearing/range offsets to fractional cells.
        float bearing_spacing_m = cell_r * cfg.beamwidth * M_PI / 180.0f;
        if (bearing_spacing_m < grid_spacing) bearing_spacing_m = grid_spacing;

        // Fractional offset in grid cells (approximate)
        float fx_off = (db * M_PI / 180.0f * cell_r) / grid_spacing;
        float fy_off = dr / grid_spacing;

        // The fractional grid position
        float fx = static_cast<float>(ix_near) + fx_off;
        float fy = static_cast<float>(iy_near) + fy_off;

        // Clamp to grid bounds
        if (fx < 0.0f) fx = 0.0f;
        if (fy < 0.0f) fy = 0.0f;
        if (fx > static_cast<float>(grid_nx - 1)) fx = static_cast<float>(grid_nx - 1);
        if (fy > static_cast<float>(grid_ny - 1)) fy = static_cast<float>(grid_ny - 1);

        // Integer lower-left and fractional part
        size_t ix0 = static_cast<size_t>(fx);
        size_t iy0 = static_cast<size_t>(fy);
        if (ix0 >= grid_nx - 1) ix0 = grid_nx - 2;
        if (iy0 >= grid_ny - 1) iy0 = grid_ny - 2;
        float wx = fx - static_cast<float>(ix0);
        float wy = fy - static_cast<float>(iy0);

        float obs_w = compute_weight(scan.bins[ibin].slant_range, alt_dist, cfg.max_alt_diff, cfg);

        add_bilinear_observation(H, ix0, iy0, wx, wy, val, obs_w, alt_dist);
        total_obs++;
      }
    }
  }

  trace::debug("  Observation operator: {} gates -> {} obs entries, {} out of bounds, alt={:.0f}m",
    total_obs, H.obs.size(), out_of_bounds, altitude);

  return H;
}
