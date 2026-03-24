#include "observation_operator.h"
#include "config.h"

#include "pch.h"
#include "cappi.h"

using namespace bom;

// Compute observation weight based on range and altitude distance.
static auto compute_weight(
      float ground_range
    , float alt_dist
    , float max_alt_diff
    , const vargrid_config& cfg
    ) -> float
{
  // Range-based weight: Gaussian decay with range_scale
  float w_range = std::exp(-0.5f * (ground_range * ground_range)
                           / (cfg.range_scale * cfg.range_scale));

  // Altitude-based weight: linear decay to zero at max_alt_diff
  float w_alt = 1.0f - std::fabs(alt_dist) / max_alt_diff;
  if (w_alt < 0.0f) w_alt = 0.0f;

  float w = w_range * w_alt;
  return std::max(w, cfg.min_weight);
}

// Build the observation operator for a single altitude layer.
// Maps radar gates to the Cartesian grid using the map projection.
// Takes the proj4 string so each thread can create its own projection context.
auto build_observation_operator(
      volume const& vol
    , string const& proj4_string
    , grid_coordinates const& coords
    , size_t grid_nx
    , size_t grid_ny
    , float altitude
    , const vargrid_config& cfg
    ) -> observation_operator
{
  observation_operator H;
  H.grid_nx = grid_nx;
  H.grid_ny = grid_ny;
  H.obs_count.resize(grid_nx * grid_ny, 0);

  // Find the top scan (excluded, same logic as CAPPI)
  size_t topscan = 0;
  for (size_t iscan = 1; iscan < vol.sweeps.size(); iscan++) {
    if (vol.sweeps[iscan].beam.elevation() > vol.sweeps[topscan].beam.elevation())
      topscan = iscan;
  }

  // Grid geometry in projected coordinates (meters).
  // Use cell edges to compute cell centers and spacing.
  auto col_edges = coords.col_edges();
  auto row_edges = coords.row_edges();

  // Cell centers
  std::vector<double> x_centers(grid_nx), y_centers(grid_ny);
  for (size_t i = 0; i < grid_nx; ++i)
    x_centers[i] = 0.5 * (col_edges[i] + col_edges[i + 1]);
  for (size_t j = 0; j < grid_ny; ++j)
    y_centers[j] = 0.5 * (row_edges[j] + row_edges[j + 1]);

  // Cell spacing (assumed uniform)
  double dx = x_centers[1] - x_centers[0];
  double dy = y_centers[1] - y_centers[0];

  auto radarloc = latlon{vol.location.lon, vol.location.lat};

  // Each call creates its own projection objects for thread safety.
  // PROJ contexts are not safe to share across threads.
  auto geo = map_projection{map_projection::default_context, "+proj=longlat +datum=WGS84"};
  auto local_proj = map_projection{map_projection::default_context, proj4_string};

  size_t total_obs = 0;
  size_t out_of_bounds = 0;

  for (size_t iscan = 0; iscan < vol.sweeps.size(); ++iscan) {
    if (iscan == topscan)
      continue;

    auto& scan = vol.sweeps[iscan];

    for (size_t ibin = 0; ibin < scan.bins.size(); ++ibin) {
      float alt_dist = scan.bins[ibin].altitude - altitude;

      if (std::fabs(alt_dist) > cfg.max_alt_diff)
        continue;

      for (size_t iray = 0; iray < scan.rays.size(); ++iray) {
        float val = scan.data[iray][ibin];

        if (std::isnan(val))
          continue;
        if (std::fabs(val - undetect) < 0.1f)
          continue;

        // Compute the lat/lon of this gate
        auto gate_ll = wgs84.bearing_range_to_latlon(
          radarloc,
          scan.rays[iray],
          scan.bins[ibin].ground_range
        );

        // Transform gate from geographic to the grid's projected CRS.
        double px = gate_ll.lon.degrees();
        double py = gate_ll.lat.degrees();
        auto sx = span<double>{&px, 1};
        auto sy = span<double>{&py, 1};
        map_projection::transform(geo, local_proj, sx, sy);

        // Convert projected coordinates to grid indices
        double fx = (px - x_centers[0]) / dx;
        double fy = (py - y_centers[0]) / dy;

        int ix = static_cast<int>(std::round(fx));
        int iy = static_cast<int>(std::round(fy));

        if (ix < 0 || ix >= static_cast<int>(grid_nx) ||
            iy < 0 || iy >= static_cast<int>(grid_ny)) {
          out_of_bounds++;
          continue;
        }

        size_t grid_idx = static_cast<size_t>(iy) * grid_nx + static_cast<size_t>(ix);

        float w = compute_weight(
          scan.bins[ibin].ground_range,
          alt_dist,
          cfg.max_alt_diff,
          cfg
        );

        H.obs.push_back({grid_idx, val, w, alt_dist});
        H.obs_count[grid_idx]++;
        total_obs++;
      }
    }
  }

  trace::debug("  Observation operator: {} obs collected, {} out of bounds, altitude={:.0f}m",
    total_obs, out_of_bounds, altitude);

  return H;
}