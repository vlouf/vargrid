#include "observation_operator.h"
#include "config.h"
#include "../types.h"

// Compute observation weight (inverse error variance) for a radar gate.
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

// Precompute projected coordinates for all radar gates (main thread only).
auto precompute_gate_projections(
      volume const& vol
    , string const& proj4_string
    ) -> gate_projections
{
  gate_projections gp;

  auto radarloc = latlon{vol.location.lon, vol.location.lat};
  auto geo = map_projection{map_projection::default_context, "+proj=longlat +datum=WGS84"};
  auto proj = map_projection{map_projection::default_context, proj4_string};

  gp.sweeps.resize(vol.sweeps.size());

  for (size_t iscan = 0; iscan < vol.sweeps.size(); ++iscan) {
    auto& scan = vol.sweeps[iscan];
    auto nrays = scan.rays.size();
    auto nbins = scan.bins.size();

    gp.sweeps[iscan].nrays = nrays;
    gp.sweeps[iscan].nbins = nbins;

    constexpr double deg2rad = M_PI / 180.0;

    std::vector<vec2d> xy(nrays * nbins);

    for (size_t iray = 0; iray < nrays; ++iray) {
      for (size_t ibin = 0; ibin < nbins; ++ibin) {
        auto gate_ll = wgs84.bearing_range_to_latlon(
          radarloc, scan.rays[iray], scan.bins[ibin].ground_range);
        size_t idx = iray * nbins + ibin;
        xy[idx] = vec2d{gate_ll.lon.degrees() * deg2rad, gate_ll.lat.degrees() * deg2rad};
      }
    }

    // Batch transform: geographic (radians) -> projected CRS
    auto xy_span = span<vec2d>{xy.data(), static_cast<ptrdiff_t>(xy.size())};
    map_projection::transform(geo, proj, xy_span);

    // Extract projected coordinates
    gp.sweeps[iscan].px.resize(nrays * nbins);
    gp.sweeps[iscan].py.resize(nrays * nbins);
    for (size_t i = 0; i < xy.size(); ++i) {
      gp.sweeps[iscan].px[i] = xy[i].x;
      gp.sweeps[iscan].py[i] = xy[i].y;
    }
  }

  return gp;
}

// Build observation operator for a single altitude layer (thread-safe).
auto build_observation_operator(
      volume const& vol
    , gate_projections const& gp
    , double x0, double y0, double dx, double dy
    , size_t grid_nx, size_t grid_ny
    , float altitude
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
    auto& sp = gp.sweeps[iscan];

    for (size_t ibin = 0; ibin < scan.bins.size(); ++ibin) {
      float alt_dist = scan.bins[ibin].altitude - altitude;
      if (std::fabs(alt_dist) > cfg.max_alt_diff) continue;

      for (size_t iray = 0; iray < scan.rays.size(); ++iray) {
        float val = scan.data[iray][ibin];
        if (std::isnan(val)) continue;
        if (std::fabs(val - undetect) < 0.1f) continue;

        size_t pidx = iray * sp.nbins + ibin;
        double fx = (sp.px[pidx] - x0) / dx;
        double fy = (sp.py[pidx] - y0) / dy;

        int ix = static_cast<int>(std::round(fx));
        int iy = static_cast<int>(std::round(fy));

        if (ix < 0 || ix >= static_cast<int>(grid_nx) ||
            iy < 0 || iy >= static_cast<int>(grid_ny)) {
          out_of_bounds++;
          continue;
        }

        size_t grid_idx = static_cast<size_t>(iy) * grid_nx + static_cast<size_t>(ix);
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