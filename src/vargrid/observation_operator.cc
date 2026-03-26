#include "observation_operator.h"
#include "config.h"

#include "pch.h"
#include "cappi.h"

using namespace bom;

// Compute observation weight (inverse error variance) for a radar gate.
//
// The weight models two physical effects:
//
// 1. Beam broadening: the radar pulse volume grows with range.
//    At slant range r, the beam cross-section is ~ (r * beamwidth)^2
//    in area. A gate at long range averages over a much larger volume
//    than a nearby gate, making it less representative of a single
//    grid cell. We model this as w_beam = (ref_range / r)^beam_power.
//    With beam_power=2, this is inverse-area scaling.
//
// 2. Altitude distance: gates further from the target CAPPI altitude
//    are less relevant. Linear decay to zero at max_alt_diff.
//
// The combined weight w = w_beam * w_alt represents R⁻¹ in the
// variational cost function: higher weight = lower observation error
// = more trust in this gate.
static auto compute_weight(
      float slant_range
    , float alt_dist
    , float max_alt_diff
    , const vargrid_config& cfg
    ) -> float
{
  // Beam volume weight: close gates get high weight, far gates get low weight.
  // Clamped so that gates closer than ref_range don't get artificially boosted.
  float effective_range = std::max(slant_range, cfg.ref_range);
  float w_beam = std::pow(cfg.ref_range / effective_range, cfg.beam_power);

  // Altitude weight: linear decay
  float w_alt = 1.0f - std::fabs(alt_dist) / max_alt_diff;
  if (w_alt < 0.0f) w_alt = 0.0f;

  float w = w_beam * w_alt;
  return std::max(w, cfg.min_weight);
}

// Precompute projected coordinates for all radar gates.
// Must be called from the main thread (single-threaded) because
// PROJ is not safe to use concurrently with shared contexts.
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

    // Build arrays of lon/lat for all gates in this sweep
    std::vector<double> lons(nrays * nbins);
    std::vector<double> lats(nrays * nbins);

    for (size_t iray = 0; iray < nrays; ++iray) {
      for (size_t ibin = 0; ibin < nbins; ++ibin) {
        auto gate_ll = wgs84.bearing_range_to_latlon(
          radarloc,
          scan.rays[iray],
          scan.bins[ibin].ground_range
        );
        size_t idx = iray * nbins + ibin;
        lons[idx] = gate_ll.lon.degrees();
        lats[idx] = gate_ll.lat.degrees();
      }
    }

    // Batch transform: geographic -> projected CRS
    auto sx = span<double>{lons.data(), lons.size()};
    auto sy = span<double>{lats.data(), lats.size()};
    map_projection::transform(geo, proj, sx, sy);

    // After transform, lons[] contains easting (x), lats[] contains northing (y)
    gp.sweeps[iscan].px = std::move(lons);
    gp.sweeps[iscan].py = std::move(lats);
  }

  return gp;
}

// Build the observation operator for a single altitude layer.
// Uses precomputed projected gate coordinates — fully thread-safe,
// no PROJ calls needed.
auto build_observation_operator(
      volume const& vol
    , gate_projections const& gp
    , double x0
    , double y0
    , double dx
    , double dy
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

  // Find the top scan (excluded)
  size_t topscan = 0;
  for (size_t iscan = 1; iscan < vol.sweeps.size(); iscan++) {
    if (vol.sweeps[iscan].beam.elevation() > vol.sweeps[topscan].beam.elevation())
      topscan = iscan;
  }

  size_t total_obs = 0;
  size_t out_of_bounds = 0;

  for (size_t iscan = 0; iscan < vol.sweeps.size(); ++iscan) {
    if (iscan == topscan)
      continue;

    auto& scan = vol.sweeps[iscan];
    auto& sp = gp.sweeps[iscan];

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

        size_t pidx = iray * sp.nbins + ibin;
        double px = sp.px[pidx];
        double py = sp.py[pidx];

        // Convert projected coordinates to grid indices
        double fx = (px - x0) / dx;
        double fy = (py - y0) / dy;

        int ix = static_cast<int>(std::round(fx));
        int iy = static_cast<int>(std::round(fy));

        if (ix < 0 || ix >= static_cast<int>(grid_nx) ||
            iy < 0 || iy >= static_cast<int>(grid_ny)) {
          out_of_bounds++;
          continue;
        }

        size_t grid_idx = static_cast<size_t>(iy) * grid_nx + static_cast<size_t>(ix);

        float w = compute_weight(
          scan.bins[ibin].slant_range,
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

  trace::debug("  Observation operator: {} obs, {} out of bounds, alt={:.0f}m",
    total_obs, out_of_bounds, altitude);

  return H;
}