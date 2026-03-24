#include "cappi.h"

auto find_ground_range_bin(array1<bin_info> const& bins, float target) -> size_t {
  // binary search
  float dr = bins[1].ground_range - bins[0].ground_range; // assume uniform spacing
  if (dr <= 0)
    return 0; // guard against invalid data
  size_t idx = std::round((target - bins[0].ground_range) / dr);
  if (idx >= bins.size())
    return bins.size() - 1;
  return idx;
}

auto find_ray(array1<angle> const& rays, angle target) -> size_t{
  // HACK - not handling 0/360 well
  auto ipos = std::upper_bound(rays.begin(), rays.end(), target);
  if (ipos == rays.end())
    return rays.size() - 1;
  if (ipos == rays.begin())
    return 0;
  if (abs(*ipos - target) < abs(*(ipos - 1) - target))
    return ipos - rays.begin();
  return ipos - rays.begin() - 1;
}

auto generate_cappi(
      volume const& vol
    , array2<latlon> const& latlons
    , float max_alt_diff
    , float idw_pwr
    , float altitude
    , const float roi
    , const float beamwidth
    ) -> array2f
{
  auto cappi = array2f{latlons.extents()};

  size_t topscan=0;
  for (size_t iscan = 1; iscan < vol.sweeps.size(); iscan++){
    if(vol.sweeps[iscan].beam.elevation() > vol.sweeps[topscan].beam.elevation())
      topscan = iscan;
  }

  // Pre-sort sweep indices by elevation (ascending).
  // This allows early termination: once we pass the target altitude,
  // we have the tightest lower/upper bounds and can stop searching.
  auto sorted_sweeps = std::vector<size_t>(vol.sweeps.size());
  for (size_t i = 0; i < sorted_sweeps.size(); ++i) sorted_sweeps[i] = i;
  std::sort(sorted_sweeps.begin(), sorted_sweeps.end(),
    [&vol](size_t a, size_t b) {
      return vol.sweeps[a].beam.elevation() < vol.sweeps[b].beam.elevation();
    });

  for (size_t y = 0; y < cappi.extents().y; ++y){
    for (size_t x = 0; x < cappi.extents().x; ++x){
      auto ll = latlons[y][x];
      auto br = wgs84.latlon_to_bearing_range(vol.location, ll);
      br.first = br.first.normalize();

      int lwr_scan = -1, upr_scan = -1;
      float lwr_val = nodata;
      float upr_val = nodata;
      float lwr_dist = nodata;
      float upr_dist = nodata;
      for (auto iscan : sorted_sweeps){
        if (iscan == topscan)
          continue;
        auto& scan = vol.sweeps[iscan];

        auto ibin = find_ground_range_bin(scan.bins, br.second);
        if (ibin >= scan.bins.size())
          continue;

        auto alt_dist = scan.bins[ibin].altitude - altitude;
        if (alt_dist > max_alt_diff)
          continue;

        auto iray = find_ray(scan.rays, br.first);
        auto val = scan.data[iray][ibin];
        if (std::isnan(val) || std::fabs(val - undetect) < .1f)
          continue;

        if (alt_dist <= 0)
        {
          // Sweeps are sorted by elevation, so each valid lower hit
          // is guaranteed to be higher altitude than the previous one.
          lwr_scan = iscan;
          lwr_val = val;
          lwr_dist = -alt_dist;
        } else {
          // First valid upper hit is the closest above.
          upr_scan = iscan;
          upr_val = val;
          upr_dist = alt_dist;
          break; // No need to check higher elevations.
        }
      }

      if (lwr_scan != -1 && upr_scan != -1)
      {
        if (lwr_dist > 0.0f)
        {
          auto idw_lwr = 1.0 / std::pow(lwr_dist, idw_pwr);
          auto idw_upr = 1.0 / std::pow(upr_dist, idw_pwr);
          auto norm = idw_lwr + idw_upr;
          cappi[y][x] = lwr_val * (idw_lwr / norm) + upr_val * (idw_upr / norm);
        }
        else // if (lwr_scan != 0)
          cappi[y][x] = lwr_val;
      }
      else if (lwr_scan != -1 && lwr_dist < roi)
        cappi[y][x] = lwr_val;
      else if (upr_scan != -1 && upr_dist < roi)
        cappi[y][x] = upr_val;
      else
        cappi[y][x] = nodata;
    }
  }

  return cappi;
}

auto precompute_grid_lookup(
      volume const& vol
    , array2<latlon> const& latlons
    ) -> grid_lookup
{
  auto ny = latlons.extents().y;
  auto nx = latlons.extents().x;
  auto nsweeps = vol.sweeps.size();

  grid_lookup lookup;
  lookup.grid_nx = nx;
  lookup.grid_ny = ny;
  lookup.nsweeps = nsweeps;
  lookup.data.resize(ny * nx * nsweeps);

  for (size_t y = 0; y < ny; ++y) {
    for (size_t x = 0; x < nx; ++x) {
      auto ll = latlons[y][x];
      auto br = wgs84.latlon_to_bearing_range(vol.location, ll);
      br.first = br.first.normalize();

      for (size_t iscan = 0; iscan < nsweeps; ++iscan) {
        auto& scan = vol.sweeps[iscan];
        auto& sl = lookup.at(y, x, iscan);
        sl.ibin = find_ground_range_bin(scan.bins, br.second);
        sl.iray = find_ray(scan.rays, br.first);
      }
    }
  }

  return lookup;
}

auto generate_cappi(
      volume const& vol
    , grid_lookup const& lookup
    , size_t grid_ny
    , size_t grid_nx
    , float max_alt_diff
    , float idw_pwr
    , float altitude
    , const float roi
    ) -> array2f
{
  auto cappi = array2f{vec2z{grid_nx, grid_ny}};

  size_t topscan = 0;
  for (size_t iscan = 1; iscan < vol.sweeps.size(); iscan++) {
    if (vol.sweeps[iscan].beam.elevation() > vol.sweeps[topscan].beam.elevation())
      topscan = iscan;
  }

  // Pre-sort sweep indices by elevation (ascending).
  auto sorted_sweeps = std::vector<size_t>(vol.sweeps.size());
  for (size_t i = 0; i < sorted_sweeps.size(); ++i) sorted_sweeps[i] = i;
  std::sort(sorted_sweeps.begin(), sorted_sweeps.end(),
    [&vol](size_t a, size_t b) {
      return vol.sweeps[a].beam.elevation() < vol.sweeps[b].beam.elevation();
    });

  for (size_t y = 0; y < grid_ny; ++y) {
    for (size_t x = 0; x < grid_nx; ++x) {

      int lwr_scan = -1, upr_scan = -1;
      float lwr_val = nodata;
      float upr_val = nodata;
      float lwr_dist = nodata;
      float upr_dist = nodata;

      for (auto iscan : sorted_sweeps) {
        if (iscan == topscan)
          continue;
        auto& scan = vol.sweeps[iscan];
        auto& sl = lookup.at(y, x, iscan);

        auto ibin = sl.ibin;
        if (ibin >= scan.bins.size())
          continue;

        auto alt_dist = scan.bins[ibin].altitude - altitude;
        if (alt_dist > max_alt_diff)
          continue;

        auto iray = sl.iray;
        auto val = scan.data[iray][ibin];
        if (std::isnan(val) || std::fabs(val - undetect) < .1f)
          continue;

        if (alt_dist <= 0)
        {
          lwr_scan = iscan;
          lwr_val = val;
          lwr_dist = -alt_dist;
        } else {
          upr_scan = iscan;
          upr_val = val;
          upr_dist = alt_dist;
          break;
        }
      }

      if (lwr_scan != -1 && upr_scan != -1)
      {
        if (lwr_dist > 0.0f)
        {
          auto idw_lwr = 1.0 / std::pow(lwr_dist, idw_pwr);
          auto idw_upr = 1.0 / std::pow(upr_dist, idw_pwr);
          auto norm = idw_lwr + idw_upr;
          cappi[y][x] = lwr_val * (idw_lwr / norm) + upr_val * (idw_upr / norm);
        }
        else
          cappi[y][x] = lwr_val;
      }
      else if (lwr_scan != -1 && lwr_dist < roi)
        cappi[y][x] = lwr_val;
      else if (upr_scan != -1 && upr_dist < roi)
        cappi[y][x] = upr_val;
      else
        cappi[y][x] = nodata;
    }
  }

  return cappi;
}