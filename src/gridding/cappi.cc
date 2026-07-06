#include "cappi.h"

static auto find_ground_range_bin(array1<bin_info> const& bins, float target) -> size_t {
  if (bins.size() < 2) return 0;
  float dr = bins[1].ground_range - bins[0].ground_range;
  if (dr <= 0) return 0;
  // Round in floating point and clamp BEFORE converting to size_t: a target
  // below the first bin used to round to -1, wrap through size_t, and clamp
  // to the LAST bin — handing near-radar cells data from maximum range.
  float idxf = std::round((target - bins[0].ground_range) / dr);
  if (idxf <= 0.0f) return 0;
  size_t idx = static_cast<size_t>(idxf);
  if (idx >= bins.size()) return bins.size() - 1;
  return idx;
}

// Angular separation in degrees, accounting for the 360 wrap.
static auto angle_sep_deg(angle a, angle b) -> double {
  double d = std::fabs((a - b).degrees());
  if (d > 180.0) d = 360.0 - d;
  return d;
}

static auto find_ray(array1<angle> const& rays, angle target) -> size_t {
  if (rays.size() < 2) return 0;
  auto ipos = std::upper_bound(rays.begin(), rays.end(), target);
  // Nearest of the two bracketing rays, treating the array as circular so a
  // bearing between the last and first ray resolves across the seam.
  size_t hi = (ipos == rays.end()) ? 0 : static_cast<size_t>(ipos - rays.begin());
  size_t lo = (hi == 0) ? rays.size() - 1 : hi - 1;
  return angle_sep_deg(rays[hi], target) < angle_sep_deg(rays[lo], target) ? hi : lo;
}

auto generate_cappi(
      volume const& vol
    , array2<latlon> const& latlons
    , float max_alt_diff
    , float idw_pwr
    , float altitude
    , const float roi
    ) -> array2f
{
  auto cappi = array2f{latlons.extents()};

  // Pre-sort sweep indices by elevation for early termination
  auto sorted_sweeps = std::vector<size_t>(vol.sweeps.size());
  for (size_t i = 0; i < sorted_sweeps.size(); ++i) sorted_sweeps[i] = i;
  std::sort(sorted_sweeps.begin(), sorted_sweeps.end(),
    [&vol](size_t a, size_t b) {
      return vol.sweeps[a].beam.elevation() < vol.sweeps[b].beam.elevation();
    });

  for (size_t y = 0; y < cappi.extents().y; ++y) {
    for (size_t x = 0; x < cappi.extents().x; ++x) {
      auto ll = latlons[y][x];
      auto br = wgs84.latlon_to_bearing_range(vol.location, ll);
      br.first = br.first.normalize();

      int lwr_scan = -1, upr_scan = -1;
      float lwr_val = nodata, upr_val = nodata;
      float lwr_dist = nodata, upr_dist = nodata;

      for (auto iscan : sorted_sweeps) {
        auto& scan = vol.sweeps[iscan];
        // Skip birdbath (vertical) scans; unlike the old "skip the highest
        // sweep" heuristic this does not discard the top real tilt.
        if (scan.beam.elevation().degrees() > 85.0) continue;

        auto ibin = find_ground_range_bin(scan.bins, br.second);
        if (ibin >= scan.bins.size()) continue;

        auto alt_dist = scan.bins[ibin].altitude - altitude;
        if (alt_dist > max_alt_diff) continue;

        auto iray = find_ray(scan.rays, br.first);
        auto val = scan.data[iray][ibin];
        if (std::isnan(val) || std::fabs(val - undetect) < .1f) continue;

        if (alt_dist <= 0) {
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

      if (lwr_scan != -1 && upr_scan != -1) {
        if (lwr_dist > 0.0f) {
          auto idw_lwr = 1.0 / std::pow(lwr_dist, idw_pwr);
          auto idw_upr = 1.0 / std::pow(upr_dist, idw_pwr);
          auto norm = idw_lwr + idw_upr;
          cappi[y][x] = lwr_val * (idw_lwr / norm) + upr_val * (idw_upr / norm);
        } else {
          cappi[y][x] = lwr_val;
        }
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
