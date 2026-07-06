#ifndef VARGRID_GRIDDING_LEROI_CORE_H
#define VARGRID_GRIDDING_LEROI_CORE_H

// Core numerics for the "leroi" gridding method (Dahl et al. 2019, as
// implemented in the Python leroi package). Two-step scheme:
//   1. Interpolate each PPI sweep horizontally onto the 2D grid.
//   2. Interpolate vertically between the stacked sweep surfaces.
//
// This header is intentionally free of BOM library dependencies so the
// numerics can be unit-tested standalone. The BOM-facing wrapper lives in
// leroi.h / leroi.cc.

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

namespace leroi {

constexpr float core_nan = std::numeric_limits<float>::quiet_NaN();
constexpr double deg2rad = 3.14159265358979323846 / 180.0;

enum class weight_type { barnes, cressman, idw, bilinear };

// Polar geometry of one sweep, decoupled from data layout.
// az_deg is sorted ascending in [0, 360); ray_index maps each sorted
// azimuth back to its original row in the sweep's data array.
struct sweep_geometry {
  float elevation_deg = 0.0f;
  std::vector<float>  az_deg;        // sorted ascending, [0, 360)
  std::vector<size_t> ray_index;     // original data row per sorted azimuth
  std::vector<float>  ground_range;  // bin centers, ascending (m)
  std::vector<float>  altitude;      // bin center altitudes (m)
  float max_az_gap_deg = 360.0f;     // largest azimuthal gap (for bilinear guard)
};

// Grid cells in polar coordinates relative to the radar.
struct cell_lookup {
  const float* bearing_deg;  // [n], in [0, 360)
  const float* range;        // [n], ground range (m)
  size_t n;
};

inline auto gate_valid(float v, float undetect_value) -> bool
{
  return !std::isnan(v) && std::fabs(v - undetect_value) > 0.1f;
}

// Build sorted geometry from raw per-sweep azimuths (any order/representation).
inline auto make_sweep_geometry(
      float elevation_deg
    , std::vector<float> const& raw_az_deg      // per data row, any range
    , std::vector<float> const& ground_range
    , std::vector<float> const& altitude
    ) -> sweep_geometry
{
  sweep_geometry g;
  g.elevation_deg = elevation_deg;
  g.ground_range = ground_range;
  g.altitude = altitude;

  size_t n = raw_az_deg.size();
  std::vector<std::pair<float, size_t>> az(n);
  for (size_t i = 0; i < n; ++i) {
    float a = std::fmod(raw_az_deg[i], 360.0f);
    if (a < 0) a += 360.0f;
    az[i] = {a, i};
  }
  std::sort(az.begin(), az.end());

  g.az_deg.resize(n);
  g.ray_index.resize(n);
  for (size_t i = 0; i < n; ++i) {
    g.az_deg[i] = az[i].first;
    g.ray_index[i] = az[i].second;
  }

  // Largest azimuthal gap including the wrap gap (used to reject bilinear
  // interpolation across sector-scan voids).
  if (n >= 2) {
    float max_gap = g.az_deg[0] + 360.0f - g.az_deg[n - 1];
    for (size_t i = 1; i < n; ++i)
      max_gap = std::max(max_gap, g.az_deg[i] - g.az_deg[i - 1]);
    g.max_az_gap_deg = max_gap;
  }
  return g;
}

// Radius of influence from azimuthal gate spacing (Dahl et al. 2019).
// Matches leroi.py get_leroy_roi: frac * (largest azimuth step) * rmax.
inline auto auto_roi(
      std::vector<sweep_geometry> const& sweeps
    , float rmax
    , float frac = 0.6f
    ) -> float
{
  float roi = 0.0f;
  for (auto& g : sweeps) {
    if (g.az_deg.size() < 2) continue;
    float max_gap = 0.0f;  // gap within sorted list, no wrap (matches python)
    for (size_t i = 1; i < g.az_deg.size(); ++i)
      max_gap = std::max(max_gap, g.az_deg[i] - g.az_deg[i - 1]);
    roi = std::max(roi, frac * static_cast<float>(max_gap * deg2rad) * rmax);
  }
  return roi;
}

// Sweep surface altitude at a given ground range (linear interp between bins).
// NaN beyond the last bin; clamped to the first bin's altitude near the radar.
inline auto height_at(sweep_geometry const& g, float r) -> float
{
  auto& gr = g.ground_range;
  if (gr.empty()) return core_nan;
  if (r <= gr.front()) return g.altitude.front();
  if (r > gr.back()) return core_nan;
  auto it = std::upper_bound(gr.begin(), gr.end(), r);
  size_t k = static_cast<size_t>(it - gr.begin()) - 1;
  if (k + 1 >= gr.size()) return g.altitude.back();
  float dr = gr[k + 1] - gr[k];
  if (dr <= 0) return g.altitude[k];
  float t = (r - gr[k]) / dr;
  return g.altitude[k] + t * (g.altitude[k + 1] - g.altitude[k]);
}

// Find up to two [begin, end) index ranges of sorted azimuths within
// +/- half_w degrees of bearing b (handles wraparound).
inline auto az_window(
      std::vector<float> const& az
    , float b
    , float half_w
    , size_t out[2][2]
    ) -> int
{
  size_t n = az.size();
  if (n == 0) return 0;
  if (half_w >= 180.0f) { out[0][0] = 0; out[0][1] = n; return 1; }

  auto rng = [&](float a0, float a1, int slot) {
    out[slot][0] = static_cast<size_t>(std::lower_bound(az.begin(), az.end(), a0) - az.begin());
    out[slot][1] = static_cast<size_t>(std::upper_bound(az.begin(), az.end(), a1) - az.begin());
  };

  float lo = b - half_w, hi = b + half_w;
  if (lo < 0.0f)         { rng(lo + 360.0f, 360.0f, 0); rng(0.0f, hi, 1); return 2; }
  if (hi >= 360.0f)      { rng(lo, 360.0f, 0); rng(0.0f, hi - 360.0f, 1); return 2; }
  rng(lo, hi, 0);
  return 1;
}

// Horizontally interpolate one sweep of one or more fields onto the grid.
// field_data[f] points to row-major [ray][bin] data in ORIGINAL ray order
// (geometry's ray_index does the remapping). out[f] must be preallocated
// to cells.n and is overwritten (NaN where no valid data).
inline void interp_sweep(
      sweep_geometry const& g
    , std::vector<const float*> const& field_data
    , size_t nbins
    , cell_lookup const& cells
    , weight_type wt
    , float roi
    , float idw_pwr
    , float undetect_value
    , std::vector<std::vector<float>>& out
    )
{
  size_t nfields = field_data.size();
  size_t nrays = g.az_deg.size();
  auto& gr = g.ground_range;

  for (size_t f = 0; f < nfields; ++f)
    std::fill(out[f].begin(), out[f].end(), core_nan);
  if (nrays == 0 || gr.empty()) return;

  std::vector<double> sw(nfields), swv(nfields);
  const float roi2 = roi * roi;
  const float barnes_denom = roi2 / 4.0f;
  const float barnes_cut = std::exp(-4.0f);

  if (wt == weight_type::bilinear) {
    // Reject interpolation across azimuthal gaps much wider than the
    // typical beam spacing (e.g. sector scan voids).
    float typical = 360.0f / static_cast<float>(nrays);
    float max_gap = std::max(4.0f * typical, 2.0f);

    for (size_t i = 0; i < cells.n; ++i) {
      float R = cells.range[i], B = cells.bearing_deg[i];
      if (R < gr.front() || R > gr.back()) continue;

      // range bracket
      auto it = std::upper_bound(gr.begin(), gr.end(), R);
      size_t k = static_cast<size_t>(it - gr.begin());
      if (k == 0) k = 1;
      if (k >= gr.size()) k = gr.size() - 1;
      size_t k0 = k - 1, k1 = k;
      float drr = gr[k1] - gr[k0];
      float fr = drr > 0 ? (R - gr[k0]) / drr : 0.0f;

      // azimuth bracket (with wraparound)
      size_t j0, j1;
      float span, fa;
      if (B < g.az_deg.front() || B >= g.az_deg.back()) {
        j0 = nrays - 1; j1 = 0;
        span = g.az_deg.front() + 360.0f - g.az_deg.back();
        float d = B - g.az_deg.back();
        if (d < 0) d += 360.0f;
        fa = span > 0 ? d / span : 0.0f;
      } else {
        auto jt = std::upper_bound(g.az_deg.begin(), g.az_deg.end(), B);
        j1 = static_cast<size_t>(jt - g.az_deg.begin());
        j0 = j1 - 1;
        span = g.az_deg[j1] - g.az_deg[j0];
        fa = span > 0 ? (B - g.az_deg[j0]) / span : 0.0f;
      }
      if (span > max_gap) continue;

      const size_t rows[2] = {g.ray_index[j0], g.ray_index[j1]};
      const size_t cols[2] = {k0, k1};
      const float wgt[2][2] = {
        {(1 - fa) * (1 - fr), (1 - fa) * fr},
        {fa * (1 - fr),       fa * fr}
      };

      for (size_t f = 0; f < nfields; ++f) {
        double s = 0, sv = 0;
        for (int a = 0; a < 2; ++a)
          for (int r = 0; r < 2; ++r) {
            float v = field_data[f][rows[a] * nbins + cols[r]];
            if (!gate_valid(v, undetect_value)) continue;
            s += wgt[a][r];
            sv += static_cast<double>(wgt[a][r]) * v;
          }
        // Require a minimal share of valid stencil weight; otherwise the
        // renormalised value would extrapolate from a far corner.
        if (s >= 0.1)
          out[f][i] = static_cast<float>(sv / s);
      }
    }
    return;
  }

  // Radius-gather methods: barnes | cressman | idw
  for (size_t i = 0; i < cells.n; ++i) {
    float R = cells.range[i], B = cells.bearing_deg[i];
    if (R - roi > gr.back() || R + roi < gr.front()) continue;

    // bin window: |r - R| <= roi
    size_t blo = static_cast<size_t>(std::lower_bound(gr.begin(), gr.end(), R - roi) - gr.begin());
    size_t bhi = static_cast<size_t>(std::upper_bound(gr.begin(), gr.end(), R + roi) - gr.begin());
    if (blo >= bhi) continue;

    // azimuth window: gates with |sin(daz)| > roi/R are all beyond roi
    float half_w = 180.0f;
    if (R > roi)
      half_w = static_cast<float>(std::asin(static_cast<double>(roi) / R) / deg2rad);
    size_t win[2][2];
    int nwin = az_window(g.az_deg, B, half_w, win);

    std::fill(sw.begin(), sw.end(), 0.0);
    std::fill(swv.begin(), swv.end(), 0.0);

    for (int wI = 0; wI < nwin; ++wI) {
      for (size_t j = win[wI][0]; j < win[wI][1]; ++j) {
        double dcos = std::cos((B - g.az_deg[j]) * deg2rad);
        size_t row = g.ray_index[j] * nbins;
        for (size_t k = blo; k < bhi; ++k) {
          float r = gr[k];
          float d2 = static_cast<float>(
            static_cast<double>(R) * R + static_cast<double>(r) * r
            - 2.0 * R * r * dcos);
          if (d2 > roi2) continue;

          float w;
          switch (wt) {
          case weight_type::barnes:
            w = std::exp(-d2 / barnes_denom);
            if (w < barnes_cut) continue;
            break;
          case weight_type::cressman:
            w = (roi2 - d2) / (roi2 + d2);
            if (w <= 0.0f) continue;
            break;
          default: {  // idw
            float d = std::sqrt(std::max(d2, 0.0f));
            d = std::max(d, 1.0f);  // avoid singularity at exact hits
            w = 1.0f / std::pow(d, idw_pwr);
            break;
          }
          }

          for (size_t f = 0; f < nfields; ++f) {
            float v = field_data[f][row + k];
            if (!gate_valid(v, undetect_value)) continue;
            sw[f] += w;
            swv[f] += static_cast<double>(w) * v;
          }
        }
      }
    }

    for (size_t f = 0; f < nfields; ++f)
      if (sw[f] > 0.0)
        out[f][i] = static_cast<float>(swv[f] / sw[f]);
  }
}

// Fill NaN heights at the bottom/top of each column with the nearest valid
// height (matches leroi.py fill_heights). Returns the number of columns
// that remain entirely NaN (grid cells outside radar range).
inline auto fill_heights(std::vector<std::vector<float>>& h, size_t ncells) -> size_t
{
  size_t ns = h.size();
  size_t unfilled = 0;
  for (size_t c = 0; c < ncells; ++c) {
    size_t first = ns, last = ns;
    for (size_t s = 0; s < ns; ++s) {
      if (!std::isnan(h[s][c])) {
        if (first == ns) first = s;
        last = s;
      }
    }
    if (first == ns) { unfilled++; continue; }
    for (size_t s = 0; s < first; ++s) h[s][c] = h[first][c];
    for (size_t s = last + 1; s < ns; ++s) h[s][c] = h[last][c];
  }
  return unfilled;
}

// Vertically interpolate stacked sweep surfaces to a constant altitude.
// heights must be non-decreasing along the sweep axis per column (after
// fill_heights). NaN surface values propagate; z outside the column's
// height span yields NaN (no extrapolation), matching leroi.py.
inline void vertical_slice(
      std::vector<std::vector<float>> const& heights
    , std::vector<std::vector<float>> const& surf
    , float z
    , float* out
    , size_t ncells
    )
{
  size_t ns = heights.size();
  for (size_t c = 0; c < ncells; ++c) {
    out[c] = core_nan;
    for (size_t s = 0; s + 1 < ns; ++s) {
      float h0 = heights[s][c], h1 = heights[s + 1][c];
      if (std::isnan(h0) || std::isnan(h1)) break;
      if (!(h0 <= z && z <= h1)) continue;
      float v0 = surf[s][c], v1 = surf[s + 1][c];
      float dh = h1 - h0;
      if (dh < 1e-3f)
        out[c] = v0;  // duplicate height (filled column edge)
      else
        out[c] = v0 + (z - h0) / dh * (v1 - v0);
      break;
    }
  }
}

}  // namespace leroi

#endif // VARGRID_GRIDDING_LEROI_CORE_H
