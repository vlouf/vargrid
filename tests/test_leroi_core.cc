// Standalone unit tests for leroi_core.h (no BOM dependencies).
#include "gridding/leroi_core.h"

#include <cstdio>
#include <cstdlib>
#include <random>

using namespace leroi;

static int failures = 0;
#define CHECK(cond, msg) do { \
  if (!(cond)) { std::printf("FAIL: %s (line %d)\n", msg, __LINE__); failures++; } \
} while (0)
#define CHECK_NEAR(a, b, tol, msg) do { \
  float _a = (a), _b = (b); \
  if (std::isnan(_a) || std::fabs(_a - _b) > (tol)) { \
    std::printf("FAIL: %s: got %g expected %g (line %d)\n", msg, _a, _b, __LINE__); failures++; } \
} while (0)

constexpr float UND = -32.0f;

// Analytic test field, linear in Cartesian coordinates.
static float f_lin(float az_deg, float r) {
  float x = r * std::sin(az_deg * deg2rad);
  float y = r * std::cos(az_deg * deg2rad);
  return 10.0f + 1e-4f * x + 2e-4f * y;
}

// Build a synthetic sweep: rays every 1 deg starting at start_az (unsorted
// wrap ordering exercises the ray_index remapping), bins every 250 m.
struct synth {
  sweep_geometry g;
  std::vector<float> data;   // [ray][bin], original ray order
  size_t nrays, nbins;
};

static synth make_synth(float elev_deg, float start_az, size_t nrays = 360, size_t nbins = 200) {
  synth s;
  s.nrays = nrays; s.nbins = nbins;
  std::vector<float> az(nrays), gr(nbins), alt(nbins);
  for (size_t i = 0; i < nrays; ++i) az[i] = std::fmod(start_az + i * (360.0f / nrays), 360.0f);
  float tane = std::tan(elev_deg * deg2rad);
  for (size_t k = 0; k < nbins; ++k) {
    gr[k] = 125.0f + k * 250.0f;
    alt[k] = gr[k] * tane;
  }
  s.g = make_sweep_geometry(elev_deg, az, gr, alt);
  s.data.resize(nrays * nbins);
  for (size_t i = 0; i < nrays; ++i)
    for (size_t k = 0; k < nbins; ++k)
      s.data[i * nbins + k] = f_lin(az[i], gr[k]);
  return s;
}

// Brute-force reference for radius-gather methods: O(all gates) scan.
static float brute_force(synth const& s, float B, float R, weight_type wt,
                         float roi, float pwr) {
  double sw = 0, swv = 0;
  for (size_t j = 0; j < s.g.az_deg.size(); ++j) {
    float a = s.g.az_deg[j];
    size_t row = s.g.ray_index[j] * s.nbins;
    for (size_t k = 0; k < s.g.ground_range.size(); ++k) {
      float r = s.g.ground_range[k];
      float d2 = (float)((double)R*R + (double)r*r - 2.0*R*r*std::cos((B - a) * deg2rad));
      if (d2 > roi * roi) continue;
      float w;
      if (wt == weight_type::barnes) {
        w = std::exp(-d2 / (roi * roi / 4.0f));
        if (w < std::exp(-4.0f)) continue;
      } else if (wt == weight_type::cressman) {
        w = (roi*roi - d2) / (roi*roi + d2);
        if (w <= 0) continue;
      } else {
        float d = std::max(std::sqrt(std::max(d2, 0.0f)), 1.0f);
        w = 1.0f / std::pow(d, pwr);
      }
      float v = s.data[row + k];
      if (!gate_valid(v, UND)) continue;
      sw += w; swv += (double)w * v;
    }
  }
  return sw > 0 ? (float)(swv / sw) : core_nan;
}

int main() {
  // --- az_window ---
  {
    std::vector<float> az(360);
    for (size_t i = 0; i < 360; ++i) az[i] = i + 0.5f;
    size_t win[2][2];
    int n = az_window(az, 180.0f, 2.0f, win);
    CHECK(n == 1 && win[0][1] - win[0][0] == 4, "az_window centered");
    // window [358.5, 2.5] endpoints are ray centers: 5 rays inclusive
    n = az_window(az, 0.5f, 2.0f, win);
    size_t cnt = (win[0][1] - win[0][0]) + (n == 2 ? win[1][1] - win[1][0] : 0);
    CHECK(n == 2 && cnt == 5, "az_window wrap low");
    n = az_window(az, 359.5f, 2.0f, win);
    cnt = (win[0][1] - win[0][0]) + (n == 2 ? win[1][1] - win[1][0] : 0);
    CHECK(n == 2 && cnt == 5, "az_window wrap high");
    n = az_window(az, 90.0f, 200.0f, win);
    CHECK(n == 1 && win[0][0] == 0 && win[0][1] == 360, "az_window full circle");
  }

  // --- height_at ---
  {
    auto s = make_synth(1.0f, 0.5f);
    float tane = std::tan(1.0f * deg2rad);
    CHECK_NEAR(height_at(s.g, 10000.0f), 10000.0f * tane, 0.5f, "height_at mid");
    CHECK_NEAR(height_at(s.g, 10.0f), 125.0f * tane, 0.5f, "height_at near clamps");
    CHECK(std::isnan(height_at(s.g, 60000.0f)), "height_at beyond range is NaN");
  }

  // --- auto_roi ---
  {
    auto s = make_synth(0.5f, 0.5f);
    std::vector<sweep_geometry> gs{s.g};
    float roi = auto_roi(gs, 100000.0f);
    CHECK_NEAR(roi, 0.6f * (float)(1.0 * deg2rad) * 100000.0f, 1.0f, "auto_roi 1deg spacing");
  }

  // --- interp_sweep: all radius methods vs brute force + analytic ---
  {
    auto s = make_synth(0.5f, 90.5f);  // start at 90.5 to exercise remapping
    float roi = 1200.0f;

    std::vector<float> cb = {45.0f, 359.8f, 0.1f, 200.0f, 123.4f};
    std::vector<float> cr = {20000.0f, 30000.0f, 15000.0f, 49000.0f, 800.0f};
    cell_lookup cells{cb.data(), cr.data(), cb.size()};

    std::vector<const float*> fd{s.data.data()};
    std::vector<std::vector<float>> out(1, std::vector<float>(cells.n));

    struct { weight_type wt; const char* name; float tol_analytic; } cases[] = {
      {weight_type::barnes,   "barnes",   0.05f},
      {weight_type::cressman, "cressman", 0.05f},
      {weight_type::idw,      "idw",      0.05f},
    };
    for (auto& c : cases) {
      interp_sweep(s.g, fd, s.nbins, cells, c.wt, roi, 2.0f, UND, out);
      for (size_t i = 0; i < cells.n; ++i) {
        float ref = brute_force(s, cb[i], cr[i], c.wt, roi, 2.0f);
        CHECK_NEAR(out[0][i], ref, 1e-3f, (std::string(c.name) + " vs brute force").c_str());
        CHECK_NEAR(out[0][i], f_lin(cb[i], cr[i]), c.tol_analytic,
                   (std::string(c.name) + " vs analytic").c_str());
      }
    }

    // beyond max range -> NaN
    std::vector<float> cb2 = {10.0f}, cr2 = {60000.0f};
    cell_lookup far{cb2.data(), cr2.data(), 1};
    std::vector<std::vector<float>> out2(1, std::vector<float>(1));
    interp_sweep(s.g, fd, s.nbins, far, weight_type::barnes, roi, 2.0f, UND, out2);
    CHECK(std::isnan(out2[0][0]), "radius: beyond range is NaN");
  }

  // --- interp_sweep: bilinear ---
  {
    auto s = make_synth(0.5f, 0.5f);
    std::vector<float> cb = {45.25f, 359.9f, 0.0f, 180.0f};
    std::vector<float> cr = {20111.0f, 30000.0f, 15000.0f, 49000.0f};
    cell_lookup cells{cb.data(), cr.data(), cb.size()};
    std::vector<const float*> fd{s.data.data()};
    std::vector<std::vector<float>> out(1, std::vector<float>(cells.n));
    interp_sweep(s.g, fd, s.nbins, cells, weight_type::bilinear, 0.0f, 2.0f, UND, out);
    for (size_t i = 0; i < cells.n; ++i)
      CHECK_NEAR(out[0][i], f_lin(cb[i], cr[i]), 0.02f, "bilinear vs analytic");

    // all 4 corners invalid -> NaN; partial corners -> renormalised value
    auto s2 = make_synth(0.5f, 0.5f);
    for (size_t k = 0; k < s2.nbins; ++k) {
      s2.data[10 * s2.nbins + k] = core_nan;   // ray az=10.5
      s2.data[11 * s2.nbins + k] = UND;        // ray az=11.5
    }
    std::vector<float> cb3 = {11.0f, 10.9f}, cr3 = {20000.0f, 20000.0f};
    cell_lookup c3{cb3.data(), cr3.data(), 2};
    std::vector<const float*> fd2{s2.data.data()};
    std::vector<std::vector<float>> out3(1, std::vector<float>(2));
    interp_sweep(s2.g, fd2, s2.nbins, c3, weight_type::bilinear, 0.0f, 2.0f, UND, out3);
    CHECK(std::isnan(out3[0][0]), "bilinear: all corners invalid -> NaN");
    CHECK(std::isnan(out3[0][1]), "bilinear: all corners invalid -> NaN (2)");

    std::vector<float> cb4 = {10.05f}, cr4 = {20000.0f};  // 95% weight on valid ray 9.5
    cell_lookup c4{cb4.data(), cr4.data(), 1};
    interp_sweep(s2.g, fd2, s2.nbins, c4, weight_type::bilinear, 0.0f, 2.0f, UND, out3);
    CHECK_NEAR(out3[0][0], f_lin(9.5f, 20000.0f), 0.05f, "bilinear: renormalised over valid");
  }

  // --- fill_heights ---
  {
    std::vector<std::vector<float>> h = {
      {core_nan, 100.0f, core_nan},
      {200.0f,   200.0f, core_nan},
      {core_nan, 300.0f, core_nan},
    };
    size_t unfilled = fill_heights(h, 3);
    CHECK(unfilled == 1, "fill_heights: one column outside range");
    CHECK_NEAR(h[0][0], 200.0f, 0.01f, "fill_heights bottom fill");
    CHECK_NEAR(h[2][0], 200.0f, 0.01f, "fill_heights top fill");
    CHECK_NEAR(h[0][1], 100.0f, 0.01f, "fill_heights untouched");
    CHECK(std::isnan(h[1][2]), "fill_heights all-NaN stays NaN");
  }

  // --- vertical_slice ---
  {
    std::vector<std::vector<float>> h = {{100.0f, 100.0f, 100.0f, 500.0f, 500.0f},
                                         {300.0f, 300.0f, 300.0f, 500.0f, 500.0f}};
    std::vector<std::vector<float>> v = {{10.0f, 10.0f, core_nan, 1.0f, 7.0f},
                                         {30.0f, 30.0f, 30.0f,    2.0f, 8.0f}};
    float out[5];
    vertical_slice(h, v, 200.0f, out, 5);
    CHECK_NEAR(out[0], 20.0f, 1e-4f, "vslice midpoint");
    CHECK(std::isnan(out[2]), "vslice NaN endpoint propagates");
    CHECK(std::isnan(out[3]), "vslice z below duplicate span");
    vertical_slice(h, v, 50.0f, out, 5);
    CHECK(std::isnan(out[0]), "vslice below lowest -> NaN");
    vertical_slice(h, v, 400.0f, out, 5);
    CHECK(std::isnan(out[1]), "vslice above highest -> NaN");
    vertical_slice(h, v, 500.0f, out, 5);
    CHECK_NEAR(out[4], 7.0f, 1e-4f, "vslice duplicate height takes lower");
    vertical_slice(h, v, 300.0f, out, 5);
    CHECK_NEAR(out[0], 30.0f, 1e-4f, "vslice exact top of bracket");
  }

  // --- randomised brute-force sweep (barnes) ---
  {
    auto s = make_synth(1.5f, 33.5f, 180, 120);  // 2 deg spacing sweep
    // poke some invalid gates
    std::mt19937 rng(42);
    std::uniform_int_distribution<size_t> ri(0, s.data.size() - 1);
    for (int i = 0; i < 2000; ++i) s.data[ri(rng)] = (i % 2) ? core_nan : UND;

    std::uniform_real_distribution<float> rb(0.0f, 360.0f), rr(500.0f, 32000.0f);
    std::vector<float> cb(50), cr(50);
    for (int i = 0; i < 50; ++i) { cb[i] = rb(rng); cr[i] = rr(rng); }
    cell_lookup cells{cb.data(), cr.data(), 50};
    std::vector<const float*> fd{s.data.data()};
    std::vector<std::vector<float>> out(1, std::vector<float>(50));

    for (auto wt : {weight_type::barnes, weight_type::cressman, weight_type::idw}) {
      interp_sweep(s.g, fd, s.nbins, cells, wt, 2500.0f, 2.0f, UND, out);
      for (int i = 0; i < 50; ++i) {
        float ref = brute_force(s, cb[i], cr[i], wt, 2500.0f, 2.0f);
        if (std::isnan(ref)) CHECK(std::isnan(out[0][i]), "random: NaN agreement");
        else CHECK_NEAR(out[0][i], ref, 1e-3f, "random vs brute force");
      }
    }
  }

  if (failures == 0) std::printf("ALL TESTS PASSED\n");
  else std::printf("%d FAILURES\n", failures);
  return failures == 0 ? 0 : 1;
}
