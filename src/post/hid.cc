#include "hid.h"
#include <cmath>
#include <algorithm>
#include <array>
#include <limits>

using namespace bom;

// Membership beta function (Dolan & Rutledge, 2009):
//   beta(x) = 1 / (1 + ((x - m) / a)^2)^b
static inline float mbf(float x, float m, float a, float b) {
  if (std::isnan(x)) return 0.0f;
  float t = (x - m) / a;
  return 1.0f / (1.0f + std::pow(t * t, b));
}

// Number of hydrometeor types
constexpr int N_TYPES = 10;

// Type indices (1-based in output, 0-based in arrays)
// 0=Drizzle, 1=Rain, 2=Ice Crystals, 3=Aggregates, 4=Wet Snow,
// 5=Vertical Ice, 6=LD Graupel, 7=HD Graupel, 8=Hail, 9=Big Drops

// Beta function parameters: {m, a, b} for each of N_TYPES species.
// S-band summer values from Dolan & Rutledge (2009) / CSU-HID.
// These are the published values for S-band radar.

struct mbf_params {
  float m, a, b;
};

struct variable_set {
  std::array<mbf_params, N_TYPES> params;
};

// S-band Reflectivity (Zh) beta function parameters [m, a, b]
// Drz    Rain   IceCr  Agg    WSnow  VIce   LDGr   HDGr   Hail   BDrop
static constexpr variable_set Zh_set = {{
  {{ 25.0f, 10.0f, 3.0f}},  // Drizzle
  {{ 35.0f, 15.0f, 3.0f}},  // Rain
  {{ 10.0f, 10.0f, 1.5f}},  // Ice Crystals
  {{ 25.0f, 15.0f, 2.0f}},  // Aggregates
  {{ 30.0f, 15.0f, 2.0f}},  // Wet Snow
  {{ 15.0f, 10.0f, 1.5f}},  // Vertical Ice
  {{ 35.0f, 15.0f, 2.0f}},  // LD Graupel
  {{ 45.0f, 10.0f, 3.0f}},  // HD Graupel
  {{ 55.0f, 15.0f, 3.0f}},  // Hail
  {{ 40.0f, 10.0f, 3.0f}},  // Big Drops
}};

// S-band Differential Reflectivity (Zdr) [m, a, b]
static constexpr variable_set Zdr_set = {{
  {{ 0.2f,  0.5f,  2.0f}},  // Drizzle
  {{ 1.5f,  1.5f,  2.0f}},  // Rain
  {{ 0.5f,  1.0f,  1.5f}},  // Ice Crystals
  {{ 0.3f,  0.5f,  2.0f}},  // Aggregates
  {{ 0.8f,  1.0f,  1.5f}},  // Wet Snow
  {{-0.2f,  0.5f,  2.0f}},  // Vertical Ice
  {{ 0.3f,  0.5f,  2.0f}},  // LD Graupel
  {{ 0.5f,  0.5f,  2.0f}},  // HD Graupel
  {{ 0.0f,  1.5f,  1.5f}},  // Hail
  {{ 2.5f,  1.0f,  2.0f}},  // Big Drops
}};

// S-band Specific Differential Phase (Kdp) [m, a, b]
static constexpr variable_set Kdp_set = {{
  {{ 0.0f,  0.3f,  2.0f}},  // Drizzle
  {{ 0.7f,  1.0f,  2.0f}},  // Rain
  {{ 0.0f,  0.3f,  1.5f}},  // Ice Crystals
  {{ 0.0f,  0.3f,  2.0f}},  // Aggregates
  {{ 0.1f,  0.5f,  1.5f}},  // Wet Snow
  {{-0.1f,  0.3f,  2.0f}},  // Vertical Ice
  {{ 0.1f,  0.5f,  2.0f}},  // LD Graupel
  {{ 0.5f,  0.5f,  2.0f}},  // HD Graupel
  {{-0.5f,  2.0f,  1.5f}},  // Hail
  {{ 0.3f,  0.5f,  2.0f}},  // Big Drops
}};

// S-band Correlation Coefficient (RhoHV) [m, a, b]
static constexpr variable_set Rho_set = {{
  {{ 0.99f, 0.01f, 3.0f}},  // Drizzle
  {{ 0.98f, 0.02f, 3.0f}},  // Rain
  {{ 0.99f, 0.02f, 2.0f}},  // Ice Crystals
  {{ 0.97f, 0.04f, 2.0f}},  // Aggregates
  {{ 0.90f, 0.05f, 2.0f}},  // Wet Snow
  {{ 0.98f, 0.02f, 2.0f}},  // Vertical Ice
  {{ 0.97f, 0.03f, 2.0f}},  // LD Graupel
  {{ 0.96f, 0.03f, 2.0f}},  // HD Graupel
  {{ 0.92f, 0.06f, 2.0f}},  // Hail
  {{ 0.95f, 0.04f, 2.0f}},  // Big Drops
}};

// S-band Temperature (T, degrees C) [m, a, b]
static constexpr variable_set T_set = {{
  {{  5.0f, 10.0f, 1.5f}},  // Drizzle
  {{  5.0f, 15.0f, 1.5f}},  // Rain
  {{-15.0f, 15.0f, 1.0f}},  // Ice Crystals
  {{-15.0f, 15.0f, 1.0f}},  // Aggregates
  {{  0.0f,  5.0f, 2.0f}},  // Wet Snow
  {{-10.0f, 15.0f, 1.0f}},  // Vertical Ice
  {{-10.0f, 15.0f, 1.0f}},  // LD Graupel
  {{-10.0f, 15.0f, 1.0f}},  // HD Graupel
  {{-10.0f, 20.0f, 1.0f}},  // Hail
  {{  5.0f, 10.0f, 1.5f}},  // Big Drops
}};

auto hid_classifier::process(
    std::map<std::string, array2f>& layer_data,
    const post_processor_context& ctx,
    const io::configuration& config) -> void
{
  size_t ny = ctx.ny;
  size_t nx = ctx.nx;
  size_t n = nx * ny;

  const float* dz  = layer_data.at("DBZH").data();
  const float* zdr = layer_data.at("ZDR").data();
  const float* rho = layer_data.at("RHOHV").data();
  const float* kdp = layer_data.at("KDP").data();

  // Temperature is optional
  const float* temp = nullptr;
  bool use_temp = (layer_data.count("TEMPERATURE") > 0);
  if (use_temp) temp = layer_data.at("TEMPERATURE").data();

  // Weights (configurable)
  float w_dz = std::stof(std::string(config.optional("hid_weight_dz", "1.5")));
  float w_dr = std::stof(std::string(config.optional("hid_weight_dr", "0.8")));
  float w_kd = std::stof(std::string(config.optional("hid_weight_kd", "1.0")));
  float w_rh = std::stof(std::string(config.optional("hid_weight_rh", "0.8")));
  float w_t  = std::stof(std::string(config.optional("hid_weight_t",  "0.4")));

  // Hybrid method (default): pol vars are weighted sum, then multiplied by Z and T.
  // weight_sum = sum of pol var weights (DR + KD + RH)
  float weight_sum = w_dr + w_kd + w_rh;

  // Output: HID class (1-10), 0 = unclassified
  auto hid = array2f{vec2z{nx, ny}};
  for (size_t i = 0; i < n; ++i)
    hid.data()[i] = 0.0f;

  // Score buffer for all types at one pixel
  std::array<float, N_TYPES> scores;

  for (size_t idx = 0; idx < n; ++idx) {
    float zh = dz[idx];
    if (std::isnan(zh)) continue;

    float zd = zdr[idx];
    float kd = kdp[idx];
    float rh = rho[idx];

    // Skip if essential pol data is missing
    if (std::isnan(zd) && std::isnan(kd) && std::isnan(rh)) continue;

    for (int c = 0; c < N_TYPES; ++c) {
      // Hybrid method: weighted sum of pol variables, then multiply by Z and T

      // Pol variable weighted sum
      float pol_sum = 0.0f;
      float pol_wt = 0.0f;

      if (!std::isnan(zd)) {
        pol_sum += w_dr * mbf(zd, Zdr_set.params[c].m, Zdr_set.params[c].a, Zdr_set.params[c].b);
        pol_wt += w_dr;
      }
      if (!std::isnan(kd)) {
        pol_sum += w_kd * mbf(kd, Kdp_set.params[c].m, Kdp_set.params[c].a, Kdp_set.params[c].b);
        pol_wt += w_kd;
      }
      if (!std::isnan(rh)) {
        pol_sum += w_rh * mbf(rh, Rho_set.params[c].m, Rho_set.params[c].a, Rho_set.params[c].b);
        pol_wt += w_rh;
      }

      float test = (pol_wt > 0.0f) ? pol_sum / pol_wt : 0.0f;

      // Multiply by Z membership
      test *= mbf(zh, Zh_set.params[c].m, Zh_set.params[c].a, Zh_set.params[c].b);

      // Multiply by T membership (if available)
      if (use_temp && !std::isnan(temp[idx])) {
        // Convert from Kelvin to Celsius if needed
        float t_c = temp[idx];
        if (t_c > 100.0f) t_c -= 273.15f;  // assume Kelvin if > 100
        test *= mbf(t_c, T_set.params[c].m, T_set.params[c].a, T_set.params[c].b);
      }

      scores[c] = test;
    }

    // Argmax → HID class (1-based)
    int best = 0;
    float best_score = scores[0];
    for (int c = 1; c < N_TYPES; ++c) {
      if (scores[c] > best_score) {
        best_score = scores[c];
        best = c;
      }
    }

    if (best_score > 0.0f)
      hid.data()[idx] = static_cast<float>(best + 1);
  }

  layer_data["HID"] = std::move(hid);
}