#include "steiner.h"
#include <cmath>
#include <algorithm>

using namespace bom;

auto steiner_classifier::process(
    std::map<std::string, array2f>& layer_data,
    const post_processor_context& ctx,
    const io::configuration& config) -> void
{
  const auto& ze_dbz = layer_data.at("DBZH");
  size_t ny = ctx.ny;
  size_t nx = ctx.nx;

  // 0: No data, 1: Stratiform, 2: Convective
  auto sclass = array2f{vec2z{nx, ny}};
  for (size_t i = 0; i < sclass.size(); ++i)
    sclass.data()[i] = 0.0f;

  // Parameters from config
  float bkg_rad_m   = std::stof(std::string(config.optional("steiner_bkg_rad", "11000.0")));
  float intense_thr = std::stof(std::string(config.optional("steiner_intense_dbz", "42.0")));
  std::string area_rel = config.optional("steiner_area_relation", "medium");
  std::string peak_rel = config.optional("steiner_peak_relation", "default");

  float bkg_rad_sq = bkg_rad_m * bkg_rad_m;
  int i_rad = static_cast<int>(std::ceil(bkg_rad_m / ctx.dx));
  int j_rad = static_cast<int>(std::ceil(bkg_rad_m / ctx.dy));

  for (size_t j = 0; j < ny; ++j) {
    for (size_t i = 0; i < nx; ++i) {
      float val = ze_dbz[j][i];
      if (std::isnan(val) || val <= undetect) continue;

      // 1. Background reflectivity (linear average in Z-space)
      double sum_linear = 0.0;
      int count = 0;

      int j0 = std::max(0, static_cast<int>(j) - j_rad);
      int j1 = std::min(static_cast<int>(ny) - 1, static_cast<int>(j) + j_rad);
      int i0 = std::max(0, static_cast<int>(i) - i_rad);
      int i1 = std::min(static_cast<int>(nx) - 1, static_cast<int>(i) + i_rad);

      for (int mj = j0; mj <= j1; ++mj) {
        for (int mi = i0; mi <= i1; ++mi) {
          float neighbor = ze_dbz[mj][mi];
          if (std::isnan(neighbor)) continue;

          float dist_sq = static_cast<float>((mi - (int)i)) * ctx.dx * static_cast<float>((mi - (int)i)) * ctx.dx
                        + static_cast<float>((mj - (int)j)) * ctx.dy * static_cast<float>((mj - (int)j)) * ctx.dy;
          if (dist_sq <= bkg_rad_sq) {
            sum_linear += std::pow(10.0, neighbor / 10.0);
            count++;
          }
        }
      }

      if (count == 0) continue;
      float ze_bkg = 10.0f * static_cast<float>(std::log10(sum_linear / count));

      // 2. Convective radius
      float conv_rad = get_conv_radius(ze_bkg, area_rel);

      // 3. Peakedness check
      float peak_req = get_peakedness(ze_bkg, peak_rel);

      if (val >= intense_thr || (val - ze_bkg) >= peak_req) {
        fill_convective(sclass, static_cast<int>(i), static_cast<int>(j), conv_rad, ctx);
      } else if (sclass[j][i] == 0.0f) {
        sclass[j][i] = 1.0f;  // Stratiform
      }
    }
  }

  layer_data["STEINER_CLASS"] = std::move(sclass);
}

auto steiner_classifier::get_conv_radius(float bkg, const std::string& rel) -> float
{
  if (rel == "small") {
    if (bkg < 30) return 1000; 
    if (bkg < 35) return 2000;
    if (bkg < 40) return 3000; 
    if (bkg < 45) return 4000; 
    return 5000;
  } else if (rel == "medium") {
    if (bkg < 25) return 1000; 
    if (bkg < 30) return 2000;
    if (bkg < 35) return 3000; 
    if (bkg < 40) return 4000; 
    return 5000;
  } else if (rel == "large") {
    if (bkg < 20) return 1000; 
    if (bkg < 25) return 2000;
    if (bkg < 30) return 3000; 
    if (bkg < 35) return 4000; 
    return 5000;
  } else if (rel == "scp") {
    if (bkg < 40) return 0; 
    if (bkg < 45) return 1000;
    if (bkg < 50) return 2000; 
    if (bkg < 55) return 6000; 
    return 8000;
  }
  throw std::runtime_error("Invalid steiner_area_relation: " + rel);
}

auto steiner_classifier::get_peakedness(float bkg, const std::string& rel) -> float
{
  if (bkg < 0)     return (rel == "sgp") ? 14.0f : 10.0f;
  if (bkg >= 42.43f) return (rel == "sgp") ? 4.0f : 0.0f;
  float base = (rel == "sgp") ? 14.0f : 10.0f;
  return base - bkg * bkg / 180.0f;
}

auto steiner_classifier::fill_convective(
    array2f& sclass, int ci, int cj, float rad,
    const post_processor_context& ctx) -> void
{
  if (rad <= 0.0f) return;

  int i_rad = static_cast<int>(std::ceil(rad / ctx.dx));
  int j_rad = static_cast<int>(std::ceil(rad / ctx.dy));
  int ny = static_cast<int>(ctx.ny);
  int nx = static_cast<int>(ctx.nx);
  float rad_sq = rad * rad;

  for (int j = std::max(0, cj - j_rad); j <= std::min(ny - 1, cj + j_rad); ++j) {
    for (int i = std::max(0, ci - i_rad); i <= std::min(nx - 1, ci + i_rad); ++i) {
      float di = static_cast<float>(i - ci) * ctx.dx;
      float dj = static_cast<float>(j - cj) * ctx.dy;
      if (di * di + dj * dj <= rad_sq)
        sclass[j][i] = 2.0f;
    }
  }
}