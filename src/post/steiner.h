#ifndef VARGRID_STEINER_H
#define VARGRID_STEINER_H

#include "post_processor.h"

// Steiner et al. (1995) convective/stratiform classification.
//
// Identifies convective cores based on reflectivity peakedness relative
// to a background mean and marks a surrounding convective radius.
//
// Required fields: DBZH
// Produced fields: STEINER_CLASS (0=no data, 1=stratiform, 2=convective)
//
// Config parameters:
//   steiner_bkg_rad       — background averaging radius in metres (default 11000)
//   steiner_intense_dbz   — intensity threshold in dBZ (default 42)
//   steiner_area_relation — "small", "medium", "large", or "scp" (default "medium")
//   steiner_peak_relation — "default" or "sgp" (default "default")
class steiner_classifier : public post_processor {
public:
  auto name() const -> std::string override { return "steiner"; }
  auto required_fields() const -> std::vector<std::string> override { return {"DBZH"}; }
  auto provided_fields() const -> std::vector<std::string> override { return {"STEINER_CLASS"}; }

  auto process(
      std::map<std::string, bom::array2f>& layer_data
    , const post_processor_context& ctx
    , const bom::io::configuration& config
    ) -> void override;

private:
  static auto get_conv_radius(float bkg, const std::string& relation) -> float;
  static auto get_peakedness(float bkg, const std::string& relation) -> float;
  static auto fill_convective(
      bom::array2f& sclass, int ci, int cj, float rad,
      const post_processor_context& ctx) -> void;
};

#endif // VARGRID_STEINER_H