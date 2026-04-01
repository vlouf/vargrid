#ifndef VARGRID_HID_H
#define VARGRID_HID_H

#include "post_processor.h"

// CSU Fuzzy-logic Hydrometeor Classification (summer, S-band).
//
// Implements the hybrid method from Dolan & Rutledge (2009) using
// membership beta functions for 10 hydrometeor types.
//
// Required fields: DBZH, ZDR, RHOHV, KDP (optional: TEMPERATURE)
// Produced fields: HID
//
// HID classes:
//   1 = Drizzle,  2 = Rain,  3 = Ice Crystals,  4 = Aggregates,
//   5 = Wet Snow, 6 = Vertical Ice, 7 = Low-Density Graupel,
//   8 = High-Density Graupel, 9 = Hail, 10 = Big Drops
//
// Config parameters:
//   hid_weight_dz  — weight for Zh     (default 1.5)
//   hid_weight_dr  — weight for Zdr    (default 0.8)
//   hid_weight_kd  — weight for Kdp    (default 1.0)
//   hid_weight_rh  — weight for RhoHV  (default 0.8)
//   hid_weight_t   — weight for T      (default 0.4)
class hid_classifier : public post_processor {
public:
  auto name() const -> std::string override { return "hid"; }
  auto required_fields() const -> std::vector<std::string> override {
    return {"DBZH", "ZDR", "RHOHV", "KDP"};
  }
  auto provided_fields() const -> std::vector<std::string> override {
    return {"HID"};
  }

  auto process(
      std::map<std::string, bom::array2f>& layer_data
    , const post_processor_context& ctx
    , const bom::io::configuration& config
    ) -> void override;
};

#endif // VARGRID_HID_H