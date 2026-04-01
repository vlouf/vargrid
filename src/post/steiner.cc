#include "post_processor.h"
#include <cmath>

namespace bom {

class steiner_classifier : public post_processor {
public:
    auto name() const -> std::string override { return "steiner"; }
    auto required_fields() const -> std::vector<std::string> override { return {"DBZH"}; }
    auto provided_fields() const -> std::vector<std::string> override { return {"STEINER_CLASS"}; }
    auto process(
        std::map<std::string, array2f>& layer_data, 
        const post_processor_context& ctx, 
        const io::configuration& config) -> void override 
    {
        const auto& ze_dbz = layer_data.at("DBZH");
        auto [ny, nx] = ze_dbz.extents();
        
        // 0: No data/Undetect, 1: Stratiform, 2: Convective
        array2f sclass{vec2z{ny, nx}};        

        // Get parameters from config
        float bkg_rad_m = std::stof(config.optional("steiner_bkg_rad", "11000.0"));
        float intense_thr = std::stof(config.optional("steiner_intense_dbz", "42.0"));
        std::string area_rel = config.optional("steiner_area_relation", "medium");
        std::string peak_rel = config.optional("steiner_peak_relation", "default");

        // Pre-calculate search windows in pixels
        int i_rad = std::ceil(bkg_rad_m / ctx.dx);
        int j_rad = std::ceil(bkg_rad_m / ctx.dy);

        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                float val = ze_dbz[j][i];
                if (std::isnan(val) || val <= undetect) continue;

                // 1. Calculate Background Reflectivity (Linear Average)
                double sum_linear = 0.0;
                int count = 0;
                
                for (int mj = std::max(0, j-j_rad); mj <= std::min((int)ny-1, j+j_rad); ++mj) {
                    for (int mi = std::max(0, i-i_rad); mi <= std::min((int)nx-1, i+i_rad); ++mi) {
                        float neighbor = ze_dbz[mj][mi];
                        if (std::isnan(neighbor)) continue;
                        
                        // Distance check
                        float dist_sq = std::pow((mi-i)*ctx.dx, 2) + std::pow((mj-j)*ctx.dy, 2);
                        if (dist_sq <= bkg_rad_m * bkg_rad_m) {
                            sum_linear += std::pow(10.0, neighbor / 10.0);
                            count++;
                        }
                    }
                }

                float ze_bkg = 10.0f * std::log10(sum_linear / count);

                // 2. Convective Radius Calculation
                float conv_rad = get_conv_radius(ze_bkg, area_rel);
                
                // 3. Peakedness Check
                float peak_req = get_peakedness(ze_bkg, peak_rel);

                if (val >= intense_thr || (val - ze_bkg) >= peak_req) {
                    // Mark as convective and fill surrounding radius
                    fill_convective(sclass, i, j, conv_rad, ctx);
                } else if (sclass[j][i] == 0) {
                    sclass[j][i] = 1; // Stratiform
                }
            }
        }
        layer_data["STEINER_CLASS"] = std::move(sclass);
    }

private:
    float get_conv_radius(float bkg, const std::string& rel) {
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
        } else {
            throw std::runtime_error("Invalid peak_relation: " + std::string(rel));
        }
    }

    float get_peakedness(float bkg, const std::string& rel) {
        if (bkg < 0) return (rel == "sgp") ? 14.0f : 10.0f;
        if (bkg >= 42.43) return (rel == "sgp") ? 4.0f : 0.0f;
        float base = (rel == "sgp") ? 14.0f : 10.0f;
        return base - std::pow(bkg, 2) / 180.0f;
    }

    void fill_convective(array2f& sclass, int ci, int cj, float rad, const post_processor_context& ctx) {
        int i_rad = std::ceil(rad / ctx.dx);
        int j_rad = std::ceil(rad / ctx.dy);
        auto [ny, nx] = sclass.extents();

        for (int j = std::max(0, cj-j_rad); j <= std::min((int)ny-1, cj+j_rad); ++j) {
            for (int i = std::max(0, ci-i_rad); i <= std::min((int)nx-1, ci+i_rad); ++i) {
                float dist_sq = std::pow((i-ci)*ctx.dx, 2) + std::pow((j-cj)*ctx.dy, 2);
                if (dist_sq <= rad * rad) {
                    sclass[j][i] = 2;
                }
            }
        }
    }
};

} // namespace bom