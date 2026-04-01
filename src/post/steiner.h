#ifndef VARGRID_STEINER_H
#define VARGRID_STEINER_H

#include "post_processor.h"

namespace bom {

/**
 * @brief Steiner et al. (1995) Convective/Stratiform Classification.
 * * This algorithm identifies convective cores based on reflectivity 
 * peakedness relative to a background mean and marks a surrounding 
 * convective radius.
 */
class steiner_classifier : public post_processor {
public:
    steiner_classifier() = default;
    virtual ~steiner_classifier() = default;

    // Interface Implementation
    auto name() const -> std::string override { return "steiner"; }
    auto required_fields() const -> std::vector<std::string> override { return {"DBZH"}; }
    auto provided_fields() const -> std::vector<std::string> override { return {"STEINER_CLASS"}; }

    /**
     * @brief Executes the Steiner classification on a 2D grid layer.
     * @param layer_data Map of available fields; must contain "DBZH".
     * @param ctx Metadata about the current grid layer (altitude, spacing).
     * @param config The global configuration object for parameter lookup.
     */
    auto process(
        std::map<std::string, array2f>& layer_data, 
        const post_processor_context& ctx, 
        const io::configuration& config
    ) -> void override;

private:
    // Internal math helpers
    auto get_conv_radius(float bkg, const std::string& relation) -> float;
    auto get_peakedness(float bkg, const std::string& relation) -> float;
    
    auto fill_convective(
        array2f& sclass, 
        int ci, 
        int cj, 
        float rad, 
        const post_processor_context& ctx
    ) -> void;
};

} // namespace bom

#endif // VARGRID_STEINER_H