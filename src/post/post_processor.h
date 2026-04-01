#ifndef VARGRID_POST_PROCESSOR_H
#define VARGRID_POST_PROCESSOR_H

#include "../types.h"
#include <map>
#include <string>
#include <vector>
#include <memory>

namespace bom {

struct post_processor_context {
    float altitude;
    float dx;
    float dy;
};

class post_processor {
public:
    virtual ~post_processor() = default;
    
    virtual auto name() const -> std::string = 0;
    virtual auto required_fields() const -> std::vector<std::string> = 0;
    virtual auto provided_fields() const -> std::vector<std::string> = 0;

    virtual auto process(
        std::map<std::string, array2f>& layer_data, 
        const post_processor_context& ctx, 
        const io::configuration& config
    ) -> void = 0;
};

// Factory function to create the pipeline based on config string
auto create_pipeline(const io::configuration& config) 
    -> std::vector<std::unique_ptr<post_processor>>;

} // namespace bom

#endif