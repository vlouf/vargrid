#ifndef VARGRID_POST_PROCESSOR_H
#define VARGRID_POST_PROCESSOR_H

#include "../types.h"
#include <map>
#include <string>
#include <vector>
#include <memory>

// Context passed to each post-processor for a single altitude layer.
struct post_processor_context {
  size_t nx;
  size_t ny;
  float altitude;
  float dx;       // grid spacing in x (m)
  float dy;       // grid spacing in y (m, positive)
};

// Abstract interface for grid-level post-processing algorithms.
// Each processor reads one or more gridded fields and produces one or more
// new fields. Processors are run after variational/CAPPI gridding and before
// writing to the output file.
class post_processor {
public:
  virtual ~post_processor() = default;

  // Human-readable name (used in logging and config).
  virtual auto name() const -> std::string = 0;

  // Fields that must exist in the layer data before this processor runs.
  virtual auto required_fields() const -> std::vector<std::string> = 0;

  // Fields that this processor adds to the layer data.
  virtual auto provided_fields() const -> std::vector<std::string> = 0;

  // Run the algorithm. Reads from layer_data and inserts new fields into it.
  virtual auto process(
      std::map<std::string, bom::array2f>& layer_data
    , const post_processor_context& ctx
    , const bom::io::configuration& config
    ) -> void = 0;
};

// Create the post-processing pipeline from a config key.
// config key: post_processors "steiner hid"
// Returns an ordered list of processors. Empty if the key is absent.
auto create_post_pipeline(const bom::io::configuration& config)
    -> std::vector<std::unique_ptr<post_processor>>;

#endif // VARGRID_POST_PROCESSOR_H