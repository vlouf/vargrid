#include "post_processor.h"
#include "steiner.h"
#include "hid.h"
#include <sstream>

using namespace bom;

// Registry of available post-processors.
// Add new processors here as they are implemented.
static auto make_processor(const std::string& name) -> std::unique_ptr<post_processor>
{
  if (name == "steiner") return std::make_unique<steiner_classifier>();
  if (name == "hid")     return std::make_unique<hid_classifier>();

  throw std::runtime_error("Unknown post-processor: '" + name + "'. "
    "Available: steiner, hid");
}

auto create_post_pipeline(const io::configuration& config)
    -> std::vector<std::unique_ptr<post_processor>>
{
  std::vector<std::unique_ptr<post_processor>> pipeline;

  std::string pp_str = config.optional("post_processors", "");
  if (pp_str.empty()) return pipeline;

  // Parse space-separated list of processor names
  std::istringstream iss(pp_str);
  std::string name;
  while (iss >> name) {
    pipeline.push_back(make_processor(name));
  }

  return pipeline;
}