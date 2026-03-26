#ifndef VARGRID_WRITER_H
#define VARGRID_WRITER_H

#include "types.h"

// Packing parameters for CF-compliant scale_factor/add_offset compression.
// value = packed_value * scale_factor + add_offset
// Stored as int16 with _FillValue = -32768.
struct packing_info {
  float valid_min;     // Minimum expected physical value
  float valid_max;     // Maximum expected physical value
  float scale_factor;  // Computed from range / 65534
  float add_offset;    // Computed from midpoint of range
  bool enabled;        // Whether this field should be packed

  // Compute scale_factor and add_offset from valid_min/valid_max.
  void compute() {
    // Map [valid_min, valid_max] to [-32767, +32767] (65534 levels)
    // Reserve -32768 for _FillValue
    scale_factor = (valid_max - valid_min) / 65534.0f;
    add_offset = valid_min + 32767.0f * scale_factor;
  }
};

// Get packing parameters for a known ODIM quantity.
// Returns enabled=false for unknown quantities (written as float32).
auto get_packing_info(const std::string& quantity) -> packing_info;

// Set CF-compliant attributes on a data variable based on the ODIM quantity name.
auto set_cf_field_attributes(io::nc::variable& var, const std::string& quantity) -> void;

// Set CF-compliant attributes on a coordinate variable.
auto set_cf_coord_attributes(io::nc::variable& var, const std::string& coord_type) -> void;

// All the information needed to write to the output NetCDF file.
struct output_context {
  io::nc::file* file;

  // Dimensions
  io::nc::dimension* dim_x;
  io::nc::dimension* dim_y;
  io::nc::dimension* dim_z;
  io::nc::dimension* dim_e;
  io::nc::dimension* dim_nrad;

  // Per-field data variables and packing info
  std::map<std::string, io::nc::variable*> data_vars;
  std::map<std::string, io::nc::variable*> nobs_vars;
  std::map<std::string, packing_info> packing;

  // Whether packing is enabled globally
  bool use_packing = false;

  // Write a gridded field, packing if appropriate.
  void write_field(const std::string& field, const array2f& data, size_t layer_index);
};

// Create the output NetCDF file with dimensions, coordinates, and data variables.
auto create_output_file(
      const std::filesystem::path& path
    , grid_coordinates const& coords
    , array1d const& y_edges
    , array2f const& lon
    , array2f const& lat
    , array1f const& altitudes
    , volume_metadata const& meta
    , const std::string& proj4_string
    , const std::string& method
    , const std::vector<std::string>& fields
    , bool output_obs_count
    , bool pack_output
    , array1d const& radar_lat
    , array1d const& radar_lon
    , array1d const& radar_alt
    ) -> std::pair<io::nc::file, output_context>;

#endif // VARGRID_WRITER_H