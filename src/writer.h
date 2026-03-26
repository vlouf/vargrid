#ifndef VARGRID_WRITER_H
#define VARGRID_WRITER_H

#include "types.h"

// Set CF-compliant attributes on a data variable based on the ODIM quantity name.
auto set_cf_field_attributes(io::nc::variable& var, const std::string& quantity) -> void;

// Set CF-compliant attributes on a coordinate variable.
auto set_cf_coord_attributes(io::nc::variable& var, const std::string& coord_type) -> void;

// All the information needed to create the output NetCDF file.
struct output_context {
  io::nc::file* file;

  // Dimensions (pointers into the file, valid for the file's lifetime)
  io::nc::dimension* dim_x;
  io::nc::dimension* dim_y;
  io::nc::dimension* dim_z;
  io::nc::dimension* dim_e;
  io::nc::dimension* dim_nrad;

  // Per-field data variables
  std::map<std::string, io::nc::variable*> data_vars;
  std::map<std::string, io::nc::variable*> nobs_vars;
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
    , array1d const& radar_lat
    , array1d const& radar_lon
    , array1d const& radar_alt
    ) -> std::pair<io::nc::file, output_context>;

#endif // VARGRID_WRITER_H
