#ifndef VARGRID_READER_H
#define VARGRID_READER_H

#include "types.h"

#include <filesystem>

// Discover all available moment names from an ODIM volume or a CF/Radial file.
// Returns field names in the order they are first encountered.
auto discover_fields(io::odim::polar_volume const& vol_odim) -> std::vector<std::string>;
auto discover_fields(std::filesystem::path const& input_path) -> std::vector<std::string>;

// Read a single moment from an ODIM volume into a volume struct.
// velocity_field: name of the velocity moment (for undetect handling).
auto read_moment_volume(
      io::odim::polar_volume const& vol_odim
    , const std::string& moment
    , const std::string& velocity_field = "VRADH"
    ) -> volume;

// Read a single moment from a CF/Radial file into a volume struct.
auto read_moment_volume(
      std::filesystem::path const& input_path
    , const std::string& moment
    , const std::string& velocity_field = "VRADH"
    ) -> volume;

// Extract volume-level metadata from an ODIM volume.
auto read_metadata(io::odim::polar_volume const& vol_odim) -> volume_metadata;

// Extract volume-level metadata from a CF/Radial file.
auto read_metadata(std::filesystem::path const& input_path) -> volume_metadata;

// Extract per-field variable metadata (units, names, descriptions) from input.
auto read_field_metadata(
      std::filesystem::path const& input_path
    , const std::vector<std::string>& fields
    ) -> std::map<std::string, variable_metadata>;

#endif // VARGRID_READER_H
