#ifndef VARGRID_READER_H
#define VARGRID_READER_H

#include "types.h"

// Discover all available moment names in an ODIM volume.
// Returns field names in the order they are first encountered.
auto discover_fields(io::odim::polar_volume const& vol_odim) -> std::vector<std::string>;

// Read a single moment from an ODIM volume into a volume struct.
// velocity_field: name of the velocity moment (for undetect handling).
auto read_moment_volume(
      io::odim::polar_volume const& vol_odim
    , const std::string& moment
    , const std::string& velocity_field = "VRADH"
    ) -> volume;

// Extract volume-level metadata from an ODIM volume.
auto read_metadata(io::odim::polar_volume const& vol_odim) -> volume_metadata;

#endif // VARGRID_READER_H
