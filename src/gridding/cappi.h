#ifndef VARGRID_GRIDDING_CAPPI_H
#define VARGRID_GRIDDING_CAPPI_H

#include "../types.h"

// Generate a CAPPI at a given altitude using inverse distance weighting.
// This is the traditional gridding method, kept as a baseline for comparison.
auto generate_cappi(
      volume const& vol
    , array2<latlon> const& latlons
    , float max_alt_diff
    , float idw_pwr
    , float altitude
    , float roi = 2500.f
    ) -> array2f;

#endif // VARGRID_GRIDDING_CAPPI_H
