#ifndef VARGRID_GRIDDING_VARIATIONAL_H
#define VARGRID_GRIDDING_VARIATIONAL_H

#include "config.h"
#include "observation_operator.h"
#include "solver.h"
#include "../types.h"

// Precompute bearing and range from radar to each grid cell (main thread).
auto precompute_grid_bearings(
      latlonalt const& radar_location
    , array2<latlon> const& latlons
    ) -> grid_bearings;

// Build observation operator for a single altitude layer (thread-safe).
auto build_observation_operator(
      volume const& vol
    , grid_bearings const& gb
    , size_t grid_nx
    , size_t grid_ny
    , float altitude
    , float grid_spacing
    , const vargrid_config& cfg
    ) -> observation_operator;

// Variational gridding from a precomputed observation operator.
auto variational_grid(
      const observation_operator& H
    , const vargrid_config& cfg
    ) -> array2f;

// Observation count per grid cell (diagnostic).
auto observation_density(const observation_operator& H) -> array2f;

#endif // VARGRID_GRIDDING_VARIATIONAL_H
