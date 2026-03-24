#ifndef VARGRID_H
#define VARGRID_H

#include "config.h"
#include "observation_operator.h"
#include "solver.h"

#include "pch.h"
using namespace bom;

// Build the observation operator for a single altitude layer.
auto build_observation_operator(
      volume const& vol
    , string const& proj4_string
    , grid_coordinates const& coords
    , size_t grid_nx
    , size_t grid_ny
    , float altitude
    , const vargrid_config& cfg
    ) -> observation_operator;

// Variational gridding of a radar volume to a single altitude layer.
auto variational_grid(
      volume const& vol
    , string const& proj4_string
    , grid_coordinates const& coords
    , size_t grid_nx
    , size_t grid_ny
    , float altitude
    , const vargrid_config& cfg
    ) -> array2f;

// Variational gridding with a precomputed observation operator.
auto variational_grid(
      const observation_operator& H
    , const vargrid_config& cfg
    ) -> array2f;

// Diagnostics: return the observation count per grid cell.
auto observation_density(const observation_operator& H) -> array2f;

#endif // VARGRID_H