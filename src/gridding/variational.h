#ifndef VARGRID_GRIDDING_VARIATIONAL_H
#define VARGRID_GRIDDING_VARIATIONAL_H

#include "config.h"
#include "observation_operator.h"
#include "solver.h"
#include "../types.h"

// Precompute projected coordinates for all radar gates (main thread only).
auto precompute_gate_projections(
      volume const& vol
    , string const& proj4_string
    ) -> gate_projections;

// Build observation operator for a single altitude layer (thread-safe).
auto build_observation_operator(
      volume const& vol
    , gate_projections const& gp
    , double x0, double y0, double dx, double dy
    , size_t grid_nx, size_t grid_ny
    , float altitude
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
