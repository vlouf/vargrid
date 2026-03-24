#ifndef VARGRID_H
#define VARGRID_H

#include "config.h"
#include "observation_operator.h"
#include "solver.h"

#include "pch.h"
using namespace bom;

// Precompute projected coordinates for all radar gates (main thread only).
auto precompute_gate_projections(
      volume const& vol
    , string const& proj4_string
    ) -> gate_projections;

// Build the observation operator for a single altitude layer (thread-safe).
auto build_observation_operator(
      volume const& vol
    , gate_projections const& gp
    , double x0    // x center of first grid cell (projected meters)
    , double y0    // y center of first grid cell (projected meters)
    , double dx    // x cell spacing (projected meters)
    , double dy    // y cell spacing (projected meters)
    , size_t grid_nx
    , size_t grid_ny
    , float altitude
    , const vargrid_config& cfg
    ) -> observation_operator;

// Variational gridding with a precomputed observation operator.
auto variational_grid(
      const observation_operator& H
    , const vargrid_config& cfg
    ) -> array2f;

// Diagnostics: return the observation count per grid cell.
auto observation_density(const observation_operator& H) -> array2f;

#endif // VARGRID_H