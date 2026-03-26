#ifndef VARGRID_GRID_H
#define VARGRID_GRID_H

#include "types.h"

// Initialise altitude levels from config (altitude_base, altitude_step, layer_count).
auto init_altitudes(io::configuration const& config) -> array1f;

// Reverse a 1D array in-place.
auto flip(array1d& data) -> void;

// Flip a 2D array vertically (top-bottom) in-place.
auto flipud(array2f& data) -> void;

#endif // VARGRID_GRID_H
