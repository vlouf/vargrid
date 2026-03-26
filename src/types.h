#ifndef VARGRID_TYPES_H
#define VARGRID_TYPES_H

// BOM library includes
#include <bom/io/configuration.h>
#include <bom/io/cf.h>
#include <bom/io/nc.h>
#include <bom/io/odim.h>
#include <bom/radar/beam_propagation.h>
#include <bom/array2.h>
#include <bom/ellipsoid.h>
#include <bom/grid_coordinates.h>
#include <bom/grid_transform.h>
#include <bom/map_projection.h>
#include <bom/trace.h>

// Standard library
#include <algorithm>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

using namespace bom;

// Sentinel values for missing data
constexpr float nodata = std::numeric_limits<float>::quiet_NaN();
constexpr float undetect = -32.0f;

// Radar gate geometry for a single range bin
struct bin_info {
  float slant_range;
  float ground_range;
  float altitude;
};

// A single radar sweep (one elevation angle)
struct sweep {
  radar::beam_propagation beam;
  array1<bin_info>        bins;  // gate geometry at bin centers
  array1<angle>           rays;  // azimuth angles at ray centers
  array2f                 data;  // [ray][bin] data values
};

// A radar volume (collection of sweeps for one moment)
struct volume {
  latlonalt     location;
  vector<sweep> sweeps;
};

// Metadata extracted from an ODIM volume
struct volume_metadata {
  array1f elevation;
  array1f nyquist;
  string source;
  string date;
  string time;
  string lowest_sweep_time;
  float beamwidth;
};

#endif // VARGRID_TYPES_H
