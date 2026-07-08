#include "writer.h"

#ifndef VARGRID_VERSION
#define VARGRID_VERSION "0.1.0"
#endif

namespace {

auto parse_cf_time_seconds(const volume_metadata& meta) -> double
{
  // Prefer lowest sweep timestamp when available (ISO 8601 style).
  if (!meta.lowest_sweep_time.empty()) {
    std::tm tm{};
    int year = 0, mon = 0, mday = 0, hour = 0, min = 0, sec = 0;
    if (std::sscanf(meta.lowest_sweep_time.c_str(), "%d-%d-%dT%d:%d:%d",
        &year, &mon, &mday, &hour, &min, &sec) == 6) {
      tm.tm_year = year - 1900;
      tm.tm_mon = mon - 1;
      tm.tm_mday = mday;
      tm.tm_hour = hour;
      tm.tm_min = min;
      tm.tm_sec = sec;
      return static_cast<double>(timegm(&tm));
    }
  }

  // ODIM fallback: date="YYYYMMDD", time="HHMMSS".
  if (meta.date.size() >= 8 && meta.time.size() >= 6) {
    std::tm tm{};
    tm.tm_year = std::stoi(meta.date.substr(0, 4)) - 1900;
    tm.tm_mon = std::stoi(meta.date.substr(4, 2)) - 1;
    tm.tm_mday = std::stoi(meta.date.substr(6, 2));
    tm.tm_hour = std::stoi(meta.time.substr(0, 2));
    tm.tm_min = std::stoi(meta.time.substr(2, 2));
    tm.tm_sec = std::stoi(meta.time.substr(4, 2));
    return static_cast<double>(timegm(&tm));
  }

  // Last resort: now.
  return static_cast<double>(time(nullptr));
}

} // namespace

// Known packing ranges for ODIM quantities.
// These are fixed ranges that cover the full physical range of each variable.
// Prefix-based lookup: try exact match first, then strip suffixes at '_' boundaries,
// then try base prefix (e.g. "PHIDP_PHIDO" → "PHIDP_PHIDO" → "PHIDP" → found).
template<typename Map>
static auto prefix_lookup(const Map& map, const std::string& key) -> typename Map::const_iterator
{
  // Exact match
  auto it = map.find(key);
  if (it != map.end()) return it;

  // Try stripping from the last '_' backward
  std::string prefix = key;
  while (true) {
    auto pos = prefix.rfind('_');
    if (pos == std::string::npos) break;
    prefix = prefix.substr(0, pos);
    it = map.find(prefix);
    if (it != map.end()) return it;
  }

  return map.end();
}

auto get_packing_info(const std::string& quantity) -> packing_info
{
  static const std::unordered_map<std::string, std::pair<float, float>> ranges = {
    // Reflectivity (dBZ)
    {"DBZH",       {-33.0f,  95.0f}},
    {"DBZV",       {-33.0f,  95.0f}},
    {"TH",         {-33.0f,  95.0f}},
    {"TV",         {-33.0f,  95.0f}},
    // Velocity (m/s)
    {"VRAD",       {-75.0f,  75.0f}},
    {"VRADH",      {-75.0f,  75.0f}},
    {"VRADDH",     {-75.0f,  75.0f}},
    {"VRADV",      {-75.0f,  75.0f}},
    // Spectrum width (m/s)
    {"WRAD",       {0.0f,    30.0f}},
    {"WRADH",      {0.0f,    30.0f}},
    {"WRADV",      {0.0f,    30.0f}},
    // Dual-pol
    {"ZDR",        {-8.0f,   12.0f}},
    {"RHOHV",      {0.0f,     1.1f}},
    {"PHIDP",      {-180.0f, 180.0f}},
    {"KDP",        {-2.0f,   15.0f}},
    // Quality
    {"SQI",        {0.0f,     1.0f}},
    {"SNR",        {-30.0f,  70.0f}},
    {"SNRH",       {-30.0f,  70.0f}},
    {"CCOR",       {-40.0f,  10.0f}},
    {"CCORH",      {-40.0f,  10.0f}},
    // Attenuation
    {"PIA",        {0.0f,    30.0f}},
    // DSD parameters
    {"NW",         {-1.0f,   10.0f}},
    {"D0",         {0.0f,     5.0f}},
    // Precipitation rate
    {"RAIN",       {0.0f,   300.0f}},
    {"SNOW",       {0.0f,   100.0f}},
    // HID (integer-like but store as float)
    {"HID",        {0.0f,    20.0f}},
    // Temperature
    {"TEMPERATURE", {-80.0f,  50.0f}},
  };

  auto it = prefix_lookup(ranges, quantity);
  if (it != ranges.end()) {
    packing_info pi;
    pi.valid_min = it->second.first;
    pi.valid_max = it->second.second;
    pi.enabled = true;
    pi.compute();
    return pi;
  }

  return {0, 0, 0, 0, false};
}

auto set_cf_field_attributes(
      io::nc::variable& var
    , const std::string& quantity
    , const variable_metadata* input_meta
    ) -> void
{
  static const std::unordered_map<std::string, std::tuple<std::string, std::string, std::string>> cf_map = {
    // Reflectivity
    {"DBZH",       {"dBZ",           "equivalent_reflectivity_factor",                          "Horizontal reflectivity"}},
    {"DBZV",       {"dBZ",           "",                                                        "Vertical reflectivity"}},
    {"TH",         {"dBZ",           "",                                                        "Total power horizontal (uncorrected Z)"}},
    {"TV",         {"dBZ",           "",                                                        "Total power vertical (uncorrected Z)"}},
    {"reflectivity", {"dBZ",         "equivalent_reflectivity_factor",                          "Reflectivity"}},
    // Velocity
    {"VRAD",       {"m s-1",         "radial_velocity_of_scatterers_away_from_instrument",      "Mean Doppler velocity"}},
    {"VRADH",      {"m s-1",         "radial_velocity_of_scatterers_away_from_instrument",      "Mean Doppler velocity (H)"}},
    {"VRADDH",     {"m s-1",         "radial_velocity_of_scatterers_away_from_instrument",      "Dealiased Doppler velocity (H)"}},
    {"VRADV",      {"m s-1",         "",                                                        "Mean Doppler velocity (V)"}},
    {"velocity",   {"m s-1",         "radial_velocity_of_scatterers_away_from_instrument",      "Mean Doppler velocity"}},
    // Spectrum width
    {"WRAD",       {"m s-1",         "doppler_spectrum_width",                                  "Doppler spectrum width"}},
    {"WRADH",      {"m s-1",         "doppler_spectrum_width",                                  "Doppler spectrum width (H)"}},
    {"WRADV",      {"m s-1",         "doppler_spectrum_width",                                  "Doppler spectrum width (V)"}},
    {"spectrum_width", {"m s-1",     "doppler_spectrum_width",                                  "Doppler spectrum width"}},
    // Dual-pol
    {"ZDR",        {"dB",            "log_differential_reflectivity_hv",                        "Differential reflectivity"}},
    {"RHOHV",      {"1",             "cross_correlation_ratio_hv",                              "Cross-correlation coefficient"}},
    {"PHIDP",      {"degrees",       "differential_phase_hv",                                   "Differential phase"}},
    {"KDP",        {"degrees km-1",  "specific_differential_phase_hv",                          "Specific differential phase"}},
    {"differential_reflectivity", {"dB", "log_differential_reflectivity_hv",                    "Differential reflectivity"}},
    {"differential_phase", {"degree", "differential_phase_hv",                                  "Differential phase"}},
    // Quality
    {"SQI",        {"1",             "",                                                        "Signal quality index"}},
    {"SNR",        {"dB",            "signal_to_noise_ratio",                                   "Signal-to-noise ratio"}},
    {"SNRH",       {"dB",            "signal_to_noise_ratio",                                   "Signal-to-noise ratio (H)"}},
    {"CCOR",       {"dB",            "",                                                        "Clutter correction"}},
    {"CCORH",      {"dB",            "",                                                        "Clutter correction (H)"}},
    {"signal_to_noise_ratio", {"dB", "signal_to_noise_ratio",                                   "Signal-to-noise ratio"}},
    {"signal_quality_index", {"1",    "",                                                       "Signal quality index"}},
    // Attenuation
    {"PIA",        {"dB",            "",                                                        "Path-integrated attenuation"}},
    // DSD parameters
    {"NW",         {"log10(m-3 mm-1)", "",                                                      "Normalised intercept parameter"}},
    {"D0",         {"mm",            "",                                                        "Median volume diameter"}},
    // Precipitation
    {"RAIN",       {"mm h-1",        "rainfall_rate",                                           "Rainfall rate"}},
    {"SNOW",       {"mm h-1",        "snowfall_rate",                                           "Snowfall rate"}},
    // Classification
    {"HID",        {"1",             "",                                                        "Hydrometeor identification"}},
    // Environment
    {"TEMPERATURE", {"C",            "air_temperature",                                         "Temperature"}},
    {"total_power", {"dBZ",          "",                                                        "Total power"}},
  };

  auto it = prefix_lookup(cf_map, quantity);
  bool has_units = false;
  if (it != cf_map.end()) {
    auto& [units, standard_name, long_name] = it->second;
    if (!units.empty()) {
      var.att_set("units", units);
      has_units = true;
    }
    if (!standard_name.empty()) var.att_set("standard_name", standard_name);
    // Use the matched long_name as a base, append the original quantity if it differs
    if (!long_name.empty()) {
      if (it->first == quantity)
        var.att_set("long_name", long_name);
      else
        var.att_set("long_name", long_name + " (" + quantity + ")");
    }
  } else {
    var.att_set("long_name", quantity);
  }

  // Preserve metadata from the input variable when available.
  if (input_meta) {
    if (!input_meta->units.empty()) {
      var.att_set("units", input_meta->units);
      has_units = true;
    }
    if (!input_meta->standard_name.empty())
      var.att_set("standard_name", input_meta->standard_name);
    if (!input_meta->long_name.empty())
      var.att_set("long_name", input_meta->long_name);
    if (!input_meta->description.empty())
      var.att_set("description", input_meta->description);
  }

  if (!has_units)
    var.att_set("units", "1");

  var.att_set("grid_mapping", "albers_conical_equal_area");
  var.att_set("coordinates", "time latitude longitude");
}

auto set_cf_coord_attributes(io::nc::variable& var, const std::string& coord_type) -> void
{
  static const std::unordered_map<std::string, std::tuple<std::string, std::string, std::string>> attrs = {
    {"latitude",        {"degrees_north", "latitude",                  ""}},
    {"longitude",       {"degrees_east",  "longitude",                 ""}},
    {"altitude",        {"m",             "altitude",                  ""}},
    {"nyquist",         {"m s-1",         "",                          "Nyquist velocity"}},
    {"radar_latitude",  {"degrees_north", "",                          "Origin latitude of the radar"}},
    {"radar_longitude", {"degrees_east",  "",                          "Origin longitude of the radar"}},
    {"radar_altitude",  {"m",             "",                          "Origin altitude of the radar"}},
  };

  auto it = attrs.find(coord_type);
  if (it != attrs.end()) {
    auto& [units, standard_name, long_name] = it->second;
    if (!units.empty())         var.att_set("units", units);
    if (!standard_name.empty()) var.att_set("standard_name", standard_name);
    if (!long_name.empty())     var.att_set("long_name", long_name);
  }
}

// Pack a float array to int16 using scale_factor/add_offset.
// NaN values are mapped to fill_value (-32768).
// Values outside [valid_min, valid_max] are clamped.
static void pack_float_to_short(
    const float* src,
    short* dst,
    size_t n,
    const packing_info& pi)
{
  constexpr short fill_value = -32768;
  float inv_scale = 1.0f / pi.scale_factor;

  for (size_t i = 0; i < n; ++i) {
    if (std::isnan(src[i])) {
      dst[i] = fill_value;
    } else {
      // Clamp to valid range
      float val = src[i];
      if (val < pi.valid_min) val = pi.valid_min;
      if (val > pi.valid_max) val = pi.valid_max;

      // Pack: packed = (value - add_offset) / scale_factor
      float packed = (val - pi.add_offset) * inv_scale;

      // Round to nearest integer, clamp to [-32767, 32767]
      int p = static_cast<int>(std::round(packed));
      if (p < -32767) p = -32767;
      if (p > 32767) p = 32767;
      dst[i] = static_cast<short>(p);
    }
  }
}

void output_context::write_field(const std::string& field, const array2f& data, size_t layer_index)
{
  auto* var = data_vars.at(field);

  if (use_packing && packing.count(field) && packing.at(field).enabled) {
    // Pack float -> int16 and write
    auto& pi = packing.at(field);
    size_t n = data.size();
    std::vector<short> packed(n);
    pack_float_to_short(data.data(), packed.data(), n, pi);

    // Write packed data using the BOM array type
    auto packed_arr = array2<short>{data.extents()};
    for (size_t i = 0; i < n; ++i)
      packed_arr.data()[i] = packed[i];

    var->write(packed_arr, {layer_index});
  } else {
    // Write as float32
    var->write(data, {layer_index});
  }
}

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
    , const std::string& altitude_reference
    , const std::vector<std::string>& fields
    , const std::map<std::string, variable_metadata>& field_metadata
    , bool output_obs_count
    , bool pack_output
    , array1d const& radar_lat
    , array1d const& radar_lon
    , array1d const& radar_alt
    ) -> std::pair<io::nc::file, output_context>
{
  auto out_file = io::nc::file{path, io_mode::create};
  output_context ctx;
  ctx.file = &out_file;
  ctx.use_packing = pack_output;

  // CF global attributes
  out_file.att_set("Conventions", "CF-1.10");
  out_file.att_set("title", "Radar volume gridded to Cartesian coordinates");
  out_file.att_set("institution", "Bureau of Meteorology");
  auto source = meta.source.empty()
    ? std::string{"BoM radar network, gridded from radar volumes"}
    : std::string{"BoM radar network ("} + meta.source + "), gridded from radar volumes";
  out_file.att_set("source", source);
  out_file.att_set("history", "Created by vargrid " VARGRID_VERSION);
  out_file.att_set("references", "UNRAVEL processing documentation and CF conventions: https://doi.org/10.1175/JTECH-D-21-0006.1");
  out_file.att_set("comment", "Georeferencing is CF-compliant via grid_mapping variable albers_conical_equal_area.");
  out_file.att_set("date_created", to_string(bom::timestamp::now()));
  out_file.att_set("projection", proj4_string);
  out_file.att_set("gridding_method", method);
  out_file.att_set("altitude_reference", altitude_reference == "radar"
    ? "height above radar station" : "height above mean sea level");

  // Dimensions
  auto& dim_x = io::cf::create_spatial_dimension(out_file, "x", "projection_x_coordinate", coords.col_units(), coords.col_edges());
  auto& dim_y = io::cf::create_spatial_dimension(out_file, "y", "projection_y_coordinate", coords.row_units(), y_edges);
  auto& dim_e = out_file.create_dimension("elevation", meta.elevation.size());
  auto& dim_a = out_file.create_dimension("z", altitudes.size());
  auto& dim_t = out_file.create_dimension("time", 1);
  auto& dim_nrad = out_file.create_dimension("nradar", 1);
  ctx.dim_x = &dim_x;
  ctx.dim_y = &dim_y;
  ctx.dim_z = &dim_a;
  ctx.dim_e = &dim_e;
  ctx.dim_nrad = &dim_nrad;

  // Coordinate variables
  auto& var_t = out_file.create_variable("time", io::nc::data_type::f64, {&dim_t});
  auto& var_e = out_file.create_variable("elevation", io::nc::data_type::f64, {&dim_e});
  auto& var_a = out_file.create_variable("z", io::nc::data_type::f64, {&dim_a});
  auto& var_lon = out_file.create_variable("longitude", io::nc::data_type::f32, {&dim_y, &dim_x}, {dim_y.size(), dim_x.size()});
  auto& var_lat = out_file.create_variable("latitude", io::nc::data_type::f32, {&dim_y, &dim_x}, {dim_y.size(), dim_x.size()});
  auto& var_nyq = out_file.create_variable("nyquist", io::nc::data_type::f32, {&dim_e}, {dim_e.size()});
  auto& var_rlat = out_file.create_variable("radar_latitude", io::nc::data_type::f64, {&dim_nrad});
  auto& var_rlon = out_file.create_variable("radar_longitude", io::nc::data_type::f64, {&dim_nrad});
  auto& var_ralt = out_file.create_variable("radar_altitude", io::nc::data_type::f64, {&dim_nrad});
  io::cf::create_grid_mapping(out_file, "albers_conical_equal_area", proj4_string, "", "m", "m");

  if (auto* var_x = out_file.find_variable("x")) {
    var_x->att_set("axis", "X");
    var_x->att_set("long_name", "x coordinate of projection");
  }
  if (auto* var_y = out_file.find_variable("y")) {
    var_y->att_set("axis", "Y");
    var_y->att_set("long_name", "y coordinate of projection");
  }

  if (auto* var_proj = out_file.find_variable("albers_conical_equal_area")) {
    var_proj->att_set("semi_major_axis", 6378137.0);
    var_proj->att_set("inverse_flattening", 298.257222101);
    var_proj->att_set("reference_ellipsoid_name", "GRS 1980");
  }

  var_t.att_set("standard_name", "time");
  var_t.att_set("units", "seconds since 1970-01-01T00:00:00Z");
  var_t.att_set("calendar", "standard");
  var_t.att_set("axis", "T");

  set_cf_coord_attributes(var_a, "altitude");
  var_a.att_set("positive", "up");
  var_a.att_set("axis", "Z");
  var_a.att_set("long_name", altitude_reference == "radar"
    ? "height above radar station" : "height above mean sea level");
  if (altitude_reference == "radar")
    var_a.att_set("standard_name", "height");  // CF: height above a reference point
  set_cf_coord_attributes(var_lon, "longitude");
  set_cf_coord_attributes(var_lat, "latitude");
  set_cf_coord_attributes(var_nyq, "nyquist");
  set_cf_coord_attributes(var_rlat, "radar_latitude");
  set_cf_coord_attributes(var_rlon, "radar_longitude");
  set_cf_coord_attributes(var_ralt, "radar_altitude");

  auto time_value = array1d{1};
  time_value[0] = parse_cf_time_seconds(meta);
  var_t.write(time_value);
  var_e.write(meta.elevation);
  var_a.write(altitudes);
  var_lon.write(lon);
  var_lat.write(lat);
  var_nyq.write(meta.nyquist);
  var_rlat.write(radar_lat);
  var_rlon.write(radar_lon);
  var_ralt.write(radar_alt);

  // Create data variables for each field
  for (auto& field : fields) {
    auto pi = get_packing_info(field);
    ctx.packing[field] = pi;
    auto it_meta = field_metadata.find(field);
    const variable_metadata* input_meta = it_meta != field_metadata.end() ? &it_meta->second : nullptr;

    if (pack_output && pi.enabled) {
      // Packed int16 variable
      auto& var = out_file.create_variable(field, io::nc::data_type::i16,
        {&dim_a, &dim_y, &dim_x}, {1, dim_y.size(), dim_x.size()});
      set_cf_field_attributes(var, field, input_meta);

      // CF packing attributes — readers use: value = packed * scale_factor + add_offset
      var.att_set("scale_factor", static_cast<double>(pi.scale_factor));
      var.att_set("add_offset", static_cast<double>(pi.add_offset));
      var.att_set("_FillValue", static_cast<short>(-32768));
      var.att_set("valid_min", static_cast<short>(-32767));
      var.att_set("valid_max", static_cast<short>(32767));

      ctx.data_vars[field] = &var;

      trace::debug("  {} packed: range [{:.1f}, {:.1f}], scale={:.6f}, offset={:.4f}",
        field, pi.valid_min, pi.valid_max, pi.scale_factor, pi.add_offset);
    } else {
      // Float32 variable (unpacked, or unknown field)
      auto& var = out_file.create_variable(field, io::nc::data_type::f32,
        {&dim_a, &dim_y, &dim_x}, {1, dim_y.size(), dim_x.size()});
      set_cf_field_attributes(var, field, input_meta);
      var.att_set("_FillValue", nodata);

      ctx.data_vars[field] = &var;

      if (pack_output && !pi.enabled)
        trace::debug("  {} unpacked (unknown field)", field);
    }

    // Observation counts are only produced by the variational method;
    // creating the variables for other methods left never-written
    // fill-value fields in the output.
    if (output_obs_count && method == "variational") {
      auto nobs_name = "nobs_" + field;
      auto& nvar = out_file.create_variable(nobs_name, io::nc::data_type::f32,
        {&dim_a, &dim_y, &dim_x}, {1, dim_y.size(), dim_x.size()});
      nvar.att_set("long_name", "Observation count for " + field);
      nvar.att_set("units", "1");
      nvar.att_set("_FillValue", nodata);
      ctx.nobs_vars[field] = &nvar;
    }
  }

  return {std::move(out_file), ctx};
}