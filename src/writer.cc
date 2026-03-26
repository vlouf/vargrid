#include "writer.h"

auto set_cf_field_attributes(io::nc::variable& var, const std::string& quantity) -> void
{
  static const std::unordered_map<std::string, std::tuple<std::string, std::string, std::string>> cf_map = {
    {"DBZH",       {"dBZ",           "equivalent_reflectivity_factor",                          "Horizontal reflectivity"}},
    {"DBZH_CLEAN", {"dBZ",           "equivalent_reflectivity_factor",                          "Filtered horizontal reflectivity"}},
    {"DBZV",       {"dBZ",           "",                                                        "Vertical reflectivity"}},
    {"TH",         {"dBZ",           "",                                                        "Total power horizontal (uncorrected Z)"}},
    {"TV",         {"dBZ",           "",                                                        "Total power vertical (uncorrected Z)"}},
    {"VRADH",      {"m s-1",         "radial_velocity_of_scatterers_away_from_instrument",      "Mean Doppler velocity (H)"}},
    {"VRADDH",     {"m s-1",         "radial_velocity_of_scatterers_away_from_instrument",      "Dealiased Doppler velocity (H)"}},
    {"VRADV",      {"m s-1",         "",                                                        "Mean Doppler velocity (V)"}},
    {"WRADH",      {"m s-1",         "doppler_spectrum_width",                                  "Doppler spectrum width (H)"}},
    {"WRADV",      {"m s-1",         "doppler_spectrum_width",                                  "Doppler spectrum width (V)"}},
    {"ZDR",        {"dB",            "log_differential_reflectivity_hv",                        "Differential reflectivity"}},
    {"RHOHV",      {"1",             "cross_correlation_ratio_hv",                              "Cross-correlation coefficient"}},
    {"PHIDP",      {"degrees",       "differential_phase_hv",                                   "Differential phase"}},
    {"KDP",        {"degrees km-1",  "specific_differential_phase_hv",                          "Specific differential phase"}},
    {"SQI",        {"1",             "",                                                        "Signal quality index"}},
    {"SNR",        {"dB",            "signal_to_noise_ratio",                                   "Signal-to-noise ratio"}},
    {"SNRH",       {"dB",            "signal_to_noise_ratio",                                   "Signal-to-noise ratio (H)"}},
    {"CCOR",       {"dB",            "",                                                        "Clutter correction"}},
    {"CCORH",      {"dB",            "",                                                        "Clutter correction (H)"}},
  };

  auto it = cf_map.find(quantity);
  if (it != cf_map.end()) {
    auto& [units, standard_name, long_name] = it->second;
    if (!units.empty())         var.att_set("units", units);
    if (!standard_name.empty()) var.att_set("standard_name", standard_name);
    if (!long_name.empty())     var.att_set("long_name", long_name);
  } else {
    var.att_set("long_name", quantity);
  }

  var.att_set("_FillValue", nodata);
  var.att_set("grid_mapping", "proj");
  var.att_set("coordinates", "latitude longitude");
}

auto set_cf_coord_attributes(io::nc::variable& var, const std::string& coord_type) -> void
{
  static const std::unordered_map<std::string, std::tuple<std::string, std::string, std::string>> attrs = {
    {"latitude",        {"degrees_north", "latitude",                  ""}},
    {"longitude",       {"degrees_east",  "longitude",                 ""}},
    {"altitude",        {"m",             "projection_z_coordinate",   ""}},
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
    , const std::vector<std::string>& fields
    , bool output_obs_count
    , array1d const& radar_lat
    , array1d const& radar_lon
    , array1d const& radar_alt
    ) -> std::pair<io::nc::file, output_context>
{
  auto out_file = io::nc::file{path, io_mode::create};
  output_context ctx;
  ctx.file = &out_file;

  // CF global attributes
  out_file.att_set("Conventions", "CF-1.10");
  out_file.att_set("title", "Radar volume gridded to Cartesian coordinates");
  out_file.att_set("institution", "Bureau of Meteorology");
  out_file.att_set("source", meta.source);
  out_file.att_set("history", "Created by vargrid 0.1.0");
  out_file.att_set("date", meta.date);
  out_file.att_set("time", meta.time);
  out_file.att_set("lowest_sweep_time", meta.lowest_sweep_time);
  out_file.att_set("date_created", to_string(bom::timestamp::now()));
  out_file.att_set("projection", proj4_string);
  out_file.att_set("gridding_method", method);

  // Dimensions
  auto& dim_x = io::cf::create_spatial_dimension(out_file, "x", "projection_x_coordinate", coords.col_units(), coords.col_edges());
  auto& dim_y = io::cf::create_spatial_dimension(out_file, "y", "projection_y_coordinate", coords.row_units(), y_edges);
  auto& dim_e = out_file.create_dimension("elevation", meta.elevation.size());
  auto& dim_a = out_file.create_dimension("z", altitudes.size());
  auto& dim_nrad = out_file.create_dimension("nradar", 1);
  ctx.dim_x = &dim_x;
  ctx.dim_y = &dim_y;
  ctx.dim_z = &dim_a;
  ctx.dim_e = &dim_e;
  ctx.dim_nrad = &dim_nrad;

  // Coordinate variables
  auto& var_e = out_file.create_variable("elevation", io::nc::data_type::f64, {&dim_e});
  auto& var_a = out_file.create_variable("z", io::nc::data_type::f64, {&dim_a});
  auto& var_lon = out_file.create_variable("longitude", io::nc::data_type::f32, {&dim_y, &dim_x}, {dim_y.size(), dim_x.size()});
  auto& var_lat = out_file.create_variable("latitude", io::nc::data_type::f32, {&dim_y, &dim_x}, {dim_y.size(), dim_x.size()});
  auto& var_nyq = out_file.create_variable("nyquist", io::nc::data_type::f32, {&dim_e}, {dim_e.size()});
  auto& var_rlat = out_file.create_variable("radar_latitude", io::nc::data_type::f64, {&dim_nrad});
  auto& var_rlon = out_file.create_variable("radar_longitude", io::nc::data_type::f64, {&dim_nrad});
  auto& var_ralt = out_file.create_variable("radar_altitude", io::nc::data_type::f64, {&dim_nrad});
  io::cf::create_grid_mapping(out_file, "proj", proj4_string, "", "m", "m");

  // Coordinate attributes
  set_cf_coord_attributes(var_a, "altitude");
  var_a.att_set("positive", "up");
  set_cf_coord_attributes(var_lon, "longitude");
  set_cf_coord_attributes(var_lat, "latitude");
  set_cf_coord_attributes(var_nyq, "nyquist");
  set_cf_coord_attributes(var_rlat, "radar_latitude");
  set_cf_coord_attributes(var_rlon, "radar_longitude");
  set_cf_coord_attributes(var_ralt, "radar_altitude");

  // Write coordinate data
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
    auto& var = out_file.create_variable(field, io::nc::data_type::f32,
      {&dim_a, &dim_y, &dim_x}, {1, dim_y.size(), dim_x.size()});
    set_cf_field_attributes(var, field);
    ctx.data_vars[field] = &var;

    if (output_obs_count) {
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
