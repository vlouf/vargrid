#include <getopt.h>
#include <atomic>
#include <chrono>
#include <functional>
#include <set>
#include "pch.h"

#include "array_operations.h"
#include "cappi.h"
#include "metadata.h"
#include "io.h"
#include "vargrid.h"

using namespace bom;

constexpr auto example_config =
R"(# vargrid configuration

# domain projection
proj4 "+proj=aea +lat_1=-32.2 +lat_2=-35.2 +lon_0=151.209 +lat_0=-33.7008 +a=6378137 +b=6356752.31414 +units=m"

# grid size
size "301 301"

# top left coordinates
left_top "-150500 150500"

# grid resolution
cell_delta "1000 -1000"

# horizontal grid units
units m

# altitude of lowest layer (m)
altitude_base 0.0

# altitude step between layers (m)
altitude_step 500.0

# number of layers
layer_count 13

# Matrix orientation
origin xy

# gridding method: "cappi" (IDW) or "variational"
gridding_method variational

# velocity field name (for undetect handling)
velocity VRADH

# Field selection (optional):
#   include_fields "DBZH VRADH ZDR"    - only grid these fields
#   exclude_fields "SQI CCOR"          - grid everything except these
# If neither is set, all fields found in the volume are gridded.
# include_fields takes precedence if both are set.

# Output observation count per field (true/false, default false)
output_obs_count false

# --- CAPPI parameters (used when gridding_method is "cappi") ---
max_alt_dist 20000
idw_pwr 2.0

# --- Variational gridding parameters ---
vargrid_alpha 1.0
vargrid_max_alt_diff 2000
vargrid_max_iterations 50
vargrid_tolerance 1e-5
vargrid_beam_power 2.0
vargrid_ref_range 10000
vargrid_min_weight 0.01
vargrid_use_nearest_init true

)";

constexpr auto try_again = "try --help for usage instructions\n";
constexpr auto usage_string =
R"(Variational gridding of radar volume to Cartesian grid

usage:
  vargrid [options] config.conf input.pvol.h5 output.nc

available options:
  -h, --help
      Show this message and exit

  -g, --generate
      Output a sample configuration file and exit

  -t, --trace=level
      Set logging level [log]
        none | status | error | warning | log | debug

Version: 0.1.0
)";

constexpr auto short_options = "hgt:";
constexpr struct option long_options[] =
{
    { "help",     no_argument,       0, 'h' }
  , { "generate", no_argument,       0, 'g' }
  , { "trace",    required_argument, 0, 't' }
  , { 0, 0, 0, 0 }
};

struct phase_timer {
  std::string name;
  std::chrono::high_resolution_clock::time_point start;

  phase_timer(std::string phase_name)
    : name(std::move(phase_name))
    , start(std::chrono::high_resolution_clock::now())
  {
    trace::log("{} ...", name);
  }

  ~phase_timer() {
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = end - start;
    trace::log("{} completed in {:.3f}s", name, dt.count());
  }
};

// Split a space-separated string into a set of tokens.
static auto split_fields(const std::string& s) -> std::set<std::string>
{
  std::set<std::string> result;
  std::istringstream iss(s);
  std::string token;
  while (iss >> token)
    result.insert(token);
  return result;
}

// Discover all available moment names in an ODIM volume.
static auto discover_fields(io::odim::polar_volume const& vol_odim) -> std::vector<std::string>
{
  std::set<std::string> seen;
  std::vector<std::string> fields;

  for (size_t iscan = 0; iscan < vol_odim.scan_count(); ++iscan) {
    auto scan_odim = vol_odim.scan_open(iscan);
    for (size_t idata = 0; idata < scan_odim.data_count(); ++idata) {
      auto data_odim = scan_odim.data_open(idata);
      auto qty = data_odim.quantity();
      if (seen.insert(qty).second)
        fields.push_back(qty);
    }
  }

  return fields;
}

// Determine which fields to grid based on include/exclude config.
static auto select_fields(
      std::vector<std::string> const& available
    , io::configuration const& config
    ) -> std::vector<std::string>
{
  std::string include_str = config.optional("include_fields", "");
  std::string exclude_str = config.optional("exclude_fields", "");

  if (!include_str.empty()) {
    auto include_set = split_fields(include_str);
    std::vector<std::string> result;
    for (auto& f : available) {
      if (include_set.count(f))
        result.push_back(f);
    }
    return result;
  }

  if (!exclude_str.empty()) {
    auto exclude_set = split_fields(exclude_str);
    std::vector<std::string> result;
    for (auto& f : available) {
      if (!exclude_set.count(f))
        result.push_back(f);
    }
    return result;
  }

  return available;
}

// Read a single moment into a volume struct.
static auto read_moment_volume(
      io::odim::polar_volume const& vol_odim
    , const std::string& moment
    , const std::string& velocity_field
    ) -> volume
{
  auto vol = volume{};
  vol.location.lat = vol_odim.latitude() * 1_deg;
  vol.location.lon = vol_odim.longitude() * 1_deg;
  vol.location.alt = vol_odim.height();

  for (size_t iscan = 0; iscan < vol_odim.scan_count(); ++iscan) {
    auto scan_odim = vol_odim.scan_open(iscan);
    if (std::fabs(scan_odim.elevation_angle() - 90) < 0.1f)
      continue;

    auto scan = sweep{};
    scan.beam = radar::beam_propagation{vol.location.alt, scan_odim.elevation_angle() * 1_deg};

    scan.bins.resize(scan_odim.bin_count());
    auto range_scale = scan_odim.range_scale();
    auto range_start = scan_odim.range_start() * 1000 + range_scale * 0.5;
    for (size_t i = 0; i < scan.bins.size(); ++i) {
      scan.bins[i].slant_range = range_start + i * range_scale;
      std::tie(scan.bins[i].ground_range, scan.bins[i].altitude) =
        scan.beam.ground_range_altitude(scan.bins[i].slant_range);
    }

    scan.rays.resize(scan_odim.ray_count());
    auto ray_scale = 360_deg / scan.rays.size();
    auto ray_start = scan_odim.ray_start() * 1_deg + ray_scale * 0.5;
    for (size_t i = 0; i < scan.rays.size(); ++i)
      scan.rays[i] = ray_start + i * ray_scale;

    for (size_t idata = 0; idata < scan_odim.data_count(); ++idata) {
      auto data_odim = scan_odim.data_open(idata);
      if (data_odim.quantity() != moment)
        continue;

      scan.data.resize(vec2z{(size_t)scan_odim.bin_count(), (size_t)scan_odim.ray_count()});
      if (moment == velocity_field)
        data_odim.read_unpack(scan.data.data(), undetect, nodata);
      else
        data_odim.read_unpack(scan.data.data(), nodata, nodata);

      vol.sweeps.push_back(std::move(scan));
      break;
    }
  }

  return vol;
}

// CF-compliant attribute mapping for ODIM quantities.
static void set_cf_attributes(io::nc::variable& var, const std::string& quantity)
{
  static const std::unordered_map<std::string, std::tuple<std::string, std::string, std::string>> cf_map = {
    {"DBZH",       {"dBZ",           "equivalent_reflectivity_factor",                          "Horizontal reflectivity"}},
    {"DBZH_CLEAN", {"dBZ",           "equivalent_reflectivity_factor",                          "Filtered horizontal reflectivity"}},
    {"DBZV",       {"dBZ",           "",                                                        "Vertical reflectivity"}},
    {"VRADH",      {"m s-1",         "radial_velocity_of_scatterers_away_from_instrument",      "Mean Doppler velocity (H)"}},
    {"VRADDH",     {"m s-1",         "radial_velocity_of_scatterers_away_from_instrument",      "Dealiased Doppler velocity (H)"}},
    {"VRADV",      {"m s-1",         "",                                                        "Mean Doppler velocity (V)"}},
    {"WRADH",      {"m s-1",         "doppler_spectrum_width",                                  "Doppler spectrum width (H)"}},
    {"ZDR",        {"dB",            "log_differential_reflectivity_hv",                        "Differential reflectivity"}},
    {"RHOHV",      {"1",             "cross_correlation_ratio_hv",                              "Cross-correlation coefficient"}},
    {"PHIDP",      {"degrees",       "differential_phase_hv",                                   "Differential phase"}},
    {"KDP",        {"degrees km-1",  "specific_differential_phase_hv",                          "Specific differential phase"}},
    {"SQI",        {"1",             "",                                                        "Signal quality index"}},
    {"SNR",        {"dB",            "signal_to_noise_ratio",                                   "Signal-to-noise ratio"}},
    {"SNRH",       {"dB",            "signal_to_noise_ratio",                                   "Signal-to-noise ratio (H)"}},
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

auto parse_vargrid_config(io::configuration const& config, float beamwidth) -> vargrid_config
{
  vargrid_config cfg;
  cfg.alpha = std::stof(config.optional("vargrid_alpha", "1.0"));
  cfg.max_alt_diff = std::stof(config.optional("vargrid_max_alt_diff", "2000"));
  cfg.max_iterations = std::stoi(config.optional("vargrid_max_iterations", "50"));
  cfg.tolerance = std::stof(config.optional("vargrid_tolerance", "1e-5"));
  cfg.beam_power = std::stof(config.optional("vargrid_beam_power", "2.0"));
  cfg.ref_range = std::stof(config.optional("vargrid_ref_range", "10000"));
  cfg.min_weight = std::stof(config.optional("vargrid_min_weight", "0.01"));
  cfg.use_nearest_init = config.optional("vargrid_use_nearest_init", "true") == "true";
  cfg.beamwidth = beamwidth;
  return cfg;
}

auto run_gridding(
      io::configuration const& config
    , std::filesystem::path const& input_path
    , std::filesystem::path const& out_path
    ) -> void
{
  trace::log("Input:  {}", input_path.string());
  trace::log("Output: {}", out_path.string());

  auto proj4_string = config["proj4"].string();
  auto proj = map_projection{map_projection::default_context, proj4_string};
  auto coords = grid_coordinates{
      config["size"]
    , config["left_top"]
    , config["cell_delta"]
    , config["units"]
    , config["units"]
    };
  auto latlons = determine_geodetic_coordinates(proj, coords);
  auto altitudes = init_altitudes(config);
  auto grid_nx = latlons.extents().x;
  auto grid_ny = latlons.extents().y;
  trace::log("Grid: {}x{}, {} altitude layers ({:.0f}m to {:.0f}m)",
    grid_nx, grid_ny, altitudes.size(),
    altitudes[0], altitudes[altitudes.size()-1]);

  // Open ODIM volume and discover fields
  io::odim::polar_volume vol_odim{input_path, io_mode::read_only};
  auto available_fields = discover_fields(vol_odim);
  auto selected_fields = select_fields(available_fields, config);

  trace::log("Available fields: {}", [&]{
    std::string s;
    for (auto& f : available_fields) { if (!s.empty()) s += " "; s += f; }
    return s;
  }());
  trace::log("Selected fields:  {}", [&]{
    std::string s;
    for (auto& f : selected_fields) { if (!s.empty()) s += " "; s += f; }
    return s;
  }());

  if (selected_fields.empty()) {
    trace::error("No fields selected for gridding.");
    return;
  }

  // Read metadata
  auto elevation = get_elevation(vol_odim);
  auto nyquist = get_nyquist(vol_odim);
  auto lowest_sweep_time = get_lowest_sweep_time(vol_odim);
  const auto& attributes = vol_odim.attributes();
  auto source = attributes["source"].get_string();
  auto date_str = attributes["date"].get_string();
  auto time_str = attributes["time"].get_string();
  auto beamwidth = attributes["beamwH"].get_real();

  auto velocity_field = config.optional("velocity", "VRADH");

  // Read all selected moments
  std::map<std::string, volume> volumes;
  {
    phase_timer t("Reading " + std::to_string(selected_fields.size()) + " moment(s)");
    for (auto& field : selected_fields) {
      volumes[field] = read_moment_volume(vol_odim, field, velocity_field);
      trace::debug("  {} : {} sweeps", field, volumes[field].sweeps.size());
    }
  }

  auto method = config.optional("gridding_method", "cappi");
  trace::log("Gridding method: {}", method);

  auto vcfg = parse_vargrid_config(config, beamwidth);
  auto output_obs_count = config.optional("output_obs_count", "false") == "true";

  // Precompute gate projections on main thread
  gate_projections gp;
  double x0 = 0, y0 = 0, dx = 1, dy = 1;
  if (method == "variational") {
    phase_timer t("Precomputing gate projections");
    gp = precompute_gate_projections(volumes.begin()->second, proj4_string);

    auto col_edges = coords.col_edges();
    auto row_edges = coords.row_edges();
    x0 = 0.5 * (col_edges[0] + col_edges[1]);
    y0 = 0.5 * (row_edges[0] + row_edges[1]);
    double x_last = 0.5 * (col_edges[grid_nx - 1] + col_edges[grid_nx]);
    double y_last = 0.5 * (row_edges[grid_ny - 1] + row_edges[grid_ny]);
    dx = (x_last - x0) / (static_cast<double>(grid_nx) - 1.0);
    dy = (y_last - y0) / (static_cast<double>(grid_ny) - 1.0);
    trace::debug("Grid origin: ({:.1f}, {:.1f}), spacing: ({:.1f}, {:.1f})", x0, y0, dx, dy);
  }

  // Prepare output coordinates
  auto y = array1d{coords.row_edges()};
  auto lon = array2f{latlons.extents()};
  auto lat = array2f{latlons.extents()};
  for (size_t i = 0; i < latlons.size(); i++) {
    lon.data()[i] = latlons.data()[i].lon.degrees();
    lat.data()[i] = latlons.data()[i].lat.degrees();
  }
  if (config["origin"].string() == "xy") {
    flip(y);
    flipud(lat);
  }

  auto radar_lat = array1d{1};
  auto radar_lon = array1d{1};
  auto radar_alt = array1d{1};
  auto& first_vol = volumes.begin()->second;
  radar_lat[0] = first_vol.location.lat.degrees();
  radar_lon[0] = first_vol.location.lon.degrees();
  radar_alt[0] = first_vol.location.alt;

  // --- Create CF-compliant output NetCDF ---
  trace::log("Creating output file: {}", out_path.string());
  auto out_file = io::nc::file{out_path, io_mode::create};

  // CF global attributes
  out_file.att_set("Conventions", "CF-1.10");
  out_file.att_set("title", "Radar volume gridded to Cartesian coordinates");
  out_file.att_set("institution", "Bureau of Meteorology");
  out_file.att_set("source", source);
  out_file.att_set("history", "Created by vargrid 0.1.0");
  out_file.att_set("date", date_str);
  out_file.att_set("time", time_str);
  out_file.att_set("lowest_sweep_time", lowest_sweep_time);
  out_file.att_set("date_created", to_string(bom::timestamp::now()));
  out_file.att_set("projection", proj4_string);
  out_file.att_set("gridding_method", method);

  // Dimensions
  auto& dim_x = io::cf::create_spatial_dimension(out_file, "x", "projection_x_coordinate", coords.col_units(), coords.col_edges());
  auto& dim_y = io::cf::create_spatial_dimension(out_file, "y", "projection_y_coordinate", coords.row_units(), y);
  auto& dim_e = out_file.create_dimension("elevation", elevation.size());
  auto& dim_a = out_file.create_dimension("z", altitudes.size());
  auto& dim_nrad = out_file.create_dimension("nradar", 1);

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

  set_nc_var_attrs(var_a, "altitude");
  var_a.att_set("positive", "up");
  set_nc_var_attrs(var_lon, "longitude");
  set_nc_var_attrs(var_lat, "latitude");
  set_nc_var_attrs(var_nyq, "nyquist");
  set_nc_var_attrs(var_rlat, "radar_latitude");
  set_nc_var_attrs(var_rlon, "radar_longitude");
  set_nc_var_attrs(var_ralt, "radar_altitude");

  var_e.write(elevation);
  var_a.write(altitudes);
  var_lon.write(lon);
  var_lat.write(lat);
  var_nyq.write(nyquist);
  var_rlat.write(radar_lat);
  var_rlon.write(radar_lon);
  var_ralt.write(radar_alt);

  // Create data variables for each field
  std::map<std::string, io::nc::variable*> data_vars;
  std::map<std::string, io::nc::variable*> nobs_vars;
  for (auto& field : selected_fields) {
    auto& var = out_file.create_variable(field, io::nc::data_type::f32,
      {&dim_a, &dim_y, &dim_x}, {1, dim_y.size(), dim_x.size()});
    set_cf_attributes(var, field);
    data_vars[field] = &var;

    if (output_obs_count) {
      auto nobs_name = "nobs_" + field;
      auto& nvar = out_file.create_variable(nobs_name, io::nc::data_type::f32,
        {&dim_a, &dim_y, &dim_x}, {1, dim_y.size(), dim_x.size()});
      nvar.att_set("long_name", "Observation count for " + field);
      nvar.att_set("units", "1");
      nvar.att_set("_FillValue", nodata);
      nobs_vars[field] = &nvar;
    }
  }

  // --- Grid all fields at all layers ---
  {
    phase_timer t("Gridding " + std::to_string(selected_fields.size()) +
                  " field(s) x " + std::to_string(altitudes.size()) + " layers");

    auto num_workers = std::min(
        static_cast<size_t>(std::thread::hardware_concurrency()),
        altitudes.size()
    );
    if (num_workers == 0) num_workers = 4;

    std::atomic<size_t> next_layer{0};
    std::mutex mut_netcdf;

    auto worker = [&]() {
      while (true) {
        auto i = next_layer.fetch_add(1);
        if (i >= altitudes.size())
          break;

        for (auto& field : selected_fields) {
          auto& vol = volumes.at(field);
          array2f gridded;

          if (method == "variational") {
            auto H = build_observation_operator(
              vol, gp, x0, y0, dx, dy,
              grid_nx, grid_ny, altitudes[i], vcfg);
            gridded = variational_grid(H, vcfg);

            if (output_obs_count) {
              auto nobs = observation_density(H);
              if (config["origin"].string() == "xy")
                flipud(nobs);
              auto lock = std::lock_guard<std::mutex>{mut_netcdf};
              nobs_vars.at(field)->write(nobs, {i});
            }
          } else {
            gridded = generate_cappi(
              vol, latlons, config["max_alt_dist"],
              config["idw_pwr"], altitudes[i]);
          }

          if (config["origin"].string() == "xy")
            flipud(gridded);

          auto lock = std::lock_guard<std::mutex>{mut_netcdf};
          data_vars.at(field)->write(gridded, {i});
        }
      }
    };

    auto threads = vector<std::thread>{};
    for (size_t w = 0; w < num_workers; ++w)
      threads.emplace_back(worker);
    for (auto& t : threads)
      t.join();
  }
}

int main(int argc, char* argv[])
{
  try
  {
    while (true)
    {
      int option_index = 0;
      int c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1)
        break;
      switch (c) {
      case 'h':
        std::cout << usage_string;
        return EXIT_SUCCESS;
      case 'g':
        std::cout << example_config;
        return EXIT_SUCCESS;
      case 't':
        trace::set_min_level(from_string<trace::level>(optarg));
        break;
      case '?':
        std::cerr << try_again;
        return EXIT_FAILURE;
      }
    }

    std::vector<std::string> positional_args;
    for (int i = optind; i < argc; ++i)
      positional_args.emplace_back(argv[i]);

    if (positional_args.size() != 3) {
      std::cerr << "Expected 3 positional arguments (config input.h5 output.nc), got "
                << positional_args.size() << "\n" << try_again;
      return EXIT_FAILURE;
    }

    auto config = io::configuration{std::ifstream{positional_args[0]}};

    auto origin = config["origin"].string();
    if (origin != "ij" && origin != "xy") {
      trace::error("Invalid origin '{}'. Must be 'ij' or 'xy'.", origin);
      return EXIT_FAILURE;
    }

    auto start = std::chrono::high_resolution_clock::now();
    run_gridding(config, positional_args[1], positional_args[2]);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt = end - start;
    trace::log("Completed in {:.3f}s: {}", dt.count(), positional_args[1]);
  }
  catch (std::exception& err)
  {
    trace::error("fatal exception: {}", format_exception(err));
    return EXIT_FAILURE;
  }
  catch (...)
  {
    trace::error("fatal exception: (unknown exception)");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}