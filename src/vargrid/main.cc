#include <getopt.h>
#include <atomic>
#include <chrono>
#include <functional>
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

# radar moment to grid
moment DBZH

# Matrix orientation
origin xy

# gridding method: "cappi" (IDW) or "variational"
gridding_method variational

# --- CAPPI parameters (used when gridding_method is "cappi") ---

# maximum distance from CAPPI altitude to use reflectivities
max_alt_dist 20000

# exponent for inverse distance weighting
idw_pwr 2.0

# --- Variational gridding parameters (used when gridding_method is "variational") ---

# Regularisation weight: smoothness vs data fit (larger = smoother)
vargrid_alpha 1.0

# Maximum altitude difference (m) for gate selection
vargrid_max_alt_diff 2000

# Maximum CG iterations
vargrid_max_iterations 50

# Convergence tolerance (relative gradient norm)
vargrid_tolerance 1e-5

# Range scale (m) for observation weight decay
vargrid_range_scale 150000

# Use nearest-gate weighted mean as initial guess
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

auto read_single_moment(
      const std::filesystem::path& path
    , const io::configuration& config
    ) -> std::pair<volume, radarset>
{
  trace::debug("Reading volume: {}", path.string());
  io::odim::polar_volume vol_odim{path, io_mode::read_only};

  radarset dset;
  const std::string& moment = config["moment"];

  dset.elevation = get_elevation(vol_odim);
  dset.nyquist = get_nyquist(vol_odim);
  dset.lowest_sweep_time = get_lowest_sweep_time(vol_odim);

  auto vol = read_moment(vol_odim, moment, config);

  const auto& attributes = vol_odim.attributes();
  dset.source = attributes["source"].get_string();
  dset.date = attributes["date"].get_string();
  dset.time = attributes["time"].get_string();
  dset.beamwidth = attributes["beamwH"].get_real();

  trace::debug("Volume read: {} (source={}, moment={}, {} sweeps)",
    path.string(), dset.source, moment, vol.sweeps.size());

  return {vol, dset};
}

auto parse_vargrid_config(io::configuration const& config) -> vargrid_config
{
  vargrid_config cfg;
  cfg.alpha = std::stof(config.optional("vargrid_alpha", "1.0"));
  cfg.max_alt_diff = std::stof(config.optional("vargrid_max_alt_diff", "2000"));
  cfg.max_iterations = std::stoi(config.optional("vargrid_max_iterations", "50"));
  cfg.tolerance = std::stof(config.optional("vargrid_tolerance", "1e-5"));
  cfg.range_scale = std::stof(config.optional("vargrid_range_scale", "150000"));
  cfg.use_nearest_init = config.optional("vargrid_use_nearest_init", "true") == "true";
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

  // Read the volume
  volume vol;
  radarset dset;
  {
    phase_timer t("Reading volume");
    auto [v, d] = read_single_moment(input_path, config);
    vol = std::move(v);
    dset = std::move(d);
  }

  auto moment_name = (string) config["moment"];
  auto method = config.optional("gridding_method", "cappi");
  trace::log("Gridding method: {}, moment: {}", method, moment_name);

  // Precompute gate projections on the main thread (PROJ is not thread-safe).
  gate_projections gp;
  double x0 = 0, y0 = 0, dx = 1, dy = 1;
  if (method == "variational") {
    phase_timer t("Precomputing gate projections");
    gp = precompute_gate_projections(vol, proj4_string);

    // Compute grid cell centers and spacing from the grid_coordinates edges.
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

  auto vcfg = parse_vargrid_config(config);

  // Get grid coordinates for output
  auto y = array1d{coords.row_edges()};
  auto lon = array2f{latlons.extents()};
  auto lat = array2f{latlons.extents()};
  for (size_t i = 0; i < latlons.size(); i++) {
    lon.data()[i] = latlons.data()[i].lon.degrees();
    lat.data()[i] = latlons.data()[i].lat.degrees();
  }
  if (config["origin"].string().compare("xy") == 0) {
    flip(y);
    flipud(lat);
  }

  // Radar site coordinates
  auto radar_lat = array1d{1};
  auto radar_lon = array1d{1};
  auto radar_alt = array1d{1};
  radar_lat[0] = vol.location.lat.degrees();
  radar_lon[0] = vol.location.lon.degrees();
  radar_alt[0] = vol.location.alt;

  // Create output file
  trace::log("Creating output file: {}", out_path.string());
  auto out_file = io::nc::file{out_path, io_mode::create};

  auto& dim_x = io::cf::create_spatial_dimension(out_file, "x", "projection_x_coordinate", coords.col_units(), coords.col_edges());
  auto& dim_y = io::cf::create_spatial_dimension(out_file, "y", "projection_y_coordinate", coords.row_units(), y);
  auto& dim_e = out_file.create_dimension("elevation", dset.elevation.size());
  auto& dim_a = out_file.create_dimension("z", altitudes.size());
  auto& dim_nrad = out_file.create_dimension("nradar", 1);

  auto& var_e = out_file.create_variable("elevation", io::nc::data_type::f64, {&dim_e});
  auto& var_a = out_file.create_variable("z", io::nc::data_type::f64, {&dim_a});
  auto& var_lon = out_file.create_variable("longitude", io::nc::data_type::f32, {&dim_y, &dim_x}, {dim_y.size(), dim_x.size()});
  auto& var_lat = out_file.create_variable("latitude", io::nc::data_type::f32, {&dim_y, &dim_x}, {dim_y.size(), dim_x.size()});
  auto& var_nyq = out_file.create_variable("nyquist", io::nc::data_type::f32, {&dim_e}, {dim_e.size()});
  auto& var_rlat = out_file.create_variable("radar_latitude", io::nc::data_type::f64, {&dim_nrad});
  auto& var_rlon = out_file.create_variable("radar_longitude", io::nc::data_type::f64, {&dim_nrad});
  auto& var_ralt = out_file.create_variable("radar_altitude", io::nc::data_type::f64, {&dim_nrad});
  io::cf::create_grid_mapping(out_file, "proj", proj4_string, "", "m", "m");

  auto& var_data = out_file.create_variable(moment_name, io::nc::data_type::f32,
    {&dim_a, &dim_y, &dim_x}, {1, dim_y.size(), dim_x.size()});

  var_e.write(dset.elevation);
  var_a.write(altitudes);
  var_lon.write(lon);
  var_lat.write(lat);
  var_nyq.write(dset.nyquist);
  var_rlat.write(radar_lat);
  var_rlon.write(radar_lon);
  var_ralt.write(radar_alt);

  set_nc_var_attrs(var_a, "altitude");
  set_nc_var_attrs(var_lon, "longitude");
  set_nc_var_attrs(var_lat, "latitude");
  set_nc_var_attrs(var_nyq, "nyquist");
  set_nc_var_attrs(var_rlat, "radar_latitude");
  set_nc_var_attrs(var_rlon, "radar_longitude");
  set_nc_var_attrs(var_ralt, "radar_altitude");
  if (moment_name == "DBZH" || moment_name == "DBZH_CLEAN")
    set_nc_var_attrs(var_data, "reflectivity");
  var_data.att_set("_FillValue", nodata);

  out_file.att_set("source", dset.source);
  out_file.att_set("date", dset.date);
  out_file.att_set("time", dset.time);
  out_file.att_set("lowest_sweep_time", dset.lowest_sweep_time);
  out_file.att_set("date_created", to_string(bom::timestamp::now()));
  out_file.att_set("projection", proj4_string);
  out_file.att_set("gridding_method", method);

  // Process each altitude layer
  {
    phase_timer t("Gridding " + std::to_string(altitudes.size()) + " layers");

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

        array2f gridded;
        if (method == "variational") {
          auto H = build_observation_operator(
            vol, gp, x0, y0, dx, dy,
            grid_nx, grid_ny, altitudes[i], vcfg);
          gridded = variational_grid(H, vcfg);
        } else {
          gridded = generate_cappi(
            vol, latlons, config["max_alt_dist"],
            config["idw_pwr"], altitudes[i]);
        }

        if (config["origin"].string().compare("xy") == 0)
          flipud(gridded);

        auto lock = std::lock_guard<std::mutex>{mut_netcdf};
        var_data.write(gridded, {i});
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