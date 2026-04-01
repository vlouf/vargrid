#include <getopt.h>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <functional>
#include <mutex>

#include "types.h"
#include "util.h"
#include "reader.h"
#include "writer.h"
#include "grid.h"
#include "gridding/variational.h"
#include "gridding/cappi.h"
#include "post/post_processor.h"

using namespace bom;

constexpr auto try_again = "try --help for usage instructions\n";
constexpr auto usage_string =
R"(Variational gridding of radar volume to Cartesian grid

usage:
  vargrid [options] config.conf input.pvol.h5 output.nc

available options:
  -h, --help
      Show this message and exit

  -g, --generate [input.pvol.h5]
      Output a sample configuration file and exit.
      If an ODIM file is provided, the projection center (lat_0, lon_0)
      is set to the radar position read from the file.

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

// Generate a configuration file, optionally reading radar position from an ODIM file.
static void generate_config(const char* odim_path)
{
  double lat0 = -33.7008;
  double lon0 = 151.209;
  double lat1 = lat0 + 2.0;
  double lat2 = lat0 - 2.0;
  std::string source_comment = "# (default radar position — pass an ODIM file to vargrid -g to auto-populate)";

  if (odim_path) {
    try {
      io::odim::polar_volume vol{std::filesystem::path(odim_path), io_mode::read_only};
      lat0 = vol.latitude();
      lon0 = vol.longitude();
      lat1 = lat0 + 5.0;  // reasonable standard parallels for AEA
      lat2 = lat0 - 5.0;

      const auto& attrs = vol.attributes();
      std::string src = attrs["source"].get_string();
      source_comment = "# Radar: " + src;
    } catch (std::exception& err) {
      std::cerr << "Warning: could not read ODIM file: " << err.what() << "\n"
                << "Using default radar position.\n";
    }
  }

  // Use a fixed-precision format for the projection values
  char proj4_line[512];
  std::snprintf(proj4_line, sizeof(proj4_line),
    "proj4 \"+proj=aea +lat_1=%.1f +lat_2=%.1f +lon_0=%.5f +lat_0=%.5f +units=m +ellps=GRS80\"",
    lat1, lat2, lon0, lat0);

  std::cout << "# vargrid configuration\n"
            << "# =====================\n"
            << source_comment << "\n"
            << "\n"
            << "# --- Grid geometry (required) ---\n"
            << "\n"
            << "# Map projection (PROJ4 string)\n"
            << proj4_line << "\n"
            << "\n"
            << "# Grid dimensions (nx ny)\n"
            << "size \"301 301\"\n"
            << "\n"
            << "# Top-left corner in projected coordinates\n"
            << "left_top \"-150500 150500\"\n"
            << "\n"
            << "# Cell size in projected coordinates (dx dy, negative dy = y increases downward)\n"
            << "cell_delta \"1000 -1000\"\n"
            << "\n"
            << "# Coordinate units\n"
            << "units m\n"
            << "\n"
            << "# Matrix orientation: \"xy\" (geographic, y-axis flipped) or \"ij\" (matrix)\n"
            << "origin xy\n"
            << "\n"
            << "# --- Altitude layers ---\n"
            << "\n"
            << "altitude_base 0.0\n"
            << "altitude_step 500.0\n"
            << "layer_count 13\n"
            << "\n"
            << "# --- Gridding method ---\n"
            << "\n"
            << "# \"variational\" (Brook et al. 2022 inspired) or \"cappi\" (IDW baseline)\n"
            << "gridding_method variational\n"
            << "\n"
            << "# --- Field selection (optional) ---\n"
            << "# include_fields \"DBZH VRADH ZDR\"\n"
            << "# exclude_fields \"SQI CCOR\"\n"
            << "\n"
            << "velocity VRADH\n"
            << "\n"
            << "# --- Output options ---\n"
            << "\n"
            << "pack_output true\n"
            << "output_obs_count false\n"
            << "\n"
            << "# --- CAPPI parameters ---\n"
            << "\n"
            << "max_alt_dist 20000\n"
            << "idw_pwr 2.0\n"
            << "\n"
            << "# --- Variational parameters ---\n"
            << "\n"
            << "vargrid_lambda_h 0.01\n"
            << "vargrid_max_alt_diff 2000\n"
            << "vargrid_max_iterations 200\n"
            << "vargrid_tolerance 1e-5\n"
            << "vargrid_beam_power 2.0\n"
            << "vargrid_ref_range 10000\n"
            << "vargrid_min_weight 0.01\n"
            << "vargrid_use_nearest_init true\n"
            << "vargrid_kappa 0\n"
            << "vargrid_range_spacing 250\n"
            << "vargrid_mask_distance 3\n";
}
static auto select_fields(
      std::vector<std::string> const& available
    , io::configuration const& config
    ) -> std::vector<std::string>
{
  std::string include_str = config.optional("include_fields", "");
  std::string exclude_str = config.optional("exclude_fields", "");

  std::set<std::string> available_set(available.begin(), available.end());

  if (!include_str.empty()) {
    auto include_set = split_fields(include_str);
    for (auto& f : include_set)
      if (!available_set.contains(f))
        trace::warning("include_fields: '{}' not found in input volume", f);
    std::vector<std::string> result;
    for (auto& f : available)
      if (include_set.contains(f)) result.push_back(f);
    return result;
  }

  if (!exclude_str.empty()) {
    auto exclude_set = split_fields(exclude_str);
    for (auto& f : exclude_set)
      if (!available_set.contains(f))
        trace::warning("exclude_fields: '{}' not found in input volume", f);
    std::vector<std::string> result;
    for (auto& f : available)
      if (!exclude_set.contains(f)) result.push_back(f);
    return result;
  }

  return available;
}

static auto parse_vargrid_config(io::configuration const& config, float beamwidth) -> vargrid_config
{
  vargrid_config cfg;
  cfg.lambda_h       = std::stof(config.optional("vargrid_lambda_h", "0.01"));
  cfg.max_alt_diff   = std::stof(config.optional("vargrid_max_alt_diff", "2000"));
  cfg.max_iterations = std::stoi(config.optional("vargrid_max_iterations", "200"));
  cfg.tolerance      = std::stof(config.optional("vargrid_tolerance", "1e-5"));
  cfg.beam_power     = std::stof(config.optional("vargrid_beam_power", "2.0"));
  cfg.ref_range      = std::stof(config.optional("vargrid_ref_range", "10000"));
  cfg.min_weight     = std::stof(config.optional("vargrid_min_weight", "0.01"));
  cfg.use_nearest_init = std::string(config.optional("vargrid_use_nearest_init", "true")) != "false";
  cfg.kappa          = std::stof(config.optional("vargrid_kappa", "0"));
  cfg.range_spacing  = std::stof(config.optional("vargrid_range_spacing", "250"));
  cfg.mask_distance_cells = std::stof(config.optional("vargrid_mask_distance", "3"));
  cfg.beamwidth      = beamwidth;
  return cfg;
}

static auto run_gridding(
      io::configuration const& config
    , std::filesystem::path const& input_path
    , std::filesystem::path const& out_path
    ) -> void
{
  trace::log("Input:  {}", input_path.string());
  trace::log("Output: {}", out_path.string());

  // --- Set up grid ---
  auto proj4_string = config["proj4"].string();
  auto proj = map_projection{map_projection::default_context, proj4_string};
  auto coords = grid_coordinates{
      config["size"], config["left_top"], config["cell_delta"],
      config["units"], config["units"]};
  auto latlons = determine_geodetic_coordinates(proj, coords);
  auto altitudes = init_altitudes(config);
  auto grid_nx = latlons.extents().x;
  auto grid_ny = latlons.extents().y;
  trace::log("Grid: {}x{}, {} altitude layers ({:.0f}m to {:.0f}m)",
    grid_nx, grid_ny, altitudes.size(), altitudes[0], altitudes[altitudes.size()-1]);

  // --- Discover and select fields ---
  io::odim::polar_volume vol_odim{input_path, io_mode::read_only};
  auto available_fields = discover_fields(vol_odim);
  auto selected_fields = select_fields(available_fields, config);

  trace::log("Available fields: {}", [&]{
    std::string s; for (auto& f : available_fields) { if (!s.empty()) s += " "; s += f; } return s;
  }());
  trace::log("Selected fields:  {}", [&]{
    std::string s; for (auto& f : selected_fields) { if (!s.empty()) s += " "; s += f; } return s;
  }());

  if (selected_fields.empty()) {
    trace::error("No fields selected for gridding.");
    return;
  }

  // --- Read metadata and volumes ---
  auto meta = read_metadata(vol_odim);
  std::string velocity_field = config.optional("velocity", "VRADH");

  std::map<std::string, volume> volumes;
  {
    phase_timer t("Reading " + std::to_string(selected_fields.size()) + " moment(s)");
    for (auto& field : selected_fields) {
      volumes[field] = read_moment_volume(vol_odim, field, velocity_field);
      trace::debug("  {} : {} sweeps", field, volumes[field].sweeps.size());
    }
  }

  std::string method = config.optional("gridding_method", "cappi");
  trace::log("Gridding method: {}", method);

  auto vcfg = parse_vargrid_config(config, meta.beamwidth);
  auto output_obs_count = std::string(config.optional("output_obs_count", "false")) == "false" ? false : true;
  auto pack_output = std::string(config.optional("pack_output", "true")) == "false" ? false : true;

  // --- Precompute grid bearings (variational only, main thread) ---
  grid_bearings gb;
  float grid_spacing = 1000.0f;  // default, will be computed from config
  if (method == "variational") {
    phase_timer t("Precomputing grid bearings");
    gb = precompute_grid_bearings(volumes.begin()->second.location, latlons);

    // Compute grid spacing from cell_delta config
    auto col_edges = coords.col_edges();
    grid_spacing = static_cast<float>(std::fabs(col_edges[1] - col_edges[0]));
    trace::debug("Grid spacing: {:.0f}m, grid bearings computed for {}x{} cells",
      grid_spacing, grid_nx, grid_ny);
  }

  // --- Prepare output coordinates ---
  auto y_edges = array1d{coords.row_edges()};
  auto lon = array2f{latlons.extents()};
  auto lat = array2f{latlons.extents()};
  for (size_t i = 0; i < latlons.size(); i++) {
    lon.data()[i] = latlons.data()[i].lon.degrees();
    lat.data()[i] = latlons.data()[i].lat.degrees();
  }
  if (config["origin"].string() == "xy") {
    flip(y_edges);
    flipud(lat);
  }

  auto radar_lat = array1d{1};
  auto radar_lon = array1d{1};
  auto radar_alt = array1d{1};
  auto& first_vol = volumes.begin()->second;
  radar_lat[0] = first_vol.location.lat.degrees();
  radar_lon[0] = first_vol.location.lon.degrees();
  radar_alt[0] = first_vol.location.alt;

  // --- Create output file ---
  trace::log("Creating output file: {}", out_path.string());
  auto [out_file, ctx] = create_output_file(
    out_path, coords, y_edges, lon, lat, altitudes, meta,
    proj4_string, method, selected_fields, output_obs_count, pack_output,
    radar_lat, radar_lon, radar_alt);

  // --- Grid all fields at all layers ---
  {
    phase_timer t("Gridding " + std::to_string(selected_fields.size()) +
                  " field(s) x " + std::to_string(altitudes.size()) + " layers");

    auto num_workers = std::min(
        static_cast<size_t>(std::thread::hardware_concurrency()),
        altitudes.size());
    if (num_workers == 0) num_workers = 4;

    std::atomic<size_t> next_layer{0};
    std::atomic<size_t> completed_tasks{0};
    size_t total_tasks = altitudes.size() * selected_fields.size();
    std::mutex mut_netcdf;

    auto worker = [&]() {
      while (true) {
        auto i = next_layer.fetch_add(1);
        if (i >= altitudes.size()) break;

        for (auto& field : selected_fields) {
          auto& vol = volumes.at(field);
          array2f gridded;

          if (method == "variational") {
            auto H = build_observation_operator(
              vol, gb, grid_nx, grid_ny, altitudes[i], grid_spacing, vcfg);
            gridded = variational_grid(H, vcfg, total_tasks, completed_tasks);

            if (output_obs_count) {
              auto nobs = observation_density(H);
              if (config["origin"].string() == "xy") flipud(nobs);
              auto lock = std::lock_guard<std::mutex>{mut_netcdf};
              ctx.nobs_vars.at(field)->write(nobs, {i});
            }
          } else {
            gridded = generate_cappi(
              vol, latlons, config["max_alt_dist"], config["idw_pwr"], altitudes[i]);
          }

          if (config["origin"].string() == "xy") flipud(gridded);

          auto lock = std::lock_guard<std::mutex>{mut_netcdf};
          ctx.write_field(field, gridded, i);
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
    bool do_generate = false;

    while (true) {
      int option_index = 0;
      int c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1) break;
      switch (c) {
      case 'h': std::cout << usage_string; return EXIT_SUCCESS;
      case 'g': do_generate = true; break;
      case 't': trace::set_min_level(from_string<trace::level>(optarg)); break;
      case '?': std::cerr << try_again; return EXIT_FAILURE;
      }
    }

    if (do_generate) {
      // Check if an ODIM file was passed as a positional argument
      const char* odim_path = (optind < argc) ? argv[optind] : nullptr;
      generate_config(odim_path);
      return EXIT_SUCCESS;
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
  catch (std::exception& err) {
    trace::error("fatal exception: {}", format_exception(err));
    return EXIT_FAILURE;
  }
  catch (...) {
    trace::error("fatal exception: (unknown exception)");
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}