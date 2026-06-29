#include "reader.h"

#include <limits>
#include <optional>

namespace {

struct cf_radial_context {
  io::nc::file file;
  std::vector<std::string> field_names;
};

auto maybe_read_double_attribute(io::nc::object const& obj, std::string const& name) -> std::optional<double>
{
  if (!obj.att_exists(name)) return std::nullopt;
  double value = 0.0;
  obj.att_get(name, value);
  return value;
}

auto maybe_read_float_attribute(io::nc::object const& obj, std::string const& name) -> std::optional<float>
{
  if (!obj.att_exists(name)) return std::nullopt;
  float value = 0.0f;
  obj.att_get(name, value);
  return value;
}

static auto make_cf_field_aliases(std::string const& field) -> std::vector<std::string>
{
  std::vector<std::string> aliases;
  if (!field.empty()) aliases.push_back(field);

  auto lower = field;
  std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });

  static const std::unordered_map<std::string, std::vector<std::string>> alias_map = {
    {"corrected_velocity", {"VRADH", "VRAD", "velocity"}},
    {"velocity", {"VRADH", "VRAD", "corrected_velocity"}},
    {"corrected_reflectivity", {"DBZH", "reflectivity"}},
    {"reflectivity", {"DBZH", "corrected_reflectivity"}},
    {"corrected_differential_reflectivity", {"ZDR", "differential_reflectivity"}},
    {"differential_reflectivity", {"ZDR", "corrected_differential_reflectivity"}},
    {"spectrum_width", {"WRADH", "WRAD"}},
    {"signal_to_noise_ratio", {"SNR"}},
    {"cross_correlation_ratio", {"RHOHV"}},
    {"corrected_specific_differential_phase", {"KDP", "differential_phase"}},
    {"differential_phase", {"PHIDP", "corrected_differential_phase"}},
    {"path_integrated_attenuation", {"PIA"}},
    {"normalized_intercept_parameter", {"NW"}},
    {"median_volume_diameter", {"D0"}},
    {"radar_estimated_rain_rate", {"RAIN"}},
    {"radar_estimated_snow_rate", {"SNOW"}},
    {"temperature", {"TEMPERATURE"}},
  };

  if (auto it = alias_map.find(lower); it != alias_map.end()) {
    for (auto const& alias : it->second) aliases.push_back(alias);
  }

  std::sort(aliases.begin(), aliases.end());
  aliases.erase(std::unique(aliases.begin(), aliases.end()), aliases.end());
  return aliases;
}

static auto read_cf_radial_variable_value(io::nc::variable const& var, std::string const& name) -> double
{
  double value = 0.0;
  if (var.att_exists(name)) {
    var.att_get(name, value);
  }
  return value;
}

static auto get_var_shape(io::nc::variable const& var) -> std::vector<size_t>
{
  std::vector<size_t> shape;
  for (auto const* dim : var.dimensions())
    shape.push_back(dim->size());
  return shape;
}

static auto normalize_scan_rays(sweep& scan) -> size_t
{
  if (scan.rays.empty())
    return 0;

  size_t min_idx = 0;
  scan.rays[0] = scan.rays[0].normalize();
  auto min_angle = scan.rays[0];

  for (size_t i = 1; i < scan.rays.size(); ++i) {
    scan.rays[i] = scan.rays[i].normalize();
    if (scan.rays[i] < min_angle) {
      min_angle = scan.rays[i];
      min_idx = i;
    }
  }

  return min_idx;
}

static void rotate_scan_data_rows(sweep& scan, size_t shift)
{
  if (shift == 0 || scan.rays.empty())
    return;

  const size_t rows = scan.rays.size();
  const size_t cols = scan.bins.size();
  std::vector<float> copy(scan.data.size());
  std::copy_n(scan.data.data(), scan.data.size(), copy.data());

  for (size_t i = 0; i < rows; ++i) {
    const size_t src = (shift + i) % rows;
    std::copy_n(copy.data() + src * cols, cols, scan.data.data() + i * cols);
  }

  std::rotate(scan.rays.begin(), scan.rays.begin() + shift, scan.rays.end());
}

static auto open_cf_radial_context(std::filesystem::path const& input_path) -> cf_radial_context
{
  return cf_radial_context{io::nc::file{input_path.string(), io_mode::read_only}, {}};
}

static auto discover_cf_radial_fields(io::nc::group const& group) -> std::vector<std::string>
{
  std::vector<std::string> fields;
  std::set<std::string> seen;

  for (auto const& var : group.variables()) {
    auto const& name = var.name();

    // Skip coordinate and auxiliary variables.
    if (name == "time" || name == "range" || name == "azimuth" || name == "elevation")
      continue;
    if (name == "latitude" || name == "longitude" || name == "altitude")
      continue;
    if (name == "fixed_angle" || name == "sweep_number" || name == "sweep_start_ray_index" || name == "sweep_end_ray_index" || name == "sweep_mode")
      continue;
    if (name == "time_coverage_start" || name == "time_coverage_end" || name == "time_reference" || name == "volume_number" || name == "instrument_type")
      continue;

    // Only keep variables that are 2-D (time, range) and look like measurements.
    auto const shape = get_var_shape(var);
    if (shape.size() != 2 || shape[0] == 0 || shape[1] == 0)
      continue;

    if (seen.insert(name).second)
      fields.push_back(name);
  }

  return fields;
}

static auto read_cf_radial_sweeps(io::nc::group const& group, volume_metadata& meta, volume& vol) -> void
{
  auto const* var_range = group.find_variable("range");
  auto const* var_azimuth = group.find_variable("azimuth");
  auto const* var_elevation = group.find_variable("elevation");
  auto const* var_fixed_angle = group.find_variable("fixed_angle");
  auto const* var_sweep_start = group.find_variable("sweep_start_ray_index");
  auto const* var_sweep_end = group.find_variable("sweep_end_ray_index");
  auto const* var_time = group.find_variable("time");

  if (!var_range || !var_azimuth || !var_elevation || !var_fixed_angle || !var_sweep_start || !var_sweep_end || !var_time) {
    throw std::runtime_error("CF/Radial file is missing required sweep coordinates");
  }

  std::vector<float> range_vals;
  {
    auto const shape = get_var_shape(*var_range);
    range_vals.resize(shape.empty() ? 0 : shape.front());
    var_range->read(range_vals);
  }

  std::vector<float> azimuth_vals;
  {
    auto const shape = get_var_shape(*var_azimuth);
    azimuth_vals.resize(shape.empty() ? 0 : shape.front());
    var_azimuth->read(azimuth_vals);
  }

  std::vector<double> elevation_vals;
  {
    auto const shape = get_var_shape(*var_elevation);
    elevation_vals.resize(shape.empty() ? 0 : shape.front());
    var_elevation->read(elevation_vals);
  }

  std::vector<float> fixed_angle_vals;
  {
    auto const shape = get_var_shape(*var_fixed_angle);
    fixed_angle_vals.resize(shape.empty() ? 0 : shape.front());
    var_fixed_angle->read(fixed_angle_vals);
  }

  std::vector<int> sweep_start_vals;
  {
    auto const shape = get_var_shape(*var_sweep_start);
    sweep_start_vals.resize(shape.empty() ? 0 : shape.front());
    var_sweep_start->read(sweep_start_vals);
  }

  std::vector<int> sweep_end_vals;
  {
    auto const shape = get_var_shape(*var_sweep_end);
    sweep_end_vals.resize(shape.empty() ? 0 : shape.front());
    var_sweep_end->read(sweep_end_vals);
  }

  std::vector<double> time_vals;
  {
    auto const shape = get_var_shape(*var_time);
    time_vals.resize(shape.empty() ? 0 : shape.front());
    var_time->read(time_vals);
  }

  if (azimuth_vals.size() != elevation_vals.size() || azimuth_vals.size() != time_vals.size())
    throw std::runtime_error("CF/Radial ray coordinate variables are inconsistent");

  if (sweep_start_vals.size() != sweep_end_vals.size() || sweep_start_vals.size() != fixed_angle_vals.size())
    throw std::runtime_error("CF/Radial sweep metadata is inconsistent");

  var_range->read(range_vals);
  var_azimuth->read(azimuth_vals);
  var_elevation->read(elevation_vals);
  var_fixed_angle->read(fixed_angle_vals);
  var_sweep_start->read(sweep_start_vals);
  var_sweep_end->read(sweep_end_vals);
  var_time->read(time_vals);

  const size_t nsweeps = fixed_angle_vals.size();
  meta.elevation = array1f{nsweeps};
  for (size_t i = 0; i < nsweeps; ++i)
    meta.elevation[i] = fixed_angle_vals[i];

  meta.nyquist = array1f{nsweeps};
  for (size_t i = 0; i < nsweeps; ++i)
    meta.nyquist[i] = -9999.f;

  auto const* var_lat = group.find_variable("latitude");
  auto const* var_lon = group.find_variable("longitude");
  auto const* var_alt = group.find_variable("altitude");
  if (var_lat && var_lon && var_alt) {
    std::vector<double> lat(1), lon(1), alt(1);
    var_lat->read(lat);
    var_lon->read(lon);
    var_alt->read(alt);
    vol.location.lat = lat.empty() ? 0.0 * 1_deg : lat.front() * 1_deg;
    vol.location.lon = lon.empty() ? 0.0 * 1_deg : lon.front() * 1_deg;
    vol.location.alt = alt.empty() ? 0.0f : static_cast<float>(alt.front());
    meta.location.lat = vol.location.lat;
    meta.location.lon = vol.location.lon;
    meta.location.alt = vol.location.alt;
  }

  auto const* var_src = group.find_variable("source");
  if (var_src) {
    std::string src;
    var_src->att_get("long_name", src);
    meta.source = src;
  } else {
    meta.source = "CF/Radial";
  }

  meta.date = "";
  meta.time = "";
  meta.lowest_sweep_time = "";
  meta.beamwidth = 1.0f;

  for (size_t isweep = 0; isweep < nsweeps; ++isweep) {
    sweep scan;
    scan.beam = radar::beam_propagation{vol.location.alt, fixed_angle_vals[isweep] * 1_deg};
    scan.bins.resize(range_vals.size());
    const auto range0 = range_vals.empty() ? 0.0f : range_vals.front();
    const auto range_step = range_vals.size() > 1 ? (range_vals[1] - range_vals[0]) : 125.0f;
    for (size_t i = 0; i < scan.bins.size(); ++i) {
      scan.bins[i].slant_range = range0 + i * range_step;
      std::tie(scan.bins[i].ground_range, scan.bins[i].altitude) =
        scan.beam.ground_range_altitude(scan.bins[i].slant_range);
    }

    const auto ray_start = static_cast<size_t>(sweep_start_vals[isweep]);
    const auto ray_end = static_cast<size_t>(sweep_end_vals[isweep]);
    const auto ray_count = ray_end >= ray_start ? (ray_end - ray_start + 1) : 0;
    scan.rays.resize(ray_count);
    for (size_t i = 0; i < scan.rays.size(); ++i) {
      const auto ray_index = ray_start + i;
      scan.rays[i] = ray_index < azimuth_vals.size() ? azimuth_vals[ray_index] * 1_deg : 0.0_deg;
    }
    scan.ray_offset = ray_start;

    scan.data = array2f{vec2z{scan.bins.size(), scan.rays.size()}};
    std::fill(scan.data.data(), scan.data.data() + scan.data.size(), nodata);
    vol.sweeps.push_back(std::move(scan));
  }
}

} // namespace

auto discover_fields(io::odim::polar_volume const& vol_odim) -> std::vector<std::string>
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

auto discover_fields(std::filesystem::path const& input_path) -> std::vector<std::string>
{
  try {
    auto file = io::nc::file{input_path.string(), io_mode::read_only};
    return discover_cf_radial_fields(file);
  } catch (std::exception const&) {
    auto vol_odim = io::odim::polar_volume{input_path, io_mode::read_only};
    return discover_fields(vol_odim);
  }
}

auto read_moment_volume(
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

    // Skip vertical-pointing scans
    if (std::fabs(scan_odim.elevation_angle() - 90) < 0.1f)
      continue;

    auto scan = sweep{};
    scan.beam = radar::beam_propagation{vol.location.alt, scan_odim.elevation_angle() * 1_deg};

    // Compute gate geometry
    scan.bins.resize(scan_odim.bin_count());
    auto range_scale = scan_odim.range_scale();
    auto range_start = scan_odim.range_start() * 1000 + range_scale * 0.5;
    for (size_t i = 0; i < scan.bins.size(); ++i) {
      scan.bins[i].slant_range = range_start + i * range_scale;
      std::tie(scan.bins[i].ground_range, scan.bins[i].altitude) =
        scan.beam.ground_range_altitude(scan.bins[i].slant_range);
    }

    // Compute ray azimuths
    scan.rays.resize(scan_odim.ray_count());
    auto ray_scale = 360_deg / scan.rays.size();
    auto ray_start = scan_odim.ray_start() * 1_deg + ray_scale * 0.5;
    for (size_t i = 0; i < scan.rays.size(); ++i)
      scan.rays[i] = ray_start + i * ray_scale;

    // Find and read the requested moment
    for (size_t idata = 0; idata < scan_odim.data_count(); ++idata) {
      auto data_odim = scan_odim.data_open(idata);
      if (data_odim.quantity() != moment)
        continue;

      scan.data.resize(vec2z{(size_t)scan_odim.bin_count(), (size_t)scan_odim.ray_count()});

      // Velocity fields: preserve undetect as a distinct value
      // Other fields: map undetect to nodata
      if (moment == velocity_field)
        data_odim.read_unpack(scan.data.data(), undetect, nodata);
      else
        data_odim.read_unpack(scan.data.data(), nodata, nodata);

      auto rotation = normalize_scan_rays(scan);
      rotate_scan_data_rows(scan, rotation);

      vol.sweeps.push_back(std::move(scan));
      break;
    }
  }

  return vol;
}

auto read_moment_volume(
      std::filesystem::path const& input_path
    , const std::string& moment
    , const std::string& velocity_field
    ) -> volume
{
  try {
    auto file = io::nc::file{input_path.string(), io_mode::read_only};
    auto vol = volume{};
    auto meta = volume_metadata{};
    read_cf_radial_sweeps(file, meta, vol);

    const std::vector<std::string> aliases = make_cf_field_aliases(moment);
    for (auto const& alias : aliases) {
      auto const* var = file.find_variable(alias);
      if (var) {
        auto shape = get_var_shape(*var);
        if (!shape.empty()) {
          const size_t nrows = shape[0];
          const size_t ncols = shape.size() > 1 ? shape[1] : 1;
          std::vector<float> values(nrows * ncols);
          var->read(values);

          const size_t nsweeps = vol.sweeps.size();
          for (size_t isweep = 0; isweep < nsweeps; ++isweep) {
            auto& scan = vol.sweeps[isweep];
            for (size_t iray = 0; iray < scan.rays.size(); ++iray) {
              const size_t ray_index = scan.ray_offset + iray;
              for (size_t igate = 0; igate < scan.bins.size(); ++igate) {
                const size_t linear = ray_index * ncols + igate;
                if (linear < values.size()) {
                  const auto value = values[linear];
                  scan.data.data()[iray * scan.bins.size() + igate] = std::isfinite(value) ? value : nodata;
                }
              }
            }
            auto rotation = normalize_scan_rays(scan);
            rotate_scan_data_rows(scan, rotation);
          }
          break;
        }
      }
    }

    return vol;
  } catch (std::exception const&) {
    auto vol_odim = io::odim::polar_volume{input_path, io_mode::read_only};
    return read_moment_volume(vol_odim, moment, velocity_field);
  }
}

auto read_metadata(io::odim::polar_volume const& vol_odim) -> volume_metadata
{
  volume_metadata meta;

  // Elevation angles
  const size_t nelev = vol_odim.scan_count();
  meta.elevation = array1f{nelev};
  for (size_t iscan = 0; iscan < nelev; ++iscan) {
    auto scan_odim = vol_odim.scan_open(iscan);
    meta.elevation[iscan] = scan_odim.elevation_angle();
  }

  // Nyquist velocities
  meta.nyquist = array1f{nelev};
  for (size_t iscan = 0; iscan < nelev; ++iscan) {
    auto scan_odim = vol_odim.scan_open(iscan);
    if (auto iatt = scan_odim.attributes().find("NI"); iatt != scan_odim.attributes().end())
      meta.nyquist[iscan] = scan_odim.attributes()["NI"].get_real();
    else
      meta.nyquist[iscan] = -9999.f;
  }

  // Lowest sweep time (midpoint of start/end)
  auto elev = 90.f;
  bom::timestamp stdate, eddate;
  for (size_t iscan = nelev - 1; iscan > 0; --iscan) {
    auto scan_odim = vol_odim.scan_open(iscan);
    if (scan_odim.elevation_angle() < elev) {
      stdate = bom::io::odim::strings_to_time(
        scan_odim.attributes()["startdate"].get_string(),
        scan_odim.attributes()["starttime"].get_string());
      eddate = bom::io::odim::strings_to_time(
        scan_odim.attributes()["enddate"].get_string(),
        scan_odim.attributes()["endtime"].get_string());
    }
    elev = scan_odim.elevation_angle();
  }
  meta.lowest_sweep_time = to_string(stdate + (eddate - stdate) / 2);

  // Global attributes
  const auto& attributes = vol_odim.attributes();
  meta.location.lat = vol_odim.latitude() * 1_deg;
  meta.location.lon = vol_odim.longitude() * 1_deg;
  meta.location.alt = vol_odim.height();
  meta.source = attributes["source"].get_string();
  meta.date = attributes["date"].get_string();
  meta.time = attributes["time"].get_string();
  meta.beamwidth = attributes["beamwH"].get_real();

  return meta;
}

auto read_metadata(std::filesystem::path const& input_path) -> volume_metadata
{
  try {
    auto file = io::nc::file{input_path.string(), io_mode::read_only};
    auto meta = volume_metadata{};
    auto vol = volume{};
    read_cf_radial_sweeps(file, meta, vol);
    return meta;
  } catch (std::exception const&) {
    auto vol_odim = io::odim::polar_volume{input_path, io_mode::read_only};
    return read_metadata(vol_odim);
  }
}
