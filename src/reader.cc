#include "reader.h"

#include <limits>
#include <optional>

namespace {

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

  // Deduplicate but preserve order: the explicitly requested field name must
  // be tried first, before any aliases.
  std::set<std::string> seen;
  std::vector<std::string> unique;
  for (auto const& a : aliases)
    if (seen.insert(a).second) unique.push_back(a);
  return unique;
}

static auto get_var_shape(io::nc::variable const& var) -> std::vector<size_t>
{
  std::vector<size_t> shape;
  for (auto const* dim : var.dimensions())
    shape.push_back(dim->size());
  return shape;
}

static auto read_optional_string_attribute(io::nc::variable const& var, const char* name) -> std::string
{
  std::string value;
  try {
    var.att_get(name, value);
  } catch (std::exception const&) {
    value.clear();
  }
  return value;
}

template <typename T>
static auto read_optional_numeric_attribute(io::nc::variable const& var, const char* name) -> std::optional<double>
{
  T value{};
  try {
    var.att_get(name, value);
    return static_cast<double>(value);
  } catch (std::exception const&) {
    return std::nullopt;
  }
}

static auto read_optional_double_attribute(io::nc::variable const& var, const char* name) -> std::optional<double>
{
  if (auto v = read_optional_numeric_attribute<double>(var, name); v.has_value()) return v;
  if (auto v = read_optional_numeric_attribute<float>(var, name); v.has_value()) return v;
  if (auto v = read_optional_numeric_attribute<long long>(var, name); v.has_value()) return v;
  if (auto v = read_optional_numeric_attribute<int>(var, name); v.has_value()) return v;
  if (auto v = read_optional_numeric_attribute<short>(var, name); v.has_value()) return v;
  if (auto v = read_optional_numeric_attribute<unsigned long long>(var, name); v.has_value()) return v;
  if (auto v = read_optional_numeric_attribute<unsigned int>(var, name); v.has_value()) return v;
  if (auto v = read_optional_numeric_attribute<unsigned short>(var, name); v.has_value()) return v;
  return std::nullopt;
}

static auto is_missing_sentinel(float value, const std::vector<float>& sentinels) -> bool
{
  for (auto s : sentinels) {
    if (!std::isfinite(s))
      continue;
    auto scale = std::max({1.0f, std::fabs(value), std::fabs(s)});
    auto tol = 1.0e-5f * scale;
    if (std::fabs(value - s) <= tol)
      return true;
  }
  return false;
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

  // Validate ray index ranges before any size_t conversion: negative or
  // out-of-range values from a malformed file must not wrap around.
  for (size_t i = 0; i < sweep_start_vals.size(); ++i) {
    if (sweep_start_vals[i] < 0 || sweep_end_vals[i] < 0
        || static_cast<size_t>(sweep_end_vals[i]) >= azimuth_vals.size())
      throw std::runtime_error("CF/Radial sweep ray indices are out of range");
  }

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
    // Skip vertical-pointing (birdbath) scans, matching the ODIM reader.
    if (std::fabs(fixed_angle_vals[isweep] - 90.0f) < 0.1f)
      continue;

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
    // libnetcdf can successfully open any HDF5 file (including ODIM polar
    // volumes), so verify this really is CF/Radial before trusting the
    // discovery result — otherwise fall through to the ODIM reader.
    if (!file.find_variable("range") || !file.find_variable("azimuth")
        || !file.find_variable("time") || !file.find_variable("sweep_start_ray_index"))
      throw std::runtime_error("not a CF/Radial file");
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
          std::vector<float> missing_sentinels;
          const auto fill_value = read_optional_double_attribute(*var, "_FillValue");
          const auto missing_value = read_optional_double_attribute(*var, "missing_value");
          const auto scale_factor = read_optional_double_attribute(*var, "scale_factor");
          const auto add_offset = read_optional_double_attribute(*var, "add_offset");
          const bool has_packing = scale_factor.has_value() && add_offset.has_value();

          auto add_missing = [&](double v) {
            auto fv = static_cast<float>(v);
            if (std::isfinite(fv))
              missing_sentinels.push_back(fv);
            if (scale_factor.has_value() && add_offset.has_value()) {
              auto unpacked = static_cast<float>(v * scale_factor.value() + add_offset.value());
              if (std::isfinite(unpacked))
                missing_sentinels.push_back(unpacked);
            }
          };

          if (fill_value.has_value()) add_missing(fill_value.value());
          if (missing_value.has_value()) add_missing(missing_value.value());

          if (shape.size() != 2)
            continue;

          const size_t dim0 = shape[0];
          const size_t dim1 = shape[1];
          std::vector<float> values(dim0 * dim1);
          var->read(values);

          size_t total_rays = 0;
          size_t max_gates = 0;
          for (auto const& s : vol.sweeps) {
            total_rays = std::max(total_rays, s.ray_offset + s.rays.size());
            max_gates = std::max(max_gates, s.bins.size());
          }

          // CF/Radial moments are usually (time, range), but some files use
          // (range, time). Detect and support both layouts.
          enum class var_layout { time_range, range_time };
          var_layout layout = var_layout::time_range;
          if (dim0 == max_gates && dim1 == total_rays)
            layout = var_layout::range_time;
          else if (dim0 == total_rays && dim1 == max_gates)
            layout = var_layout::time_range;
          else if (dim0 == max_gates)
            layout = var_layout::range_time;

          const size_t time_len = (layout == var_layout::time_range) ? dim0 : dim1;
          const size_t range_len = (layout == var_layout::time_range) ? dim1 : dim0;

          const size_t nsweeps = vol.sweeps.size();
          for (size_t isweep = 0; isweep < nsweeps; ++isweep) {
            auto& scan = vol.sweeps[isweep];
            for (size_t iray = 0; iray < scan.rays.size(); ++iray) {
              const size_t ray_index = scan.ray_offset + iray;
              if (ray_index >= time_len)
                continue;
              for (size_t igate = 0; igate < scan.bins.size(); ++igate) {
                if (igate >= range_len)
                  continue;
                const size_t linear = (layout == var_layout::time_range)
                  ? (ray_index * dim1 + igate)
                  : (igate * dim1 + ray_index);
                if (linear < values.size()) {
                  const auto value = values[linear];
                  const auto unpacked = has_packing
                    ? static_cast<float>(value * scale_factor.value() + add_offset.value())
                    : value;
                  scan.data.data()[iray * scan.bins.size() + igate] =
                    (std::isfinite(value) && !is_missing_sentinel(value, missing_sentinels))
                    ? unpacked
                    : nodata;
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

  // Lowest sweep time (midpoint of start/end of the lowest-elevation scan).
  // The previous loop skipped scan 0 entirely — for ascending volumes that
  // is precisely the lowest sweep — and underflowed when nelev was 0.
  meta.lowest_sweep_time = "";
  if (nelev > 0) {
    size_t lowest = 0;
    for (size_t i = 1; i < nelev; ++i)
      if (meta.elevation[i] < meta.elevation[lowest])
        lowest = i;
    try {
      auto scan_odim = vol_odim.scan_open(lowest);
      auto stdate = bom::io::odim::strings_to_time(
        scan_odim.attributes()["startdate"].get_string(),
        scan_odim.attributes()["starttime"].get_string());
      auto eddate = bom::io::odim::strings_to_time(
        scan_odim.attributes()["enddate"].get_string(),
        scan_odim.attributes()["endtime"].get_string());
      meta.lowest_sweep_time = to_string(stdate + (eddate - stdate) / 2);
    } catch (std::exception const&) {
      trace::warning("Could not read start/end time of the lowest sweep");
    }
  }

  // Global attributes. Optional ones fall back to defaults instead of
  // aborting the run (matches the CF/Radial reader's behaviour).
  const auto& attributes = vol_odim.attributes();
  meta.location.lat = vol_odim.latitude() * 1_deg;
  meta.location.lon = vol_odim.longitude() * 1_deg;
  meta.location.alt = vol_odim.height();

  auto get_str = [&](const char* name) -> std::string {
    if (auto it = attributes.find(name); it != attributes.end())
      return attributes[name].get_string();
    trace::debug("ODIM attribute '{}' missing, using empty string", name);
    return "";
  };
  meta.source = get_str("source");
  meta.date = get_str("date");
  meta.time = get_str("time");

  if (auto it = attributes.find("beamwH"); it != attributes.end()) {
    meta.beamwidth = attributes["beamwH"].get_real();
  } else {
    meta.beamwidth = 1.0f;
    trace::debug("ODIM attribute 'beamwH' missing, assuming 1.0 deg beamwidth");
  }

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

auto read_field_metadata(
      std::filesystem::path const& input_path
    , const std::vector<std::string>& fields
    ) -> std::map<std::string, variable_metadata>
{
  std::map<std::string, variable_metadata> out;

  try {
    auto file = io::nc::file{input_path.string(), io_mode::read_only};

    // Only attempt CF/Radial variable metadata extraction.
    if (!file.find_variable("range") || !file.find_variable("azimuth")
        || !file.find_variable("time") || !file.find_variable("sweep_start_ray_index"))
      throw std::runtime_error("not a CF/Radial file");

    for (auto const& field : fields) {
      const std::vector<std::string> aliases = make_cf_field_aliases(field);
      for (auto const& alias : aliases) {
        auto const* var = file.find_variable(alias);
        if (!var)
          continue;

        variable_metadata meta;
        meta.units = read_optional_string_attribute(*var, "units");
        meta.standard_name = read_optional_string_attribute(*var, "standard_name");
        meta.long_name = read_optional_string_attribute(*var, "long_name");
        meta.description = read_optional_string_attribute(*var, "description");

        if (meta.description.empty())
          meta.description = read_optional_string_attribute(*var, "comment");

        out[field] = std::move(meta);
        break;
      }
    }
  } catch (std::exception const&) {
    // ODIM and unsupported formats typically do not provide CF variable
    // metadata in a directly mappable form; leave per-field metadata empty.
  }

  return out;
}
