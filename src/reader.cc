#include "reader.h"

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

      vol.sweeps.push_back(std::move(scan));
      break;
    }
  }

  return vol;
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
  meta.source = attributes["source"].get_string();
  meta.date = attributes["date"].get_string();
  meta.time = attributes["time"].get_string();
  meta.beamwidth = attributes["beamwH"].get_real();

  return meta;
}
