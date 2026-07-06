#include "leroi.h"

#include <atomic>
#include <stdexcept>
#include <thread>

auto parse_leroi_config(io::configuration const& config) -> leroi_config
{
  leroi_config cfg;

  std::string w = config.optional("leroi_weight", "barnes");
  if      (w == "barnes")   cfg.weight = leroi::weight_type::barnes;
  else if (w == "cressman") cfg.weight = leroi::weight_type::cressman;
  else if (w == "idw")      cfg.weight = leroi::weight_type::idw;
  else if (w == "bilinear") cfg.weight = leroi::weight_type::bilinear;
  else
    throw std::runtime_error(
      "invalid leroi_weight '" + w + "' (expected barnes, cressman, idw or bilinear)");

  cfg.roi              = std::stof(config.optional("leroi_roi", "0"));
  cfg.idw_pwr          = std::stof(config.optional("leroi_idw_pwr", "2.0"));
  cfg.ground_elevation = std::stof(config.optional("leroi_ground_elevation", "-999"));
  return cfg;
}

// Sweep indices sorted by elevation, excluding birdbath (vertical) scans.
static auto used_sweep_indices(volume const& vol) -> std::vector<size_t>
{
  std::vector<size_t> idx;
  for (size_t i = 0; i < vol.sweeps.size(); ++i)
    if (vol.sweeps[i].beam.elevation().degrees() < 85.0)
      idx.push_back(i);
  std::sort(idx.begin(), idx.end(), [&vol](size_t a, size_t b) {
    return vol.sweeps[a].beam.elevation() < vol.sweeps[b].beam.elevation();
  });
  return idx;
}

static auto make_geometry(sweep const& s) -> leroi::sweep_geometry
{
  std::vector<float> az(s.rays.size());
  for (size_t i = 0; i < s.rays.size(); ++i)
    az[i] = static_cast<float>(s.rays[i].degrees());

  std::vector<float> gr(s.bins.size()), alt(s.bins.size());
  for (size_t i = 0; i < s.bins.size(); ++i) {
    gr[i] = s.bins[i].ground_range;
    alt[i] = s.bins[i].altitude;
  }

  return leroi::make_sweep_geometry(
    static_cast<float>(s.beam.elevation().degrees()), az, gr, alt);
}

auto leroi_precompute(
      std::map<std::string, volume> const& volumes
    , grid_bearings const& gb
    , leroi_config const& cfg
    , float max_altitude
    ) -> leroi_grids
{
  leroi_grids pre;
  pre.nx = gb.nx;
  pre.ny = gb.ny;
  size_t ncells = gb.nx * gb.ny;

  auto const& ref_vol = volumes.begin()->second;
  auto sweep_idx = used_sweep_indices(ref_vol);
  size_t ns = sweep_idx.size();
  if (ns < 2)
    throw std::runtime_error("leroi requires at least 2 non-vertical sweeps");

  // Field list and per-field sanity check (all volumes come from the same
  // ODIM file, so sweep structure should match the reference volume).
  std::vector<std::string> field_names;
  for (auto& [name, vol] : volumes) {
    if (vol.sweeps.size() != ref_vol.sweeps.size()) {
      trace::warning("leroi: field {} has {} sweeps (expected {}), skipping",
        name, vol.sweeps.size(), ref_vol.sweeps.size());
      continue;
    }
    field_names.push_back(name);
  }
  size_t nfields = field_names.size();

  // Geometry is shared across fields (same file), so build it once per sweep.
  std::vector<leroi::sweep_geometry> geoms(ns);
  for (size_t s = 0; s < ns; ++s)
    geoms[s] = make_geometry(ref_vol.sweeps[sweep_idx[s]]);

  // Radius of influence (unused by the bilinear stencil)
  float roi = cfg.roi;
  if (cfg.weight != leroi::weight_type::bilinear && roi <= 0.0f) {
    float rmax = std::sqrt(gb.max_range * gb.max_range + max_altitude * max_altitude);
    roi = leroi::auto_roi(geoms, rmax);
    trace::log("leroi: radius of influence set to {:.0f} m", roi);
  }

  leroi::cell_lookup cells{gb.cell_bearing_deg.data(), gb.cell_range.data(), ncells};

  // Allocate outputs (all map keys inserted before threads start; workers
  // only touch preallocated slots through raw pointers).
  pre.heights.assign(ns, std::vector<float>(ncells, leroi::core_nan));
  std::vector<std::vector<std::vector<float>>*> surface_ptrs(nfields);
  for (size_t f = 0; f < nfields; ++f) {
    pre.surfaces[field_names[f]].assign(ns, std::vector<float>(ncells, leroi::core_nan));
    surface_ptrs[f] = &pre.surfaces[field_names[f]];
  }

  // Interpolate each sweep (all fields at once); parallel over sweeps.
  std::atomic<size_t> next{0};
  auto worker = [&]() {
    std::vector<std::vector<float>> out(nfields, std::vector<float>(ncells));
    while (true) {
      size_t s = next.fetch_add(1);
      if (s >= ns) break;

      std::vector<const float*> field_data(nfields);
      size_t nbins = ref_vol.sweeps[sweep_idx[s]].bins.size();
      for (size_t f = 0; f < nfields; ++f)
        field_data[f] = volumes.at(field_names[f]).sweeps[sweep_idx[s]].data.data();

      leroi::interp_sweep(geoms[s], field_data, nbins, cells,
        cfg.weight, roi, cfg.idw_pwr, undetect, out);

      for (size_t f = 0; f < nfields; ++f)
        (*surface_ptrs[f])[s] = out[f];

      // Sweep surface heights (geometry only)
      auto& h = pre.heights[s];
      if (s == 0 && geoms[s].elevation_deg <= cfg.ground_elevation) {
        // Anchor the lowest sweep at ground level where it has coverage.
        // Note vargrid altitudes are ASL, so "ground" is the radar station
        // altitude (leroi.py works radar-relative, where this would be 0).
        float ground_alt = static_cast<float>(ref_vol.location.alt);
        for (size_t c = 0; c < ncells; ++c)
          if (!std::isnan(leroi::height_at(geoms[s], cells.range[c])))
            h[c] = ground_alt;
      } else {
        for (size_t c = 0; c < ncells; ++c)
          h[c] = leroi::height_at(geoms[s], cells.range[c]);
      }
    }
  };

  size_t nthreads = std::min(static_cast<size_t>(std::thread::hardware_concurrency()), ns);
  if (nthreads == 0) nthreads = 1;
  std::vector<std::thread> threads;
  for (size_t t = 0; t < nthreads; ++t)
    threads.emplace_back(worker);
  for (auto& t : threads)
    t.join();

  auto outside = leroi::fill_heights(pre.heights, ncells);
  if (outside > 0)
    trace::debug("leroi: {} grid columns outside radar range", outside);

  return pre;
}

auto leroi_slice(
      leroi_grids const& pre
    , std::string const& field
    , float altitude
    ) -> array2f
{
  auto out = array2f{vec2z{pre.nx, pre.ny}};
  auto it = pre.surfaces.find(field);
  if (it == pre.surfaces.end()) {
    std::fill(out.data(), out.data() + out.size(), nodata);
    return out;
  }
  leroi::vertical_slice(pre.heights, it->second, altitude, out.data(), pre.nx * pre.ny);
  return out;
}
