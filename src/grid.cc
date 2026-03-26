#include "grid.h"

auto init_altitudes(io::configuration const& config) -> array1f
{
  auto alts = array1f{config["layer_count"]};
  auto base = float(config["altitude_base"]);
  auto step = float(config["altitude_step"]);
  for (size_t i = 0; i < alts.size(); ++i)
    alts[i] = base + step * i;
  return alts;
}

auto flip(array1d& data) -> void
{
  auto n = data.size();
  for (size_t i = 0; i < n / 2; i++)
    std::swap(data[i], data[n - 1 - i]);
}

auto flipud(array2f& data) -> void
{
  auto ny = data.extents().y;
  auto nx = data.extents().x;
  for (size_t y = 0; y < ny / 2; ++y)
    for (size_t x = 0; x < nx; ++x)
      std::swap(data[y][x], data[ny - 1 - y][x]);
}
