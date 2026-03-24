#include "array_operations.h"

// auto argmin(array1f const& x) -> int{  
//   int pos = std::distance(x.begin(), std::min_element(x.begin(), x.end()));
//   return pos;
// }

auto copy_mask(array2f const& ref, array2f& dest) -> void
{
  for (size_t y = 0; y < ref.extents().y; ++y)
  {
    for (size_t x = 0; x < ref.extents().x; ++x)
    {
      if(std::isnan(ref[y][x]) || std::fabs(ref[y][x] - undetect) < 0.001f)
        dest[y][x] = nodata;
      if(std::fabs(dest[y][x] - undetect) < 0.001f)
        dest[y][x] = nodata;
    }
  }
}

auto flip(array1d& data) -> void
{
  auto n = data.size();
  for(size_t i = 0; i < n / 2; i++)
  {
    std::swap(data[i], data[n - 1 - i]);
  }
}

auto flipud(array2f& data) -> void
{
  auto ny = data.extents().y;
  auto nx = data.extents().x;
  for (size_t y = 0; y < ny / 2; ++y)
  {
    for (size_t x = 0; x < nx; ++x)
    {
      std::swap(data[y][x], data[ny - 1 - y][x]);
    }
  }
}

auto fliplr(array2f& data) -> void
{
  auto ny = data.extents().y;
  auto nx = data.extents().x;
  for (size_t y = 0; y < ny; ++y)
  {
    for (size_t x = 0; x < nx / 2; ++x)
    {
      std::swap(data[y][x], data[y][nx - 1 - x]);
    }
  }
}

auto mean(const array1f& x) -> float {
  auto sum = 0.0f;
  auto count = 0;

  for (const auto& value : x) {
    if (std::isnan(value))
      continue;

    sum += value;
    count++;
  }

  if (count == 0)
    return nodata;

  return sum / count;
}