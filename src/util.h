#ifndef VARGRID_UTIL_H
#define VARGRID_UTIL_H

#include "types.h"
#include <chrono>
#include <set>

// RAII timer that logs phase name on entry and elapsed time on exit.
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

  // Non-copyable
  phase_timer(const phase_timer&) = delete;
  phase_timer& operator=(const phase_timer&) = delete;
};

// Split a space-separated string into a set of tokens.
inline auto split_fields(const std::string& s) -> std::set<std::string>
{
  std::set<std::string> result;
  std::istringstream iss(s);
  std::string token;
  while (iss >> token)
    result.insert(token);
  return result;
}

#endif // VARGRID_UTIL_H
