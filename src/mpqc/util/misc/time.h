#ifndef MPQC_UTILITY_TIME_H
#define MPQC_UTILITY_TIME_H

#include <chrono>

#include <madness/world/MADworld.h>

namespace mpqc {

using time_point = std::chrono::high_resolution_clock::time_point;

inline time_point now() { return std::chrono::high_resolution_clock::now(); }

inline std::chrono::system_clock::time_point system_now() {
  return std::chrono::system_clock::now();
}

inline double duration_in_s(time_point const &t0, time_point const &t1) {
  return std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0)
      .count();
}

inline time_point fenced_now(madness::World &world) {
  world.gop.fence();
  return now();
}
/**
 *  return fence() time based on bool fence
 */
inline time_point now(madness::World &world, bool fence) {
  if (fence) {
    world.gop.fence();
  }
  return now();
}

}  // namespace mpqc

#endif  // MPQC_UTILITY_TIME_H
