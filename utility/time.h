#pragma once
#ifndef MPQC_UTILITY_TIME_H
#define MPQC_UTILITY_TIME_H

#include "../common/typedefs.h"

#include <chrono>

namespace mpqc {
namespace utility {
namespace time {

using t_point = std::chrono::high_resolution_clock::time_point;

inline t_point now() { return std::chrono::high_resolution_clock::now(); }

inline double duration_in_s(t_point const &t0, t_point const &t1) {
    return std::chrono::duration_cast<std::chrono::duration<double>>(t1 - t0)
        .count();
}

inline t_point fenced_now(madness::World &world){
    world.gop.fence();
    return now();
}

} // namespace time


} // namespace utility
} // namespace mpqc

namespace mpqc_time = mpqc::utility::time;

#endif // MPQC_UTILITY_TIME_H
