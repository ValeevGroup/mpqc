#pragma once
#ifndef MPQC_UTILITY_PARALLELBREAKPOINT_H
#define MPQC_UTILITY_PARALLELBREAKPOINT_H

#include <tiledarray.h>

namespace mpqc {
namespace utility {

inline void parallel_break_point(madness::World &world, volatile int debug) {
    if (0 != debug) {
        char hostname[256];
        gethostname(hostname, sizeof(hostname));
        printf("PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        if (world.rank() == 1) {
            while (0 != debug) sleep(5);
        }
    }
    world.gop.fence();
}

} // namespace utility
} // namespace namespace mpqc

#endif // MPQC_UTILITY_PARALLELBREAKPOINT_H
