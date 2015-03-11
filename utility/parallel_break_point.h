#pragma once
#ifndef TCC_UTILITY_PARALLELBREAKPOINT_H
#define TCC_UTILITY_PARALLELBREAKPOINT_H

#include "../include/tiledarray.h"

namespace tcc {
namespace utility {

void parallal_break_point(madness::World &world, volatile int debug) {
    if (0 != debug) {
        char hostname[256];
        gethostname(hostname, sizeof(hostname));
        printf("PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        if (world.rank() == 0) {
            while (0 != debug) sleep(5);
        }
    }
    world.gop.fence();
}

} // namespace utility
} // namespace namespace tcc

#endif // TCC_UTILITY_PARALLELBREAKPOINT_H
