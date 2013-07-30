
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#ifndef _chemistry_qc_lmp2_dgemminfo_h
#define _chemistry_qc_lmp2_dgemminfo_h

#ifndef USE_MICROTIME
#include <sys/time.h>
#include <time.h>
#endif

namespace sc {

#ifdef USE_MICROTIME

  static inline void microtime(unsigned *lo, unsigned *hi) {
    __asm __volatile (
		      ".byte 0x0f; .byte 0x31; movl    %%edx,%0; movl    %%eax,%1"
		      : "=g" (*hi), "=g" (*lo) :: "eax", "edx");
  }

  static inline double cpu_walltime(void) {
    unsigned lo, hi;
    microtime(&lo, &hi);
    unsigned long long ticks = lo + (((unsigned long long)hi)<<32);
    double secs = ticks*(1.0/(1496.612*1.0e6));
    return secs;
  }

#else

  /// Returns the time in seconds.
  static inline double cpu_walltime(void) {
    struct timeval tod;
    gettimeofday(&tod,0);
    return tod.tv_sec + 0.000001 * tod.tv_usec;
  }

#endif

  /// Records information about the time take to perform a DGEMM operation.
  extern void count_dgemm(int n, int l, int m, double t);
  
}

#endif
