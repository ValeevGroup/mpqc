/*
 * timer.h
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

#ifndef _util_misc_timer_h
#define _util_misc_timer_h

#include <mpqc_config.h>

#define DEPRECATE_TIMER 1

#if DEPRECATE_TIMER
# ifndef _util_misc_regtime_cc
#  warning "util/misc/timer.h is deprecated"
# endif
# define TIMER_DEPRECATED DEPRECATED
#else
# define TIMER_DEPRECATED
#endif

namespace sc {

  // These functions are all deprecated.  Please see the "Exceptions
  // and Region Timers" section in the MPQC manual
  // (http://www.mpqc.org/mpqc-html/develop.html#scexcepttimer) for
  // information on how to do timing calls in an exception-safe manner.
  void tim_enter(const char *) TIMER_DEPRECATED;
  void tim_exit(const char *) TIMER_DEPRECATED;
  void tim_change(const char *) TIMER_DEPRECATED;
  void tim_set_default(const char *) TIMER_DEPRECATED;
  void tim_enter_default() TIMER_DEPRECATED;
  void tim_exit_default() TIMER_DEPRECATED;
  void tim_print(int) TIMER_DEPRECATED;

}

#endif /* _util_misc_timer_h */
