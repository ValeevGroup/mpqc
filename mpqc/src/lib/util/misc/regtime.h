//
// regtime.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_misc_regtime_h
#define _util_misc_regtime_h

#include <iostream>
#include <util/class/class.h>

class TimedRegion;

DescribedClass_REF_fwddec(RegionTimer);
REF_fwddec(KeyVal);

class RegionTimer: public DescribedClass {
#  define CLASSNAME RegionTimer
#  define HAVE_KEYVAL_CTOR
#  include <util/class/classd.h>
  protected:
    int wall_time_;
    int cpu_time_;
    int flops_;

    TimedRegion *top_;
    TimedRegion *current_;
    TimedRegion *default_;

  public:
    RegionTimer(const char *topname = "total",
                int cpu_time = 0, int wall_time = 1);
    RegionTimer(const RefKeyVal &);
    ~RegionTimer();
    void enter(const char * = 0);
    void change(const char *newname, const char * oldname = 0);
    void exit(const char * = 0);
    void set_default(const char *);
    void unset_default();
    void enter_default();
    void exit_default();
    virtual void print(std::ostream& = ExEnv::out()) const;

    void update_top() const;

    int nregion() const;
    void get_region_names(const char *names[]) const;
    void get_wall_times(double *) const;
    void get_cpu_times(double *) const;
    void get_flops(double *) const;
    void get_depth(int *) const;

    double get_wall_time() const;
    double get_cpu_time() const;
    double get_flops() const;

    void add_wall_time(const char *, double);
    void add_cpu_time(const char *, double);
    void add_flops(const char *, double);

    static RegionTimer *default_regiontimer();
    static void set_default_regiontimer(const RefRegionTimer &);
};
DescribedClass_REF_dec(RegionTimer);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
