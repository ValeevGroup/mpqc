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
#include <string>
#include <util/class/class.h>

namespace sc {

class TimedRegion {
  private:
    char *name_;
    TimedRegion *up_;
    TimedRegion *subregions_;
    TimedRegion *next_;
    TimedRegion *prev_;
    double cpu_time_;
    double wall_time_;
    double cpu_enter_;
    double wall_enter_;
    double flops_;
    double flops_enter_;

    TimedRegion *insert_after(const char *name);
    TimedRegion *insert_before(const char *name);
  public:
    TimedRegion(const char *name);
    ~TimedRegion();
    const char *name() const { return name_; }
    TimedRegion *findinsubregion(const char *);
    void cpu_enter(double);
    void wall_enter(double);
    void flops_enter(double);
    void cpu_exit(double);
    void wall_exit(double);
    void flops_exit(double);
    void cpu_add(double t) { cpu_time_ += t; }
    void wall_add(double t) { wall_time_ += t; }
    void flops_add(double t) { flops_ += t; }
    TimedRegion *up() const { return up_; }
    TimedRegion *subregions() const { return subregions_; }
    TimedRegion *next() const { return next_; }
    TimedRegion *prev() const { return prev_; }

    /// Include the regions in r in this object's regions.
    void merge(const TimedRegion* r);

    int nregion();
    void get_region_names(const char *names[]);
    void get_wall_times(double *);
    void get_cpu_times(double *);
    void get_flops(double *);
    void get_depth(int *, int depth = 0);
};

/** The RegionTimer class is used to record the time spent in a section of
code.  During the run of a code, enter and exit members are called to begin
and end timed sections.  The print member is used to display the obtained
times.  Multiple enter calls for a region with the same name aggregate the
timings. Nested regions are supported. */
class RegionTimer: public DescribedClass {
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
    RegionTimer(const Ref<KeyVal> &);
    ~RegionTimer();
    void enter(const char * = 0);
    void change(const char *newname, const char * oldname = 0);
    void exit(const char * = 0, bool do_not_throw = false);
    void set_default(const char *);
    void unset_default();
    void enter_default();
    void exit_default();
    virtual void print(std::ostream& = ExEnv::out0()) const;

    /// Include the regions in r in this object's regions.
    void merge(const Ref<RegionTimer> &r);

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
    static void set_default_regiontimer(const Ref<RegionTimer> &);
};

/** The Timer class uses RegionTimer to time intervals in an exception safe
manner.  It will automatically call RegionTimer::enter when its constructor
is called and RegionTimer::exit when its destructor is called.  The reset
member can also result in RegionTimer's enter and exit routines being
called.  The programmer is responsible for making sure that timers are
exited in the reverse of the order that they are entered.  */
class Timer {
    Ref<RegionTimer> timer_;
    std::string name_;
    bool active_;
  public:
    /** Start timing a region using the default RegionTimer and activate
        the timer.  If a null name pointer is given, then the
        timer will not be activated. */
    Timer(const char *name);
    /** Start timing a region using the given RegionTimer. If a null name
        pointer is given, then the timer will not be activated. */
    Timer(const Ref<RegionTimer> &, const char *name);
    /** Stop timing a region, if active. */
    ~Timer();
    /** Stop timing the current region, if active.  If a new region name is
        passed in, start timing with that name.  If no region name is
        given, the Timer will be deactivated.  */
    void reset(const char * = 0);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
