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

#ifndef _util_misc_regtime_h
#define _util_misc_regtime_h

#include <iostream>
#include <string>
#include <list>
#ifdef MPQC_NEW_FEATURES
#include <chrono>
#endif
#include <mpqc_config.h>
#include <util/class/class.h>

namespace sc {

#if MPQC_NEW_FEATURES
class MultiThreadTimer;
#endif

/** TimedRegion is a helper class for RegionTimer. */
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
    void acquire_subregion(TimedRegion*);
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

    /** reports the current values of cpu_time, wall_time, and flops  */
    //@{
    double cpu_time() const { return cpu_time_; }
    double wall_time() const { return wall_time_; }
    double flops() const { return flops_; }
    //@}

    /// Include the regions in r in this object's regions.
    void merge(const TimedRegion* r);

    int nregion();
    void get_region_names(const char *names[]);
    void get_wall_times(double *);
    void get_cpu_times(double *);
    void get_flops(double *);
    void get_depth(int *, int depth = 0);

#if MPQC_NEW_FEATURES
    friend class MultiThreadTimer;
#endif
    friend class RegionTimer;
};

/** The RegionTimer class is used to record the time spent in a section of
code.  Except for the creation of an initial RegionTimer, this class should
usually not be used directly.  Instead use the Timer class to control the
RegionTimer in an exception safe manner. */
class RegionTimer: public DescribedClass {
  protected:
    int wall_time_;
    int cpu_time_;
    int flops_;

    TimedRegion *top_;
    TimedRegion *current_;
    TimedRegion *default_;
    std::list<TimedRegion*> defaults_;

  public:
    RegionTimer(const char *topname = "total",
                int cpu_time = 0, int wall_time = 1);
    RegionTimer(const Ref<KeyVal> &);
    ~RegionTimer();
    void enter(const char * = 0);
    void change(const char *newname, const char * oldname = 0);
    void exit(const char * = 0, bool do_not_throw = false);
    double cpu_time(const char * name = 0) const;
    double wall_time(const char * name = 0) const;
    double flops(const char * name = 0) const;
    void set_default(const char *);
    void unset_default();
    void enter_default();
    void exit_default();
    virtual void print(std::ostream& = ExEnv::out0()) const;
    void reset();

    /// Include the regions in r in this object's regions.
    void merge(const Ref<RegionTimer> &r);

    void update_top() const;

    int nregion() const;
    void get_region_names(const char *names[]) const;
    void get_wall_times(double *) const;
    void get_cpu_times(double *) const;
    void get_flops(double *) const;
    void get_depth(int *) const;

    /// @return the time reported by the system clock (the number of seconds since Epoch)
    /// @note precision is about a microsecond or less (see documentation for \c gettimeofday ).
    static double get_wall_time();
    /// @return the CPU time (in seconds) used by this process, as reported by \c getrusage
    static double get_cpu_time();
    static double get_flops();

#ifdef MPQC_NEW_FEATURES
    /// @return the time_point object for this moment in time
    static std::chrono::time_point<std::chrono::high_resolution_clock>
    get_time_point();
#endif

    void add_wall_time(const char *, double);
    void add_cpu_time(const char *, double);
    void add_flops(const char *, double);

    /// Add r as a subregion of current_.  After invocation, r is owned by this RegionTimer
    void acquire_timed_region(TimedRegion* r);

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
    int depth_;
    int default_depth_;
    bool default_entered_;
  public:
    /** Start timing a region using the default RegionTimer and activate
        the timer.  If a null name pointer is given, then the
        timer will not be activated. */
    Timer(const char *name);
    /** Start timing a region using the default RegionTimer. */
    Timer(const std::string &name);
    /** Start timing a region using the given RegionTimer. If a null name
        pointer is given, then the timer will not be activated. */
    Timer(const Ref<RegionTimer> &, const char *name);
    /** Start timing a region using the given RegionTimer. */
    Timer(const Ref<RegionTimer> &, const std::string &name);
    /** Construct a Timer object using the giving RegionTimer, but do
        not begin timing a region. */
    Timer(const Ref<RegionTimer> &);
    /** Construct a Timer object using the default RegionTimer and do not
     * begin timing a region. */
    Timer();
    /** Stop timing a region, if active. */
    ~Timer();
    /** Print the timings held by this object's RegionTimer. */
    void print(std::ostream& = ExEnv::out0()) const;
    /** Stop timing the current region, if active.  If a new region name is
        passed in, start timing with that name.  If no region name is
        given, the Timer will be deactivated.  This member is
        deprecated. */
    DEPRECATED void reset(const char * = 0);

    /** Begin timing, using the given timing region name.  Nested regions
        are supported.  That is, after a region is entered, another may be
        entered before the first is exited.  The time for the nested region
        is included in the time for the containing region.  */
    //@{
    void enter(const char *region);
    void enter(const std::string &region) { enter(region.c_str()); }
    //@}
    /** Change the current timing region to the one specified by the given
        name. */
    //@{
    void change(const char *region);
    void change(const std::string &region) { change(region.c_str()); }
    //@}
    /** Exit the current timing region.  The name optionally can be given
        for a consistency check. */
    //@{
    void exit(const char *region = 0);
    void exit(const std::string &region) { this->exit(region.c_str()); }
    //@}
    /** Default timing regions are provided as a performance optimization.
        The set_default member is used to specify the default region.
        Timing is begun with enter_default, which does not need to lookup
        the region as the enter member must do.  Timing is stopped with
        exit_default.  Default regions can be nested, and unset_default
        restores the default in effect before the previous set_default
        call.
    */
    //@{
    void set_default(const char *region);
    void set_default(const std::string &r) { set_default(r.c_str()); }
    void unset_default();
    void enter_default();
    void exit_default();
    //@}
    /** Query the timers and flop counter for the region */
    //@{
    double cpu_time(const char *region) const;
    double cpu_time(const std::string &region) const { return cpu_time(region.c_str()); }
    double wall_time(const char *region) const;
    double wall_time(const std::string &region) const { return wall_time(region.c_str()); }
    double flops(const char *region) const;
    double flops(const std::string &region) const { return flops(region.c_str()); }
    //@}

#if MPQC_NEW_FEATURES
    void insert(const MultiThreadTimer& timer);
#endif
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
