//
// regtime.cc
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

#ifdef HAVE_CONFIG_H
#  include <mpqc_config.h>
#endif

#include <stdexcept>

#include <math.h>
#include <iostream>
#include <iomanip>

// getrusage and gettimeofday don't exit under SUNMOS
// so if NX is being used call dclock() instead.
#ifdef HAVE_NX
#include <nx.h>
#define HAVE_WALL_TIME 1
#define HAVE_CPU_TIME 0
#else //HAVE_NX
#include <time.h>
#include <sys/types.h>
#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif
#ifdef HAVE_SYS_TIMES_H
#  include <sys/times.h>
#endif
#ifdef HAVE_SYS_RESOURCE_H
#  include <sys/resource.h>
#endif
#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif
#define HAVE_WALL_TIME 1
#define HAVE_CPU_TIME 1
#endif //HAVE_NX

#ifdef HAVE_PERF
#  define HAVE_FLOPS 1
#else
#  define HAVE_FLOPS 0
#endif

#if HAVE_FLOPS
extern "C" {
#  include <perf.h>
}
#endif

#ifdef MPQC_NEW_FEATURES
#include <chrono>
#endif

// AIX 3.2 has broken include files, likewise SunOS
#if defined(_AIX32) || defined(__sun)
extern "C" {
int getrusage (
  int Who,
  struct rusage *RUsage); }
#endif

#include <util/keyval/keyval.h>
#include <util/misc/regtime.h>
#define _util_misc_regtime_cc
#include <util/misc/timer.h>
#include <util/misc/scexception.h>
#ifdef MPQC_NEW_FEATURES
#include <util/misc/thread_timer.h>
#endif

using namespace std;
using namespace sc;

namespace sc {

//////////////////////////////////////////////////////////////////////

TimedRegion::TimedRegion(const char *name)
{
  name_ = strcpy(new char[strlen(name)+1], name);
  flops_ = wall_time_ = cpu_time_ = 0.0;
  up_ = 0;
  subregions_ = 0;
  next_ = prev_ = 0;
}

TimedRegion::~TimedRegion()
{
  delete[] name_;
  if (subregions_) while (subregions_->prev_) subregions_ = subregions_->prev_;
  delete subregions_;
  delete next_;
}

int
TimedRegion::nregion()
{
  int n = 1;
  if (subregions_) while (subregions_->prev_) subregions_ = subregions_->prev_;
  for (TimedRegion *i = subregions_; i!=0; i=i->next_) {
      n += i->nregion();
    }
  return n;
}

void
TimedRegion::get_region_names(const char *names[])
{
  names[0] = name();
  int n = 1;
  if (subregions_) while (subregions_->prev_) subregions_ = subregions_->prev_;
  for (TimedRegion *i = subregions_; i!=0; i=i->next_) {
      i->get_region_names(names + n);
      n += i->nregion();
    }
}

void
TimedRegion::get_depth(int *depth, int current_depth)
{
  depth[0] = current_depth;
  int n = 1;
  if (subregions_) while (subregions_->prev_) subregions_ = subregions_->prev_;
  for (TimedRegion *i = subregions_; i!=0; i=i->next_) {
      i->get_depth(depth + n, current_depth + 1);
      n += i->nregion();
    }
}

void
TimedRegion::get_wall_times(double *t)
{
  t[0] = wall_time_;
  int n = 1;
  if (subregions_) while (subregions_->prev_) subregions_ = subregions_->prev_;
  for (TimedRegion *i = subregions_; i!=0; i=i->next_) {
      i->get_wall_times(t + n);
      n += i->nregion();
    }
}

void
TimedRegion::get_cpu_times(double *t)
{
  t[0] = cpu_time_;
  int n = 1;
  if (subregions_) while (subregions_->prev_) subregions_ = subregions_->prev_;
  for (TimedRegion *i = subregions_; i!=0; i=i->next_) {
      i->get_cpu_times(t + n);
      n += i->nregion();
    }
}

void
TimedRegion::get_flops(double *t)
{
  t[0] = flops_;
  int n = 1;
  if (subregions_) while (subregions_->prev_) subregions_ = subregions_->prev_;
  for (TimedRegion *i = subregions_; i!=0; i=i->next_) {
      i->get_flops(t + n);
      n += i->nregion();
    }
}

TimedRegion *
TimedRegion::findinsubregion(const char *soughtname)
{
  if (!subregions_) {
      subregions_ = new TimedRegion(soughtname);
      subregions_->up_ = this;
      return subregions_;
    }
  int cmp = strcmp(subregions_->name_, soughtname);
  if (cmp < 0) {
      do {
          if (!subregions_->next_) {
              return subregions_->insert_after(soughtname);
            }
          subregions_ = subregions_->next_;
        } while ((cmp = strcmp(subregions_->name_, soughtname)) < 0);
      if (cmp == 0) return subregions_;
      subregions_ = subregions_->insert_before(soughtname);
    }
  else if (cmp > 0) {
      do {
          if (!subregions_->prev_) {
              return subregions_->insert_before(soughtname);
            }
          subregions_ = subregions_->prev_;
        } while ((cmp = strcmp(subregions_->name_, soughtname)) > 0);
      if (cmp == 0) return subregions_;
      subregions_ = subregions_->insert_after(soughtname);
    }
  return subregions_;
}

TimedRegion *
TimedRegion::insert_after(const char *name)
{
  TimedRegion *res = new TimedRegion(name);
  res->prev_ = this;
  res->next_ = this->next_;
  if (res->next_) res->next_->prev_ = res;
  res->up_ = up_;
  this->next_ = res;
  return res;
}

TimedRegion *
TimedRegion::insert_before(const char *name)
{
  TimedRegion *res = new TimedRegion(name);
  res->next_ = this;
  res->prev_ = this->prev_;
  if (res->prev_) res->prev_->next_ = res;
  res->up_ = up_;
  this->prev_ = res;
  return res;
}

void
TimedRegion::acquire_subregion(TimedRegion* reg)
{
  //TimedRegion* subreg = subregions_;
  //TimedRegion* new_reg = findinsubregion(reg->name_);
  //new_reg->merge(reg);
}

void
TimedRegion::cpu_enter(double t)
{
  cpu_enter_ = t;
}

void
TimedRegion::wall_enter(double t)
{
  wall_enter_ = t;
}

void
TimedRegion::flops_enter(double f)
{
  flops_enter_ = f;
}

void
TimedRegion::cpu_exit(double t)
{
  cpu_time_ += t - cpu_enter_;
  cpu_enter_ = t;
}

void
TimedRegion::wall_exit(double t)
{
  wall_time_ += t - wall_enter_;
  wall_enter_ = t;
}

void
TimedRegion::flops_exit(double f)
{
  flops_ += f - flops_enter_;
  flops_enter_ = f;
}

void
TimedRegion::merge(const TimedRegion* r)
{
  if (!r) return;

  const TimedRegion *start = r;
  while (start->prev_) start = start->prev_;
  for (const TimedRegion *riter = start;
       riter; riter = riter->next_) {
      TimedRegion *subr = findinsubregion(riter->name());
      subr->cpu_time_  += riter->cpu_time_;
      subr->wall_time_ += riter->wall_time_;
      subr->flops_     += riter->flops_;
      subr->merge(riter->subregions_);
    }
}

//////////////////////////////////////////////////////////////////////

static ClassDesc RegionTimer_cd(
    typeid(RegionTimer),"RegionTimer",1,"public DescribedClass");

RegionTimer::RegionTimer(const Ref<KeyVal> &keyval)
{
  KeyValValueboolean yes(1);
  KeyValValueboolean no(0);
  KeyValValuestring defname("total");

  wall_time_ = keyval->booleanvalue("wall_time",yes);
  cpu_time_ = keyval->booleanvalue("cpu_time",yes);
  flops_ = keyval->booleanvalue("flops",no);

#if !HAVE_CPU_TIME
  cpu_time_ = 0;
#endif
#if !HAVE_WALL_TIME
  wall_time_ = 0;
#endif
#if !HAVE_FLOPS
  flops_ = 0;
#endif

#if HAVE_FLOPS
  if (flops_) {
      if (perf_reset() || perf_set_config(0, PERF_FLOPS) || perf_start())
          flops_ = 0;
    }
#endif

  std::string topname = keyval->stringvalue("name", defname);
  top_ = new TimedRegion(topname.c_str());
  if (cpu_time_) top_->cpu_enter(get_cpu_time());
  if (wall_time_) top_->wall_enter(get_wall_time());
  if (flops_) top_->flops_enter(get_flops());
  current_ = top_;
}

RegionTimer::RegionTimer(const char *topname, int cpu_time, int wall_time):
  wall_time_(0),
  cpu_time_(0),
  flops_(0),
  default_(0)
{
#if HAVE_CPU_TIME
  cpu_time_ = cpu_time;
#endif
#if HAVE_WALL_TIME
  wall_time_ = wall_time;
#endif
  top_ = new TimedRegion(topname);
  if (cpu_time_) top_->cpu_enter(get_cpu_time());
  if (wall_time_) top_->wall_enter(get_wall_time());
  if (flops_) top_->flops_enter(get_flops());
  current_ = top_;
}

RegionTimer::~RegionTimer()
{
  delete top_;
}

void
RegionTimer::reset()
{
  std::string topname(top_->name());
  delete top_;
  top_ = new TimedRegion(topname.c_str());
  if (cpu_time_) top_->cpu_enter(get_cpu_time());
  if (wall_time_) top_->wall_enter(get_wall_time());
  if (flops_) top_->flops_enter(get_flops());
  current_ = top_;
}

double
RegionTimer::get_cpu_time()
{
#if defined(HAVE_NX)
  return 0.0;
#endif
  double res;
  struct rusage r;
  getrusage(RUSAGE_SELF,&r);
  res = r.ru_utime.tv_sec + r.ru_stime.tv_sec;
  res += 0.000001 * ( r.ru_utime.tv_usec + r.ru_stime.tv_usec );
  return res;
}

double
RegionTimer::get_wall_time()
{
#if defined(HAVE_NX)
  return dclock();
#endif
  struct timeval tod;
  gettimeofday(&tod,0);
  return tod.tv_sec + 0.000001 * tod.tv_usec;
}

double
RegionTimer::get_flops()
{
#if !HAVE_FLOPS
  return 0.0;
#else
  unsigned long long counter;
  perf_read(0,&counter);
  return (double)counter;
#endif
}

#ifdef MPQC_NEW_FEATURES
std::chrono::time_point<std::chrono::high_resolution_clock>
RegionTimer::get_time_point()
{
  return std::chrono::high_resolution_clock::now();
}
#endif

void
RegionTimer::enter(const char *name)
{
  current_ = current_->findinsubregion(name);
  if (cpu_time_) current_->cpu_enter(get_cpu_time());
  if (wall_time_) current_->wall_enter(get_wall_time());
  if (flops_) current_->flops_enter(get_flops());
}

void
RegionTimer::exit(const char *name, bool do_not_throw)
{
  if (!current_ || (name && strcmp(name, current_->name()))) {
      if (do_not_throw) {
          // we have an error but cannot throw.  ignore this call
          return;
        }
      else {
          throw ProgrammingError("region mismatch",
                                 __FILE__, __LINE__, this->class_desc());
        }
    }
  if (cpu_time_) current_->cpu_exit(get_cpu_time());
  if (wall_time_) current_->wall_exit(get_wall_time());
  if (flops_) current_->flops_exit(get_flops());
  if (! current_->up()) {
      if (do_not_throw) {
          // we have an error but cannot throw.  ignore this call
          return;
        }
      else {
          throw ProgrammingError("tried to exit top level",
                                 __FILE__, __LINE__, this->class_desc());
        }
    }
  current_ = current_->up();
}

double
RegionTimer::wall_time(const char *name) const
{
  double result;
  if (!current_ || !name) {
    throw ProgrammingError("no region name given, but default does not exists",
                           __FILE__, __LINE__, this->class_desc());
  }
  if (name && !strcmp(name, current_->name()))
    result = current_->wall_time();
  else
    result = current_->findinsubregion(name)->wall_time();

  return result;
}

double
RegionTimer::cpu_time(const char *name) const
{
  double result;
  if (!current_ || !name) {
    throw ProgrammingError("no region name given, but default does not exists",
                           __FILE__, __LINE__, this->class_desc());
  }
  if (name && !strcmp(name, current_->name()))
    result = current_->cpu_time();
  else
    result = current_->findinsubregion(name)->cpu_time();

  return result;
}

double
RegionTimer::flops(const char *name) const
{
  double result;
  if (!current_ || !name) {
    throw ProgrammingError("no region name given, but default does not exists",
                           __FILE__, __LINE__, this->class_desc());
  }
  if (name && !strcmp(name, current_->name()))
    result = current_->flops();
  else
    result = current_->findinsubregion(name)->flops();

  return result;
}

void
RegionTimer::add_wall_time(const char *name, double t)
{
  if (wall_time_) {
    current_ = current_->findinsubregion(name);
    current_->wall_add(t);
    current_ = current_->up();
    }
}

void
RegionTimer::add_cpu_time(const char *name, double t)
{
  if (cpu_time_) {
    current_ = current_->findinsubregion(name);
    current_->cpu_add(t);
    current_ = current_->up();
    }
}

void
RegionTimer::add_flops(const char *name, double t)
{
  if (flops_) {
    current_ = current_->findinsubregion(name);
    current_->flops_add(t);
    current_ = current_->up();
    }
}


void
RegionTimer::enter_default()
{
  if (cpu_time_) default_->cpu_enter(get_cpu_time());
  if (wall_time_) default_->wall_enter(get_wall_time());
  if (flops_) default_->flops_enter(get_flops());
}

void
RegionTimer::exit_default()
{
  if (cpu_time_) default_->cpu_exit(get_cpu_time());
  if (wall_time_) default_->wall_exit(get_wall_time());
  if (flops_) default_->flops_exit(get_flops());
}

void
RegionTimer::set_default(const char *name)
{
  if (default_) { defaults_.push_front(default_); }

  default_ = current_->findinsubregion(name);
}

void
RegionTimer::unset_default()
{
  if (defaults_.begin() != defaults_.end()) {
      default_ = defaults_.front();
      defaults_.pop_front();
    }
  else {
      default_ = 0;
    }
}

void
RegionTimer::change(const char *newname, const char *oldname)
{
  if (!current_ || (oldname && strcmp(oldname, current_->name()))) {
      ExEnv::errn() << "RegionTimer::change("
           << "\"" << newname << "\","
           << "\"" << oldname << "\""
           << "):"
           << " current region"
           << " (\"" << current_->name() << "\")"
           << " doesn't match name"
           << endl;
      abort();
    }
  double cpu=0.0, wall=0.0, flops=0.0;
  if (cpu_time_) current_->cpu_exit(cpu = get_cpu_time());
  if (wall_time_) current_->wall_exit(wall = get_wall_time());
  if (flops_) current_->flops_exit(flops = get_flops());
  if (! current_->up()) {
      ExEnv::errn() << "RegionTimer::change: already at top level" << endl;
      abort();
    }
  current_ = current_->up();
  current_ = current_->findinsubregion(newname);
  if (cpu_time_) current_->cpu_enter(cpu);
  if (wall_time_) current_->wall_enter(wall);
  if (flops_) current_->flops_enter(flops);
}

int
RegionTimer::nregion() const
{
  return top_->nregion();
}

void
RegionTimer::get_region_names(const char *region_names[]) const
{
  top_->get_region_names(region_names);
}

void
RegionTimer::get_cpu_times(double *cpu_time) const
{
  top_->get_cpu_times(cpu_time);
}

void
RegionTimer::get_wall_times(double *wall_time) const
{
  top_->get_wall_times(wall_time);
}

void
RegionTimer::get_flops(double *flops) const
{
  top_->get_flops(flops);
}

void
RegionTimer::get_depth(int *depth) const
{
  top_->get_depth(depth);
}

void
RegionTimer::update_top() const
{
  if (cpu_time_) top_->cpu_exit(get_cpu_time());
  if (wall_time_) top_->wall_exit(get_wall_time());
  if (flops_) top_->flops_exit(get_flops());
}

void
RegionTimer::print(ostream& o) const
{
  update_top();

  int n = nregion();
  double *cpu_time = 0;
  double *wall_time = 0;
  double *flops = 0;
  const char *flops_name = 0;
  if (cpu_time_) {
      cpu_time = new double[n];
      get_cpu_times(cpu_time);
    }
  if (wall_time_) {
      wall_time = new double[n];
      get_wall_times(wall_time);
    }
  if (flops_) {
      flops = new double[n];
      get_flops(flops);
      if (cpu_time_) {
        for (int i=0; i<n; i++) {
          if (fabs(cpu_time[i]) > 1.0e-10) flops[i] /= cpu_time[i]*1000000.;
          else flops[i] = 0.0;
          }
        flops_name = "MFLOP/S";
        }
      else if (wall_time_) {
        for (int i=0; i<n; i++) {
          if (fabs(wall_time[i]) > 1.0e-10) flops[i] /= wall_time[i]*1000000.;
          else flops[i] = 0.0;
          }
        flops_name = "MFLOP/WS";
        }
      else {
        for (int i=0; i<n; i++) {
          flops[i] /= 1000000.;
          }
        flops_name = "mflops";
        }
    }
  const char **names = new const char*[n];
  get_region_names(names);
  int *depth = new int[n];
  get_depth(depth);

  int i,j;
  int maxwidth = 0;
  double maxcputime = 0.0;
  double maxwalltime = 0.0;
  double maxflops = 0.0;
  for (i=0; i<n; i++) {
      int width = strlen(names[i]) + 2 * depth[i] + 2;
      if (width > maxwidth) maxwidth = width;
      if (cpu_time_ && cpu_time[i] > maxcputime) maxcputime = cpu_time[i];
      if (wall_time_ && wall_time[i] > maxwalltime) maxwalltime = wall_time[i];
      if (flops_ && flops[i] > maxflops) maxflops = flops[i];
    }

  size_t maxwallwidth = 4;
  while (maxwalltime >= 10.0) { maxwalltime/=10.0; maxwallwidth++; }

  size_t maxcpuwidth = 4;
  while (maxcputime >= 10.0) { maxcputime/=10.0; maxcpuwidth++; }

  size_t maxflopswidth = 4;
  if (flops_) {
    while (maxflops >= 10.0) { maxflops/=10.0; maxflopswidth++; }
    if (maxflopswidth < strlen(flops_name)) maxflopswidth = strlen(flops_name);
    }

  o.setf(ios::right);
  for (i=0; i<maxwidth; i++) o << " ";
  if (cpu_time_) o << " " << setw(maxcpuwidth) << "CPU";
  if (wall_time_) o << " " << setw(maxwallwidth) << "Wall";
  if (flops_) o << " " << setw(maxflopswidth) << flops_name;
  o << endl;

  o.setf(ios::fixed);
  streamsize oldprecision = o.precision();
  o.precision(2);
  for (i=0; i<n; i++) {
      int width = strlen(names[i]) + 2 * depth[i] + 2;
      for (j=0; j<depth[i]; j++) o << "  ";
      o << names[i] << ": ";
      for (j=width; j<maxwidth; j++) o << " ";
      if (cpu_time_) {
          o << " " << setw(maxcpuwidth) << cpu_time[i];
        }
      if (wall_time_) {
          o << " " << setw(maxwallwidth) << wall_time[i];
        }
      if (flops_) {
          o << " " << setw(maxflopswidth) << flops[i];
        }
      o << endl;
    }
  o.precision(oldprecision);
  o.unsetf(ios::fixed);

  delete[] cpu_time;
  delete[] wall_time;
  delete[] flops;
  delete[] names;
  delete[] depth;
}

void
RegionTimer::merge(const Ref<RegionTimer> &r)
{
  r->update_top();
  current_->merge(r->top_);
}

static Ref<RegionTimer> default_regtimer;

RegionTimer *
RegionTimer::default_regiontimer()
{
  return default_regtimer.pointer();
}

void
RegionTimer::set_default_regiontimer(const Ref<RegionTimer>& t)
{
  default_regtimer = t;
}

void
RegionTimer::acquire_timed_region(TimedRegion* reg)
{
  current_->merge(reg);

}

//////////////////////////////////////////////////////////////////////
// Timer functions

Timer::Timer(const char *name):
  timer_(),
  depth_(0),
  default_depth_(0),
  default_entered_(false)
{
  timer_ = RegionTimer::default_regiontimer();
  if (timer_ && name != 0) {
      enter(name);
    }
}

Timer::Timer(const Ref<RegionTimer>&t, const char *name):
  timer_(t),
  depth_(0),
  default_depth_(0),
  default_entered_(false)
{
  if (timer_ && name != 0) {
      enter(name);
    }
}

Timer::Timer(const std::string &name):
  timer_(),
  depth_(0),
  default_depth_(0),
  default_entered_(false)
{
  timer_ = RegionTimer::default_regiontimer();
  if (timer_) {
      enter(name);
    }
}

Timer::Timer(const Ref<RegionTimer>&t, const std::string &name):
      timer_(t),
      depth_(0),
      default_depth_(0),
      default_entered_(false)
{
  if (timer_) {
      enter(name);
    }
}

Timer::Timer(const Ref<RegionTimer>&t):
          timer_(t),
          depth_(0),
          default_depth_(0),
          default_entered_(false)
{
}

Timer::Timer():
          timer_(),
          depth_(0),
          default_depth_(0),
          default_entered_(false)
{
  timer_ = RegionTimer::default_regiontimer();
}

void
Timer::print(std::ostream &o) const
{
  if (timer_) timer_->print(o);
}

Timer::~Timer()
{
  // NB: This unwind code could be incorrect for complicated
  // uses of default timing regions.
  if (timer_) {
      while (depth_ > 0) {
          this->exit();
        }
      if (default_entered_) exit_default();
      while (default_depth_ > 0) {
          unset_default();
        }
    }
}

void
Timer::enter(const char *name)
{
  if (timer_) {
      timer_->enter(name);
      depth_++;
    }
}

void
Timer::change(const char *name)
{
  if (timer_) {
      timer_->change(name);
    }
}

void
Timer::exit(const char *name)
{
  if (timer_) {
      timer_->exit(name,true);
      depth_--;
    }
}

void
Timer::set_default(const char *name)
{
  if (timer_) {
      timer_->set_default(name);
      default_depth_++;
    }
}

void
Timer::unset_default()
{
  if (timer_) {
      timer_->unset_default();
      default_depth_--;
    }
}

void
Timer::enter_default()
{
  if (timer_) {
      timer_->enter_default();
      default_entered_ = true;
    }
}

void
Timer::exit_default()
{
  if (timer_) {
      timer_->exit_default();
      default_entered_ = false;
    }
}

void
Timer::reset(const char *name)
{
  if (depth_ > 0) {
      if (name) {
          change(name);
        }
      else {
          this->exit();
        }
    }
  else {
      if (name) {
          enter(name);
        }
    }
}

double
Timer::cpu_time(const char* region) const {
  return timer_->cpu_time(region);
}

double Timer::wall_time(const char* region) const {
  return timer_->wall_time(region);
}

double Timer::flops(const char* region) const {
  return timer_->flops(region);
}

#ifdef MPQC_NEW_FEATURES
void Timer::insert(const MultiThreadTimer& timer) {
  TimedRegion* reg = timer.make_timed_region();
  timer_->acquire_timed_region(reg);
  delete reg;
  return;
}
#endif

//////////////////////////////////////////////////////////////////////
// Shorthand to manipulate the global region timer

void
tim_enter(const char *name) {
  if (default_regtimer) default_regtimer->enter(name);
}

void
tim_exit(const char *name)
{
  if (default_regtimer) default_regtimer->exit(name);
}

void
tim_set_default(const char *name)
{
  if (default_regtimer) default_regtimer->set_default(name);
}

void
tim_enter_default()
{
  if (default_regtimer) default_regtimer->enter_default();
}

void
tim_exit_default()
{
  if (default_regtimer) default_regtimer->exit_default();
}

void
tim_change(const char *name)
{
  if (default_regtimer) default_regtimer->change(name);
}

void
tim_print(int)
{
  if (default_regtimer) default_regtimer->print();
}

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
