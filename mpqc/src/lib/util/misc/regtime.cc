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

#ifdef __GNUC__
#pragma implementation
#endif

#ifdef HAVE_CONFIG_H
#  include <scconfig.h>
#endif

#include <iostream.h>
#include <iomanip.h>

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

// AIX 3.2 has broken include files, likewise SunOS
#if defined(_AIX32) || defined(__sun)
extern "C" {
int getrusage (
  int Who,
  struct rusage *RUsage); }
#endif

#include <util/misc/regtime.h>
#include <util/misc/timer.h>

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

    TimedRegion *insert_after(const char *name);
    TimedRegion *insert_before(const char *name);
  public:
    TimedRegion(const char *name);
    ~TimedRegion();
    const char *name() const { return name_; }
    TimedRegion *findinsubregion(const char *);
    void cpu_enter(double);
    void wall_enter(double);
    void cpu_exit(double);
    void wall_exit(double);
    TimedRegion *up() const { return up_; }

    int nregion();
    void get_region_names(const char *names[]);
    void get_wall_times(double *);
    void get_cpu_times(double *);
    void get_depth(int *, int depth = 0);
};

//////////////////////////////////////////////////////////////////////

TimedRegion::TimedRegion(const char *name)
{
  name_ = strcpy(new char[strlen(name)+1], name);
  wall_time_ = cpu_time_ = 0.0;
  up_ = 0;
  subregions_ = 0;
  next_ = prev_ = 0;
}

TimedRegion::~TimedRegion()
{
  delete name_;
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

//////////////////////////////////////////////////////////////////////

RegionTimer::RegionTimer(const char *topname, int cpu_time, int wall_time):
  cpu_time_(0),
  wall_time_(0),
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
  current_ = top_;
}

RegionTimer::~RegionTimer()
{
  delete top_;
}

double
RegionTimer::get_cpu_time()
{
#ifdef HAVE_NX
  return 0.0;
#else
  double res;
  struct rusage r;
  getrusage(RUSAGE_SELF,&r);
  res = r.ru_utime.tv_sec + r.ru_stime.tv_sec;
  res += 0.000001 * ( r.ru_utime.tv_usec + r.ru_stime.tv_usec );
  return res;
#endif
}

double
RegionTimer::get_wall_time()
{
#ifdef HAVE_NX
  return dclock();
#else
  struct timeval tod;
  gettimeofday(&tod,0);
  return tod.tv_sec + 0.000001 * tod.tv_usec;
#endif
}

void
RegionTimer::enter(const char *name)
{
  current_ = current_->findinsubregion(name);
  if (cpu_time_) current_->cpu_enter(get_cpu_time());
  if (wall_time_) current_->wall_enter(get_wall_time());
}

void
RegionTimer::exit(const char *name)
{
  if (!current_ || (name && strcmp(name, current_->name()))) {
      cerr << "TimeRegion::exit(\"" << name << "\"):"
           << " current region"
           << " (\"" << current_->name() << "\")"
           << " doesn't match name"
           << endl;
      abort();
    }
  if (cpu_time_) current_->cpu_exit(get_cpu_time());
  if (wall_time_) current_->wall_exit(get_wall_time());
  if (! current_->up()) {
      cerr << "RegionTimer::exit: already at top level" << endl;
      abort();
    }
  current_ = current_->up();
}

void
RegionTimer::enter_default()
{
  if (cpu_time_) default_->cpu_enter(get_cpu_time());
  if (wall_time_) default_->wall_enter(get_wall_time());
}

void
RegionTimer::exit_default()
{
  if (cpu_time_) default_->cpu_exit(get_cpu_time());
  if (wall_time_) default_->wall_exit(get_wall_time());
}

void
RegionTimer::set_default(const char *name)
{
  default_ = current_->findinsubregion(name);
}

void
RegionTimer::unset_default()
{
  default_ = 0;
}

void
RegionTimer::change(const char *newname, const char *oldname)
{
  if (!current_ || (oldname && strcmp(oldname, current_->name()))) {
      cerr << "RegionTimer::change("
           << "\"" << newname << "\","
           << "\"" << oldname << "\""
           << "):"
           << " current region"
           << " (\"" << current_->name() << "\")"
           << " doesn't match name"
           << endl;
      abort();
    }
  double cpu, wall;
  if (cpu_time_) current_->cpu_exit(cpu = get_cpu_time());
  if (wall_time_) current_->wall_exit(wall = get_wall_time());
  if (! current_->up()) {
      cerr << "RegionTimer::change: already at top level" << endl;
      abort();
    }
  current_ = current_->up();
  current_ = current_->findinsubregion(newname);
  if (cpu_time_) current_->cpu_enter(cpu);
  if (wall_time_) current_->wall_enter(wall);
}

int
RegionTimer::nregion()
{
  return top_->nregion();
}

void
RegionTimer::get_region_names(const char *region_names[])
{
  top_->get_region_names(region_names);
}

void
RegionTimer::get_cpu_times(double *cpu_time)
{
  top_->get_cpu_times(cpu_time);
}

void
RegionTimer::get_wall_times(double *wall_time)
{
  top_->get_wall_times(wall_time);
}

void
RegionTimer::get_depth(int *depth)
{
  top_->get_depth(depth);
}

void
RegionTimer::update_top()
{
  if (cpu_time_) top_->cpu_exit(get_cpu_time());
  if (wall_time_) top_->wall_exit(get_wall_time());
}

void
RegionTimer::print(ostream& o)
{
  update_top();

  int n = nregion();
  double *cpu_time = 0;
  double *wall_time = 0;
  if (cpu_time_) {
      cpu_time = new double[n];
      get_cpu_times(cpu_time);
    }
  if (wall_time_) {
      wall_time = new double[n];
      get_wall_times(wall_time);
    }
  const char **names = new const char*[n];
  get_region_names(names);
  int *depth = new int[n];
  get_depth(depth);

  int i,j;
  int maxwidth = 0;
  double maxcputime = 0.0;
  double maxwalltime = 0.0;
  for (i=0; i<n; i++) {
      int width = strlen(names[i]) + 2 * depth[i] + 2;
      if (width > maxwidth) maxwidth = width;
      if (cpu_time_ && cpu_time[i] > maxcputime) maxcputime = cpu_time[i];
      if (wall_time_ && wall_time[i] > maxwalltime) maxwalltime = wall_time[i];
    }

  int maxwallwidth = 4;
  while (maxwalltime > 1.0) { maxwalltime/=10.0; maxwallwidth++; }

  int maxcpuwidth = 4;
  while (maxcputime > 1.0) { maxcputime/=10.0; maxcpuwidth++; }

  o.setf(ios::right);
  for (i=0; i<maxwidth; i++) o << " ";
  if (cpu_time_) o << setw(maxcpuwidth+1) << " CPU";
  if (wall_time_) o << setw(maxwallwidth+1) << " Wall";
  o << endl;

  o.setf(ios::fixed);
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
      o << endl;
    }

  delete[] cpu_time;
  delete[] wall_time;
  delete[] names;
  delete[] depth;
}

static RefRegionTimer default_regtimer;

RegionTimer *
RegionTimer::default_regiontimer()
{
  return default_regtimer.pointer();
}

void
RegionTimer::set_default_regiontimer(const RefRegionTimer& t)
{
  default_regtimer = t;
}

//////////////////////////////////////////////////////////////////////
// C compatibility routines

void
tim_enter(const char *name) {
  if (default_regtimer.nonnull()) default_regtimer->enter(name);
}

void
tim_exit(const char *name)
{
  if (default_regtimer.nonnull()) default_regtimer->exit(name);
}

void
tim_set_default(const char *name)
{
  if (default_regtimer.nonnull()) default_regtimer->set_default(name);
}

void
tim_enter_default()
{
  if (default_regtimer.nonnull()) default_regtimer->enter_default();
}

void
tim_exit_default()
{
  if (default_regtimer.nonnull()) default_regtimer->exit_default();
}

void
tim_change(const char *name)
{
  if (default_regtimer.nonnull()) default_regtimer->change(name);
}

void
tim_print(int)
{
  if (default_regtimer.nonnull()) default_regtimer->print();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
