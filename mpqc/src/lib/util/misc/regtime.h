
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_misc_regtime_h
#define _util_misc_regtime_h

#include <iostream.h>
#include <util/class/class.h>

class TimedRegion;

class RefRegionTimer;

class RegionTimer: public VRefCount {
  protected:
    int wall_time_;
    int cpu_time_;

    TimedRegion *top_;
    TimedRegion *current_;

  public:
    RegionTimer(const char *topname = "total",
                int cpu_time = 0, int wall_time = 1);
    ~RegionTimer();
    void enter(const char *);
    void change(const char *newname, const char * oldname = 0);
    void exit(const char * = 0);
    virtual void print(ostream& = cout);

    void update_top();

    int nregion();
    void get_region_names(const char *names[]);
    void get_wall_times(double *);
    void get_cpu_times(double *);
    void get_depth(int *);

    double get_wall_time();
    double get_cpu_time();

    static RegionTimer *default_regiontimer();
    static void set_default_regiontimer(const RefRegionTimer &);
};
REF_dec(RegionTimer);

#endif
