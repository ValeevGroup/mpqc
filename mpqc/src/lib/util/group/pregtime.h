
#ifdef __GNUC__
#pragma interface
#endif

#ifndef _util_group_pregtime_h
#define _util_group_pregtime_h

#include <iostream.h>
#include <util/misc/regtime.h>
#include <util/group/message.h>

class ParallelRegionTimer: public RegionTimer {
  protected:
    RefMessageGrp msg_;
  public:
    ParallelRegionTimer(const RefMessageGrp&,
                        const char *topname = "total",
                        int cpu_time = 0, int wall_time = 1);
    ~ParallelRegionTimer();

    void print(ostream& = cout);
};

#endif
