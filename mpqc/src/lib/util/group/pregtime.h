//
// pregtime.h
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

#ifndef _util_group_pregtime_h
#define _util_group_pregtime_h

#include <iostream>
#include <util/misc/regtime.h>
#include <util/group/message.h>

namespace sc {

/// This is a parallel-away derivative of RegionTimer. It will compute
/// region timings by summing them across nodes. Timings will be printed by
/// node zero only.
class ParallelRegionTimer: public RegionTimer {
  protected:
    Ref<MessageGrp> msg_;

    void send_subregions(int node, const TimedRegion *r) const;
    void recv_subregions(int node, TimedRegion *r) const;
    void all_reduce_regions() const;

  public:
    ParallelRegionTimer(const Ref<MessageGrp>&,
                        const char *topname = "total",
                        int cpu_time = 0, int wall_time = 1);
    ParallelRegionTimer(const Ref<KeyVal>&);
    ~ParallelRegionTimer();

    void print(std::ostream& = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
