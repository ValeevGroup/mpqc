//
// thpuma.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#ifndef _util_group_thpuma_h
#define _util_group_thpuma_h

#include <pthread.h>
#include <util/group/thread.h>

namespace sc {
  
/** The PumaThreadGrp class privides a concrete thread group
    appropriate for the intel teraflops machine.
*/
class PumaThreadGrp: public ThreadGrp {
  private:
    volatile unsigned int flag_;
    
  public:
    PumaThreadGrp();
    PumaThreadGrp(const Ref<KeyVal>&);
    ~PumaThreadGrp();

    int start_threads();
    int wait_threads();
    Ref<ThreadLock> new_lock();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
