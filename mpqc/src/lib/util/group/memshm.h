//
// memshm.h
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

#ifndef _util_group_memshm_h
#define _util_group_memshm_h

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <iostream>

#include <util/group/globcnt.h>
#include <util/group/memmsg.h>


/** The ShmMemoryGrp concrete class provides an implementation of
MsgMemoryGrp.  It uses SYSV IPC to provided shared memory in a system that
provide shared memory in hardware.  It is very fast and reliable. */
class ShmMemoryGrp: public MsgMemoryGrp {
#define CLASSNAME ShmMemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    int nregion_;
    int *shmid_;
    void **attach_address_;
    GlobalCounter lock_;
    GlobalCounter *update_;
    void *data_;
    void *memory_;
    Pool *pool_;
    RangeLock *rangelock_; // the locks_ member of the base class is ignored

    void clear_release_count();
    void wait_for_release();
    void note_release();
    void obtain_lock();
    void release_lock();

    void cleanup();
    int attach_memory(void *ataddress, int size);
    void detach_memory();
  public:
    ShmMemoryGrp(const RefMessageGrp& msg);
    ShmMemoryGrp(const RefKeyVal&);
    ~ShmMemoryGrp();

    void set_localsize(int);
    void *localdata();

    void *obtain_readwrite(distsize_t offset, int size);
    void *obtain_readonly(distsize_t offset, int size);
    void release_read(void *data, distsize_t offset, int size);
    void release_write(void *data, distsize_t offset, int size);

    virtual void sum_reduction(double *data, distsize_t doffset, int dsize);

    void print(std::ostream &o = ExEnv::out()) const;
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
