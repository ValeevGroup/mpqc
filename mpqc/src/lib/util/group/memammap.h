//
// memammap.h --- definition of memory group which uses DEC mmap
//
// Copyright (C) 1998 Limit Point Systems, Inc.
//
// Author: Edward T. Seidl <seidl@janed.com>
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

#ifndef _util_group_memammap_h
#define _util_group_memammap_h

#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <iostream.h>

#include <util/group/globcnt.h>
#include <util/group/memmsg.h>

//. The \clsnm{AlphaMMapMessageGrp} class is an implementation of
//. \clsnmref{MessageGrp} that allows multiple process to be
//. started that communication with shared memory.
class AlphaMMapMemoryGrp: public MsgMemoryGrp {
#define CLASSNAME AlphaMMapMemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
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
  public:
    AlphaMMapMemoryGrp(const RefMessageGrp& msg);
    AlphaMMapMemoryGrp(const RefKeyVal&);
    ~AlphaMMapMemoryGrp();

    void set_localsize(int);

    void *obtain_readwrite(int offset, int size);
    void *obtain_readonly(int offset, int size);
    void release_read(void *data, int offset, int size);
    void release_write(void *data, int offset, int size);

    virtual void sum_reduction(double *data, int doffset, int dsize);

    void print(ostream &o = cout);
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
