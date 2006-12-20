//
// memarmci.h
// based on memshm.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: SNL
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

#ifndef _util_group_memarmci_h
#define _util_group_memarmci_h

#include <iostream>

#include <util/group/memrdma.h>

namespace sc {

/** The ARMCIMemoryGrp concrete class provides an implementation of
    MsgMemoryGrp.  It uses the ARMCI interface. */
class ARMCIMemoryGrp: public RDMAMemoryGrp {
  private:
    void **all_data_;
    void init();
    void finalize();
    Ref<ThreadLock> armci_lock_;
  public:
    ARMCIMemoryGrp(const Ref<MessageGrp>& msg);
    ARMCIMemoryGrp(const Ref<KeyVal>&);
    ~ARMCIMemoryGrp();

    void set_localsize(size_t);

    void retrieve_data(void *, int node, long offset, long size, int lock);
    void replace_data(void *, int node, long offset, long size, int unlock);
    void sum_data(double *data, int node, long doffset, long dsize);

    void sync();
    void deactivate();

    void* malloc_local(size_t nbyte);
    void free_local(void *data);

    void print(std::ostream &o = ExEnv::out0()) const;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
