//
// memproc.h
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

#ifndef _util_group_memproc_h
#define _util_group_memproc_h

#include <sys/types.h>

#include <util/group/memmsg.h>

class ProcMemoryGrp: public MemoryGrp {
#define CLASSNAME ProcMemoryGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  private:
    char *data_;
  public:
    ProcMemoryGrp();
    ProcMemoryGrp(const RefKeyVal&);
    ~ProcMemoryGrp();

    void set_localsize(int);

    void *obtain_readwrite(int offset, int size);
    void *obtain_readonly(int offset, int size);
    void release_read(void *data, int offset, int size);
    void release_write(void *data, int offset, int size);

    void sync();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
