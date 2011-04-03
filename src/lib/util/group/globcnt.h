//
// globcnt.h
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

#ifndef _util_group_globcnt_h
#define _util_group_globcnt_h

namespace sc {

/// The GlobalCounter class allows processes on the same SMP
/// node to share a counter using SysV IPC semaphores.
/// A process can create a GlobalCounter using the void CTOR.
/// This process can share the string representation of the
/// counter with other processes.  They can then use the const
/// char * CTOR to create global counters that reference the
/// same global counter.
class GlobalCounter {
  private:
    int semid_;
    int controls_release_;

    void cleanup();
    
  public:
    GlobalCounter();
    void initialize();
    void initialize(const char *stringrep);
    ~GlobalCounter();

    char *stringrep();

    void wait_for_zero();
    void operator += (int);
    void operator ++();
    void operator --();
    void operator ++(int) { operator++(); }
    void operator --(int) { operator--(); }
    void operator = (int);

    int val();
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
