//
// messshm.h
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

#ifndef _util_group_messshm_h
#define _util_group_messshm_h

#include <util/group/message.h>

/** The ShmMessageGrp class is an implementation of MessageGrp that allows
multiple process to be started that communicate with shared memory.  This
only provides improved performance if you have multiple CPU's in a
symmetric multiprocessor configuration.  Nonetheless, it is quite useful on
a single CPU for tracking down bugs.

The ShmMessageGrp KeyVal constructor takes a single keyword that specifies
the number of processors.  Here is an example of a ParsedKeyVal input that
creates a ShmMessageGrp that runs on four processors.

\begin{alltt}
message<ShmMessageGrp>: n = 4
\end{alltt}

*/
class ShmMessageGrp: public intMessageGrp {
#define CLASSNAME ShmMessageGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  protected:
    void basic_send(int target, int type, void* data, int nbyte);
    void basic_recv(int type, void* data, int nbyte);
    int basic_probe(int type);
    void initialize(int nprocs);
    void initialize();

    // Information about the last message received or probed.
    int last_type_;
    int last_source_;
    int last_size_; // the size in bytes

    void set_last_type(int a) { last_type_ = a; }
    void set_last_source(int a) { last_source_ = a; }
    void set_last_size(int a) { last_size_ = a; }
  public:
    ShmMessageGrp(); // read nprocs from environmental variable NUMPROC
    ShmMessageGrp(const RefKeyVal&);
    ShmMessageGrp(int nprocs);
    ~ShmMessageGrp();
    void sync();
 
    int last_source();
    int last_size();
    int last_type();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
