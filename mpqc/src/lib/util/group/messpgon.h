//
// messpgon.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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

#ifndef _util_group_messpgon_h
#define _util_group_messpgon_h

#include <util/group/message.h>

class ParagonMessageGrp: public intMessageGrp {
#define CLASSNAME ParagonMessageGrp
#define HAVE_KEYVAL_CTOR
#include <util/class/classd.h>
  protected:
    void basic_send(int target, int type, void* data, int nbyte);
    void basic_recv(int type, void* data, int nbyte);
    int basic_probe(int type);
    void initialize();
  public:
    ParagonMessageGrp();
    ParagonMessageGrp(const RefKeyVal&);
    ~ParagonMessageGrp();
    void sync();
 
    int last_source();
    int last_size();
    int last_type();

    void reduce(double*, int n, GrpReduce<double>&,
                double*scratch = 0, int target = -1);
    void reduce(int*, int n, GrpReduce<int>&,
                int*scratch = 0, int target = -1);
    void reduce(char*, int n, GrpReduce<char>&,
                char*scratch = 0, int target = -1);
    void reduce(unsigned char*, int n, GrpReduce<unsigned char>&,
                unsigned char*scratch = 0, int target = -1);
    void reduce(short*, int n, GrpReduce<short>&,
                short*scratch = 0, int target = -1);
    void reduce(float*, int n, GrpReduce<float>&,
                float*scratch = 0, int target = -1);
    void reduce(long*, int n, GrpReduce<long>&,
                long*scratch = 0, int target = -1);

    void raw_bcast(void* data, int nbyte, int from);

    void raw_collect(void *whole, const int *lengths,
                     const void *part, int bytes_per_datum=1);
};

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
