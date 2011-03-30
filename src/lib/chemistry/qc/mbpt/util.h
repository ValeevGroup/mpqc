//
// util.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifndef _chemistry_qc_mbpt_util_h
#define _chemistry_qc_mbpt_util_h

#include <util/group/message.h>

namespace sc {

class BiggestContribs {
  private:
    int nindex_;
    int **indices_;
    double *vals_;
    int ncontrib_;
    int maxcontrib_;
  public:
    BiggestContribs(int nindex, int maxcontrib);
    ~BiggestContribs();
    double val(int i) { return vals_[i]; }
    int ncontrib() { return ncontrib_; }
    const int *indices(int i) { return indices_[i]; }
    void insert(double val, const int *);
    void insert(double val, int i0, int i1);
    void insert(double val, int i0, int i1, int i2, int i3);
    void insert(double val, int i0, int i1, int i2, int i3, int i4);
    void combine(const Ref<MessageGrp> &grp);
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
