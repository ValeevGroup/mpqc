//
// int2e.cc
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

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/utils.h>

using namespace sc;

Int2eV3::Int2eV3(Integral *integral,
                 const Ref<GaussianBasisSet>& b1,
                 const Ref<GaussianBasisSet>& b2,
                 const Ref<GaussianBasisSet>& b3,
                 const Ref<GaussianBasisSet>& b4,
                 int order, size_t storage) :
  integral_(integral),
  grp_(integral->messagegrp()),
  store(0),
  int_Qvec(0),
  int_Rvec(0)
{
  bs1_ = b1;
  bs2_ = b2;
  bs3_ = b3;
  bs4_ = b4;

  int_unit2 = bs2_ == 0;
  int_unit4 = bs4_ == 0;

  transform_init();
  int_initialize_offsets2();
  int_initialize_erep(storage,order,bs1_,bs2_,bs3_,bs4_);

  // bounds in IntV3 can only be computed under certain circumstances
  const bool compute_bounds = !(int_unit2 || int_unit4) &&
                              ((bs1_ == bs2_)&&(bs1_ == bs3_)&&(bs1_ == bs4_));
  if (compute_bounds) {
    if (order==0) {
      init_bounds();
      }
    else if (order==1) {
      init_bounds_1der();
    }
  }
}

Int2eV3::~Int2eV3()
{
  transform_done();
  int_done_offsets2();
  int_done_erep();
  if (int_integral_storage) {
    done_storage();
  }
  done_bounds();
  done_bounds_1der();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
