//
// int1e.cc
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

#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/utils.h>

using namespace sc;

Int1eV3::Int1eV3(Integral *integral,
                 const Ref<GaussianBasisSet>&b1,
                 const Ref<GaussianBasisSet>&b2,
                 int order)
{
  integral_ = integral;

  exponent_weighted = -1;
  scale_shell_result = 0;
  result_scale_factor = 1.0;
  three_center = 0;
  init_order = -1;
  buff = 0;
  cartesianbuffer = 0;
  cartesianbuffer_scratch = 0;

  bs1_ = b1;
  bs2_ = b2;

  transform_init();
  int_initialize_offsets1();
  int_initialize_1e(0,order);
}

Int1eV3::~Int1eV3()
{
  transform_done();
  int_done_1e();
  int_done_offsets1();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
