//
// r12ia.cc
//
// Copyright (C) 2002 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
// Maintainer: EV
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
#pragma implementation
#endif

#include <stdlib.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/mbptr12/r12ia.h>

using namespace std;
using namespace sc;

/*--------------------------------
  R12IntsAcc
 --------------------------------*/

R12IntsAcc::R12IntsAcc(int num_te_types, int nbasis1, int nbasis2, int nocc, int nfzc) :
  num_te_types_(num_te_types), nbasis1_(nbasis1), nbasis2_(nbasis2), nocc_(nocc), nfzc_(nfzc)
{
  nocc_act_ = nocc_ - nfzc_;
  nbasis__2_ = nbasis1_*nbasis2_;
  blksize_ = nbasis__2_*sizeof(double);
  blocksize_ = blksize_*num_te_types_;
  committed_ = false;
}

R12IntsAcc::~R12IntsAcc()
{
}


///////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
