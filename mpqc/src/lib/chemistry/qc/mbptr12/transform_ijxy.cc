//
// transform_ijxy.cc
//
// Copyright (C) 2004 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <stdexcept>

#include <util/misc/formio.h>
#include <util/state/state_bin.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/transform_ijxy.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------
  TwoBodyMOIntsTransform_ijxy
 -----------*/
static ClassDesc TwoBodyMOIntsTransform_ijxy_cd(
  typeid(TwoBodyMOIntsTransform_ijxy),"TwoBodyMOIntsTransform_ijxy",1,"public TwoBodyMOIntsTransform",
  0, 0, create<TwoBodyMOIntsTransform_ijxy>);

TwoBodyMOIntsTransform_ijxy::TwoBodyMOIntsTransform_ijxy(const Ref<MOIntsTransformFactory>& factory,
                                                         const Ref<MOIndexSpace>& space1, const Ref<MOIndexSpace>& space2,
                                                         const Ref<MOIndexSpace>& space3, const Ref<MOIndexSpace>& space4) :
  TwoBodyMOIntsTransform(factory,space1,space2,space3,space4)
{
}

TwoBodyMOIntsTransform_ijxy::TwoBodyMOIntsTransform_ijxy(StateIn& si) : TwoBodyMOIntsTransform(si)
{
}

TwoBodyMOIntsTransform_ijxy::~TwoBodyMOIntsTransform_ijxy()
{
}

void
TwoBodyMOIntsTransform_ijxy::save_data_state(StateOut& so)
{
  TwoBodyMOIntsTransform::save_data_state(so);
}

//////////////////////////////////////////////////////
// Compute required (dynamic) memory
// for a given batch size of the transformation
//
// Only arrays allocated before exiting the loop over
// i-batches are included here  - only these arrays
// affect the batch size.
//////////////////////////////////////////////////////
distsize_t
TwoBodyMOIntsTransform_ijxy::compute_transform_dynamic_memory_(int ni) const
{
  int nproc = msg_->n();
  int nthread = thr_->nthread();

  int rank2 = space2_->rank();
  int nbasis2 = space2_->basis()->nbasis();
  int nfuncmax3 = space3_->basis()->max_nfunction_in_shell();
  int nfuncmax4 = space4_->basis()->max_nfunction_in_shell();
  int rank3 = space3_->rank();
  int nbasis4 = space4_->basis()->nbasis();

  // compute nij as nij on node 0, since nij on node 0 is >= nij on other nodes
  int index = 0;
  int nij = 0;
  for (int i=0; i<ni; i++) {
    for (int j=0; j<rank2; j++) {
      if (index++ % nproc == 0) nij++;
    }
  }

  distsize_t memsize = sizeof(double)*(num_te_types_*((distsize_t)nthread * ni * nbasis2 * nfuncmax3 * nfuncmax4 // iqrs
						     + (distsize_t)ni * rank2 * nfuncmax3 * nfuncmax4  // ijrs
						     + (distsize_t)nij * rank3 * nbasis4 // ijxs - buffer of 3 q.t. and higher
						     // transformed integrals
						     )
				       + (distsize_t)rank3 * nbasis4 // xs or xy
				       );

  return memsize;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
