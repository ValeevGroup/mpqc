//
// transform_iRjS.cc
//
// Copyright (C) 2008 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
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

#include <chemistry/qc/lcao/transform_iRjS.h>
#include <chemistry/qc/lcao/transform_13inds.h>
#include <math/distarray4/distarray4_memgrp.h>
#include <math/distarray4/distarray4_node0file.h>
#ifdef HAVE_MPIIO
#  include <math/distarray4/distarray4_mpiiofile.h>
#endif
#include <util/group/memory.h>
#include <util/group/memregion.h>


using namespace std;
using namespace sc;

/*-----------
  TwoBodyMOIntsTransform_iRjS
 -----------*/
static ClassDesc TwoBodyMOIntsTransform_iRjS_cd(
  typeid(TwoBodyMOIntsTransform_iRjS),"TwoBodyMOIntsTransform_iRjS",1,"public TwoBodyMOIntsTransform",
  0, 0, create<TwoBodyMOIntsTransform_iRjS>);

TwoBodyMOIntsTransform_iRjS::TwoBodyMOIntsTransform_iRjS(const std::string& name, const Ref<MOIntsTransformFactory>& factory,
                                                         const Ref<TwoBodyIntDescr>& tbintdescr,
                                                         const Ref<OrbitalSpace>& space1,
                                                         const Ref<OrbitalSpace>& space2,
                                                         const Ref<OrbitalSpace>& space3,
                                                         const Ref<OrbitalSpace>& space4) :
  TwoBodyMOIntsTransform(name,factory,tbintdescr,space1,space2,space3,space4)
{
  init_vars();
}

TwoBodyMOIntsTransform_iRjS::TwoBodyMOIntsTransform_iRjS(StateIn& si) : TwoBodyMOIntsTransform(si)
{
  init_vars();
}

TwoBodyMOIntsTransform_iRjS::~TwoBodyMOIntsTransform_iRjS()
{
}

void
TwoBodyMOIntsTransform_iRjS::save_data_state(StateOut& so)
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
TwoBodyMOIntsTransform_iRjS::compute_transform_dynamic_memory_(int ni) const
{
  // dynamic contribution from 1+2 QT
  const int nthread = thr_->nthread();
  const distsize_t memsize12 = (distsize_t) nthread *
                         TwoBodyMOIntsTransform_13Inds::compute_required_dynamic_memory(*this,ni);

  // integrals held by MemoryGrp
  // compute nij as nij on node 0, since nij on node 0 is >= nij on other nodes
  const int nproc = msg()->n();
  const int rank3 = space3()->rank();
  int nij = compute_nij(ni, rank3, nproc, 0);
  const distsize_t memsize_memgrp = num_te_types() * nij * (distsize_t) memgrp_blksize();

  // determine the peak memory requirements
  return memsize_memgrp + memsize12;
}

size_t
TwoBodyMOIntsTransform_iRjS::memgrp_blksize() const
{
  const int nbasis2 = space2()->basis()->nbasis();
  const int nbasis4 = space4()->basis()->nbasis();
  return nbasis2*nbasis4*sizeof(double);
}

void
TwoBodyMOIntsTransform_iRjS::init_acc()
{
  if (ints_acc_.nonnull())
    return;

  const int nij = compute_nij(batchsize_, space4()->rank(), msg_->n(), msg_->me());
  const size_t localmem = num_te_types() * nij * memgrp_blksize();

  //
  // NOTE: results come out stored as YX for the sake maximizing the size of messages
  //       in TwoBodyMOIntsTransform_13Inds !!!
  //
  switch (ints_method_) {

  case MOIntsTransform::StoreMethod::mem_only:
    if (npass_ > 1)
      throw std::runtime_error("TwoBodyMOIntsTransform_iRjS::init_acc() -- cannot use MemoryGrp-based accumulator in multi-pass transformations");
    {
      // use a subset of a MemoryGrp provided by TransformFactory
      set_memgrp(new MemoryGrpRegion(mem(),localmem));
      ints_acc_ = new DistArray4_MemoryGrp(mem(), num_te_types(),
                                           space1()->rank(), space3()->rank(),
                                           space2()->rank(), space4()->rank(),
                                           memgrp_blksize(),
                                           DistArray4Storage_YX);
    }
    break;

  case MOIntsTransform::StoreMethod::mem_posix:
    // if can do in one pass, use the factory hints about how data will be used
    if (npass_ == 1 && !factory()->hints().data_persistent()) {
      // use a subset of a MemoryGrp provided by TransformFactory
      set_memgrp(new MemoryGrpRegion(mem(),localmem));
      ints_acc_ = new DistArray4_MemoryGrp(mem(), num_te_types(),
                                           space1()->rank(), space3()->rank(),
                                           space2()->rank(), space4()->rank(),
                                           memgrp_blksize(),
                                           DistArray4Storage_YX);
      break;
    }
    // else use the next case

  case MOIntsTransform::StoreMethod::posix:
    ints_acc_ = new DistArray4_Node0File((file_prefix_+"."+name_).c_str(), num_te_types(),
                                         space1()->rank(), space3()->rank(),
                                         space2()->rank(), space4()->rank(),
                                         DistArray4Storage_YX);
    break;

#ifdef HAVE_MPIIO
  case MOIntsTransform::StoreMethod::mem_mpi:
    // if can do in one pass, use the factory hints about how data will be used
    if (npass_ == 1 && !factory()->hints().data_persistent()) {
      // use a subset of a MemoryGrp provided by TransformFactory
      set_memgrp(new MemoryGrpRegion(mem(),localmem));
      ints_acc_ = new DistArray4_MemoryGrp(mem(), num_te_types(),
                                           space1()->rank(), space3()->rank(),
                                           space2()->rank(), space4()->rank(),
                                           memgrp_blksize(),
                                           DistArray4Storage_YX);
      break;
    }
    // else use the next case

  case MOIntsTransform::StoreMethod::mpi:
    ints_acc_ = new DistArray4_MPIIOFile_Ind((file_prefix_+"."+name_).c_str(), num_te_types(),
                                             space1()->rank(), space3()->rank(),
                                             space2()->rank(), space4()->rank(),
                                             DistArray4Storage_YX);
    break;
#endif

  default:
    throw std::runtime_error("TwoBodyMOIntsTransform_iRjS::init_acc() -- invalid integrals store method");
  }

  // now safe to use memorygrp
}

void
TwoBodyMOIntsTransform_iRjS::check_int_symm(double threshold) throw (ProgrammingError)
{
  // this is a partial transform, hence the integrals have not been symmetrized yet
  // assume all is well
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
