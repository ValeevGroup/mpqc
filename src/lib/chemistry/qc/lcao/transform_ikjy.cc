//
// transform_ikjy.cc
//
// Copyright (C) 2004 Edward Valeev
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

#include <stdexcept>
#include <cassert>

#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/group/memregion.h>
#include <util/state/state_bin.h>
#include <util/ref/ref.h>
#include <math/scmat/local.h>
#include <chemistry/qc/lcao/transform_ikjy.h>

// set to 1 when finished rewriting DistArray4_MemoryGrp
#define HAVE_R12IA_MEMGRP 1
#if HAVE_R12IA_MEMGRP
  #include <math/distarray4/distarray4_memgrp.h>
#endif
#include <cassert>

#include <math/distarray4/distarray4_node0file.h>

// set to 1 when finished rewriting DistArray4_MPIIO
#define HAVE_R12IA_MPIIO 1
#ifdef HAVE_MPIIO
#  if HAVE_R12IA_MPIIO
  #include <math/distarray4/distarray4_mpiiofile.h>
#  endif
#endif
#include <chemistry/qc/lcao/transform_123inds.h>

using namespace std;
using namespace sc;

inline int max(int a,int b) { return (a > b) ? a : b;}

/*-----------
  TwoBodyMOIntsTransform_ikjy
 -----------*/
static ClassDesc TwoBodyMOIntsTransform_ikjy_cd(
  typeid(TwoBodyMOIntsTransform_ikjy),"TwoBodyMOIntsTransform_ikjy",1,"public TwoBodyMOIntsTransform",
  0, 0, create<TwoBodyMOIntsTransform_ikjy>);

TwoBodyMOIntsTransform_ikjy::TwoBodyMOIntsTransform_ikjy(const std::string& name, const Ref<MOIntsTransformFactory>& factory,
                                                         const Ref<TwoBodyIntDescr>& tbintdescr,
                                                         const Ref<OrbitalSpace>& space1, const Ref<OrbitalSpace>& space2,
                                                         const Ref<OrbitalSpace>& space3, const Ref<OrbitalSpace>& space4) :
  TwoBodyMOIntsTransform(name,factory,tbintdescr,space1,space2,space3,space4)
{
  init_vars();
}

TwoBodyMOIntsTransform_ikjy::TwoBodyMOIntsTransform_ikjy(StateIn& si) : TwoBodyMOIntsTransform(si)
{
  init_vars();
}

TwoBodyMOIntsTransform_ikjy::~TwoBodyMOIntsTransform_ikjy()
{
}

void
TwoBodyMOIntsTransform_ikjy::save_data_state(StateOut& so)
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
TwoBodyMOIntsTransform_ikjy::compute_transform_dynamic_memory_(int ni) const
{
  // dynamic contribution from 1+2+3 QT
  const int nthread = thr_->nthread();
  const distsize_t memsize123 = (distsize_t) nthread *
                         TwoBodyMOIntsTransform_123Inds::compute_required_dynamic_memory(*this,ni);

  // dynamic contributon from 4 QT
  const int rank2 = space2()->rank();
  const int rank4 = space4()->rank();
  const int nbasis4 = space4()->basis()->nbasis();
  const distsize_t memsize4 = (distsize_t)sizeof(double)*
                                 (rank4 * nbasis4 + // coefs4
                                  rank2 * rank4     // xy
                                 );

  // integrals held by MemoryGrp
  // compute nij as nij on node 0, since nij on node 0 is >= nij on other nodes
  const int nproc = msg()->n();
  const int rank3 = space3()->rank();
  const int nij = compute_nij(ni, rank3, nproc, 0);
  const distsize_t memsize_memgrp = num_te_types() * nij * (distsize_t) memgrp_blksize();

  // determine the peak memory requirements
  return memsize_memgrp + std::max(memsize123,memsize4);
}

size_t
TwoBodyMOIntsTransform_ikjy::memgrp_blksize() const
{
  const int nbasis4 = space4_->basis()->nbasis();
  const int rank4 = space4_->rank();
  const int dim4 = (nbasis4 > rank4) ? nbasis4 : rank4;
  return space2_->rank()*dim4*sizeof(double);
}

void
TwoBodyMOIntsTransform_ikjy::init_acc()
{
  if (ints_acc_.nonnull())
    return;

  const int nij = compute_nij(batchsize_, space3_->rank(), msg_->n(), msg_->me());
  const size_t localmem = num_te_types() * nij * memgrp_blksize();

  switch (ints_method_) {

  case MOIntsTransform::StoreMethod::mem_only:
    {
      if (npass_ > 1)
        throw AlgorithmException("TwoBodyMOIntsTransform_ikjy::init_acc() -- cannot use MemoryGrp-based accumulator in multi-pass transformations; add more memory or swtch to other store methods",
                               __FILE__, __LINE__);
#if HAVE_R12IA_MEMGRP
      {
        // use a subset of a MemoryGrp provided by TransformFactory
        set_memgrp(new MemoryGrpRegion(mem(),localmem));
        ints_acc_ = new DistArray4_MemoryGrp(mem(), num_te_types(),
                                             space1_->rank(), space3_->rank(), space2_->rank(), space4_->rank(),
                                             memgrp_blksize());
      }
#else
      MPQC_ASSERT(false);
#endif
      break;
    }

  case MOIntsTransform::StoreMethod::mem_posix:
    try {
      // if can do in one pass, use the factory hints about how data will be used
      if (npass_ == 1 && !factory()->hints().data_persistent()) {
#if HAVE_R12IA_MEMGRP
        // use a subset of a MemoryGrp provided by TransformFactory
        set_memgrp(new MemoryGrpRegion(mem(),localmem));
        ints_acc_ = new DistArray4_MemoryGrp(mem(), num_te_types(),
                                             space1_->rank(), space3_->rank(), space2_->rank(), space4_->rank(),
                                             memgrp_blksize());
#else
        MPQC_ASSERT(false);
#endif
        break;
      }
    }
    catch (...) {}
    // else use the next case

  case MOIntsTransform::StoreMethod::posix:
    ints_acc_ = new DistArray4_Node0File((file_prefix_+"."+name_).c_str(), num_te_types(),
                                         space1_->rank(), space3_->rank(), space2_->rank(), space4_->rank());
    break;

#ifdef HAVE_MPIIO
  case MOIntsTransform::StoreMethod::mem_mpi:
    try {
      // if can do in one pass, use the factory hints about how data will be used
      if (npass_ == 1 && !factory()->hints().data_persistent()) {
#if HAVE_R12IA_MEMGRP
        // use a subset of a MemoryGrp provided by TransformFactory
        set_memgrp(new MemoryGrpRegion(mem(),localmem));
        ints_acc_ = new DistArray4_MemoryGrp(mem(), num_te_types(),
                                             space1_->rank(), space3_->rank(), space2_->rank(), space4_->rank(),
                                             memgrp_blksize());
#else
        MPQC_ASSERT(false);
#endif
        break;
      }
    }
    catch (...) {}
    // else use the next case

  case MOIntsTransform::StoreMethod::mpi:
#if HAVE_R12IA_MPIIO
    ints_acc_ = new DistArray4_MPIIOFile_Ind((file_prefix_+"."+name_).c_str(), num_te_types(),
                                             space1_->rank(), space3_->rank(), space2_->rank(), space4_->rank());
#else
    MPQC_ASSERT(false);
#endif
    break;
#endif

  default:
    throw InputError("TwoBodyMOIntsTransform_ikjy::init_acc() -- invalid integrals store method", __FILE__, __LINE__);
  }

  // now safe to use MemoryGrp
}

void
TwoBodyMOIntsTransform_ikjy::check_int_symm(double threshold) throw (ProgrammingError)
{
  Ref<DistArray4> iacc = ints_distarray4();
  iacc->activate();
  int num_te_types = iacc->num_te_types();
  int ni = iacc->ni();
  int nj = iacc->nj();
  int nk = iacc->nx();
  int ny = iacc->ny();
  vector<unsigned int> isyms = space1_->orbsym();
  vector<unsigned int> jsyms = space3_->orbsym();
  vector<unsigned int> ksyms = space2_->orbsym();
  vector<unsigned int> ysyms = space4_->orbsym();

  int me = msg_->me();
  vector<int> twi_map;
  int ntasks_with_ints = iacc->tasks_with_access(twi_map);
  if (!iacc->has_access(me))
    return;

  int ij=0;
  for(int i=0; i<ni; i++) {
    int isym = isyms[i];
    for(int j=0; j<nj; j++, ij++) {
      int jsym = jsyms[j];
      if (ij%ntasks_with_ints != twi_map[me])
        continue;

      for(int t=0; t<num_te_types; t++) {
        const double* ints = iacc->retrieve_pair_block(i,j,static_cast<DistArray4::tbint_type>(t));
        int ky=0;
        for(int k=0; k<nk; k++) {
          int ksym = ksyms[k];
          for(int y=0; y<ny; y++, ky++) {
            int ysym = ysyms[y];
            if ( (isym^jsym^ksym^ysym) != 0 && fabs(ints[ky]) > threshold) {
              ExEnv::outn() << scprintf("Integral type=%d i=%d k=%d j=%d y=%d should be zero\n",t,i,k,j,y);
              throw ProgrammingError("TwoBodyMOIntsTransform_ikjy::check_int_symm() -- nonzero nonsymmetric integrals are detected",
                                     __FILE__, __LINE__);
            }
          }
        }
        iacc->release_pair_block(i,j,static_cast<DistArray4::tbint_type>(t));
      }
    }
  }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
