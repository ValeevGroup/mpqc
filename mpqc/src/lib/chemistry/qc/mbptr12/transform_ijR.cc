//
// transform_ijR.cc
//
// Copyright (C) 2009 Edward Valeev
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

#ifdef __GNUG__
#pragma implementation
#endif

// includes go here
#include<chemistry/qc/mbptr12/transform_ijR.h>
#include <chemistry/qc/mbptr12/r12ia_memgrp.h>
#include <chemistry/qc/mbptr12/r12ia_node0file.h>
#ifdef HAVE_MPIIO
#  include <chemistry/qc/mbptr12/r12ia_mpiiofile.h>
#endif
#include <util/group/memory.h>
#include <util/group/memregion.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/print.h>

using namespace sc;

static ClassDesc TwoBodyThreeCenterMOIntsTransform_ijR_cd(
  typeid(TwoBodyThreeCenterMOIntsTransform_ijR),"TwoBodyThreeCenterMOIntsTransform_ijR",1,
  "public TwoBodyThreeCenterMOIntsTransform",
  0, 0, create<TwoBodyThreeCenterMOIntsTransform_ijR>);

TwoBodyThreeCenterMOIntsTransform_ijR::~TwoBodyThreeCenterMOIntsTransform_ijR() {}

TwoBodyThreeCenterMOIntsTransform_ijR::TwoBodyThreeCenterMOIntsTransform_ijR(const std::string& name,
  const Ref<MOIntsTransformFactory>& factory,
  const Ref<TwoBodyThreeCenterIntDescr>& tbintdescr,
  const Ref<OrbitalSpace>& space1,
  const Ref<OrbitalSpace>& space2,
  const Ref<OrbitalSpace>& space3) :
    TwoBodyThreeCenterMOIntsTransform(name, factory, tbintdescr,
                                      space1, space2, space3)
  {
    init_vars();
  }

TwoBodyThreeCenterMOIntsTransform_ijR::TwoBodyThreeCenterMOIntsTransform_ijR(StateIn& si) :
  TwoBodyThreeCenterMOIntsTransform(si) {
  assert(false);
  init_vars();
}

void
TwoBodyThreeCenterMOIntsTransform_ijR::save_data_state(StateOut& so) {
  assert(false);
}

int
TwoBodyThreeCenterMOIntsTransform_ijR::compute_transform_batchsize(size_t mem_static, int rank_R)
{
  // Check is have enough for even static objects
  size_t mem_dyn = 0;
  const size_t max_memory = factory_->memory();
  if (max_memory <= mem_static)
    return 0;
  else
    mem_dyn = max_memory - mem_static;

  // Determine if calculation is possible at all (i.e., if nR=nmaxfunR possible)
  const int nmaxfunR = space3()->basis()->max_nfunction_in_shell();
  distsize_t maxdyn = compute_transform_dynamic_memory(nmaxfunR);
  if (maxdyn > mem_dyn) {
    return 0;
  }

  int nR = nmaxfunR+1;
  while (nR<=rank_R) {
    maxdyn = compute_transform_dynamic_memory(nR);
    if (maxdyn >= mem_dyn) {
      nR--;
      break;
    }
    nR++;
  }
  if (nR > rank_R) nR = rank_R;

  return nR;
}

distsize_t
TwoBodyThreeCenterMOIntsTransform_ijR::compute_transform_dynamic_memory(int batchsize) const
{
  const unsigned int rank1 = space1()->rank();
  const unsigned int rank2 = space2()->rank();
  const unsigned int rank3 = space3()->rank();
  if (batchsize == -1)
    batchsize = rank3;
  // can only do 1-pass transformation
  assert(batchsize == rank3);

  return num_te_types() * rank1 * static_cast<distsize_t>(rank2 * rank3 * sizeof(double));
}

void
TwoBodyThreeCenterMOIntsTransform_ijR::extra_memory_report(std::ostream& os) const
{
  const int rank_R_restart = space3()->rank() - restart_orbital_;
  const int nmaxfunR = space3()->basis()->max_nfunction_in_shell();

  os << indent
     << "Number of passes:               " << (rank_R_restart+batchsize_-1)/batchsize_
     << endl;
  os << indent
     << "Memory required for one pass:   "
     << compute_transform_dynamic_memory(rank_R_restart)+static_memory_
     << " Bytes"
     << endl;
  os << indent
     << "Minimum memory required:        "
     << compute_transform_dynamic_memory(nmaxfunR)+static_memory_
     << " Bytes"
   << endl;
}

void
TwoBodyThreeCenterMOIntsTransform_ijR::init_acc() {
  if (ints_acc_.nonnull())
    return;

  const int nproc = mem_->n();
  const size_t blksize = space2()->rank() * space3()->rank() * sizeof(double);
  const size_t localmem = (num_te_types() * space1()->rank() * blksize  + nproc - 1) / nproc;

  switch (ints_method_) {

  case MOIntsTransform::StoreMethod::mem_only:
    {
      // use a subset of a MemoryGrp provided by TransformFactory
      Ref<MemoryGrp> mem = new MemoryGrpRegion(factory()->mem(),localmem);
      ints_acc_ = new R12IntsAcc_MemoryGrp(mem, num_te_types(),
                                           1, space1()->rank(),
                                           space2()->rank(), space3()->rank(),
                                           blksize);
    }
    break;

  case MOIntsTransform::StoreMethod::mem_posix:
    // if can do in one pass, use the factory hints about how data will be used
    if (!factory()->hints().data_persistent()) {
      // use a subset of a MemoryGrp provided by TransformFactory
      Ref<MemoryGrp> mem = new MemoryGrpRegion(factory()->mem(),localmem);
      ints_acc_ = new R12IntsAcc_MemoryGrp(mem, num_te_types(),
                                           1, space1()->rank(),
                                           space2()->rank(), space3()->rank(),
                                           blksize);
      break;
    }
    // else use the next case

  case MOIntsTransform::StoreMethod::posix:
    ints_acc_ = new R12IntsAcc_Node0File((file_prefix_+"."+name_).c_str(), num_te_types(),
                                         1, space1()->rank(), space2()->rank(), space3()->rank());
    break;

#if HAVE_MPIIO
  case MOIntsTransform::StoreMethod::mem_mpi:
    // if can do in one pass, use the factory hints about how data will be used
    if (!factory()->hints().data_persistent()) {
      // use a subset of a MemoryGrp provided by TransformFactory
      Ref<MemoryGrp> mem = new MemoryGrpRegion(factory()->mem(),localmem);
      ints_acc_ = new R12IntsAcc_MemoryGrp(mem, num_te_types(),
                                           1, space1()->rank(),
                                           space2()->rank(), space3()->rank(),
                                           blksize);
      break;
    }
    // else use the next case

  case MOIntsTransform::StoreMethod::mpi:
    ints_acc_ = new R12IntsAcc_MPIIOFile_Ind((file_prefix_+"."+name_).c_str(), num_te_types(),
                                             1, space1()->rank(), space2()->rank(), space3()->rank());
    break;
#endif

  default:
    throw std::runtime_error("TwoBodyThreeCenterMOIntsTransform_ijR::init_acc() -- invalid integrals store method");
  }

}

void
TwoBodyThreeCenterMOIntsTransform_ijR::compute() {

  // if all integrals are already computed, do nothing
  if (space3()->rank() == restart_orbital_)
    return;

  std::string tim_label("tbint_tform_ijR ");
  tim_label += this->name();
  Timer tim(tim_label);

  const int nproc = mem_->n();
  const int me = mem_->me();

  const Ref<GaussianBasisSet>& b1 = this->space1()->basis();
  const Ref<GaussianBasisSet>& b2 = this->space2()->basis();
  const Ref<GaussianBasisSet>& b3 = this->space3()->basis();

  const int num_te_types = this->num_te_types();
  const int n1 = this->space1()->rank();
  const int n2 = this->space2()->rank();
  const int n3 = b3->nbasis();
  const int n23 = n2*n3;
  const distsize_t ijR_globalsize = (((static_cast<distsize_t>(n1))*n23)*num_te_types)*sizeof(double);
  const int ni_local = (n1 + nproc - 1)/ nproc;
  const size_t memgrp_blocksize = (static_cast<size_t>(n23))*sizeof(double);
  const size_t ijR_localsize = num_te_types * ni_local * memgrp_blocksize;
  this->alloc_mem(ijR_localsize);
  memset(mem_->localdata(), 0, ijR_localsize);

  const int nbasis1 = b1->nbasis();
  const int nbasis2 = b2->nbasis();
  int nfuncmax3 = b3->max_nfunction_in_shell();

  // get scratch storage
  double** pq_ints = new double*[num_te_types];
  for(int te_type=0; te_type<num_te_types; te_type++)
    pq_ints[te_type] = new double[nbasis1 * nbasis2 * nfuncmax3];
  double* iq_ints = new double[n1 * nbasis2 * nfuncmax3];
  memset(iq_ints, 0, n1 * nbasis2 * nfuncmax3 * sizeof(double));
  double* jR_ints = new double[n2 * nfuncmax3];
  memset(jR_ints, 0, n2 * nfuncmax3 * sizeof(double));
  RefSCMatrix mocoefs1 = space1_->coefs();   // (pi)
  double* pi = new double[nbasis1 * n1];
  mocoefs1.convert(pi);  mocoefs1 = 0;
  RefSCMatrix mocoefs2 = space2_->coefs();   // (pi)
  double* qj = new double[nbasis2 * n2];
  mocoefs2.convert(qj);  mocoefs2 = 0;

  // get the integral evaluator
  Ref<Integral> integral = tbintdescr_->factory();
  integral->set_basis(b1, b2, b3);
  Ref<TwoBodyThreeCenterInt> inteval = tbintdescr_->inteval();
  const Ref<TwoBodyOperSetDescr>& descr = inteval->descr();
  const double **buffer = new const double*[num_te_types];
  for(int te_type=0; te_type<num_te_types; te_type++)
    buffer[te_type] = inteval->buffer( descr->opertype(te_type) );

  // distribute work by basis3
  // TODO 1) use DistShell 2) use threads
  for (int s3 = 0; s3 < b3->nshell(); ++s3) {
    if (s3 % nproc != me) continue;

    const int s3offset = b3->shell_to_function(s3);
    const int nf3 = b3->shell(s3).nfunction();
    const int nb2f3 = nf3 * nbasis2;

    // compute the entire (p1 p2| block
    for (int s1 = 0; s1 < b1->nshell(); ++s1) {
      const int s1offset = b1->shell_to_function(s1);
      const int nf1 = b1->shell(s1).nfunction();

      for (int s2 = 0; s2 < b2->nshell(); ++s2) {
        const int s2offset = b2->shell_to_function(s2);
        const int nf2 = b2->shell(s2).nfunction();
        const int nf23 = nf2 * nf3;

        // compute shell triplet
        inteval->compute_shell(s1, s2, s3);

        // copy buffer into pq_ints
        // for each f1 copy s2 s3 block to ((f1+s1offset) * nbasis2 + s2offset) * nf3
        for(int te_type=0; te_type < num_te_types; ++te_type) {
          const double* f1s2s3_src = buffer[te_type];
          double* f1s2s3_dst = pq_ints[te_type] + (s1offset * nbasis2 + s2offset) * nf3;
          for(int f1=0; f1<nf1; ++f1, f1s2s3_src += nf23, f1s2s3_dst += nb2f3) {
            std::copy(f1s2s3_src, f1s2s3_src + nf23, f1s2s3_dst);
          }
        }

      }
    }

#if 0
    if (s3 == 1) {
    for(int i=0; i<nbasis1; ++i) {
      for(int j=0; j<nbasis2; ++j) {
        const double value = pq_ints[(i*nbasis2 + j)*nf3];
        ExEnv::outn() << "i = " << i << " j = " << j << " R = 1  value = " << value << endl;
      }
    }
    }
#endif

    // **** transform (pq|R) to (ij|R) ****
    const char notransp = 'n';
    const char transp = 't';
    const double one = 1.0;
    const double zero = 0.0;
    for(int te_type=0; te_type < num_te_types; ++te_type) {
      // (iq|R) = (ip) * (pq|R)
      F77_DGEMM(&notransp,&transp,&nb2f3,&n1,&nbasis1,&one,pq_ints[te_type],&nb2f3,
                pi,&n1,&zero,iq_ints,&nb2f3);
      // for each i: (ij|R) = (jq) * (iq|R)
      const double* qR_ints = iq_ints;
      for(int i=0; i<n1; ++i) {
        F77_DGEMM(&notransp,&transp,&nf3,&n2,&nbasis2,&one,qR_ints,&nf3,
                  qj,&n2,&zero,jR_ints,&nf3);

        // accumulate to each memory location
        const int i_proc = i % nproc;
        const int i_local = i / nproc;
        const size_t i_offset = (i_local * num_te_types + te_type) * static_cast<size_t>(n23) ;
        size_t ijR_offset = i_offset + s3offset;
        double* R_ptr = jR_ints;
        for(int j=0; j<n2; ++j) {
#if 0
    if (s3 == 1) {
      const double value = *R_ptr;
      ExEnv::outn() << "i = " << i << " j = " << j << " R = 1  i_proc = " << i_proc
                    << " ijR_offset = " << ijR_offset << "  value = " << value << endl;
    }
#endif
          mem_->sum_reduction_on_node(R_ptr, ijR_offset , nf3, i_proc);
          ijR_offset += n3;
          R_ptr += nf3;
        }

        qR_ints += nb2f3;
      } // end of i loop
    } // end of te_type loop

  } // end of loop over R shells

  mem_->sync();

  // dump integrals from MemoryGrp to R12IntsAcc
  ints_acc_->activate();
  detail::store_memorygrp(ints_acc_,mem_,0,1,memgrp_blocksize);
  if (ints_acc_->data_persistent()) ints_acc_->deactivate();

  restart_orbital_ = space3()->rank();
  tim.exit();
  ExEnv::out0() << indent << "Built TwoBodyMOIntsTransform_ijR: name = " << this->name() << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
