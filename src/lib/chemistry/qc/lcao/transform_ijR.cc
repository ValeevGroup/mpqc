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

#include <cassert>
#include<chemistry/qc/lcao/transform_ijR.h>
#include <math/distarray4/distarray4_memgrp.h>
#include <math/distarray4/distarray4_node0file.h>
#ifdef HAVE_MPIIO
#  include <math/distarray4/distarray4_mpiiofile.h>
#endif
#include <util/group/memory.h>
#include <util/group/memregion.h>
#include <util/misc/consumableresources.h>
#include <math/scmat/blas.h>
#include <util/misc/print.h>

using namespace std;
using namespace sc;

ClassDesc TwoBodyThreeCenterMOIntsTransform_ijR::class_desc_(
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
  MPQC_ASSERT(false);
  init_vars();
}

void
TwoBodyThreeCenterMOIntsTransform_ijR::save_data_state(StateOut& so) {
  TwoBodyThreeCenterMOIntsTransform::save_data_state(so);
}

int
TwoBodyThreeCenterMOIntsTransform_ijR::compute_transform_batchsize(size_t mem_static, int rank_R)
{
  // Check is have enough for even static objects
  size_t mem_dyn = 0;
  const size_t max_memory = ConsumableResources::get_default_instance()->memory();
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
  MPQC_ASSERT(batchsize == rank3);

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
  if (ints_acc_)
    return;

  const int nproc = mem_->n();
  const size_t blksize = space2()->rank() * space3()->rank() * sizeof(double);
  const size_t n1_local = (space1()->rank() + nproc - 1) / nproc ;
  const size_t localmem =  num_te_types() * n1_local * blksize;

  switch (ints_method_) {

  case MOIntsTransform::StoreMethod::mem_only:
    {
      // use a subset of a MemoryGrp provided by TransformFactory
      set_memgrp(new MemoryGrpRegion(mem(),localmem));
      ints_acc_ = new DistArray4_MemoryGrp(mem(), num_te_types(),
                                           1, space1()->rank(),
                                           space2()->rank(), space3()->rank(),
                                           blksize);
    }
    break;

  case MOIntsTransform::StoreMethod::mem_posix:
    // if can do in one pass, use the factory hints about how data will be used
    if (!factory()->hints().data_persistent()) {
      try {
        // use a subset of a MemoryGrp provided by TransformFactory
        set_memgrp(new MemoryGrpRegion(mem(),localmem));
        ints_acc_ = new DistArray4_MemoryGrp(mem(), num_te_types(),
                                             1, space1()->rank(),
                                             space2()->rank(), space3()->rank(),
                                             blksize);
        break;
      }
      catch(...) {}
    }
    // else use the next case

  case MOIntsTransform::StoreMethod::posix:
    ints_acc_ = new DistArray4_Node0File((file_prefix_+"."+name_).c_str(), num_te_types(),
                                         1, space1()->rank(), space2()->rank(), space3()->rank());
    break;

#ifdef HAVE_MPIIO
  case MOIntsTransform::StoreMethod::mem_mpi:
    // if can do in one pass, use the factory hints about how data will be used
    if (!factory()->hints().data_persistent()) {
      try {
        // use a subset of a MemoryGrp provided by TransformFactory
        set_memgrp(new MemoryGrpRegion(mem(),localmem));
        ints_acc_ = new DistArray4_MemoryGrp(mem(), num_te_types(),
                                             1, space1()->rank(),
                                             space2()->rank(), space3()->rank(),
                                             blksize);
        break;
      }
      catch (...) {}
    }
    // else use the next case

  case MOIntsTransform::StoreMethod::mpi:
    ints_acc_ = new DistArray4_MPIIOFile_Ind((file_prefix_+"."+name_).c_str(), num_te_types(),
                                             1, space1()->rank(), space2()->rank(), space3()->rank());
    break;
#endif

  default:
    throw InputError("TwoBodyThreeCenterMOIntsTransform_ijR::init_acc() -- invalid integrals store method",
                     __FILE__, __LINE__);
  }

}

void
TwoBodyThreeCenterMOIntsTransform_ijR::compute() {

  // if all integrals are already computed, do nothing
  if (space3()->rank() == restart_orbital_)
    return;

  // determine whether space1, space2, and space3 are AO spaces
  Ref<AOSpaceRegistry> aoidxreg = this->factory()->ao_registry();
  const bool space1_is_ao = aoidxreg->value_exists( this->space1() );
  const bool space2_is_ao = aoidxreg->value_exists( this->space2() );
  const bool space3_is_ao = aoidxreg->value_exists( this->space3() );

  if (space1_is_ao)
    compute_pjR();
  else
    compute_ijR();

  restart_orbital_ = space3()->rank();
}

void
TwoBodyThreeCenterMOIntsTransform_ijR::compute_ijR() {

  // determine whether space1, space2, and space3 are AO spaces
  Ref<AOSpaceRegistry> aoidxreg = this->factory()->ao_registry();
  const bool space1_is_ao = false;
  const bool space2_is_ao = aoidxreg->value_exists( this->space2() );
  const bool space3_is_ao = aoidxreg->value_exists( this->space3() );

  std::string tim_label("tbint_tform_ijR ");
  tim_label += this->name();
  Timer tim(tim_label);

  const int nproc = mem_->n();
  const int me = mem_->me();

  const Ref<GaussianBasisSet>& b1 = this->space1()->basis();
  const Ref<GaussianBasisSet>& b2 = this->space2()->basis();
  const Ref<GaussianBasisSet>& b3 = this->space3()->basis();
  const bool b1_equiv_b2 = b1->equiv(b2);

  const int num_te_types = this->num_te_types();
  const blasint n1 = this->space1()->rank();
  const blasint n2 = this->space2()->rank();
  const blasint n3 = b3->nbasis();
  const blasint n23 = n2*n3;
  const distsize_t ijR_globalsize = (((static_cast<distsize_t>(n1))*n23)*num_te_types)*sizeof(double);
  const int ni_local = (n1 + nproc - 1)/ nproc;
  const size_t memgrp_blocksize = (static_cast<size_t>(n23))*sizeof(double);
  const size_t ijR_localsize = num_te_types * ni_local * memgrp_blocksize;
  this->alloc_mem(ijR_localsize);
  memset(mem_->localdata(), 0, ijR_localsize);

  const blasint nbasis1 = b1->nbasis();
  const blasint nbasis2 = b2->nbasis();
  int nfuncmax3 = b3->max_nfunction_in_shell();

  // get scratch storage
  double** pq_ints = new double*[num_te_types];
  const size_t pq_ints_size = nbasis1 * nbasis2 * nfuncmax3;
  for(int te_type=0; te_type<num_te_types; te_type++) {
    pq_ints[te_type] = new double[pq_ints_size];
  }

  double* iq_ints;
  double* pi;
  if (!space1_is_ao) {
    iq_ints = new double[n1 * nbasis2 * nfuncmax3];
    memset(iq_ints, 0, n1 * nbasis2 * nfuncmax3 * sizeof(double));
    RefSCMatrix mocoefs1 = space1_->coefs();   // (pi)
    pi = new double[nbasis1 * n1];
    mocoefs1.convert(pi);  mocoefs1 = 0;
  }

  double* jR_ints;
  double* qj;
  if (!space2_is_ao) {
    jR_ints = new double[n2 * nfuncmax3];
    memset(jR_ints, 0, n2 * nfuncmax3 * sizeof(double));
    RefSCMatrix mocoefs2 = space2_->coefs();   // (pi)
    qj = new double[nbasis2 * n2];
    mocoefs2.convert(qj);  mocoefs2 = 0;
  }

  // get the integral evaluator
  Ref<Integral> integral = tbintdescr_->factory();
  integral->set_basis(b1, b2, b3);
  Ref<TwoBodyThreeCenterInt> inteval = tbintdescr_->inteval();
  const Ref<TwoBodyOperSetDescr>& descr = inteval->descr();
  const double **buffer = new const double*[num_te_types];
  for(int te_type=0; te_type<num_te_types; te_type++)
    buffer[te_type] = inteval->buffer( descr->opertype(te_type) );

  // distribute work by basis3
  // TODO 1) use DistShell -- not too important now
  // TODO 2) use threads
  for (int s3 = 0; s3 < b3->nshell(); ++s3) {
    if (s3 % nproc != me) continue;

    const int s3offset = b3->shell_to_function(s3);
    const blasint nf3 = b3->shell(s3).nfunction();
    const blasint nb2f3 = nf3 * nbasis2;
    const blasint nb1nb2f3 = nbasis1 * nb2f3;
    for(int te_type=0; te_type < num_te_types; ++te_type) {
      memset(pq_ints[te_type], 0, nb1nb2f3*sizeof(double));
    }

    // compute the entire (p1 p2| block
    for (int s1 = 0; s1 < b1->nshell(); ++s1) {
      const int s1offset = b1->shell_to_function(s1);
      const int nf1 = b1->shell(s1).nfunction();

      const int s2_fence = b1_equiv_b2 ? s1+1 : b2->nshell();
      for (int s2 = 0; s2 < s2_fence; ++s2) {

        const bool copy_to_s2_s1 = b1_equiv_b2 && s1 != s2;

        const int s2offset = b2->shell_to_function(s2);
        const int nf2 = b2->shell(s2).nfunction();
        const int nf23 = nf2 * nf3;

        // skip the insignificant integrals
        const double log2_cauchy_bound = inteval->log2_shell_bound(s1, s2, s3);
        if (log2_cauchy_bound < this->log2_epsilon()) {
          continue;
        }

        // compute shell triplet
        inteval->compute_shell(s1, s2, s3);

        // compare the actual bound with the estimated bound
        if (0) {
          const double actual_max = * std::max_element(buffer[0], buffer[0]+nf23*nf1, abs_less<double>());
          ExEnv::outn() << scprintf("s1=%d s2=%d s3=%d log2_cauchy=%10.5e log2_actual_max=%10.5e\n",
                                    s1, s2, s3, log2_cauchy_bound, log(abs(actual_max))/log(2.0));
        }

        // copy buffer into pq_ints
        // for each f1 copy s2 s3 block to ((f1+s1offset) * nbasis2 + s2offset) * nf3
        for(int te_type=0; te_type < num_te_types; ++te_type) {
          const double* f1s2s3_src = buffer[te_type];
          double* f1s2s3_dst = pq_ints[te_type] + (s1offset * nbasis2 + s2offset) * nf3;
          for(int f1=0; f1<nf1; ++f1, f1s2s3_src += nf23, f1s2s3_dst += nb2f3) {
            std::copy(f1s2s3_src, f1s2s3_src + nf23, f1s2s3_dst);
          }
        }
        if (copy_to_s2_s1) { // if also need s2 s1
          // for each f1f2 copy s3 block to ((f2+s2offset) * nbasis1 + f1+s1offset) * nf3, note that nbasis1 equiv nbasis2
          for(int te_type=0; te_type < num_te_types; ++te_type) {
            for(int f2=0; f2<nf2; ++f2) {
              double* f2f1s3_dst = pq_ints[te_type] + ((f2+s2offset) * nbasis1 + s1offset) * nf3;
              const double* f1f2s3_src = buffer[te_type] + f2*nf3;
              for(int f1=0; f1<nf1; ++f1, f1f2s3_src += nf23, f2f1s3_dst += nf3) {
                std::copy(f1f2s3_src, f1f2s3_src + nf3, f2f1s3_dst);
              }
            }
          }
        }

      }
    }

#if 0
    if (s3 == 1) {
    for(int i=0; i<nbasis1; ++i) {
      for(int j=0; j<nbasis2; ++j) {
        const double value = pq_ints[0][(i*nbasis2 + j)*nf3];
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
      if (!space1_is_ao) {
        F77_DGEMM(&notransp,&transp,
                  &nb2f3,&n1,&nbasis1,
                  &one,pq_ints[te_type],&nb2f3,
                  pi,&n1,
                  &zero,iq_ints,&nb2f3);
      }
      else
        iq_ints = pq_ints[te_type];

      // for each i: (ij|R) = (jq) * (iq|R)
      const double* qR_ints = iq_ints;
      for(int i=0; i<n1; ++i) {

        if (!space2_is_ao)
          F77_DGEMM(&notransp,&transp,
                    &nf3,&n2,&nbasis2,
                    &one,qR_ints,&nf3,
                    qj,&n2,
                    &zero,jR_ints,&nf3);
        else
          jR_ints = const_cast<double*>(qR_ints);

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

  // transform R now, if necessary
  if (not space3_is_ao) {
    throw FeatureNotImplemented("TwoBodyThreeCenterMOIntsTransform_ijR: non-AO R space in <ij|R> is not yet supported",
                                __FILE__, __LINE__, this->class_desc());
  }

  // dump integrals from MemoryGrp to DistArray4
  ints_acc_->activate();
  detail::store_memorygrp(ints_acc_,mem_,0,1,memgrp_blocksize);
  if (ints_acc_->data_persistent()) ints_acc_->deactivate();

  // cleanup
  delete[] pq_ints[0];
  delete[] pq_ints;
  if (!space1_is_ao) {
    delete[] iq_ints;
    delete[] pi;
  }
  if (!space2_is_ao) {
    delete[] jR_ints;
    delete[] qj;
  }
  delete[] buffer;

  tim.exit();
  ExEnv::out0() << indent << "Built TwoBodyMOIntsTransform_ijR: name = " << this->name() << std::endl;
}

void
TwoBodyThreeCenterMOIntsTransform_ijR::compute_pjR() {

  // foreach p shell
  //   compute qS integrals
  //   transform qS -> jR
  //   store jR to ints_acc

  long int shellset_counter = 0;;

  // determine whether space1, space2, and space3 are AO spaces
  Ref<AOSpaceRegistry> aoidxreg = this->factory()->ao_registry();
  const bool space1_is_ao = true;
  const bool space2_is_ao = aoidxreg->value_exists( this->space2() );
  const bool space3_is_ao = aoidxreg->value_exists( this->space3() );

  std::string tim_label("tbint_tform_ijR ");
  tim_label += this->name();
  Timer tim(tim_label);

  const int nproc = mem_->n();
  const int me = mem_->me();

  const Ref<GaussianBasisSet>& b1 = this->space1()->basis();
  const Ref<GaussianBasisSet>& b2 = this->space2()->basis();
  const Ref<GaussianBasisSet>& b3 = this->space3()->basis();
  const bool b1_equiv_b2 = b1->equiv(b2);

  const int num_te_types = this->num_te_types();
  const blasint n1 = this->space1()->rank();
  const blasint n2 = this->space2()->rank();
  const blasint n3 = this->space3()->rank();
  const blasint n23 = n2*n3;
  const blasint nbasis2 = b2->nbasis();
  const blasint nbasis3 = b3->nbasis();
  const blasint nb23 = nbasis2 * nbasis3;
  int nfuncmax1 = b1->max_nfunction_in_shell();

  Ref<DistArray4_MemoryGrp> ints_acc_cast; ints_acc_cast << ints_acc_;
  const bool need_memgrp = ints_acc_cast;
  if (need_memgrp) {
    const distsize_t ijR_globalsize = (((static_cast<distsize_t>(n1))*n23)*num_te_types)*sizeof(double);
    const int ni_local = (n1 + nproc - 1)/ nproc;
    const size_t memgrp_blocksize = (static_cast<size_t>(n23))*sizeof(double);
    const size_t ijR_localsize = num_te_types * ni_local * memgrp_blocksize;
    this->alloc_mem(ijR_localsize);
    memset(mem_->localdata(), 0, ijR_localsize);
  }

  // get scratch storage
  const size_t qS_ints_size = nfuncmax1 * nbasis2 * nbasis3;
  double* qS_ints_ptr = allocate<double>(num_te_types * qS_ints_size);
  std::vector< double* > qS_ints(num_te_types);
  qS_ints[0] = qS_ints_ptr;
  for(int te_type=1; te_type<num_te_types; te_type++) {
    qS_ints[te_type] = qS_ints[te_type - 1] + qS_ints_size;
  }

  double* jS_ints;
  double* qj;
  if (!space2_is_ao) {
    jS_ints = allocate<double>(n2 * nbasis3);
    RefSCMatrix mocoefs2 = space2_->coefs();   // (pi)
    qj = allocate<double>(nbasis2 * n2);
    mocoefs2.convert(qj);  mocoefs2 = 0;
  }

  double* jR_ints;
  double* SR;
  if (!space3_is_ao) {
    jR_ints = allocate<double>(n23);
    RefSCMatrix mocoefs3 = space3_->coefs();   // (SR)
    SR = allocate<double>(nbasis3 * n3);
    mocoefs3.convert(SR);  mocoefs3 = 0;
  }

  // get the integral evaluator
  Ref<Integral> integral = tbintdescr_->factory();
  integral->set_basis(b1, b2, b3);
  Ref<TwoBodyThreeCenterInt> inteval = tbintdescr_->inteval();
  const Ref<TwoBodyOperSetDescr>& descr = inteval->descr();
  std::vector<const double*> buffer(num_te_types);
  for(int te_type=0; te_type<num_te_types; te_type++)
    buffer[te_type] = inteval->buffer( descr->opertype(te_type) );

  ints_acc_->activate();

  std::vector<int> tasks_with_access;
  const int ntasks_with_access = ints_acc_->tasks_with_access(tasks_with_access);

  // distribute work by basis1
  // TODO 1) use DistShell -- not too important now
  // TODO 2) use threads
  if (tasks_with_access[me] != -1) { // work, if can write to ints_acc
    for (int s1 = 0; s1 < b1->nshell(); ++s1) {
      const int s1offset = b1->shell_to_function(s1);
      const int nf1 = b1->shell(s1).nfunction();

      // static round-robin load distribution
      if (s1 % ntasks_with_access != tasks_with_access[me])
        continue;

      //////////
      // compute the entire (... q| S) block where ... is shell s1
      // multiple threads will cooperate here
      //////////

      const size_t f1nb2nb3 = nf1 * nbasis2 * nbasis3;
      for (int te_type = 0; te_type < num_te_types; ++te_type) {
        memset(qS_ints[te_type], 0, f1nb2nb3 * sizeof(double));
      }

      for (int s2 = 0; s2 < b2->nshell(); ++s2) {

        const int s2offset = b2->shell_to_function(s2);
        const int nf2 = b2->shell(s2).nfunction();

        for (int s3 = 0; s3 < b3->nshell(); ++s3) {

          // skip the insignificant integrals
          const double log2_cauchy_bound = inteval->log2_shell_bound(s1, s2, s3);
          if (log2_cauchy_bound < this->log2_epsilon()) {
            continue;
          }
          ++shellset_counter;

          const int s3offset = b3->shell_to_function(s3);
          const blasint nf3 = b3->shell(s3).nfunction();
          const blasint nf23 = nf2 * nf3;

          // compute shell triplet
          inteval->compute_shell(s1, s2, s3);

          // compare the actual bound with the estimated bound
          if (0) {
            const double actual_max = *std::max_element(buffer[0],
                                                        buffer[0] + nf1 * nf23,
                                                        abs_less<double>());
            ExEnv::outn()
                << scprintf(
                    "s1=%d s2=%d s3=%d log2_cauchy=%10.5e log2_actual_max=%10.5e\n",
                    s1, s2, s3, log2_cauchy_bound,
                    log(abs(actual_max)) / log(2.0));
          }

          // copy buffer into qS_ints
          // for each f1 copy s2 s3 block to the isometric block with left top corner at
          // (f1 * nbasis2 + s2offset) * nbasis3 + s3offset
          for (int te_type = 0; te_type < num_te_types; ++te_type) {
            const double* f1s2s3_src = buffer[te_type];
            double* f1s2s3_dst = qS_ints[te_type] + s2offset * nbasis3
                + s3offset;
            for (int f1 = 0; f1 < nf1; ++f1, f1s2s3_dst += nb23) {
              for (int f2 = 0; f2 < nf2; ++f2, f1s2s3_src += nf3) {
                std::copy(f1s2s3_src, f1s2s3_src + nf3,
                          f1s2s3_dst + f2 * nbasis3);
              }
            }
          }

        }
      }

      // **** transform (pq|S) to (pj|S) ****
      const char notransp = 'n';
      const char transp = 't';
      const double one = 1.0;
      const double zero = 0.0;
      for (int te_type = 0; te_type < num_te_types; ++te_type) {

        // for each p: (pq|S) = (jq) * (pq|S)
        const double* pqS_ints = qS_ints[te_type];
        for (int p = 0; p < nf1; ++p, pqS_ints += nb23) {

          if (!space2_is_ao)
            C_DGEMM(transp, notransp, nbasis2, n2, nbasis3, one, qj, n2,
                    pqS_ints, nbasis3, zero, jS_ints, nbasis3);
          else
            jS_ints = const_cast<double*>(pqS_ints);

          // (pj|R) = (pj|S) * (S|R)
          if (!space3_is_ao)
            C_DGEMM(notransp, notransp, n2, nbasis3, n3, one, jS_ints, nbasis3,
                    SR, n3, zero, jR_ints, n3);
          else
            jR_ints = const_cast<double*>(jS_ints);

          // store jR block to ints_acc
          ints_acc_->store_pair_block(0, s1offset + p, te_type, jR_ints);

        } // end of p loop

      } // end of te_type loop

    } // end of loop over p shells
  } // tasks with access

  GrpSumReduce<long int> sumred;
  MessageGrp::get_default_messagegrp()->reduce(&shellset_counter, 1, sumred);
  mem_->sync();

  if (debug() > 0 && mem_->me() == 0) {
    ExEnv::out0() << indent << "total    # of shell sets = " << b1->nshell() * b2->nshell() * b3->nshell() << endl;
    ExEnv::out0() << indent << "computed # of shell sets = " << shellset_counter << endl;
  }

  if (ints_acc_->data_persistent()) ints_acc_->deactivate();

  // cleanup
  deallocate(qS_ints_ptr);
  if (not space2_is_ao) {
    deallocate(qj);
    deallocate(jS_ints);
  }
  if (not space3_is_ao) {
    deallocate(SR);
    deallocate(jR_ints);
  }

  tim.exit();
  ExEnv::out0() << indent << "Built TwoBodyMOIntsTransform_ijR: name = " << this->name() << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR::class_desc_(
  typeid(TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR),"TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR",1,
  "public TwoBodyThreeCenterMOIntsTransform_ijR",
  0, 0, create<TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR>);

TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR::~TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR() {}

TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR::TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR(const std::string& name,
  const Ref<TwoBodyThreeCenterMOIntsTransform_ijR>& iqR_tform, const Ref<OrbitalSpace>& space2) :
    TwoBodyThreeCenterMOIntsTransform_ijR(name, iqR_tform->factory(), iqR_tform->intdescr(),
                                          iqR_tform->space1(), space2, iqR_tform->space3()),
                                          iqR_tform_(iqR_tform)
{
  Ref<AOSpaceRegistry> aoidxreg = this->factory()->ao_registry();
  // determine whether space1, space2, and space3 are AO spaces
  const bool space1_is_ao = aoidxreg->value_exists( this->space1() );
  const bool space2_is_ao = aoidxreg->value_exists( this->space2() );
  const bool space3_is_ao = aoidxreg->value_exists( this->space3() );
  MPQC_ASSERT(!space2_is_ao);
  MPQC_ASSERT(space3_is_ao);
  MPQC_ASSERT( aoidxreg->value_exists(iqR_tform->space2()) );
}

TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR::TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR(StateIn& si) :
  TwoBodyThreeCenterMOIntsTransform_ijR(si) {
  MPQC_ASSERT(false);
  init_vars();
}

void
TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR::save_data_state(StateOut& so) {
  TwoBodyThreeCenterMOIntsTransform::save_data_state(so);
}

int
TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR::compute_transform_batchsize(size_t mem_static, int rank_R)
{
  // assuming can hold qR and jR matrices in-core
  return rank_R;
}

distsize_t
TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR::compute_transform_dynamic_memory(int batchsize) const
{
  const bool using_memgrp = (ints_method_ == MOIntsTransform::StoreMethod::mem_only) ||
    (!factory()->hints().data_persistent() && (ints_method_ == MOIntsTransform::StoreMethod::mem_posix ||
                                               ints_method_ == MOIntsTransform::StoreMethod::mem_mpi));
  const unsigned int rank1 = space1()->rank();
  const unsigned int rank2 = space2()->rank();
  const unsigned int rank3 = space3()->rank();
  if (batchsize == -1)
    batchsize = rank3;
  // can only do 1-pass transformation
  MPQC_ASSERT(batchsize == rank3);

  const distsize_t result = using_memgrp ?
                              num_te_types() * rank1 * static_cast<distsize_t>(rank2 * rank3 * sizeof(double)) :
                              0;
  return result;
}

void
TwoBodyThreeCenterMOIntsTransform_ijR_using_iqR::compute() {

  // if all integrals are already computed, do nothing
  if (space3()->rank() == restart_orbital_)
    return;

  std::string tim_label("tbint_tform_ijR_using_iqR ");
  tim_label += this->name();
  Timer tim(tim_label);

  const int nproc = mem_->n();
  const int me = mem_->me();

  const int num_te_types = this->num_te_types();
  const int n1 = this->space1()->rank();
  const int n2 = this->space2()->rank();
  const int n3 = this->space3()->rank();
  const int n23 = n2*n3;
  const int nbasis2 = this->space2()->basis()->nbasis();

  Ref<DistArray4_MemoryGrp> ints_acc_cast; ints_acc_cast << ints_acc_;
  const bool need_memgrp = ints_acc_cast;
  if (need_memgrp) {
    const distsize_t ijR_globalsize = (((static_cast<distsize_t>(n1))*n23)*num_te_types)*sizeof(double);
    const int ni_local = (n1 + nproc - 1)/ nproc;
    const size_t memgrp_blocksize = (static_cast<size_t>(n23))*sizeof(double);
    const size_t ijR_localsize = num_te_types * ni_local * memgrp_blocksize;
    this->alloc_mem(ijR_localsize);
    memset(mem_->localdata(), 0, ijR_localsize);
  }

  iqR_tform_->compute();
  Ref<DistArray4> iqR_data = iqR_tform_->ints_acc();
  Ref<DistArray4> ijR_data = ints_acc_;
  iqR_data->activate();
  ijR_data->activate();
  {
    std::vector<int> writers;
    const int nwriters = ijR_data->tasks_with_access(writers);

    if (ijR_data->has_access(me)) {

      // scratch for holding transformed vectors
      double* C_jR = new double[n2 * n3];

      // AO->MO coefficents, rows are AOs
      double* tform = new double[nbasis2 * n2];
      space2()->coefs().convert(tform);

      int task_count = 0;
      for(int te_type=0; te_type<num_te_types; ++te_type) {
        for(int i = 0; i < n1; ++i, ++task_count) {

          // distribute work in round robin
          if (not iqR_data->is_local(0,i))
            continue;

          const double* C_qR = iqR_data->retrieve_pair_block(0, i, te_type);

          // transform
          C_DGEMM('t', 'n',
                  n2, n3, nbasis2,
                  1.0, tform, n2,
                  C_qR, n3,
                  0.0, C_jR, n3);

          // write
          ijR_data->store_pair_block(0, i, te_type, C_jR);

          // release this block
          iqR_data->release_pair_block(0, i, te_type);

        }
      }

      delete[] C_jR;
      delete[] tform;
    }

  }
  if (iqR_data->data_persistent()) iqR_data->deactivate();
  if (ijR_data->data_persistent()) ijR_data->deactivate();

  mem_->sync();

  restart_orbital_ = space3()->rank();
  tim.exit();
  ExEnv::out0() << indent << "Built TwoBodyMOIntsTransform_ijR_from_iqR: name = " << this->name() << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
