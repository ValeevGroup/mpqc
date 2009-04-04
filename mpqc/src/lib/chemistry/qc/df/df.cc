//
// df.cc
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
#include <chemistry/qc/df/df.h>
#include <util/state/statein.h>
#include <chemistry/qc/mbptr12/orbitalspace.h>
#include <chemistry/qc/mbptr12/moints_runtime.h>
#include <chemistry/qc/mbptr12/distarray4.h>
#include <chemistry/qc/mbptr12/distarray4_memgrp.h>
#include <chemistry/qc/mbptr12/transform_ijR.h>
#include <chemistry/qc/mbptr12/svd.h>
#include <chemistry/qc/mbptr12/fockbuilder.h>

using namespace sc;

static ClassDesc DensityFitting_cd(
  typeid(DensityFitting),"DensityFitting",1,
  "virtual public SavableState",
  0, 0, create<DensityFitting>);

DensityFitting::~DensityFitting() {}

DensityFitting::DensityFitting(const Ref<MOIntsRuntime>& mointsruntime,
                               const std::string& kernel_key,
                               const Ref<OrbitalSpace>& space1,
                               const Ref<OrbitalSpace>& space2,
                               const Ref<GaussianBasisSet>& fitting_basis) :
                                 runtime_(mointsruntime),
                                 space1_(space1),
                                 space2_(space2),
                                 fbasis_(fitting_basis)
                                 {
  // only Coulomb fitting is supported at the moment
  if (kernel_key != std::string("1/r_{12}"))
    throw FeatureNotImplemented("Non-Coulomb fitting kernels are not supported",__FILE__,__LINE__);

  // fitting dimension
  RefSCDimension fdim = new SCDimension(fbasis_->nbasis(), "");
  Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
  kernel_ = kit->symmmatrix(fdim);
}

DensityFitting::DensityFitting(StateIn& si) :
  SavableState(si) {
  runtime_ << SavableState::restore_state(si);
  fbasis_ << SavableState::restore_state(si);
  space1_ << SavableState::restore_state(si);
  space2_ << SavableState::restore_state(si);
  cC_ << SavableState::restore_state(si);
  C_ << SavableState::restore_state(si);

  int count;
  detail::FromStateIn<RefSymmSCMatrix>::get(kernel_,si,count);
}

void
DensityFitting::save_data_state(StateOut& so) {
  SavableState::save_state(runtime_.pointer(),so);
  SavableState::save_state(fbasis_.pointer(),so);
  SavableState::save_state(space1_.pointer(),so);
  SavableState::save_state(space2_.pointer(),so);
  SavableState::save_state(cC_.pointer(),so);
  SavableState::save_state(C_.pointer(),so);

  int count;
  detail::ToStateOut<RefSymmSCMatrix>::put(kernel_,so,count);
}

RefSCDimension
DensityFitting::product_dimension() const {
  return new SCDimension(space1_->rank() * space2_->rank(), "");
}

void
DensityFitting::compute()
{
  if (cC_.nonnull() && kernel_ && C_.nonnull()) // nothing to compute then
    return;
  const std::string name = ParsedDensityFittingKey::key(space1_->id(),
                                                        space2_->id(),
                                                        AOSpaceRegistry::instance()->value(fbasis_)->id());
  std::string tim_label("DensityFitting ");
  tim_label += name;
  Timer tim(tim_label);

  const Ref<Integral>& integral = this->integral();

  // compute cC_ first
  const Ref<AOSpaceRegistry>& aoidxreg = AOSpaceRegistry::instance();
  // TODO need non-Coulomb fitting kernels?
  const std::string cC_key = ParsedTwoBodyThreeCenterIntKey::key(space1_->id(),
                                                                 aoidxreg->value(fbasis_)->id(),
                                                                 space2_->id(),
                                                                 "ERI", "");
  Ref<TwoBodyThreeCenterMOIntsTransform> tform = runtime()->runtime_3c()->get(cC_key);
  tform->compute();
  cC_ = tform->ints_acc();

#if 1
  {
    cC_->activate();
    const int ni = space1_->rank();
    const int nj = space2_->rank();
    const int nR = fbasis_->nbasis();
    RefSCMatrix cC = kernel_.kit()->matrix(new SCDimension(nR),
                                           this->product_dimension());
    cC.assign(0.0);

    std::vector<int> readers;
    const int nreaders = cC_->tasks_with_access(readers);
    const int me = runtime()->factory()->msg()->me();

    if (cC_->has_access(me)) {
      RefSCMatrix cC_jR = kernel_.kit()->matrix(new SCDimension(nj),
                                                new SCDimension(nR));
      for (int i = 0; i < ni; ++i) {

        if (i % nreaders != readers[me])
          continue;

        const double* cC_jR_buf = cC_->retrieve_pair_block(0, i,
                                                           TwoBodyOper::eri);
        cC_jR.assign(cC_jR_buf);
        cC.assign_subblock(cC_jR.t(), 0, nR - 1, i * nj, (i + 1) * nj - 1);
      }
    }

    // broadcast to all tasks

    // print out the result
    //cC.print("DensityFitting: cC");

    runtime()->factory()->mem()->sync();
  }
#endif

  // compute the kernel second
  {
    const Ref<AOSpaceRegistry>& aoidxreg = AOSpaceRegistry::instance();
    const std::string kernel_key = ParsedTwoBodyTwoCenterIntKey::key(aoidxreg->value(fbasis_)->id(),
                                                                     aoidxreg->value(fbasis_)->id(),
                                                                     "ERI", "");
    RefSCMatrix kernel_rect = runtime()->runtime_2c()->get(kernel_key);
    const int nfbs = kernel_.dim().n();
    kernel_.assign_subblock(kernel_rect, 0, nfbs-1, 0, nfbs-1);
  }
  //kernel_.print("DensityFitting::kernel");

  // solve the linear system C_ * kernel = cC_
  //
  // create C_ as a clone cC_
  // parallelize over all tasks that can read from cC_
  // loop over i1
  //   decide if I will handle this i
  //   fetch jR block
  //   solve the system
  //   save jR block of the solution to C_
  // end of loop over i
  C_ = cC_->clone();
  C_->activate();
  cC_->activate();
  {
    const int me = runtime()->factory()->msg()->me();
    std::vector<int> writers;
    const int nwriters = C_->tasks_with_access(writers);

    if (C_->has_access(me)) {

      const int n1 = space1_->rank();
      const int n2 = space2_->rank();
      const int n3 = fbasis_->nbasis();
      // scratch for holding solution vectors
      double* C_jR = new double[n2 * n3];
      // convert kernel_ to a packed upper-triangle form
      double* kernel_packed = new double[n3 * (n3 + 1) / 2];
      kernel_->convert(kernel_packed);

      for (int i = 0; i < n1; ++i) {

        // distribute work in round robin
        if (i % nwriters != writers[me])
          continue;

        const double* cC_jR = cC_->retrieve_pair_block(0, i, TwoBodyOper::eri);

        // solve the linear system
        sc::exp::lapack_linsolv_symmnondef(kernel_packed, n3, C_jR, cC_jR, n2);

        // write
        C_->store_pair_block(0, i, TwoBodyOper::eri, C_jR);

        // release this block
        cC_->release_pair_block(0, i, TwoBodyOper::eri);

      }

      delete[] C_jR;
      delete[] kernel_packed;
    }

  }
  if (C_->data_persistent()) C_->deactivate();
  if (cC_->data_persistent()) cC_->deactivate();

  runtime()->factory()->mem()->sync();

  tim.exit();
  ExEnv::out0() << indent << "Built DensityFitting: name = " << name << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

#if 0
static ClassDesc testBasisProductDecomposition_cd(
  typeid(test::BasisProductDecomposition),"BasisProductDecomposition",1,
  "virtual public SavableState",
  0, 0, 0);
#endif

test::BasisProductDecomposition::~BasisProductDecomposition() {}

test::BasisProductDecomposition::BasisProductDecomposition(const Ref<Integral>& integral,
                                                     const Ref<GaussianBasisSet>& basis1,
                                                     const Ref<GaussianBasisSet>& basis2) :
                                                       integral_(integral) {
  basis_[0] = basis1;
  basis_[1] = basis2;

  // compute product dim
  const int nbasis1 = basis1->nbasis();
  const int nbasis2 = basis2->nbasis();
  pdim_ = new SCDimension(nbasis1 * nbasis2, "");

}

test::BasisProductDecomposition::BasisProductDecomposition(StateIn& si) :
  SavableState(si) {
  integral_ << SavableState::restore_state(si);
  basis_[0] << SavableState::restore_state(si);
  basis_[1] << SavableState::restore_state(si);
  pdim_ << SavableState::restore_state(si);
}

void
test::BasisProductDecomposition::save_data_state(StateOut& so) {
  SavableState::save_state(integral_.pointer(),so);
  SavableState::save_state(basis_[0].pointer(),so);
  SavableState::save_state(basis_[1].pointer(),so);
  SavableState::save_state(pdim_.pointer(),so);
}

/////////////////////////////////////////////////////////////////////////////

static ClassDesc testDensityFitting_cd(
  typeid(test::DensityFitting),"test::DensityFitting",1,
  "public BasisProductDecomposition",
  0, 0, create<test::DensityFitting>);

test::DensityFitting::~DensityFitting() {}

test::DensityFitting::DensityFitting(const Ref<Integral>& integral,
                               const Ref<GaussianBasisSet>& basis1,
                               const Ref<GaussianBasisSet>& basis2,
                               const Ref<GaussianBasisSet>& fitting_basis) :
                                 BasisProductDecomposition(integral,basis1,basis2),
                                 fbasis_(fitting_basis) {

  // fitting dimension
  RefSCDimension fdim = new SCDimension(fbasis_->nbasis(), "");

  Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
  C_ = kit->matrix(fdim, this->product_dimension());
  cC_ = kit->matrix(fdim, this->product_dimension());
  kernel_ = kit->symmmatrix(fdim);
}

test::DensityFitting::DensityFitting(StateIn& si) :
  BasisProductDecomposition(si) {
  fbasis_ << SavableState::restore_state(si);

  int count;
  detail::FromStateIn<RefSCMatrix>::get(C_,si,count);
  detail::FromStateIn<RefSymmSCMatrix>::get(kernel_,si,count);
  detail::FromStateIn<RefSCMatrix>::get(cC_,si,count);
}

void
test::DensityFitting::save_data_state(StateOut& so) {
  SavableState::save_state(fbasis_.pointer(),so);

  int count;
  detail::ToStateOut<RefSCMatrix>::put(C_,so,count);
  detail::ToStateOut<RefSymmSCMatrix>::put(kernel_,so,count);
  detail::ToStateOut<RefSCMatrix>::put(cC_,so,count);
}

void
test::DensityFitting::compute()
{
  const Ref<Integral>& integral = this->integral();

  // compute the kernel first
  {
    const Ref<GaussianBasisSet>& b = fbasis_;
    integral->set_basis(b, b);
    Ref<TwoBodyTwoCenterInt> coulomb2int = integral->electron_repulsion2();
    const double* buffer = coulomb2int->buffer();
    for (int s1 = 0; s1 < b->nshell(); ++s1) {
      const int s1offset = b->shell_to_function(s1);
      const int nf1 = b->shell(s1).nfunction();
      for (int s2 = 0; s2 <= s1; ++s2) {
        const int s2offset = b->shell_to_function(s2);
        const int nf2 = b->shell(s2).nfunction();

        // compute shell doublet
        coulomb2int->compute_shell(s1, s2);

        // copy buffer into kernel_
        const double* bufptr = buffer;
        for(int f1=0; f1<nf1; ++f1) {
          for(int f2=0; f2<nf2; ++f2, ++bufptr) {
            kernel_.set_element(f1+s1offset, f2+s2offset, *bufptr);
          }
        }

      }
    }
  }

  kernel_.print("test::DensityFitting::kernel");

  // compute the conjugate coefficient matrix
  {
    const Ref<GaussianBasisSet>& b1 = this->basis(0);
    const Ref<GaussianBasisSet>& b2 = this->basis(1);
    const Ref<GaussianBasisSet>& b3 = fbasis_;
    const int nbasis2 = b2->nbasis();
    integral->set_basis(b1, b2, b3);
    Ref<TwoBodyThreeCenterInt> coulomb3int = integral->electron_repulsion3();
    const double* buffer = coulomb3int->buffer();
    for (int s1 = 0; s1 < b1->nshell(); ++s1) {
      const int s1offset = b1->shell_to_function(s1);
      const int nf1 = b1->shell(s1).nfunction();
      for (int s2 = 0; s2 < b2->nshell(); ++s2) {
        const int s2offset = b2->shell_to_function(s2);
        const int nf2 = b2->shell(s2).nfunction();

        for (int s3 = 0; s3 < b3->nshell(); ++s3) {
          const int s3offset = b3->shell_to_function(s3);
          const int nf3 = b3->shell(s3).nfunction();

          // compute shell triplet
          coulomb3int->compute_shell(s1, s2, s3);

          // copy buffer into kernel_
          const double* bufptr = buffer;
          for(int f1=0; f1<nf1; ++f1) {
            const int s12offset = (s1offset+f1) * nbasis2 + s2offset;
            for(int f2=0; f2<nf2; ++f2) {
              for(int f3=0; f3<nf3; ++f3, ++bufptr) {
                cC_.set_element(f3+s3offset, f2+s12offset, *bufptr);
              }
            }
          }

        }
      }
    }
  }

#if 1
  cC_->print("test::DensityFitting: cC");
#endif

  // solve the linear system
  // TODO parallelize. Can do trivially by dividing RHS vectors into subsets
  sc::exp::lapack_linsolv_symmnondef(kernel_, C_, cC_);

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
