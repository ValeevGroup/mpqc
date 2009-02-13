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
#include <chemistry/qc/mbptr12/svd.h>

using namespace sc;
using namespace sc::test;

static ClassDesc BasisProductDecomposition_cd(
  typeid(BasisProductDecomposition),"BasisProductDecomposition",1,
  "virtual public SavableState",
  0, 0, 0);

BasisProductDecomposition::~BasisProductDecomposition() {}

BasisProductDecomposition::BasisProductDecomposition(const Ref<Integral>& integral,
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

BasisProductDecomposition::BasisProductDecomposition(StateIn& si) :
  SavableState(si) {
  integral_ << SavableState::restore_state(si);
  basis_[0] << SavableState::restore_state(si);
  basis_[1] << SavableState::restore_state(si);
  pdim_ << SavableState::restore_state(si);
}

void
BasisProductDecomposition::save_data_state(StateOut& so) {
  SavableState::save_state(integral_.pointer(),so);
  SavableState::save_state(basis_[0].pointer(),so);
  SavableState::save_state(basis_[1].pointer(),so);
  SavableState::save_state(pdim_.pointer(),so);
}

/////////////////////////////////////////////////////////////////////////////

static ClassDesc DensityFitting_cd(
  typeid(DensityFitting),"DensityFitting",1,
  "public BasisProductDecomposition",
  0, 0, create<DensityFitting>);

DensityFitting::~DensityFitting() {}

DensityFitting::DensityFitting(const Ref<Integral>& integral,
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

DensityFitting::DensityFitting(StateIn& si) :
  BasisProductDecomposition(si) {
  fbasis_ << SavableState::restore_state(si);

  int count;
  detail::FromStateIn<RefSCMatrix>::get(C_,si,count);
  detail::FromStateIn<RefSymmSCMatrix>::get(kernel_,si,count);
  detail::FromStateIn<RefSCMatrix>::get(cC_,si,count);
}

void
DensityFitting::save_data_state(StateOut& so) {
  SavableState::save_state(fbasis_.pointer(),so);

  int count;
  detail::ToStateOut<RefSCMatrix>::put(C_,so,count);
  detail::ToStateOut<RefSymmSCMatrix>::put(kernel_,so,count);
  detail::ToStateOut<RefSCMatrix>::put(cC_,so,count);
}

void
DensityFitting::compute()
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

  // solve the linear system
  // TODO parallelize. Can do trivially by dividing RHS vectors into subsets
  sc::exp::lapack_linsolv_symmnondef(kernel_, C_, cC_);

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
