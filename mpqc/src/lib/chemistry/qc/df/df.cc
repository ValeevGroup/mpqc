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
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/fockbuilder.h>

#define USE_KERNEL_INVERSE 0

using namespace sc;

ClassDesc DensityFitting::class_desc_(
  typeid(DensityFitting),"DensityFitting",1,
  "virtual public SavableState",
  0, 0, create<DensityFitting>);

DensityFitting::~DensityFitting() {}

DensityFitting::DensityFitting(const Ref<MOIntsRuntime>& mointsruntime,
                               const std::string& kernel_key,
                               SolveMethod solver,
                               const Ref<OrbitalSpace>& space1,
                               const Ref<OrbitalSpace>& space2,
                               const Ref<GaussianBasisSet>& fitting_basis) :
                                 runtime_(mointsruntime),
                                 space1_(space1),
                                 space2_(space2),
                                 fbasis_(fitting_basis),
                                 kernel_key_(kernel_key),
                                 solver_(solver)
                                 {
  // only Coulomb fitting is supported at the moment
  if (kernel_key_ != std::string("1/r_{12}"))
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
  si.get(kernel_key_);
  int solvemethod; si.get(solvemethod); solver_ = static_cast<SolveMethod>(solvemethod);
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
  so.put(kernel_key_);
  so.put((int)solver_);
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
  const Ref<AOSpaceRegistry>& aoidxreg = this->runtime()->factory()->ao_registry();
  const std::string name = ParsedDensityFittingKey::key(space1_->id(),
                                                        space2_->id(),
                                                        aoidxreg->value(fbasis_)->id());
  std::string tim_label("DensityFitting ");
  tim_label += name;
  Timer tim(tim_label);

  const Ref<Integral>& integral = this->integral();

  // compute cC_ first
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

    if (cC_->data_persistent()) cC_->deactivate();
    runtime()->factory()->mem()->sync();
  }
#endif

  // compute the kernel second
  {
    const std::string kernel_key = ParsedTwoBodyTwoCenterIntKey::key(aoidxreg->value(fbasis_)->id(),
                                                                     aoidxreg->value(fbasis_)->id(),
                                                                     "ERI", "");
    RefSCMatrix kernel_rect = runtime()->runtime_2c()->get(kernel_key);
    const int nfbs = kernel_.dim().n();
    kernel_.assign_subblock(kernel_rect, 0, nfbs-1, 0, nfbs-1);
  }
  //kernel_.print("DensityFitting::kernel");

  // solve the linear system C_ * kernel = cC_
  // directly
  // OR
  // as C_ = cC_ * kernel^-1
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
      std::vector<double> C_jR(n2 * n3);

      std::vector<double> kernel_i; // only needed for inverse method
      std::vector<double> kernel_packed;  // only needed for factorized methods
      std::vector<double> kernel_factorized;
      std::vector<int> ipiv;

      // factorize or invert kernel
      switch (solver_) {
        case SolveMethod_InverseBunchKaufman:
        case SolveMethod_InverseCholesky:
        {
          RefSymmSCMatrix kernel_i_mat = kernel_.clone();
          kernel_i_mat.assign(kernel_);
          if (solver_ == SolveMethod_InverseBunchKaufman)
            exp::lapack_invert_symmnondef(kernel_i_mat, 1e10);
          if (solver_ == SolveMethod_InverseCholesky)
            exp::lapack_invert_symmposdef(kernel_i_mat, 1e10);

          // convert kernel_i to dense rectangular form
          kernel_i.resize(n3 * n3);
          int rc = 0;
          for (int r = 0; r < n3; ++r) {
            for (int c = 0; c < n3; ++c, ++rc) {
              kernel_i[rc] = kernel_i_mat.get_element(r, c);
            }
          }
          kernel_i_mat = 0;
        }
        break;

        case SolveMethod_Cholesky:
        case SolveMethod_RefinedCholesky:
        {
          // convert kernel_ to a packed upper-triangle form
          kernel_packed.resize(n3 * (n3 + 1) / 2);
          kernel_->convert(&(kernel_packed[0]));
          // factorize kernel_ using diagonal pivoting from LAPACK's DSPTRF
          kernel_factorized.resize(n3 * (n3 + 1) / 2);
          sc::exp::lapack_cholesky_symmposdef(kernel_,
                                              &(kernel_factorized[0]),
                                              1e10);
        }
        break;

        case SolveMethod_BunchKaufman:
        case SolveMethod_RefinedBunchKaufman:
        {
          // convert kernel_ to a packed upper-triangle form
          kernel_packed.resize(n3 * (n3 + 1) / 2);
          kernel_->convert(&(kernel_packed[0]));
          // factorize kernel_ using diagonal pivoting from LAPACK's DSPTRF
          kernel_factorized.resize(n3 * (n3 + 1) / 2);
          ipiv.resize(n3);
          sc::exp::lapack_dpf_symmnondef(kernel_, &(kernel_factorized[0]),
                                         &(ipiv[0]), 1e10);
        }
        break;

        default:
          throw ProgrammingError("unknown solve method", __FILE__, __LINE__, class_desc());
      }

      for (int i = 0; i < n1; ++i) {

        // distribute work in round robin
        if (i % nwriters != writers[me])
          continue;

        const double* cC_jR = cC_->retrieve_pair_block(0, i, TwoBodyOper::eri);

        bool refine_solution = true;
        // solve the linear system
        switch (solver_) {
          case SolveMethod_InverseCholesky:
          case SolveMethod_InverseBunchKaufman:
          {
            C_DGEMM('n', 'n', n2, n3, n3, 1.0, cC_jR, n3, &(kernel_i[0]), n3,
                    0.0, &(C_jR[0]), n3);
          }
          break;

          case SolveMethod_Cholesky:
            refine_solution = false;
          case SolveMethod_RefinedCholesky:
          {
            sc::exp::lapack_linsolv_cholesky_symmposdef(&(kernel_packed[0]), n3,
                                                        &(kernel_factorized[0]),
                                                        &(C_jR[0]), cC_jR, n2,
                                                        refine_solution);
          }
          break;

          case SolveMethod_BunchKaufman:
            refine_solution = false;
          case SolveMethod_RefinedBunchKaufman:
          {
            //sc::exp::lapack_linsolv_symmnondef(&(kernel_packed[0]), n3, &(C_jR[0]), cC_jR, n2);
            sc::exp::lapack_linsolv_dpf_symmnondef(&(kernel_packed[0]), n3,
                                                   &(kernel_factorized[0]),
                                                   &(ipiv[0]), &(C_jR[0]), cC_jR, n2,
                                                   refine_solution);
          }
          break;

          default:
            throw ProgrammingError("unknown solve method", __FILE__, __LINE__, class_desc());
        }

        // write
        C_->store_pair_block(0, i, 0, &(C_jR[0]));

        // release this block
        cC_->release_pair_block(0, i, TwoBodyOper::eri);

      }

    }

  }
  if (C_->data_persistent()) C_->deactivate();
  if (cC_->data_persistent()) cC_->deactivate();

  runtime()->factory()->mem()->sync();

  tim.exit();
  ExEnv::out0() << indent << "Built DensityFitting: name = " << name << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc TransformedDensityFitting::class_desc_(
  typeid(TransformedDensityFitting),"TransformedDensityFitting",1,
  "public DensityFitting",
  0, 0, create<TransformedDensityFitting>);

TransformedDensityFitting::~TransformedDensityFitting() {}

TransformedDensityFitting::TransformedDensityFitting(const Ref<MOIntsRuntime>& rtime,
                                                     const std::string& kernel_key,
                                                     DensityFitting::SolveMethod solver,
                                                     const Ref<OrbitalSpace>& space1,
                                                     const Ref<OrbitalSpace>& space2,
                                                     const Ref<GaussianBasisSet>& fitting_basis,
                                                     const Ref<DistArray4>& mo1_ao2_df) :
                       DensityFitting(rtime, kernel_key, solver, space1, space2, fitting_basis),
                       mo1_ao2_df_(mo1_ao2_df)
{
  // make sure that mo1_ao2_df is a valid input
  assert(mo1_ao2_df->ni() == 1 &&
         mo1_ao2_df->nj() == space1->rank() &&
         mo1_ao2_df->nx() == space2->basis()->nbasis() &&
         mo1_ao2_df->ny() == fitting_basis->nbasis());
}

TransformedDensityFitting::TransformedDensityFitting(StateIn& si) :
  DensityFitting(si) {
  mo1_ao2_df_ << SavableState::restore_state(si);
}

void
TransformedDensityFitting::save_data_state(StateOut& so) {
  SavableState::save_state(mo1_ao2_df_.pointer(),so);
}

void
TransformedDensityFitting::compute()
{
  if (C_.nonnull()) // nothing to compute then
    return;
  Ref<AOSpaceRegistry> aoidxreg = this->runtime()->factory()->ao_registry();
  const std::string name = ParsedDensityFittingKey::key(this->space1()->id(),
                                                        this->space2()->id(),
                                                        aoidxreg->value(this->fbasis())->id());
  std::string tim_label("TransformedDensityFitting ");
  tim_label += name;
  Timer tim(tim_label);

  const int n1 = space1()->rank();
  const int n2 = space2()->rank();
  const int n2_ao = space2()->basis()->nbasis();
  const int n3 = fbasis()->nbasis();

  Ref<DistArray4> C_ao = mo1_ao2_df_;
  DistArray4Dimensions Cdims(1,
                             1, n1, n2, n3,
                             DistArray4Storage_XY);
  C_ = C_ao->clone(Cdims);
  C_->activate();
  C_ao->activate();
  {
    const int me = runtime()->factory()->msg()->me();
    std::vector<int> writers;
    const int nwriters = C_->tasks_with_access(writers);

    if (C_->has_access(me)) {

      // scratch for holding transformed vectors
      double* C_jR = new double[n2 * n3];

      // AO->MO coefficents, rows are AOs
      double* tform = new double[n2_ao * n2];
      space2()->coefs().convert(tform);

      for (int i = 0; i < n1; ++i) {

        // distribute work in round robin
        if (i % nwriters != writers[me])
          continue;

        const double* C_qR = C_ao->retrieve_pair_block(0, i, 0);

        // transform
        C_DGEMM('t', 'n',
                n2, n3, n2_ao,
                1.0, tform, n2,
                C_qR, n3,
                0.0, C_jR, n3);

        // write
        C_->store_pair_block(0, i, 0, C_jR);

        // release this block
        C_ao->release_pair_block(0, i, 0);

      }

      delete[] C_jR;
      delete[] tform;
    }

  }
  if (C_->data_persistent()) C_->deactivate();
  if (C_ao->data_persistent()) C_ao->deactivate();

  runtime()->factory()->mem()->sync();

  tim.exit();
  ExEnv::out0() << indent << "Built TransformedDensityFitting: name = " << name << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc PermutedDensityFitting::class_desc_(
  typeid(PermutedDensityFitting),"PermutedDensityFitting",1,
  "public DensityFitting",
  0, 0, create<PermutedDensityFitting>);

PermutedDensityFitting::~PermutedDensityFitting() {}

PermutedDensityFitting::PermutedDensityFitting(const Ref<MOIntsRuntime>& rtime,
                                               const std::string& kernel_key,
                                               DensityFitting::SolveMethod solver,
                                               const Ref<OrbitalSpace>& space1,
                                               const Ref<OrbitalSpace>& space2,
                                               const Ref<GaussianBasisSet>& fitting_basis,
                                               const Ref<DistArray4>& df21) :
                 DensityFitting(rtime, kernel_key, solver, space1, space2, fitting_basis),
                 df21_(df21)
{
  // make sure that mo1_ao2_df is a valid input
  assert(df21->ni() == 1 &&
         df21->nj() == space2->rank() &&
         df21->nx() == space1->rank() &&
         df21->ny() == fitting_basis->nbasis());
}

PermutedDensityFitting::PermutedDensityFitting(const Ref<DensityFitting>& df21) :
                       DensityFitting(df21->runtime(), df21->kernel_key(),
                                      df21->solver(),
                                      df21->space2(), df21->space1(),   // swapped!
                                      df21->fbasis()),
                       df21_(df21->C())
{}

PermutedDensityFitting::PermutedDensityFitting(StateIn& si) :
  DensityFitting(si) {
  df21_ << SavableState::restore_state(si);
}

void
PermutedDensityFitting::save_data_state(StateOut& so) {
  SavableState::save_state(df21_.pointer(),so);
}

void
PermutedDensityFitting::compute()
{
  if (C_.nonnull()) // nothing to compute then
    return;

  Ref<AOSpaceRegistry> aoidxreg = this->runtime()->factory()->ao_registry();
  const std::string name = ParsedDensityFittingKey::key(this->space1()->id(),
                                                        this->space2()->id(),
                                                        aoidxreg->value(this->fbasis())->id());
  std::string tim_label("PermutedDensityFitting ");
  tim_label += name;
  Timer tim(tim_label);

  Ref<DistArray4> C21 = df21_;
  C21->activate();
  C_ = permute23(C21, 1000000000);
  if (C_->data_persistent()) C_->deactivate();
  if (C21->data_persistent()) C21->deactivate();

  runtime()->factory()->mem()->sync();

  tim.exit();
  ExEnv::out0() << indent << "Built PermutedDensityFitting: name = " << name << std::endl;
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
