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

// includes go here
#include <chemistry/qc/lcao/df.h>
#include <util/state/statein.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <chemistry/qc/lcao/moints_runtime.h>
#include <math/distarray4/distarray4.h>
#include <math/distarray4/distarray4_memgrp.h>
#include <chemistry/qc/lcao/transform_ijR.h>
#include <math/scmat/svd.h>
#include <math/scmat/blas.h>
#include <math/optimize/gaussianfit.h>
#include <math/optimize/gaussianfit.timpl.h>
#include <chemistry/qc/lcao/fockbuilder.h>
#include <Eigen/Dense>
#include <util/misc/sharedptr.h>

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
  if (not DensityFittingParams::valid_kernel(kernel_key_))
    throw FeatureNotImplemented("Invalid or unsupported fitting kernel",__FILE__,__LINE__);

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
  FromStateIn(kernel_,si,count);
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
  ToStateOut(kernel_,so,count);
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
  const std::string fbasis_space_id = aoidxreg->value(fbasis_)->id();
  const std::string name = ParsedDensityFittingKey::key(space1_->id(),
                                                        space2_->id(),
                                                        fbasis_space_id,
                                                        kernel_key_);

  std::string tim_label("DensityFitting ");
  tim_label += name;
  Timer tim(tim_label);

  const Ref<Integral>& integral = this->integral();

  // this assumes there is only one nontrivial parameter
  std::string params_key = runtime()->factory()->df_info()->params()->intparams_key();
  TwoBodyOper::type kernel_oper;
  TwoBodyOperSet::type operset;
  if (kernel_key_ == "coulomb") {
    operset = TwoBodyOperSet::ERI;
    kernel_oper = TwoBodyOper::eri;
    params_key = "";
  }
  else if (kernel_key_ == "delta") {
    operset = TwoBodyOperSet::DeltaFunction;
    kernel_oper = TwoBodyOper::delta;
    params_key = "";
  }
  else if (kernel_key_.find("exp") != std::string::npos) {
    operset = TwoBodyOperSet::G12NC;
    kernel_oper = TwoBodyOper::r12_0_g12;
  }
  const std::string operset_key = TwoBodyOperSetDescr::instance(operset)->key();
  const unsigned int ints_type_idx = TwoBodyOperSetDescr::instance(operset)->opertype(kernel_oper);

  // compute cC_ first
  std::string cC_key;
  cC_key = ParsedTwoBodyThreeCenterIntKey::key(space1_->id(),
                                               aoidxreg->value(fbasis_)->id(),
                                               space2_->id(),
                                               operset_key, params_key);
  Ref<TwoBodyThreeCenterMOIntsTransform> tform = runtime()->runtime_3c()->get(cC_key);
  tform->compute();
  cC_ = tform->ints_acc();

  // compute the kernel second
  {
    std::string kernel_ints_key;

    kernel_ints_key = ParsedTwoBodyTwoCenterIntKey::key(aoidxreg->value(fbasis_)->id(),
                                                        aoidxreg->value(fbasis_)->id(),
                                                        operset_key, params_key);
    RefSCMatrix kernel_rect = runtime()->runtime_2c()->get(kernel_ints_key);
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
  DistArray4Dimensions Cdims(1,
                             1, cC_->nj(), cC_->nx(), cC_->ny(),
                             DistArray4Storage_XY);
  C_ = cC_->clone(Cdims);
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
      std::vector<blasint> ipiv;
#if 0
      Eigen::Map<Eigen::MatrixXd> C_jR_map(&(C_jR[0]), n2, n3);
      typedef Eigen::HouseholderQR<Eigen::MatrixXd> HQR_type;
      typedef Eigen::ColPivHouseholderQR<Eigen::MatrixXd> CPHQR_type;
      typedef Eigen::FullPivHouseholderQR<Eigen::MatrixXd> FPHQR_type;
      std::unique_ptr<HQR_type> solverHQR = 0;
      std::unique_ptr<CPHQR_type> solverCPHQR = 0;
      std::unique_ptr<FPHQR_type> solverFPHQR = 0;
      std::unique_ptr<Eigen::MatrixXd> kernel_mat = 0;
      if(solver_ > 100) { // Eigen based solvers
        // TODO call lapack directly and stop using Eigen as an intermediary
        // TODO use a map to wrap pointer if matrix is local or replicated
        assert(kernel_.n()==n3);
        kernel_mat = std::unique_ptr<Eigen::MatrixXd>(new Eigen::MatrixXd(n3, n3));
        kernel_.convert2RefSCMat().convert(kernel_mat->data());
      }
#endif

      // factorize or invert kernel
      switch (solver_) {
        case SolveMethod_InverseBunchKaufman:
        case SolveMethod_InverseCholesky:
        {
          RefSymmSCMatrix kernel_i_mat;
          if (not runtime_->runtime_2c_inv()->key_exists(fbasis_space_id)) {
            RefSymmSCMatrix kernel_i_mat = kernel_.copy();
            if (solver_ == SolveMethod_InverseBunchKaufman)
              lapack_invert_symmnondef(kernel_i_mat, 1e10);
            if (solver_ == SolveMethod_InverseCholesky)
              lapack_invert_symmposdef(kernel_i_mat, 1e10);
            runtime_->runtime_2c_inv()->add(fbasis_space_id, kernel_i_mat);
          }
          kernel_i_mat = runtime_->runtime_2c_inv()->value(fbasis_space_id);

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
          sc::lapack_cholesky_symmposdef(kernel_,
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
          sc::lapack_dpf_symmnondef(kernel_, &(kernel_factorized[0]),
                                         &(ipiv[0]), 1e10);
        }
        break;
#if 0

        case SolveMethod_HouseholderQR:
        {
          solverHQR = std::unique_ptr<HQR_type>(new HQR_type(*kernel_mat));
        }
        break;

        case SolveMethod_ColPivHouseholderQR:
        {
          // Make an Eigen matrix
          solverCPHQR = std::unique_ptr<CPHQR_type>(new CPHQR_type(*kernel_mat));
        }
        break;

        case SolveMethod_FullPivHouseholderQR:
        {
          solverFPHQR = std::unique_ptr<FPHQR_type>(new FPHQR_type(*kernel_mat));
        }
        break;
#endif

        default:
          throw ProgrammingError("unknown solve method", __FILE__, __LINE__, class_desc());
      }

      for (int i = 0; i < n1; ++i) {

        // work on local blocks
        if (not cC_->is_local(0, i))
          continue;

        const double* cC_jR = cC_->retrieve_pair_block(0, i, ints_type_idx);
#if 0
        Eigen::Map<const Eigen::MatrixXd> cC_jR_map(cC_jR, n2, n3);
#endif

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
            sc::lapack_linsolv_cholesky_symmposdef(&(kernel_packed[0]), n3,
                                                        &(kernel_factorized[0]),
                                                        &(C_jR[0]), cC_jR, n2,
                                                        refine_solution);
          }
          break;

          case SolveMethod_BunchKaufman:
            refine_solution = false;
          case SolveMethod_RefinedBunchKaufman:
          {
            //sc::lapack_linsolv_symmnondef(&(kernel_packed[0]), n3, &(C_jR[0]), cC_jR, n2);
            sc::lapack_linsolv_dpf_symmnondef(&(kernel_packed[0]), n3,
                                                   &(kernel_factorized[0]),
                                                   &(ipiv[0]), &(C_jR[0]), cC_jR, n2,
                                                   refine_solution);
          }
          break;

#if 0 // Doesn't work
          case SolveMethod_HouseholderQR:
          {
            C_jR_map.transpose() = solverHQR->solve(cC_jR_map.transpose());
          }
          break;

          case SolveMethod_ColPivHouseholderQR:
          {
            C_jR_map.transpose() = solverCPHQR->solve(cC_jR_map.transpose());
          }
          break;

          case SolveMethod_FullPivHouseholderQR:
          {
            C_jR_map.transpose() = solverFPHQR->solve(cC_jR_map.transpose());
          }
          break;
#endif

          default:
            throw ProgrammingError("unknown solve method", __FILE__, __LINE__, class_desc());
        }

        // write
        C_->store_pair_block(0, i, 0, &(C_jR[0]));

        // release this block
        cC_->release_pair_block(0, i, ints_type_idx);

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
                                                        aoidxreg->value(this->fbasis())->id(),
                                                        kernel_key());
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
                                                        aoidxreg->value(this->fbasis())->id(),
                                                        kernel_key());
  std::string tim_label("PermutedDensityFitting ");
  tim_label += name;
  Timer tim(tim_label);

  Ref<DistArray4> C21 = df21_;
  C21->activate();
  C_ = permute23(C21);
  if (C_->data_persistent()) C_->deactivate();
  if (C21->data_persistent()) C21->deactivate();

  runtime()->factory()->mem()->sync();

  tim.exit();
  ExEnv::out0() << indent << "Built PermutedDensityFitting: name = " << name << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
