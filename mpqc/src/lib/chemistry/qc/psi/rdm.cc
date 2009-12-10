//
// rdm.cc
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

#include <chemistry/qc/psi/rdm.h>

using namespace sc;

ClassDesc
PsiRDMTwo::class_desc_(typeid(PsiRDMTwo),
                     "PsiRDMTwo",
                     1,               // version
                     "public RDM<Two>", // must match parent
                     0,               // change to create<PsiRDMTwo> if this class is DefaultConstructible
                     create<PsiRDMTwo>, // change to 0 if this class is not KeyValConstructible
                     create<PsiRDMTwo>  // change to 0 if this class is not StateInConstructible
                     );

PsiRDMTwo::PsiRDMTwo(const Ref<KeyVal>& kv) : RDM<Two>(kv) {
  psiwfn_ = require_dynamic_cast<PsiWavefunction*>(
        kv->describedclassvalue("wfn").pointer(),
        "PsiRDMTwo::PsiRDMTwo\n"
        );
}

PsiRDMTwo::PsiRDMTwo(StateIn& si) : RDM<Two>(si) {
  psiwfn_ << SavableState::restore_state(si);
  if (psiwfn_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

PsiRDMTwo::~PsiRDMTwo() {
}

void
PsiRDMTwo::save_data_state(StateOut& so) {
  SavableState::save_state(psiwfn_.pointer(), so);
}

Ref<PsiRDMTwo::cumulant_type> sc::PsiRDMTwo::cumulant() const
{
  return new PsiRDMCumulantTwo(const_cast<PsiRDMTwo*>(this));
}

void sc::PsiRDMTwo::release_block(SpinCase2 spin, size_t bra, double *blk) const
{
  throw "not yet implemented";
}

const double *sc::PsiRDMTwo::obtain_block(SpinCase2 spin, size_t bra) const
{
  throw "not yet implemented";
}

void sc::PsiRDMTwo::compute()
{
  psiwfn_->compute();
}

size_t sc::PsiRDMTwo::ndim(SpinCase2 spincase) const
{
  const int nmo = psiwfn_->oso_dimension().n();
  switch (spincase) {
    case AlphaAlpha:
    case BetaBeta:
      return nmo * (nmo-1) / 2;
    case AlphaBeta:
      return nmo * nmo;
  }
  assert(false);  // unreachable
}

RefSymmSCMatrix sc::PsiRDMTwo::scmat(SpinCase2 spin) const {

  if (scmat_[spin].nonnull()) return scmat_[spin];
  if (psiwfn_->spin_polarized() && spin == BetaBeta)
    return scmat(AlphaAlpha);

  Ref<PsiCorrWavefunction> psicorrwfn_;
  psicorrwfn_ << psiwfn_;
  if (psicorrwfn_.nonnull()) {
    scmat_[spin] = psicorrwfn_->twopdm_dirac(spin);
    return scmat_[spin];
  }

  Ref<PsiSCF> psiscfwfn_;
  psiscfwfn_ << psiwfn_;
  if (psiscfwfn_.nonnull()) {
    Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
    const SpinCase1 spin1 = case1(spin);
    const SpinCase1 spin2 = case2(spin);
    const RefSymmSCMatrix opdm1 = psiscfwfn_->mo_density(spin1);
    const RefSymmSCMatrix opdm2 = psiscfwfn_->mo_density(spin2);
    const int n = opdm1.n();

    RefSCDimension dim = new SCDimension( spin == AlphaBeta ? n*n : n*(n-1)/2 );
    RefSymmSCMatrix tpdm = kit->symmmatrix(dim); tpdm.assign(0.0);

    int b12 = 0;
    for(int b1=0; b1<n; ++b1) {
      const int b2fence = (spin == AlphaBeta) ? n : b1;
      for(int b2=0; b2< b2fence; ++b2, ++b12) {

        int k12 = 0;
        for(int k1=0; k1<n; ++k1) {
          const double gamma_b1_k1 = opdm1.get_element(b1,k1);
          const double gamma_b2_k1 = (spin != AlphaBeta) ? opdm1.get_element(b2,k1) : 0.0;

          const int k2fence = (spin == AlphaBeta) ? n : k1;
          for(int k2=0; k2< k2fence; ++k2, ++k12) {
            double value = gamma_b1_k1 * opdm2.get_element(b2,k2);
            if (spin != AlphaBeta) value -= opdm1.get_element(b1,k2) * gamma_b2_k1;
            tpdm.accumulate_element(b12,k12,value);
          }
        }
      }
    }
    scmat_[spin] = tpdm;
    return scmat_[spin];
  } // PsiSCF wavefunction

  throw FeatureNotImplemented("PsiRDMTwo::scmat() is not implemented for this PsiWavefunction",__FILE__,__LINE__);
}



/////////////////////

ClassDesc
PsiRDMCumulantTwo::class_desc_(typeid(PsiRDMCumulantTwo),
			       "PsiRDMCumulantTwo",
			       1,               // version
			       "public RDMCumulant<Two>", // must match parent
			       0,               // change to create<PsiRDMCumulantTwo> if this class is DefaultConstructible
			       0, // change to 0 if this class is not KeyValConstructible
			       create<PsiRDMCumulantTwo>  // change to 0 if this class is not StateInConstructible
    );

PsiRDMCumulantTwo::PsiRDMCumulantTwo(const Ref<PsiRDMTwo>& density) : density_(density), RDMCumulant<Two>(density) {
}

PsiRDMCumulantTwo::PsiRDMCumulantTwo(StateIn& si) : RDMCumulant<Two>(si) {
  density_ << SavableState::restore_state(si);
  if (density_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

PsiRDMCumulantTwo::~PsiRDMCumulantTwo() {
}

void
PsiRDMCumulantTwo::save_data_state(StateOut& so) {
  SavableState::save_state(density_.pointer(), so);
}

void sc::PsiRDMCumulantTwo::release_block(SpinCase2 spin, size_t bra, double *blk) const
{
  throw "not yet implemented";
}

const double *sc::PsiRDMCumulantTwo::obtain_block(SpinCase2 spin, size_t bra) const
{
  throw "not yet implemented";
}

void sc::PsiRDMCumulantTwo::compute()
{
  density_->compute();
}

RefSymmSCMatrix sc::PsiRDMCumulantTwo::scmat(SpinCase2 spin) const {

  if (scmat_[spin].nonnull()) return scmat_[spin];
  if (!density_->psiwfn()->spin_polarized() && spin == BetaBeta)
    return scmat(AlphaAlpha);

  const SpinCase1 spin1 = case1(spin);
  const SpinCase1 spin2 = case2(spin);
  const RefSymmSCMatrix opdm1 = density_->psiwfn()->mo_density(spin1);
  const RefSymmSCMatrix opdm2 = density_->psiwfn()->mo_density(spin2);
  const RefSymmSCMatrix tpdm = density_->scmat(spin);
  const int n = opdm1.n();
  RefSymmSCMatrix lambda = tpdm.copy();

  int b12 = 0;
  for(int b1=0; b1<n; ++b1) {
    const int b2fence = (spin == AlphaBeta) ? n : b1;
    for(int b2=0; b2< b2fence; ++b2, ++b12) {

      int k12 = 0;
      for(int k1=0; k1<n; ++k1) {
        const double gamma_b1_k1 = opdm1.get_element(b1,k1);
        const double gamma_b2_k1 = (spin != AlphaBeta) ? opdm1.get_element(b2,k1) : 0.0;

        const int k2fence = (spin == AlphaBeta) ? n : k1;
        for(int k2=0; k2< k2fence; ++k2, ++k12) {
          double value = - (gamma_b1_k1 * opdm2.get_element(b2,k2));
          if (spin != AlphaBeta) value += opdm1.get_element(b1,k2) * gamma_b2_k1;
          lambda.accumulate_element(b12,k12,value);
        }
      }
    }
  }

  scmat_[spin] = lambda;
  return scmat_[spin];
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
