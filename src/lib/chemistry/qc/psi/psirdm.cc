//
// psirdm.cc
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
#include <chemistry/qc/psi/psirdm.h>
#include <chemistry/qc/psi/psici.h>

using namespace sc;

namespace {
  Ref<OrbitalSpace> orbs_from_psiwfn(const Ref<PsiWavefunction>& wfn, SpinCase1 spin) {
    Ref<PsiSCF> psiscfwfn_;
    psiscfwfn_ << wfn;
    if (psiscfwfn_.nonnull()) {
      return psiscfwfn_->orbs_sb(spin);
    }

    Ref<PsiRASCI> psiciwfn_;
    psiciwfn_ << wfn;
    if (psiciwfn_.nonnull()) {
      return psiciwfn_->occ(spin);
    }

    Ref<PsiCorrWavefunction> psicorrwfn_;
    psicorrwfn_ << wfn;
    if (psicorrwfn_.nonnull()) {
      return psicorrwfn_->orbs_sb(spin);
    }

    throw std::logic_error("don't know how to compute OrbitalSpace from this PsiWavefunction");
  }
}

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
  wfn_ = require_dynamic_cast<PsiWavefunction*>(
        kv->describedclassvalue("wfn").pointer(),
        "PsiRDMTwo::PsiRDMTwo\n"
        );
}

PsiRDMTwo::PsiRDMTwo(StateIn& si) : RDM<Two>(si) {
  wfn_ << SavableState::restore_state(si);
  if (wfn_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

PsiRDMTwo::~PsiRDMTwo() {
}

void
PsiRDMTwo::save_data_state(StateOut& so) {
  SavableState::save_state(wfn_.pointer(), so);
}

Ref<PsiRDMTwo::cumulant_type> sc::PsiRDMTwo::cumulant() const
{
  return new RDMCumulant<Two>(const_cast<PsiRDMTwo*>(this));
}

Ref< RDM<One> > sc::PsiRDMTwo::rdm_m_1() const {
  return new PsiRDMOne(wfn_);
}

Ref<OrbitalSpace>
sc::PsiRDMTwo::orbs(SpinCase1 spin) const {
  try {
    return orbs_from_psiwfn(wfn_, spin);
  }
  catch (...) {
    throw ProgrammingError("PsiRDMTwo::orbs() not defined for this PsiWavefunction",
                           __FILE__,
                           __LINE__);
  }
}

RefSymmSCMatrix sc::PsiRDMTwo::scmat(SpinCase2 spin) const {

  if (scmat_[spin].nonnull()) return scmat_[spin];
  if (!wfn_->spin_polarized() && spin == BetaBeta)
    return scmat(AlphaAlpha);

  Ref<PsiRASCI> psiciwfn_;
  psiciwfn_ << wfn_;
  if (psiciwfn_.nonnull()) {
    return psiciwfn_->twopdm_occ(spin);
  }

  Ref<PsiCorrWavefunction> psicorrwfn_;
  psicorrwfn_ << wfn_;
  if (psicorrwfn_.nonnull()) {
    scmat_[spin] = psicorrwfn_->twopdm_dirac(spin);
    return scmat_[spin];
  }

  Ref<PsiSCF> psiscfwfn_;
  psiscfwfn_ << wfn_;
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

#if 0
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

RefSymmSCMatrix sc::PsiRDMCumulantTwo::scmat(SpinCase2 spin) const {

  if (scmat_[spin].nonnull()) return scmat_[spin];
  if (!density_->wfn()->spin_polarized() && spin == BetaBeta)
    return scmat(AlphaAlpha);

  const SpinCase1 spin1 = case1(spin);
  const SpinCase1 spin2 = case2(spin);
  const RefSymmSCMatrix opdm1 = density_->wfn()->mo_density(spin1);
  const RefSymmSCMatrix opdm2 = density_->wfn()->mo_density(spin2);
  const RefSymmSCMatrix tpdm = density_->scmat(spin);
  const int n = opdm1.n();
  RefSymmSCMatrix lambda = tpdm.copy();

  int b12 = 0;
  for(int b1=0; b1<n; ++b1) {
    const int b2fence = (spin == AlphaBeta) ? n : b1;
    for(int b2=0; b2< b2fence; ++b2, ++b12) {

      int k12 = 0;
      for(int k1=0; k1<n && k12<=b12; ++k1) {
        const double gamma_b1_k1 = opdm1.get_element(b1,k1);
        const double gamma_b2_k1 = (spin != AlphaBeta) ? opdm1.get_element(b2,k1) : 0.0;

        const int k2fence = (spin == AlphaBeta) ? n : k1;
        for(int k2=0; k2<k2fence && k12<=b12; ++k2, ++k12) {
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
#endif

////////////////////

ClassDesc
PsiRDMOne::class_desc_(typeid(PsiRDMOne),
                     "PsiRDMOne",
                     1,               // version
                     "public RDM<One>", // must match parent
                     0,               // change to create<PsiRDMTwo> if this class is DefaultConstructible
                     create<PsiRDMOne>, // change to 0 if this class is not KeyValConstructible
                     create<PsiRDMOne>  // change to 0 if this class is not StateInConstructible
                     );

PsiRDMOne::PsiRDMOne(const Ref<KeyVal>& kv) : RDM<One>(kv) {
  wfn_ = require_dynamic_cast<PsiWavefunction*>(
        kv->describedclassvalue("wfn").pointer(),
        "PsiRDMOne::PsiRDMOne\n"
        );
}

PsiRDMOne::PsiRDMOne(StateIn& si) : RDM<One>(si) {
  wfn_ << SavableState::restore_state(si);
  if (wfn_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

PsiRDMOne::PsiRDMOne(const Ref<PsiWavefunction>& wfn) : RDM<One>(wfn.pointer()), wfn_(wfn) {
}

PsiRDMOne::~PsiRDMOne() {
}

void
PsiRDMOne::save_data_state(StateOut& so) {
  SavableState::save_state(wfn_.pointer(), so);
}

Ref<OrbitalSpace>
PsiRDMOne::orbs(SpinCase1 spin) const {
  try {
    return orbs_from_psiwfn(wfn_, spin);
  }
  catch (...) {
    throw ProgrammingError("PsiRDMTwo::orbs() not defined for this PsiWavefunction",
                           __FILE__,
                           __LINE__);
  }
}

RefSymmSCMatrix
PsiRDMOne::scmat(SpinCase1 spin) const {
  Ref<PsiRASCI> psiciwfn;
  psiciwfn << wfn_;
  if (psiciwfn.nonnull()) {
    return psiciwfn->onepdm_occ(spin);
  }

  return wfn_->mo_density(spin);
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc
PsiSpinFreeRDMTwo::class_desc_(typeid(PsiSpinFreeRDMTwo),
                     "PsiSpinFreeRDMTwo",
                     1,               // version
                     "public SpinFreeRDM<Two>", // must match parent
                     0,               // change to create<PsiRDMTwo> if this class is DefaultConstructible
                     create<PsiSpinFreeRDMTwo>, // change to 0 if this class is not KeyValConstructible
                     create<PsiSpinFreeRDMTwo>  // change to 0 if this class is not StateInConstructible
                     );

PsiSpinFreeRDMTwo::PsiSpinFreeRDMTwo(const Ref<KeyVal>& kv) : SpinFreeRDM<Two>(kv) {
  wfn_ = require_dynamic_cast<PsiWavefunction*>(
        kv->describedclassvalue("wfn").pointer(),
        "PsiSpinFreeRDMTwo::PsiSpinFreeRDMTwo\n"
        );
}

PsiSpinFreeRDMTwo::PsiSpinFreeRDMTwo(StateIn& si) : SpinFreeRDM<Two>(si) {
  wfn_ << SavableState::restore_state(si);
  if (wfn_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

PsiSpinFreeRDMTwo::~PsiSpinFreeRDMTwo() {
}

void
PsiSpinFreeRDMTwo::save_data_state(StateOut& so) {
  SavableState::save_state(wfn_.pointer(), so);
}

Ref< SpinFreeRDM<One> > sc::PsiSpinFreeRDMTwo::rdm_m_1() const {
  return new PsiSpinFreeRDMOne(wfn_);
}

Ref<OrbitalSpace>
sc::PsiSpinFreeRDMTwo::orbs() const {
  try {
    return orbs_from_psiwfn(wfn_, Alpha);
  }
  catch (...) {
    throw ProgrammingError("PsiSpinFreeRDMTwo::orbs() not defined for this PsiWavefunction",
                           __FILE__,
                           __LINE__);
  }
}

RefSymmSCMatrix sc::PsiSpinFreeRDMTwo::scmat() const {

  if (scmat_.nonnull()) return scmat_;

  Ref<PsiRASCI> psiciwfn_;
  psiciwfn_ << wfn_;
  if (psiciwfn_.nonnull()) {
    return psiciwfn_->twopdm_occ();
  }

  Ref<PsiCorrWavefunction> psicorrwfn_;
  psicorrwfn_ << wfn_;
  if (psicorrwfn_.nonnull()) {
    scmat_ = psicorrwfn_->twopdm_dirac();
    return scmat_;
  }

  Ref<PsiSCF> psiscfwfn_;
  psiscfwfn_ << wfn_;
  if (psiscfwfn_.nonnull()) {
    Ref<SCMatrixKit> kit = SCMatrixKit::default_matrixkit();
    const RefSymmSCMatrix opdm_a = psiscfwfn_->mo_density(Alpha);
    const RefSymmSCMatrix opdm_b = psiscfwfn_->mo_density(Beta);
    MPQC_ASSERT(opdm_a.n() == opdm_b.n());
    const int n = opdm_a.n();

    RefSCDimension dim = new SCDimension(n*n);
    RefSymmSCMatrix tpdm = kit->symmmatrix(dim); tpdm.assign(0.0);

    int b12 = 0;
    for(int b1=0; b1<n; ++b1) {
      for(int b2=0; b2<n; ++b2, ++b12) {

        int k12 = 0;
        for(int k1=0; k1<n; ++k1) {
          const double gamma_a_b1_k1 = opdm_a.get_element(b1,k1);
          const double gamma_a_b2_k1 = opdm_a.get_element(b2,k1);
          const double gamma_b_b1_k1 = opdm_b.get_element(b1,k1);
          const double gamma_b_b2_k1 = opdm_b.get_element(b2,k1);

          for(int k2=0; k2<n; ++k2, ++k12) {

            const double gamma_a_b2_k2 = opdm_a.get_element(b2,k2);
            const double gamma_a_b1_k2 = opdm_a.get_element(b1,k2);
            const double gamma_b_b2_k2 = opdm_b.get_element(b2,k2);
            const double gamma_b_b1_k2 = opdm_b.get_element(b1,k2);

            double value = (gamma_a_b1_k1 * gamma_a_b2_k2 - gamma_a_b2_k1 * gamma_a_b1_k2)
                + (gamma_b_b1_k1 * gamma_b_b2_k2 - gamma_b_b2_k1 * gamma_b_b1_k2)
                + gamma_a_b1_k1 * gamma_b_b2_k2
                + gamma_b_b1_k1 * gamma_a_b2_k2;

            tpdm.accumulate_element(b12,k12,value);
          }
        }
      }
    }
    scmat_ = tpdm;
    return scmat_;
  } // PsiSCF wavefunction

  throw FeatureNotImplemented("PsiSpinFreeRDMTwo::scmat() is not implemented for this PsiWavefunction",__FILE__,__LINE__);
}

const Ref<DistArray4>&
PsiSpinFreeRDMTwo::da4() const {
  if (da4_.null()) {

    Ref<PsiRASCI> psiciwfn_;
    psiciwfn_ << wfn_;
    if (psiciwfn_.nonnull()) {
      RefSymmSCMatrix rdm2_scmat = psiciwfn_->twopdm_occ();
      const int n = psiciwfn_->occ(Alpha)->rank();
      MPQC_ASSERT(n == int(sqrt(rdm2_scmat.n())));
      da4_ = make_distarray4(1, n, n, n, n);
      da4_->activate();

      std::vector<double> k1k2_buf(n*n);
      int b12 = 0;
      for(int b1=0; b1<n; ++b1) {
        for(int b2=0; b2<n; ++b2, ++b12) {

          RefSCVector b12_row = rdm2_scmat.get_row(b12);
          b12_row.convert(&k1k2_buf[0]);

          da4_->store_pair_block(b1, b2, 0, &(k1k2_buf[0]));
        }
      }

      if (da4_->data_persistent()) da4_->deactivate();

      return da4_;

    }

    Ref<PsiCorrWavefunction> psicorrwfn_;
    psicorrwfn_ << wfn_;
    if (psicorrwfn_.nonnull()) {
      throw ProgrammingError("PsiSpinFreeRDMTwo::da4() -- DistArray4-based storage is not implemented for PsiCorrWavefunction",
                             __FILE__, __LINE__);
    }

    Ref<PsiSCF> psiscfwfn_;
    psiscfwfn_ << wfn_;
    if (psiscfwfn_.nonnull()) {
      const RefSymmSCMatrix opdm_a = psiscfwfn_->mo_density(Alpha);
      const RefSymmSCMatrix opdm_b = psiscfwfn_->mo_density(Beta);
      MPQC_ASSERT(opdm_a.n() == opdm_b.n());
      const int n = opdm_a.n();

      da4_ = make_distarray4(1, n, n, n, n);
      da4_->activate();

      std::vector<double> k1k2_buf(n*n);
      int b12 = 0;
      for(int b1=0; b1<n; ++b1) {
        for(int b2=0; b2<n; ++b2, ++b12) {

          int k12 = 0;
          for(int k1=0; k1<n; ++k1) {
            const double gamma_a_b1_k1 = opdm_a.get_element(b1,k1);
            const double gamma_a_b2_k1 = opdm_a.get_element(b2,k1);
            const double gamma_b_b1_k1 = opdm_b.get_element(b1,k1);
            const double gamma_b_b2_k1 = opdm_b.get_element(b2,k1);

            for(int k2=0; k2<n; ++k2, ++k12) {

              const double gamma_a_b2_k2 = opdm_a.get_element(b2,k2);
              const double gamma_a_b1_k2 = opdm_a.get_element(b1,k2);
              const double gamma_b_b2_k2 = opdm_b.get_element(b2,k2);
              const double gamma_b_b1_k2 = opdm_b.get_element(b1,k2);

              double value = (gamma_a_b1_k1 * gamma_a_b2_k2 - gamma_a_b2_k1 * gamma_a_b1_k2)
                  + (gamma_b_b1_k1 * gamma_b_b2_k2 - gamma_b_b2_k1 * gamma_b_b1_k2)
                  + gamma_a_b1_k1 * gamma_b_b2_k2
                  + gamma_b_b1_k1 * gamma_a_b2_k2;

              k1k2_buf[k12] = (b12 == k12 ? 1.0 : 2.0) * value;
            }
          }

          da4_->store_pair_block(b1, b2, 0, &(k1k2_buf[0]));
        }
      }

      if (da4_->data_persistent()) da4_->deactivate();

      return da4_;
    }

    throw FeatureNotImplemented("PsiSpinFreeRDMTwo::da4() is not implemented for this PsiWavefunction",__FILE__,__LINE__);
  }
  // da4_ is nonzero
  return da4_;
}

/////////////////////////////////////////////////////////////////////////////

ClassDesc
PsiSpinFreeRDMOne::class_desc_(typeid(PsiSpinFreeRDMOne),
                     "PsiSpinFreeRDMOne",
                     1,               // version
                     "public SpinFreeRDM<One>", // must match parent
                     0,               // change to create<PsiRDMTwo> if this class is DefaultConstructible
                     create<PsiSpinFreeRDMOne>, // change to 0 if this class is not KeyValConstructible
                     create<PsiSpinFreeRDMOne>  // change to 0 if this class is not StateInConstructible
                     );

PsiSpinFreeRDMOne::PsiSpinFreeRDMOne(const Ref<KeyVal>& kv) : SpinFreeRDM<One>(kv) {
  wfn_ = require_dynamic_cast<PsiWavefunction*>(
        kv->describedclassvalue("wfn").pointer(),
        "PsiSpinFreeRDMOne::PsiRDMOne\n"
        );
}

PsiSpinFreeRDMOne::PsiSpinFreeRDMOne(StateIn& si) : SpinFreeRDM<One>(si) {
  wfn_ << SavableState::restore_state(si);
  if (wfn_.null())
    throw ProgrammingError("failed constructor",__FILE__,__LINE__);
}

PsiSpinFreeRDMOne::PsiSpinFreeRDMOne(const Ref<PsiWavefunction>& wfn) : SpinFreeRDM<One>(wfn.pointer()), wfn_(wfn) {
}

PsiSpinFreeRDMOne::~PsiSpinFreeRDMOne() {
}

void
PsiSpinFreeRDMOne::save_data_state(StateOut& so) {
  SavableState::save_state(wfn_.pointer(), so);
}

Ref<OrbitalSpace>
PsiSpinFreeRDMOne::orbs() const {
  try {
    return orbs_from_psiwfn(wfn_, Alpha);
  }
  catch (...) {
    throw ProgrammingError("PsiSpinFreeRDMOne::orbs() not defined for this PsiWavefunction",
                           __FILE__,
                           __LINE__);
  }
}

RefSymmSCMatrix
PsiSpinFreeRDMOne::scmat() const {
  Ref<PsiRASCI> psiciwfn;
  psiciwfn << wfn_;
  if (psiciwfn.nonnull()) {
    return psiciwfn->onepdm_occ(Alpha) + psiciwfn->onepdm_occ(Beta);
  }
  else
    return wfn_->mo_density(Alpha) + wfn_->mo_density(Beta);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
