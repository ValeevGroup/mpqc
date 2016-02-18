//
// psicc_pt2r12.cc
//
// Copyright (C) 2002 Edward Valeev
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

#include <cmath>
#include <cassert>
#include <ccfiles.h>
#include <math/scmat/local.h>
#include <chemistry/qc/lcao/utils.h>
#include <util/misc/print.h>
#include <chemistry/qc/lcao/utils.impl.h>
#include <math/mmisc/pairiter.impl.h>
#include <chemistry/qc/psi/psicc_pt2r12.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>

#define TEST_ViaT1 0

using namespace sc;
using namespace sc::fastpairiter;

namespace {
  void _print(SpinCase2 spin,
              const Ref<DistArray4>& mat,
              const char* label);
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCC_PT2R12_cd(typeid(PsiCC_PT2R12), "PsiCC_PT2R12", 2,
                                   "public PsiCC", 0, 0, 0);

PsiCC_PT2R12::PsiCC_PT2R12(const Ref<KeyVal>&keyval) :
  PsiCC(keyval), cabs_singles_energy_(0.0) {
  if (!replace_Lambda_with_T_)
    throw FeatureNotImplemented("PsiCC_PT2R12::PsiCC_PT2R12() -- cannot properly use Lambdas yet",__FILE__,__LINE__);

  // if world not given, make this the center of a new World
  Ref<WavefunctionWorld> world; world << keyval->describedclassvalue("world", KeyValValueRefDescribedClass(0));
  if (world.null())
    world = new WavefunctionWorld(keyval);
  if (world.null())
    throw InputError("PsiCC_PT2R12 requires a WavefunctionWorld; input did not specify it, neither could it be constructed",
                     __FILE__, __LINE__, "world");
  if (world->wfn() == 0) world->set_wfn(this);

  Ref<RefWavefunction> refinfo = RefWavefunctionFactory::make(world, this->reference(), false,
                                                              this->nfzc(), this->nfzv());
  r12world_ = new R12WavefunctionWorld(keyval, refinfo);
  Ref<R12Technology> r12tech = r12world_->r12tech();
  cabs_singles_ = keyval->booleanvalue("cabs_singles",KeyValValueboolean(true));
  const bool openshell = this->reference()->spin_polarized();
  spinadapted_ = keyval->booleanvalue("spinadapted", KeyValValueboolean(openshell ? 0 : 1));

  pccsd_alpha_ = keyval->doublevalue("pccsd_alpha", KeyValValuedouble(1.0));
  pccsd_beta_ = keyval->doublevalue("pccsd_beta", KeyValValuedouble(1.0));
  pccsd_gamma_ = keyval->doublevalue("pccsd_gamma", KeyValValuedouble(1.0));

  r12eval_ = 0;
  mp2r12_energy_ = 0;

  // cannot do gbc = false yet
  if (!r12tech->gbc())
    throw FeatureNotImplemented("PsiCC_PT2R12::PsiCC_PT2R12() -- gbc = false is not yet implemented",__FILE__,__LINE__);
  // cannot do coupling=true either
  if (r12tech->coupling())
    throw FeatureNotImplemented("PsiCC_PT2R12::PsiCC_PT2R12() -- coupling = true is not yet implemented",__FILE__,__LINE__);
}

PsiCC_PT2R12::~PsiCC_PT2R12() {
}

PsiCC_PT2R12::PsiCC_PT2R12(StateIn&s) :
  PsiCC(s) {
  if (this->class_version() < 2)
    throw ProgrammingError("cannot use archives older than version 2", __FILE__, __LINE__, class_desc());

  r12eval_ << SavableState::restore_state(s);
  r12world_ << SavableState::restore_state(s);
  mp2r12_energy_ << SavableState::restore_state(s);

  if((r12world()->r12tech()->ansatz()->orbital_product_GG()==R12Technology::OrbProdGG_pq) ||
     (r12world()->r12tech()->ansatz()->orbital_product_gg()==R12Technology::OrbProdgg_pq)) {
    throw InputError("PsiCC_PT2R12::PsiCC_PT2R12 -- pq Ansatz not allowed",__FILE__,__LINE__);
  }

  int spinadapted; s.get(spinadapted); spinadapted_ = (bool)spinadapted;
  s.get(cabs_singles_);
  s.get(cabs_singles_energy_);
}

void PsiCC_PT2R12::save_data_state(StateOut&s) {
  PsiCC::save_data_state(s);
  SavableState::save_state(r12eval_.pointer(),s);
  SavableState::save_state(r12world_.pointer(),s);
  SavableState::save_state(mp2r12_energy_.pointer(),s);

  s.put((int)spinadapted_);
  s.put(cabs_singles_);
  s.put(cabs_singles_energy_);
}

void PsiCC_PT2R12::write_basic_input(int convergence) {
  Ref<PsiInput> input = get_psi_input();
  input->write_keyword("ccenergy:convergence", convergence);
  input->write_keyword("ccenergy:maxiter", maxiter_);
  input->write_keyword("ccenergy:pccsd_alpha", pccsd_alpha_);
  input->write_keyword("ccenergy:pccsd_beta", pccsd_beta_);
  input->write_keyword("ccenergy:pccsd_gamma", pccsd_gamma_);
}

void PsiCC_PT2R12::compute_ept2r12() {

  r12world()->initialize();

  // compute Psi3 CC wave function
  PsiCorrWavefunction::compute();

  // to compute intermediates make sure r12eval_ is ready
  if (r12eval_.null()) {
    r12eval_ = new R12IntEval(r12world());
    r12eval_->debug(debug_);
  }
  Ref<R12Technology> r12tech = r12world_->r12tech();
  Ref<R12EnergyIntermediates> r12intermediates = new R12EnergyIntermediates(r12eval_,r12tech->stdapprox());

  //
  // Obtain CC amplitudes from Psi and copy into local matrices
  //
  RefSCMatrix T1[NSpinCases1];
  Ref<DistArray4> T2[NSpinCases2];
  Ref<DistArray4> L2[NSpinCases2];
  // T1
  const int nspincases1 = r12eval()->nspincases1();
  const int nspincases2 = (r12eval()->spin_polarized() ? 3 : 2);
  for(int s=0; s<nspincases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);
    RefSCMatrix T1_psi = this->T1(spin);
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    T1[s] = localkit->matrix(T1_psi.rowdim(), T1_psi.coldim());
    T1[s]->convert(T1_psi);
    if (debug() >= DefaultPrintThresholds::mostN2) {
      T1[spin].print(prepend_spincase(spin,"CCSD T1 amplitudes:").c_str());
    }
  }
  if (nspincases1 == 1) {
    T1[Beta] = T1[Alpha];
  }
  // T2
  for(int s=0; s<nspincases2; ++s) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;

    T2[s] = this->T2_da4(spincase2, "t");

    if (debug() >= DefaultPrintThresholds::mostO2N2) {
      _print(spincase2, T2[spincase2], prepend_spincase(spincase2,"CCSD T2 amplitudes:").c_str());
    }
  }
  // L2
  if (need_lambda_) {
  for(int s=0; s<nspincases2; ++s) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;

    L2[s] = this->Lambda2_da4(spincase2);

    if (debug() >= DefaultPrintThresholds::mostO2N2) {
      _print(spincase2, L2[spincase2], prepend_spincase(spincase2,"CCSD L2 amplitudes:").c_str());
    }
  }
  }

  // compute (2)_R12 energy as MP2-R12 energy with dressed V intermediate

  // include cabs singles energy?
  cabs_singles_energy_ = 0.0;
  if (cabs_singles_) {
    cabs_singles_energy_ = r12eval()->emp2_cabs_singles(T1[Alpha],T1[Beta]);
    //cabs_singles_energy_ = r12eval()->emp2_cabs_singles(); // test: use MP2 CABS Singles
  }

  const bool diag = r12eval()->r12world()->r12tech()->ansatz()->diag();
  if (diag == true) { // -> pass T1 and T2 amplitudes to the diagonal MP2-R12 energy evaluator

    // Pass T1 to r12intermediates
    for(int s=0; s<NSpinCases1; s++) {
      const SpinCase1 spin = static_cast<SpinCase1>(s);
      r12intermediates->assign_T1_cc(spin,T1[spin]);
    }
    //_print(AlphaBeta, T2[0], prepend_spincase(AlphaBeta,"CCSD T2[0] amplitudes:").c_str());
    // Pass T2 to r12intermediates
    for(int s=0; s<nspincases2; ++s) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      if (r12eval()->dim_oo(spincase2).n() == 0)
        continue;
      r12intermediates->assign_T2_cc(spincase2,T2[s]);
    }
    // Pass L2 to r12intermediates, if needed
    if (need_lambda_) {
    for(int s=0; s<nspincases2; ++s) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      if (r12eval()->dim_oo(spincase2).n() == 0)
        continue;
      r12intermediates->assign_L2_cc(spincase2,L2[s]);
    }
    }

  // Import the Psi CCSD one-particle density
  if (compute_1rdm_) {
    // Obtain CC one-particle density from Psi and copy into local matrices
    //
    RefSCMatrix D[NSpinCases1];
    for(int s = 0; s < nspincases1; ++s) {
      const SpinCase1 spin = static_cast<SpinCase1>(s);
      RefSCMatrix D_psicc = this->Onerdm(spin);
      Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
      D[spin] = localkit->matrix(D_psicc.rowdim(), D_psicc.coldim());
      D[spin]->convert(D_psicc);

      if (debug() >= DefaultPrintThresholds::mostN2) {
      D[spin].print(prepend_spincase(spin,"CCSD one-particle density:").c_str());
      }
    }

    if (nspincases1 == 1) {
      D[Beta] = D[Alpha];
    }

      // -> pass D to the diagonal MP2-R12 energy evaluator
      for(int s = 0; s < NSpinCases1; s++) {
        const SpinCase1 spin = static_cast<SpinCase1>(s);
        r12intermediates->assign_1rdm_cc(spin,D[spin]);
      }

      if (!strncmp(dertype_, "first", 10)) {
        RefSCMatrix D_orbs[NSpinCases1];
        // T1 & cabs I1 need to be ready for CABS_Singles orbital relaxation Z-vector
        compute_onerdm_relax(r12intermediates, D_orbs[Alpha], D_orbs[Beta]);
      if (debug() >= DefaultPrintThresholds::mostN2) {
          D_orbs[Alpha].print(prepend_spincase(Alpha,"CCSD_F12 one-particle density from relaxation:").c_str());
          if (nspincases1 != 1) {
            D_orbs[Beta].print(prepend_spincase(Beta,"CCSD_F12 one-particle density from relaxation:").c_str());
          }
      }
        // pass orbital relaxation 1rdm to the diagonal MP2-R12 energy evaluator
        for(int s = 0; s < NSpinCases1; s++) {
          const SpinCase1 spin = static_cast<SpinCase1>(s);
          r12intermediates->assign_1rdm_relax(spin,D_orbs[spin]);
        }

        // tests:
//        for(int s = 0; s < nspincases1; ++s) {
//          const SpinCase1 spin = static_cast<SpinCase1>(s);
//
//        // test: print the Z-vector from PSI3
//        RefSCMatrix X_psi = this->Onerdm_relax_X(spin);
//        X_psi.print(prepend_spincase(spin,"PSI3 Z-vector X:").c_str());
//        // test: print the relaxation effect from PSI3
//        RefSCMatrix Dorbs_psi = this->Onerdm_relax_D(Alpha);
//        Dorbs_psi.print(prepend_spincase(Alpha,"PSI3 Dorbs:").c_str());
//
//        // test: print the Z-vector from F12 contribution
//        RefSCMatrix Xf12 = Onerdm_X_F12(spin, r12eval_, debug());
//        Xf12.print(prepend_spincase(spin,"F12 Z-vector X:").c_str());
//
//          if (cabs_singles_) {
//            // test: print the Z-vector from CABS Singles contribution
//            RefSCMatrix X_cabs = Onerdm_X_CABS_Singles(spin, r12eval_, r12intermediates, debug());
//            X_cabs.print(prepend_spincase(spin,"CABS_Singles Z-vector X:").c_str());
//          }
//        }

      } // end of dertype_

  } // end of compute_1rdm_


  } // end of diag = true clause
  else { // diag = false:
         // compute dressed V intermediate of CC(2)_R12 method here explicitly if non-diagonal ansatz is used
         // (diagonal MP2-R12 energy evaluator computes all intermediates)
#if 0
    ExEnv::out0() << indent << "Trying out MPQC3-based R12 code" << std::endl;
    r12eval()->compute_ccr12_1rdm(T1[Alpha], T2);
#endif

    Ref<DistArray4> Vab[NSpinCases2];
    Ref<DistArray4> Via[NSpinCases2];
    Ref<DistArray4> Vai[NSpinCases2];
    RefSCMatrix Vij[NSpinCases2];
    RefSymmSCMatrix B[NSpinCases2];
    RefSymmSCMatrix X[NSpinCases2];
    Ref<DistArray4> A[NSpinCases2];

    for (int s=0; s<nspincases2; s++) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const SpinCase1 spin1 = case1(spincase2);
      const SpinCase1 spin2 = case2(spincase2);
      if (r12eval()->dim_oo(spincase2).n() == 0)
        continue;

      const Ref<OrbitalSpace>& p1 = r12eval()->orbs(spin1);
      const Ref<OrbitalSpace>& p2 = r12eval()->orbs(spin2);
      const Ref<OrbitalSpace>& x1 = r12eval()->GGspace(spin1);
      const Ref<OrbitalSpace>& x2 = r12eval()->GGspace(spin2);
      const Ref<OrbitalSpace>& v1 = r12eval()->vir_act(spin1);
      const Ref<OrbitalSpace>& v2 = r12eval()->vir_act(spin2);
      const Ref<OrbitalSpace>& o1 = r12eval()->occ_act(spin1);
      const Ref<OrbitalSpace>& o2 = r12eval()->occ_act(spin2);

      std::vector< Ref<DistArray4> > Vpq_vec = r12eval()->V_distarray4(spincase2, p1, p2);
      MPQC_ASSERT(Vpq_vec.size() == 1);
      Ref<DistArray4> Vpq = Vpq_vec[0];
#define SKIP_R12INTEVAL_COMPUTE 0
#if !SKIP_R12INTEVAL_COMPUTE
      Vij[s] = r12eval()->V(spincase2);
      X[s] = r12eval()->X(spincase2);
      B[s] = r12eval()->B(spincase2);
#endif

      // extract Vab and Via from Vpq
      map(Vpq, x1, x2, p1, p2, Vab[s], x1, x2, v1, v2);
      map(Vpq, x1, x2, p1, p2, Via[s], x1, x2, o1, v2);
      if (r12tech->ebc() == false) {
        std::vector< Ref<DistArray4> > Avec = A_distarray4(spincase2, r12eval_);
        A[s] = Avec[0];
      }

      // Vai[AlphaBeta] is also needed if spin-polarized reference is used
      if (spincase2 == AlphaBeta && r12world()->refwfn()->spin_polarized())
        map(Vpq, x1, x2, p1, p2, Vai[s], x1, x2, v1, o2);

      if (debug() >= DefaultPrintThresholds::mostO2N2) {
        _print(spincase2, Vpq, prepend_spincase(spincase2,"Vpq matrix").c_str());
        _print(spincase2, Vab[s], prepend_spincase(spincase2,"Vab matrix").c_str());
        _print(spincase2, Via[s], prepend_spincase(spincase2,"Via matrix").c_str());
        if (Vai[s])
          _print(spincase2, Vai[s], prepend_spincase(spincase2,"Vai matrix").c_str());
        if (A[s])
          _print(spincase2, A[s], prepend_spincase(spincase2,"A matrix").c_str());
      }
      if (debug() >= DefaultPrintThresholds::mostO4) {
        Vij[s].print(prepend_spincase(spincase2,"Vij matrix").c_str());
      }
    } // end of spincase2 loop
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

    // Compute Hamiltonian matrix elements
    RefSCMatrix H1_R0[NSpinCases2];
    for (int s=0; s<nspincases2; s++) {
      const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
      const SpinCase1 spin1 = case1(spincase2);
      const SpinCase1 spin2 = case2(spincase2);
      if (r12eval()->dim_oo(spincase2).n() == 0)
        continue;

      const Ref<OrbitalSpace>& x1 = r12eval()->GGspace(spin1);
      const Ref<OrbitalSpace>& x2 = r12eval()->GGspace(spin2);
      const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
      const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
      const Ref<OrbitalSpace>& vir1_act = r12eval()->vir_act(spin1);
      const Ref<OrbitalSpace>& vir2_act = r12eval()->vir_act(spin2);
      const bool p1_equiv_p2 = (occ1_act == occ2_act) && (vir1_act == vir2_act);


      // Vij is the leading term in H1_R0
      H1_R0[s] = Vij[s].clone();
      H1_R0[s].assign(Vij[s]);

      // the rest of terms (Vbar) are bigger than o^4, hence will be computed in DistArray4
      {
        DistArray4Dimensions dims(1, x1->rank(), x2->rank(), occ1_act->rank(), occ2_act->rank());
        Ref<DistArray4> HT = Vab[s]->clone(dims);

        // the leading term in <R|(HT)|0> is T2.Vab
        // if not assuming EBC then also include coupling matrix term
        if (r12tech->ebc() == false) { // Vab += A
          axpy(A[s], 1.0, Vab[s]);
        }
        contract34( HT, 1.0,
                    Vab[s], 0,
                    T2[s], 0 );

        if (debug() >= DefaultPrintThresholds::allO4)
          _print(spincase2, HT, prepend_spincase(spincase2,"<R|(H*T2)|0>").c_str());

        // the next term is T1.Vai. It's third-order if BC hold, second-order otherwise
        if ( completeness_order_for_intermediates_ >= 3 ||
            (completeness_order_for_intermediates_ >= 2 && !r12eval()->bc()) ) {

          // Via . T1
          Ref<DistArray4> VT1 = HT->clone();
          contract4(Via[s],T1[spin2].t(),VT1);
          if (p1_equiv_p2) {
            // NOTE: V_{xy}^{ia} T_a^j + V_{xy}^{aj} T_a^i = 2 * symm(V_{xy}^{ia} T_a^j)
            symmetrize(VT1);
            if (debug() >= DefaultPrintThresholds::allO4)
              _print(spincase2, VT1, prepend_spincase(spincase2,"1/2<R|(H*T1)|0>").c_str());
            axpy(VT1, 2.0, HT);
          }
          else {
            axpy(VT1, 1.0, HT);
            // Vai^t . T1
            contract3(Vai[s],T1[spin1].t(),VT1);
            if (debug() >= DefaultPrintThresholds::allO4)
              _print(spincase2, VT1, prepend_spincase(spincase2,"<R|(H*T1)|0>").c_str());
            axpy(VT1, 1.0, HT);
          }

        }

        if (debug() >= DefaultPrintThresholds::mostO4)
          _print(spincase2, HT, prepend_spincase(spincase2,"<R|(H*T)|0>").c_str());

        RefSCMatrix HT_scmat = H1_R0[s].clone();
        HT_scmat << HT;
        H1_R0[s].accumulate(HT_scmat);
      }

      if (debug() >= DefaultPrintThresholds::O4) {
        H1_R0[s].print(prepend_spincase(spincase2,"<R|Hb|0>").c_str());
      }
    }
    // Make H1_R0[AlphaAlpha] if this is a closed-shell case (it isn't computed then)
    if (nspincases2 == 2) {
      H1_R0[AlphaAlpha] = r12eval()->V(AlphaAlpha).clone();
      antisymmetrize(H1_R0[AlphaAlpha], H1_R0[AlphaBeta], r12eval()->GGspace(Alpha),
                     r12eval()->occ_act(Alpha));
    }

    // compute the second-order correction: E2 = - H1_0R . H0_RR^{-1} . H1_R0 = C . H1_R0
    // H0_RR is the usual B of the standard MP-R12 theory
    // if using screening approximations and replacing lambdas with Ts H1_0R = transpose(H1_R0)
    // E2 = - Vbar . B^{-1} . Vbar, where Vbar = V + VT
    for(int i=0; i<NSpinCases2; i++) {
      SpinCase2 spincase2i=static_cast<SpinCase2>(i);
      if (r12eval()->dim_oo(spincase2i).n() == 0)
        continue;
      r12intermediates->assign_V(spincase2i,H1_R0[spincase2i]);
    }
  } // end of diag == false clause

  Ref<MP2R12Energy> r12energy = construct_MP2R12Energy(r12intermediates,
                                                       false,
                                                       debug(),
                                                       diag);
  std::vector<double> E2(NSpinCases2,0.0);
  const int num_unique_spincases2 = (reference_->spin_polarized() ? 3 : 2);
  r12energy->compute();
  for (int s=0; s<num_unique_spincases2; s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;

    E2[s]=r12energy->ef12tot(spincase2);
  }

  double e2 = cabs_singles_energy();
  ExEnv::out0() << indent << "E2(AB)        = "<< scprintf("%20.15lf",E2[AlphaBeta])
                << std::endl;
  if (num_unique_spincases2 > 2) {
    e2 += E2[AlphaBeta] + E2[AlphaAlpha] + E2[BetaBeta];
    ExEnv::out0() << indent << "E2(BB)        = "<< scprintf("%20.15lf",E2[BetaBeta])
                  << std::endl;
    ExEnv::out0() << indent << "E2(AA)        = "<< scprintf("%20.15lf",E2[AlphaAlpha])
                  << std::endl;
  }
  else {
    e2 += E2[AlphaBeta] + 2.0*E2[AlphaAlpha];
    ExEnv::out0() << indent << "E2(AA)        = "<< scprintf("%20.15lf",2.0*E2[AlphaAlpha])
                  << std::endl;
    ExEnv::out0() << indent << "E2(s)         = "<< scprintf("%20.15lf",E2[AlphaBeta] - E2[AlphaAlpha])
                  << std::endl;
    ExEnv::out0() << indent << "E2(t)         = "<< scprintf("%20.15lf",3.0*E2[AlphaAlpha])
                  << std::endl;
  }
  if (cabs_singles_) {
    ExEnv::out0() << indent << "E2(CABS)      = "<< scprintf("%20.15lf",cabs_singles_energy_)
                  << std::endl;
  }
  ExEnv::out0() << indent << "E2            = "<< scprintf("%20.15lf",e2)
                << std::endl;

  set_energy(reference_energy() + e2);
}

double
PsiCC_PT2R12::cabs_singles_energy()
{
  return cabs_singles_energy_;
}

// Compute orbital relaxation contribution for
// CCSD_F12 one-electron density
void PsiCC_PT2R12::compute_onerdm_relax(const Ref<R12EnergyIntermediates>& r12intermediates,
                                        RefSCMatrix& Dorbs_alpha,
                                        RefSCMatrix& Dorbs_beta)
{
  // grab orbital info
  psi::PSIO& psio = exenv()->psio();

  const SpinCase1 spin1 = Alpha;
  // get # of occupied and unoccupied orbitals of spin S per irrep
  const std::vector<unsigned int>& occ1pi = reference()->occpi(spin1);
  const std::vector<unsigned int>& uocc1pi = reference()->uoccpi(spin1);

  // obtain # of orbitals per irrep
  std::vector<unsigned int> occ1pioff(nirrep_);
  std::vector<unsigned int> uocc1pioff(nirrep_);
  occ1pioff[0] = 0;
  uocc1pioff[0] = 0;

  unsigned int nocc1 = occ1pi[0];
  unsigned int nuocc1 = uocc1pi[0];
  for (unsigned int irrep = 1; irrep < nirrep_; ++irrep) {
    occ1pioff[irrep] = occ1pioff[irrep-1] + occ1pi[irrep-1];
    uocc1pioff[irrep] = uocc1pioff[irrep-1] + uocc1pi[irrep-1];
    nocc1 += occ1pi[irrep];
    nuocc1 += uocc1pi[irrep];
  }

  unsigned int na1i1_dpd = 0;
  for (unsigned int h = 0; h < nirrep_; ++h)
    na1i1_dpd += uocc1pi[h] * occ1pi[h];

  // DPD of orbital product spaces
  unsigned int nemai = 0;
  std::vector<size_t> a1i1_pi(nirrep_);
  for (unsigned int h = 0; h < nirrep_; ++h) {
    size_t nai = 0;
    for (unsigned int g = 0; g < nirrep_; ++g) {
      nai += (size_t)uocc1pi[g] * occ1pi[h ^ g];
    }
    a1i1_pi[h] = nai;
    nemai += a1i1_pi[h] ;
  }

  // compute the orbital Z vector contribution from F12
  RefSCMatrix X_F12[NSpinCases1];
  RefSCMatrix Xf12_alpha = Onerdm_X_F12(spin1, r12eval_, debug());
                           //Onerdm_X_CABS_Singles(spin1, r12eval_, r12intermediates, debug());
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
  X_F12[spin1] = localkit->matrix(Xf12_alpha.rowdim(), Xf12_alpha.coldim());
  X_F12[spin1]->convert(Xf12_alpha);
  //X_F12[spin1].assign(0.0);

  psio.open(CC_OEI, PSIO_OPEN_OLD);
  psio.open(CC_MISC, PSIO_OPEN_OLD);
  if (reference()->reftype() == PsiSCF::rhf) {

    double* Xai_ccsd = new double[na1i1_dpd];
    std::fill_n(Xai_ccsd, na1i1_dpd, 0.0);
    //ExEnv::out0() << std::endl << "X nai_dpd: " << na1i1_dpd << std::endl;

    psio.read_entry(CC_OEI, "XAI",
                    reinterpret_cast<char*>(Xai_ccsd), na1i1_dpd*sizeof(double));

    // add X_F12 to X
    double* Xai = new double[na1i1_dpd];
    std::fill_n(Xai, na1i1_dpd, 0.0);

    double* iter_Xai = Xai;
    const double* iter_Xai_ccsd = Xai_ccsd;
    for (unsigned int h = 0; h < nirrep_; ++h) {
      const unsigned int a_offset = uocc1pioff[h];
      const unsigned int i_offset = occ1pioff[h];

      for (int a = 0; a < uocc1pi[h]; ++a)
        for (int i = 0; i< occ1pi[h]; ++i, ++iter_Xai, ++iter_Xai_ccsd)
         *iter_Xai = - *iter_Xai_ccsd + X_F12[spin1].get_element(a+a_offset, i+i_offset);
    }
    delete[] Xai_ccsd;

    // Grab only irrep 0 of the orbital Hessian
//      unsigned int nai_dpd_h0 = uocc1pi[0] * occ1pi[0];
    unsigned int nai_dpd_h0 = a1i1_pi[0];
    //ExEnv::out0() << std::endl << "nai_dpd of irrep 0: " << nai_dpd_h0 << std::endl;
    double** A = block_matrix(nai_dpd_h0, nai_dpd_h0);

    psio_address A_address = PSIO_ZERO;
    for(int ai = 0; ai < nai_dpd_h0; ai++) {
      psio.read(CC_MISC, "A(EM,AI)",
                reinterpret_cast<char*> (A[ai]), nai_dpd_h0*sizeof(double),
                A_address, &A_address);
     }
//    // test: print of A
//    ExEnv::out0() << std::endl << "A(EM,AI)" << std::endl;
//    for (int em = 0; em < nai_dpd_h0; em++) {
//      ExEnv::out0() << em << ":  ";
//      for(int ai = 0; ai < nai_dpd_h0; ai++) {
//          ExEnv::out0() << A[em][ai] << " ";
//      }
//      ExEnv::out0() << std::endl;
//    }

//    // test: print of A in all irreps
//    {
//      double** A = block_matrix(nemai, nemai);
//
//      psio_address A_address = PSIO_ZERO;
//      for(int ai = 0; ai < nemai; ai++) {
//        psio.read(CC_MISC, "A(EM,AI)",
//                  reinterpret_cast<char*> (A[ai]), nemai*sizeof(double),
//                  A_address, &A_address);
//       }
//      ExEnv::out0() << std::endl << "A(EM,AI) with all irreps" << std::endl;
//      for (int em = 0; em < nemai; em++) {
//        ExEnv::out0() << em << ":  ";
//        for(int ai = 0; ai < nemai; ai++) {
//            ExEnv::out0() << A[em][ai] << " ";
//        }
//        ExEnv::out0() << std::endl;
//      }
//    }

     FILE* outfile = tmpfile ();
     pople(A, Xai, nai_dpd_h0, 1, 1e-12, outfile, 0);
     fclose (outfile);

     psio.close(CC_OEI, 1);
     psio.close(CC_MISC, 1);

     RefSCDimension rowdim = new SCDimension(nuocc1);
     RefSCDimension coldim = new SCDimension(nocc1);
     Dorbs_alpha = matrixkit()->matrix(rowdim, coldim);
     Dorbs_alpha.assign(0.0);

     unsigned int ai = 0;
     for (unsigned int h = 0; h < nirrep_; ++h) {
       const unsigned int a_offset = uocc1pioff[h];
       const unsigned int i_offset = occ1pioff[h];

       for (int a = 0; a < uocc1pi[h]; ++a)
         for (int i = 0; i< occ1pi[h]; ++i, ++ai)
           Dorbs_alpha.set_element(a+a_offset, i+i_offset,Xai[ai]);
     }

     delete[] Xai;

     Dorbs_beta = Dorbs_alpha;

     // test: print the relaxation effect from PSI3
//       RefSCMatrix Dorbs_psi = this->Onerdm_relax_D(Alpha);
//       Dorbs_psi.print(prepend_spincase(Alpha,"PSI3 Dorbs:").c_str());
  } else if (reference()->reftype() == PsiSCF::uhf) {

      const SpinCase1 spin2 = Beta;
      const std::vector<unsigned int>& occ2pi = reference()->occpi(spin2);
      const std::vector<unsigned int>& uocc2pi = reference()->uoccpi(spin2);

      RefSCMatrix Xf12_beta = Onerdm_X_F12(spin2, r12eval_, debug());
      X_F12[spin2] = localkit->matrix(Xf12_beta.rowdim(), Xf12_beta.coldim());
      X_F12[spin2]->convert(Xf12_beta);
//      X_F12[spin2].assign(0.0);

      std::vector<unsigned int> occ2pioff(nirrep_);
      std::vector<unsigned int> uocc2pioff(nirrep_);
      occ2pioff[0] = 0;
      uocc2pioff[0] = 0;

      unsigned int nocc2 = occ2pi[0];
      unsigned int nuocc2 = uocc2pi[0];
      for (unsigned int irrep = 1; irrep < nirrep_; ++irrep) {
        occ2pioff[irrep] = occ2pioff[irrep-1] + occ2pi[irrep-1];
        uocc2pioff[irrep] = uocc2pioff[irrep-1] + uocc2pi[irrep-1];
        nocc2 += occ2pi[irrep];
        nuocc2 += uocc2pi[irrep];
      }

      unsigned int na2i2_dpd = 0;
      for (unsigned int h = 0; h < nirrep_; ++h)
        na2i2_dpd += uocc2pi[h] * occ2pi[h];

      double* const Xai_ccsd = new double[na1i1_dpd + na2i2_dpd];
      std::fill_n(Xai_ccsd, na1i1_dpd + na2i2_dpd, 0.0);
      //ExEnv::out0() << std::endl << "X na1i1_dpd: " << na1i1_dpd << std::endl;
      //ExEnv::out0() << std::endl << "X na2i2_dpd: " << na2i2_dpd << std::endl;

      double* iter_Xai_ccsd = Xai_ccsd;
      psio.read_entry(CC_OEI, "XAI",
                      reinterpret_cast<char*>(iter_Xai_ccsd), na1i1_dpd*sizeof(double));

      iter_Xai_ccsd += na1i1_dpd;
      psio.read_entry(CC_OEI, "Xai",
                      reinterpret_cast<char*>(iter_Xai_ccsd), na1i1_dpd*sizeof(double));

      // add Alpha and Beta X_F12 to X
      double* const Xai = new double[na1i1_dpd + na2i2_dpd];
      std::fill_n(Xai, na1i1_dpd + na2i2_dpd, 0.0);

      double* iter_X = Xai;
      const double* iter_X_ccsd = Xai_ccsd;
      for (unsigned int h = 0; h < nirrep_; ++h) {
        const unsigned int a_offset = uocc1pioff[h];
        const unsigned int i_offset = occ1pioff[h];

        for (int a = 0; a < uocc1pi[h]; ++a)
          for (int i = 0; i< occ1pi[h]; ++i, ++iter_X, ++iter_X_ccsd)
           *iter_X = - *iter_X_ccsd + X_F12[spin1].get_element(a+a_offset, i+i_offset);
      }
      for (unsigned int h = 0; h < nirrep_; ++h) {
        const unsigned int a_offset = uocc2pioff[h];
        const unsigned int i_offset = occ2pioff[h];

        for (int a = 0; a < uocc2pi[h]; ++a)
          for (int i = 0; i< occ2pi[h]; ++i, ++iter_X, ++iter_X_ccsd)
           *iter_X = - *iter_X_ccsd + X_F12[spin2].get_element(a+a_offset, i+i_offset);
      }
      delete[] Xai_ccsd;

      std::vector<size_t> a2i2_pi(nirrep_);
      for (unsigned int h = 0; h < nirrep_; ++h) {
        size_t nai = 0;
        for (unsigned int g = 0; g < nirrep_; ++g) {
          nai += (size_t)uocc2pi[g] * occ2pi[h ^ g];
        }
        a2i2_pi[h] = nai;
      }

      // Grab only irrep 0 of the orbital Hessian
      unsigned int na1i1_dpd_h0 = a1i1_pi[0];
      unsigned int na2i2_dpd_h0 = a2i2_pi[0];
      unsigned int nai_dpd_h0  = na1i1_dpd_h0 + na2i2_dpd_h0;
      //ExEnv::out0() << std::endl << "na1i1_dpd of irrep 0: " << na1i1_dpd_h0 << std::endl;
      //ExEnv::out0() << std::endl << "na2i2_dpd of irrep 0: " << na2i2_dpd_h0 << std::endl;

      double** A = block_matrix(nai_dpd_h0, nai_dpd_h0);

      psio_address A1_address = PSIO_ZERO;
      for(int a1i1 = 0; a1i1 < na1i1_dpd_h0; a1i1++) {
        psio.read(CC_MISC, "A(AI,BJ)",
                  reinterpret_cast<char*> (A[a1i1]), na1i1_dpd_h0*sizeof(double),
                  A1_address, &A1_address);
       }

      psio_address A2_address = PSIO_ZERO;
      for(int a2i2 = na1i1_dpd_h0; a2i2 < nai_dpd_h0; a2i2++) {
        double* iter_A = A[a2i2] + na1i1_dpd_h0;
        psio.read(CC_MISC, "A(ai,bj)",
                  reinterpret_cast<char*> (iter_A), na2i2_dpd_h0*sizeof(double),
                  A2_address, &A2_address);
      }

      psio_address A12_address = PSIO_ZERO;
      for(int a1i1 = 0; a1i1 < na1i1_dpd_h0; a1i1++) {
        double* iter_A = A[a1i1] + na1i1_dpd_h0;
        psio.read(CC_MISC, "A(AI,bj)",
                  reinterpret_cast<char*> (iter_A), na2i2_dpd_h0*sizeof(double),
                  A12_address, &A12_address);

        for(int a2i2 = na1i1_dpd_h0; a2i2 < nai_dpd_h0; a2i2++) {
            A[a2i2][a1i1] = A[a1i1][a2i2];
        }
      }

      FILE* outfile = tmpfile ();
      pople(A, Xai, nai_dpd_h0, 1, 1e-12, outfile, 0);
      fclose (outfile);

      psio.close(CC_OEI, 1);
      psio.close(CC_MISC, 1);

      RefSCDimension rowdim1 = new SCDimension(nuocc1);
      RefSCDimension coldim1 = new SCDimension(nocc1);
      Dorbs_alpha = matrixkit()->matrix(rowdim1, coldim1);
      Dorbs_alpha.assign(0.0);

      unsigned int ai = 0;
      for (unsigned int h = 0; h < nirrep_; ++h) {
        const unsigned int a_offset = uocc1pioff[h];
        const unsigned int i_offset = occ1pioff[h];

        for (int a = 0; a < uocc1pi[h]; ++a)
          for (int i = 0; i< occ1pi[h]; ++i, ++ai)
            Dorbs_alpha.set_element(a+a_offset, i+i_offset, Xai[ai]);
      }

      RefSCDimension rowdim2 = new SCDimension(nuocc2);
      RefSCDimension coldim2 = new SCDimension(nocc2);
      Dorbs_beta = matrixkit()->matrix(rowdim2, coldim2);
      Dorbs_beta.assign(0.0);

      for (unsigned int h = 0; h < nirrep_; ++h) {
        const unsigned int a_offset = uocc2pioff[h];
        const unsigned int i_offset = occ2pioff[h];

        for (int a = 0; a < uocc2pi[h]; ++a)
          for (int i = 0; i< occ2pi[h]; ++i, ++ai)
            Dorbs_beta.set_element(a+a_offset, i+i_offset, Xai[ai]);
      }
      delete[] Xai;

      // test: print the relaxation effect from PSI3
//        RefSCMatrix Dorbs_psi_alpha = this->Onerdm_relax_D(Alpha);
//        Dorbs_psi_alpha.print(prepend_spincase(Alpha,"PSI3 Dorbs:").c_str());
//        RefSCMatrix Dorbs_psi_beta = this->Onerdm_relax_D(Beta);
//        Dorbs_psi_beta.print(prepend_spincase(Beta,"PSI3 Dorbs:").c_str());
  }

}
// end of one-electron relaxation

void PsiCC_PT2R12::print(std::ostream&o) const {
  o << indent << "PsiCC_PT2R12:" << std::endl;
  o << incindent;
  o << indent << "Spin-adapted algorithm: " << (spinadapted_ ? "true" : "false") << std::endl;
  o << indent << "Include CABS singles? : " << (cabs_singles_ ? "true" : "false") << std::endl;
  if (pccsd_alpha_ != 1.0 || pccsd_beta_ != 1.0 || pccsd_gamma_ != 1.0) {
    o << indent << "pCCSD(alpha,beta, gamma)     : (" << pccsd_alpha_ << "," << pccsd_beta_ << "," << pccsd_gamma_ << ")" << std::endl;
  }
  if (cabs_singles_) {
    o << indent << "  E(CABS singles) = " << scprintf("%25.15lf", cabs_singles_energy_)
                                          << std::endl;
  }

  r12world()->print(o);
  PsiCC::print(o);
  o << decindent;
}

void
PsiCC_PT2R12::obsolete() {
  r12eval_ = 0;
  cabs_singles_energy_ = 0.0;
  r12world_->world()->obsolete();
  r12world_->obsolete();
  PsiCC::obsolete();
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_PT2R12_cd(typeid(PsiCCSD_PT2R12), "PsiCCSD_PT2R12", 1,
                                    "public PsiCC", 0, create<PsiCCSD_PT2R12>,
                                    create<PsiCCSD_PT2R12>);

PsiCCSD_PT2R12::PsiCCSD_PT2R12(const Ref<KeyVal>&keyval) :
  PsiCC_PT2R12(keyval), eccsd_(NAN) {
  if (!replace_Lambda_with_T_)
    throw FeatureNotImplemented("PsiCCSD_PT2R12::PsiCCSD_PT2R12() -- cannot properly use Lambdas yet",__FILE__,__LINE__);
}

PsiCCSD_PT2R12::~PsiCCSD_PT2R12() {
}

PsiCCSD_PT2R12::PsiCCSD_PT2R12(StateIn&s) :
  PsiCC_PT2R12(s) {
  s.get(eccsd_);
}

void PsiCCSD_PT2R12::save_data_state(StateOut&s) {
  PsiCC_PT2R12::save_data_state(s);
  s.put(eccsd_);
}

void PsiCCSD_PT2R12::write_input(int convergence) {
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiCorrWavefunction::write_input(convergence);
  input->write_keyword("psi:wfn", "ccsd");
  if (need_lambda_) input->write_keyword("psi:jobtype", "oeprop");
  // make sure Psi uses semicanonical orbitals for ROHF-CCSD (normally it would use spin-restricted orbitals)
  const bool openshell_ref = this->reference()->spin_polarized();
  if (openshell_ref)
    input->write_keyword("psi:semicanonical", "true");
  PsiCC_PT2R12::write_basic_input(convergence);
  input->close();
}

void PsiCCSD_PT2R12::compute() {

  PsiCC_PT2R12::compute_ept2r12();

  // read Psi3 CCSD energy
  {
    psi::PSIO& psio = exenv()->psio();
    psio.open(CC_INFO, PSIO_OPEN_OLD);
    psio.read_entry(CC_INFO, "CCSD Energy", reinterpret_cast<char*>(&eccsd_),
                    sizeof(double));
    psio.close(CC_INFO, 1);
  }

  const double e2 = value() - reference_energy();
  ExEnv::out0() << indent << "ECCSD         = "<< scprintf("%20.15lf",eccsd_)
                << std::endl;
  ExEnv::out0() << indent << "ECCSD_PT2R12  = "<< scprintf("%20.15lf",e2 + eccsd_)
                << std::endl;
  ExEnv::out0() << indent << "ECCSD_PT2R12+REF  = "<< scprintf("%20.15lf", reference_energy() + e2 + eccsd_)
                << std::endl;

  set_energy(eccsd_ + value());
}

void PsiCCSD_PT2R12::print(std::ostream&o) const {
  o << indent << "PsiCCSD_PT2R12:" << std::endl;
  o << incindent;
  PsiCC_PT2R12::print(o);
  o << decindent;
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_PT2R12T_cd(typeid(PsiCCSD_PT2R12T), "PsiCCSD_PT2R12T", 1,
                                    "public PsiCC_PT2R12", 0, create<PsiCCSD_PT2R12T>,
                                    create<PsiCCSD_PT2R12T>);

PsiCCSD_PT2R12T::PsiCCSD_PT2R12T(const Ref<KeyVal>&keyval) :
  PsiCC_PT2R12(keyval), eccsd_(NAN), e_t_(NAN) {
  if (!replace_Lambda_with_T_)
    throw FeatureNotImplemented("PsiCCSD_PT2R12T::PsiCCSD_PT2R12T() -- cannot properly use Lambdas yet",__FILE__,__LINE__);
}

PsiCCSD_PT2R12T::~PsiCCSD_PT2R12T() {
}

PsiCCSD_PT2R12T::PsiCCSD_PT2R12T(StateIn&s) :
  PsiCC_PT2R12(s) {
  s.get(eccsd_);
  s.get(e_t_);
}

void PsiCCSD_PT2R12T::save_data_state(StateOut&s) {
  PsiCC_PT2R12::save_data_state(s);
  s.put(eccsd_);
  s.put(e_t_);
}

void PsiCCSD_PT2R12T::write_input(int convergence) {
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiCorrWavefunction::write_input(convergence);
  input->write_keyword("psi:wfn", "ccsd_t");
  PsiCC_PT2R12::write_basic_input(convergence);
  input->close();
}

void PsiCCSD_PT2R12T::compute() {
  PsiCC_PT2R12::compute_ept2r12();
  const double e2 = value() - reference_energy();
  // read Psi3 CCSD energy
  {
    psi::PSIO& psio = exenv()->psio();
    psio.open(CC_INFO, PSIO_OPEN_OLD);
    psio.read_entry(CC_INFO, "CCSD Energy", reinterpret_cast<char*>(&eccsd_),
                    sizeof(double));
    psio.close(CC_INFO, 1);
    e_t_ = exenv()->chkpt().rd_e_t();
  }

  ExEnv::out0() << indent << "ECCSD         = "<< scprintf("%20.15lf",eccsd_)
                << std::endl;
  ExEnv::out0() << indent << "ECCSD_PT2R12  = "<< scprintf("%20.15lf",e2 + eccsd_)
                << std::endl;
  ExEnv::out0() << indent << "E(T)          = " << scprintf("%20.15lf",e_t_) << std::endl;
  ExEnv::out0() << indent << "ECCSD_PT2R12T = " << scprintf("%20.15lf",e_t_ + e2 + eccsd_)
                << std::endl;
  ExEnv::out0() << indent << "ECCSD_PT2R12T+REF = "<< scprintf("%20.15lf", reference_energy() + e_t_ + e2 + eccsd_)
                << std::endl;
  set_energy(e_t_ + e2 + eccsd_ + reference_energy());
}

void PsiCCSD_PT2R12T::print(std::ostream&o) const {
  o << indent << "PsiCCSD_PT2R12T:" << std::endl;
  o << incindent;
  PsiCC_PT2R12::print(o);
  o << decindent;
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCC3_PT2R12_cd(typeid(PsiCC3_PT2R12), "PsiCC3_PT2R12", 1,
                                  "public PsiCC", 0, create<PsiCC3_PT2R12>,
                                  create<PsiCC3_PT2R12>);

PsiCC3_PT2R12::PsiCC3_PT2R12(const Ref<KeyVal>&keyval) :
  PsiCC_PT2R12(keyval), ecc3_(NAN) {
  if (!replace_Lambda_with_T_)
    throw FeatureNotImplemented("PsiCC3_PT2R12::PsiCC3_PT2R12() -- cannot properly use Lambdas yet",__FILE__,__LINE__);
}

PsiCC3_PT2R12::~PsiCC3_PT2R12() {
}

PsiCC3_PT2R12::PsiCC3_PT2R12(StateIn&s) :
  PsiCC_PT2R12(s) {
  s.get(ecc3_);
}

void PsiCC3_PT2R12::save_data_state(StateOut&s) {
  PsiCC_PT2R12::save_data_state(s);
  s.put(ecc3_);
}

void PsiCC3_PT2R12::write_input(int convergence) {
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiCorrWavefunction::write_input(convergence);
  input->write_keyword("psi:wfn", "cc3");
  PsiCC_PT2R12::write_basic_input(convergence);
  input->close();
}

void PsiCC3_PT2R12::compute() {
  PsiCC_PT2R12::compute_ept2r12();
  const double e2 = value() - reference_energy();
  // read Psi3 CC3 energy
  {
    psi::PSIO& psio = exenv()->psio();
    psio.open(CC_INFO, PSIO_OPEN_OLD);
    psio.read_entry(CC_INFO, "CC3 Energy", reinterpret_cast<char*>(&ecc3_),
                    sizeof(double));
    psio.close(CC_INFO, 1);
  }

  ExEnv::out0() << indent << "ECC3_PT2R12   = " << scprintf("%20.15lf",e2 + ecc3_)
                << std::endl;
  ExEnv::out0() << indent << "ECC3_PT2R12+REF = "<< scprintf("%20.15lf", ecc3_ + value())
                << std::endl;
  set_energy(ecc3_ + value());
}

void PsiCC3_PT2R12::print(std::ostream&o) const {
  o << indent << "PsiCC3_PT2R12:" << std::endl;
  o << incindent;
  PsiCC_PT2R12::print(o);
  o << decindent;
}

////////////////////////

namespace {
  void _print(SpinCase2 spin,
             const Ref<DistArray4>& mat,
             const char* label) {
    if (mat->msg()->me() == 0) {
      const size_t nij = (spin != AlphaBeta && mat->ni() == mat->nj()) ? mat->ni() * (mat->ni()-1) / 2 : mat->ni() * mat->nj();
      const size_t nxy = (spin != AlphaBeta && mat->nx() == mat->ny()) ? mat->nx() * (mat->nx()-1) / 2 : mat->nx() * mat->ny();
      RefSCMatrix scmat = SCMatrixKit::default_matrixkit()->matrix(new SCDimension(nij), new SCDimension(nxy));
      scmat << mat;
      scmat.print(label);
    }
  }

} // anonymous namespace

