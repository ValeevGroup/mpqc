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
#include <ccfiles.h>
#include <math/scmat/local.h>
#include <chemistry/qc/lcao/utils.h>
#include <util/misc/print.h>
#include <chemistry/qc/lcao/utils.impl.h>
#include <math/mmisc/pairiter.impl.h>
#include <chemistry/qc/psi/psicc_pt2r12.h>

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

    T2[s] = this->T2_distarray4(spincase2);

    if (debug() >= DefaultPrintThresholds::mostO2N2) {
      _print(spincase2, T2[spincase2], prepend_spincase(spincase2,"CCSD T2 amplitudes:").c_str());
    }
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

  } // end of diag = true clause
  else { // diag = false:
         // compute dressed V intermediate of CC(2)_R12 method here explicitly if non-diagonal ansatz is used
         // (diagonal MP2-R12 energy evaluator computes all intermediates)

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
      const Ref<OrbitalSpace>& x1 = r12eval()->xspace(spin1);
      const Ref<OrbitalSpace>& x2 = r12eval()->xspace(spin2);
      const Ref<OrbitalSpace>& v1 = r12eval()->vir_act(spin1);
      const Ref<OrbitalSpace>& v2 = r12eval()->vir_act(spin2);
      const Ref<OrbitalSpace>& o1 = r12eval()->occ_act(spin1);
      const Ref<OrbitalSpace>& o2 = r12eval()->occ_act(spin2);

      std::vector< Ref<DistArray4> > Vpq_vec = r12eval()->V_distarray4(spincase2, p1, p2);
      assert(Vpq_vec.size() == 1);
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

#if 0
      for (unsigned int xy=0; xy<nxy; ++xy) {
        unsigned int ab = 0;
        for (unsigned int a=0; a<nv1; ++a) {
          const unsigned int p = v1_to_p1[a];
          for (unsigned int b=0; b<nv1; ++b, ++ab) {
            const unsigned int q = v2_to_p2[b];
            const unsigned int pq = p*np2 + q;
            const double elem = Vpq[s].get_element(xy, pq);
            Vab[s].set_element(xy, ab, elem);
          }
        }
      }
#endif

      // Vai[AlphaBeta] is also needed if spin-polarized reference is used
      if (spincase2 == AlphaBeta && r12world()->refwfn()->spin_polarized())
        map(Vpq, x1, x2, p1, p2, Vai[s], x1, x2, v1, o2);

#if 0
      for (unsigned int xy=0; xy<nxy; ++xy) {
        unsigned int ia = 0;
        for (unsigned int i=0; i<no1; ++i) {
          const unsigned int p = o1_to_p1[i];
          for (unsigned int a=0; a<nv2; ++a, ++ia) {
            const unsigned int q = v1_to_p1[a];
            const unsigned int pq = p*np2 + q;
            const double elem = Vpq[s].get_element(xy, pq);
            Via[s].set_element(xy, ia, elem);
          }
        }
      }
#endif

      if (debug() >= DefaultPrintThresholds::mostO2N2) {
        _print(spincase2, Vpq, prepend_spincase(spincase2,"Vpq matrix").c_str());
        _print(spincase2, Vab[s], prepend_spincase(spincase2,"Vab matrix").c_str());
        _print(spincase2, Via[s], prepend_spincase(spincase2,"Via matrix").c_str());
        if (Vai[s].nonnull())
          _print(spincase2, Vai[s], prepend_spincase(spincase2,"Vai matrix").c_str());
        if (A[s].nonnull())
          _print(spincase2, A[s], prepend_spincase(spincase2,"A matrix").c_str());
      }
      if (debug() >= DefaultPrintThresholds::mostO4) {
        Vij[s].print(prepend_spincase(spincase2,"Vij matrix").c_str());
      }
#if TEST_V
      RefSCMatrix Vab_test = r12eval()->V(spincase2,vir1_act,vir2_act);
      Vab_test.print("Vab matrix (test)");
      (Vab[s] - Vab_test).print("Vab - Vab (test): should be 0");
#endif
#if TEST_V
      RefSCMatrix Via_test = r12eval()->V(spincase2,occ1_act,vir2_act);
      Via_test.print("Via matrix (test)");
      (Via[s] - Via_test).print("Via - Via (test): should be 0");
#endif
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

      const Ref<OrbitalSpace>& x1 = r12eval()->xspace(spin1);
      const Ref<OrbitalSpace>& x2 = r12eval()->xspace(spin2);
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
      antisymmetrize(H1_R0[AlphaAlpha], H1_R0[AlphaBeta], r12eval()->xspace(Alpha),
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

  // compute (2)_R12 energy as MP2-R12 energy with dressed V intermediate
  Ref<MP2R12Energy> r12energy = construct_MP2R12Energy(r12intermediates,debug(),diag);
  std::vector<double> E2(NSpinCases2,0.0);
  const int num_unique_spincases2 = (reference_->spin_polarized() ? 3 : 2);
  r12energy->compute();
  for (int s=0; s<num_unique_spincases2; s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;

    E2[s]=r12energy->ef12tot(spincase2);
  }

  // include cabs singles energy?
  cabs_singles_energy_ = cabs_singles_ ? r12eval()->emp2_cabs_singles() : 0.0;

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

