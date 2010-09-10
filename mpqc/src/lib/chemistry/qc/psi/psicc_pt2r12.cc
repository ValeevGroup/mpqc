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

#ifdef __GNUC__
#pragma implementation
#endif

#include <cmath>
#include <ccfiles.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/mbptr12/utils.impl.h>
#include <chemistry/qc/mbptr12/pairiter.impl.h>
#include <chemistry/qc/psi/psicc_pt2r12.h>

#define TEST_ViaT1 0

using namespace sc;
using namespace sc::fastpairiter;

namespace {
  /// convert <xy|pq> to <xy|ab> using sparse maps a->p and b->q
  template <sc::fastpairiter::PairSymm PSymm_pq, sc::fastpairiter::PairSymm PSymm_ab>
  void xypq_to_xyab(const RefSCMatrix& xypq,
                    RefSCMatrix& xyab,
                    const MOIndexMap& a_to_p,
                    const MOIndexMap& b_to_q,
                    const int np,
                    const int nq);
  RefSCMatrix contract_Via_T1ja(bool part2, const RefSCMatrix& V, const RefSCMatrix& T1,
                                const Ref<OrbitalSpace>& i, const Ref<OrbitalSpace>& a, const Ref<OrbitalSpace>& j);
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_PT2R12_cd(typeid(PsiCCSD_PT2R12), "PsiCCSD_PT2R12", 2,
                                   "public PsiCC", 0, create<PsiCCSD_PT2R12>,
                                   create<PsiCCSD_PT2R12>);

PsiCCSD_PT2R12::PsiCCSD_PT2R12(const Ref<KeyVal>&keyval) :
  PsiCC(keyval), eccsd_(NAN), cabs_singles_energy_(0.0) {
  if (!replace_Lambda_with_T_)
    throw FeatureNotImplemented("PsiCCSD_PT2R12::PsiCCSD_PT2R12() -- cannot properly use Lambdas yet",__FILE__,__LINE__);

  Ref<WavefunctionWorld> world = new WavefunctionWorld(keyval, this);
  Ref<RefWavefunction> refinfo = RefWavefunctionFactory::make(world, this->reference(), false,
                                                              this->nfzc(), this->nfzv());
  r12world_ = new R12WavefunctionWorld(keyval, refinfo);
  Ref<R12Technology> r12tech = r12world_->r12tech();
  cabs_singles_ = keyval->booleanvalue("cabs_singles",KeyValValueboolean((int)false));
  const bool openshell = this->reference()->spin_polarized();
  spinadapted_ = keyval->booleanvalue("spinadapted", KeyValValueboolean(openshell ? 0 : 1));

  r12eval_ = 0;
  mp2r12_energy_ = 0;

  // cannot do gbc = false yet
  if (!r12tech->gbc())
    throw FeatureNotImplemented("PsiCCSD_PT2R12::PsiCCSD_PT2R12() -- gbc = false is not yet implemented",__FILE__,__LINE__);
  // cannot do coupling=true either
  if (r12tech->coupling())
    throw FeatureNotImplemented("PsiCCSD_PT2R12::PsiCCSD_PT2R12() -- coupling = true is not yet implemented",__FILE__,__LINE__);

  set_desired_value_accuracy(desired_value_accuracy());
}

PsiCCSD_PT2R12::~PsiCCSD_PT2R12() {
  // need to manually break up a cycle of smart pointers
  R12WavefunctionWorld* r12world_ptr = r12world_.pointer();
  r12world_.clear();
  r12world_ptr->dereference();
  r12world_ptr->~R12WavefunctionWorld();
}

PsiCCSD_PT2R12::PsiCCSD_PT2R12(StateIn&s) :
  PsiCC(s) {
  if (this->class_version() < 2)
    throw ProgrammingError("cannot use archives older than version 2", __FILE__, __LINE__, class_desc());

  r12eval_ << SavableState::restore_state(s);
  r12world_ << SavableState::restore_state(s);
  mp2r12_energy_ << SavableState::restore_state(s);

  if((r12world()->r12tech()->ansatz()->orbital_product_GG()==R12Technology::OrbProdGG_pq) ||
     (r12world()->r12tech()->ansatz()->orbital_product_gg()==R12Technology::OrbProdgg_pq)) {
    throw InputError("PsiCCSD_PT2R12::PsiCCSD_PT2R12 -- pq Ansatz not allowed",__FILE__,__LINE__);
  }

  int spinadapted; s.get(spinadapted); spinadapted_ = (bool)spinadapted;
  s.get(cabs_singles_);
  s.get(cabs_singles_energy_);
  s.get(eccsd_);
}

int PsiCCSD_PT2R12::gradient_implemented() const {
  return 0;
}

void PsiCCSD_PT2R12::save_data_state(StateOut&s) {
  PsiCC::save_data_state(s);
  SavableState::save_state(r12eval_.pointer(),s);
  SavableState::save_state(r12world_.pointer(),s);
  SavableState::save_state(mp2r12_energy_.pointer(),s);

  s.put((int)spinadapted_);
  s.put(cabs_singles_);
  s.put(cabs_singles_energy_);
  s.put(eccsd_);
}

void PsiCCSD_PT2R12::write_input(int convergence) {
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiCorrWavefunction::write_input(convergence);
  input->write_keyword("psi:wfn", "ccsd");
  input->write_keyword("ccenergy:convergence", convergence);
  input->write_keyword("ccenergy:maxiter", maxiter_);
  input->close();
}

void PsiCCSD_PT2R12::compute() {
  // compute Psi3 CCSD wave function
  PsiCorrWavefunction::compute();
  // read Psi3 CCSD energy
  {
    psi::PSIO& psio = exenv()->psio();
    psio.open(CC_INFO, PSIO_OPEN_OLD);
    psio.read_entry(CC_INFO, "CCSD Energy", reinterpret_cast<char*>(&eccsd_),
                    sizeof(double));
    psio.close(CC_INFO, 1);
  }

  // to compute intermediates make sure r12eval_ is ready
  if (r12eval_.null()) {
    r12world()->world()->memory(this->memory_);
    r12eval_ = new R12IntEval(r12world());
    r12eval_->debug(debug_);
  }

  RefSCMatrix Vpq[NSpinCases2];
  RefSCMatrix Vab[NSpinCases2];
  RefSCMatrix Via[NSpinCases2];
  RefSCMatrix Vai[NSpinCases2];
  RefSCMatrix Vij[NSpinCases2];
  RefSymmSCMatrix B[NSpinCases2];
  RefSymmSCMatrix X[NSpinCases2];
  RefSCMatrix A[NSpinCases2];

  Ref<R12Technology> r12tech = r12world_->r12tech();

  const int nspincases2 = r12eval()->nspincases2();
  for (int s=0; s<nspincases2; s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;

    const Ref<OrbitalSpace>& p1 = r12eval()->orbs(spin1);
    const Ref<OrbitalSpace>& p2 = r12eval()->orbs(spin2);
    const unsigned int np1 = p1->rank();
    const unsigned int np2 = p2->rank();

    const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
    const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
    const Ref<OrbitalSpace>& vir1_act = r12eval()->vir_act(spin1);
    const Ref<OrbitalSpace>& vir2_act = r12eval()->vir_act(spin2);

    Vpq[s] = r12eval()->V(spincase2, p1, p2);
    Vij[s] = r12eval()->V(spincase2);
    X[s] = r12eval()->X(spincase2);
    B[s] = r12eval()->B(spincase2);

    const Ref<OrbitalSpace>& v1 = r12eval()->vir_act(spin1);
    const Ref<OrbitalSpace>& v2 = r12eval()->vir_act(spin2);
    const unsigned int nv1 = v1->rank();
    const unsigned int nv2 = v2->rank();
    const Ref<OrbitalSpace>& o1 = r12eval()->occ_act(spin1);
    const Ref<OrbitalSpace>& o2 = r12eval()->occ_act(spin2);
    const unsigned int no1 = o1->rank();
    const unsigned int no2 = o2->rank();
    const RefSCDimension v1v2dim = (spincase2 != AlphaBeta) ? pairdim<AntiSymm>(nv1,nv2) : pairdim<ASymm>(nv1,nv2);
    const RefSCDimension o1v2dim = pairdim<ASymm>(no1,nv2);

    // extract Vab and Via from Vpq
    typedef MOIndexMap OrbMap;
    OrbMap v1_to_p1(*p1<<*v1);
    OrbMap v2_to_p2(*p2<<*v2);
    OrbMap o1_to_p1(*p1<<*o1);
    OrbMap o2_to_p2(*p2<<*o2);
    Vab[s] = Vpq[s].kit()->matrix(Vpq[s].rowdim(), v1v2dim);
    Via[s] = Vpq[s].kit()->matrix(Vpq[s].rowdim(), o1v2dim);
    {
      if (spincase2 != AlphaBeta)
        xypq_to_xyab<AntiSymm,AntiSymm>(Vpq[s],Vab[s],v1_to_p1,v2_to_p2,np1,np2);
      else
        xypq_to_xyab<ASymm,ASymm>(Vpq[s],Vab[s],v1_to_p1,v2_to_p2,np1,np2);
      if (r12tech->ebc() == false) {
        A[s] = r12eval()->A(spincase2);
      }
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

    if (spincase2 != AlphaBeta)
      xypq_to_xyab<AntiSymm,ASymm>(Vpq[s],Via[s],o1_to_p1,v2_to_p2,np1,np2);
    else
      xypq_to_xyab<ASymm,ASymm>(Vpq[s],Via[s],o1_to_p1,v2_to_p2,np1,np2);
    // Vai[AlphaBeta] is also needed if spin-polarized reference is used
    if (spincase2 == AlphaBeta && r12world()->ref()->spin_polarized()) {
      const RefSCDimension v1o2dim = pairdim<ASymm>(nv1,no2);
      Vai[s] = Vpq[s].kit()->matrix(Vpq[s].rowdim(), v1o2dim);
      xypq_to_xyab<ASymm,ASymm>(Vpq[s],Vai[s],v1_to_p1,o2_to_p2,np1,np2);
      Vai[s].print("Vai");
    }

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
      Vpq[s].print("Vpq matrix");
      Vab[s].print("Vab matrix");
      Via[s].print("Via matrix");
      if (Vai[s].nonnull())
        Vai[s].print("Vai matrix");
    }
    if (debug() >= DefaultPrintThresholds::mostO4) {
      Vij[s].print("Vij matrix");
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

  //
  // Obtain CC amplitudes from Psi and copy into local matrices
  //
  RefSCMatrix T1[NSpinCases1];
  RefSCMatrix T2[NSpinCases2];

  // get T1
  const int nspincases1 = r12eval()->nspincases1();
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

  // get T2
  for(int s=0; s<nspincases2; ++s) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;

    RefSCMatrix T2_psi = this->T2(spincase2);
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;
    T2[s] = localkit->matrix(T2_psi.rowdim(), T2_psi.coldim());
    T2[s]->convert(T2_psi);

    if (debug() >= DefaultPrintThresholds::mostO2N2) {
      T2[spincase2].print(prepend_spincase(spincase2,"CCSD T2 amplitudes:").c_str());
    }
  }

  // Compute Hamiltonian matrix elements
  RefSCMatrix H1_R0[NSpinCases2];
  for (int s=0; s<nspincases2; s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;

    const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
    const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
    const Ref<OrbitalSpace>& vir1_act = r12eval()->vir_act(spin1);
    const Ref<OrbitalSpace>& vir2_act = r12eval()->vir_act(spin2);
    const bool p1_equiv_p2 = (occ1_act == occ2_act) && (vir1_act == vir2_act);


    // Vij is the leading term in H1_R0
    H1_R0[s] = Vij[s].clone();
    H1_R0[s].assign(Vij[s]);

    {

      // the leading term in <R|(HT)|0> is T2.Vab
      // if not assuming EBC then also include coupling matrix term
      RefSCMatrix VT2 = r12tech->ebc() ? Vab[s] * T2[s].t() : (Vab[s] + A[s]) * T2[s].t();
      // E = Vbar B^-1 Vbar, where Vbar = V + VT
      RefSCMatrix HT = VT2;  VT2 = 0;

      // the next term is T1.Vai. It's third-order if BC hold, second-order otherwise
      if ( completeness_order_for_intermediates_ >= 3 ||
          (completeness_order_for_intermediates_ >= 2 && !r12eval()->bc()) ) {

        // Via . T1
        RefSCMatrix VT1_2 = contract_Via_T1ja(true,Via[s],T1[spin2],occ1_act,vir2_act,occ2_act);
        RefSCMatrix VT1 = VT1_2.clone(); VT1.assign(VT1_2);
        //VT1_2.print("VT1_2");
        if (p1_equiv_p2) {
          // NOTE: V_{xy}^{ia} T_a^j + V_{xy}^{aj} T_a^i = 2 * symm(V_{xy}^{ia} T_a^j
          VT1_2.scale(2.0);
          const Ref<OrbitalSpace>& xspace1 = r12eval()->xspace(spin1);
          const Ref<OrbitalSpace>& xspace2 = r12eval()->xspace(spin2);
          if (spincase2 == AlphaBeta)
            symmetrize12<false,ASymm,ASymm>(VT1, VT1_2, xspace1, xspace2, occ1_act, occ2_act);
          else
            symmetrize12<false,AntiSymm,ASymm>(VT1, VT1_2, xspace1, xspace2, occ1_act, occ2_act);
        }
        else {
          // Vai^t . T1
          RefSCMatrix VT1_1 = contract_Via_T1ja(false,Vai[s],T1[spin1],occ2_act,vir1_act,occ1_act);
          //VT1_1.print("VT1_1");
          VT1.accumulate(VT1_1);
        }
        // If needed, antisymmetrize VT1
        if (spincase2 != AlphaBeta) {
          RefSCMatrix VT1Anti = HT.clone();
          const Ref<OrbitalSpace>& xspace1 = r12eval()->xspace(spin1);
          const Ref<OrbitalSpace>& xspace2 = r12eval()->xspace(spin2);
          symmetrize<false,AntiSymm,ASymm,AntiSymm,AntiSymm>(VT1Anti,VT1,xspace1,xspace2,occ1_act,occ2_act);
          VT1 = VT1Anti;
        }
        if (debug() >= DefaultPrintThresholds::mostO2N2) {
          VT1.print("Via.T1 + T1.Vai");
        }

        HT.accumulate(VT1);

      }

      if (debug() >= DefaultPrintThresholds::mostO2N2) {
        (HT).print(prepend_spincase(spincase2,"<R|(H*T)|0>").c_str());
      }
      H1_R0[s].accumulate(HT);
    }

    if (debug() >= DefaultPrintThresholds::O2N2) {
      H1_R0[s].print(prepend_spincase(spincase2,"<R|Hb|0>").c_str());
    }
  }
  // Make H1_R0[AlphaAlpha] if this is a closed-shell case (it isn't computed then)
  if (nspincases2 == 1) {
    H1_R0[AlphaAlpha] = r12eval()->V(AlphaAlpha).clone();
    antisymmetrize(H1_R0[AlphaAlpha], H1_R0[AlphaBeta], r12eval()->xspace(Alpha),
                   r12eval()->occ_act(Alpha));
  }

  // compute the second-order correction: E2 = - H1_0R . H0_RR^{-1} . H1_R0 = C . H1_R0
  // H0_RR is the usual B of the standard MP-R12 theory
  // if using screening approximations and replacing lambdas with Ts H1_0R = transpose(H1_R0)
  const bool diag = r12eval()->r12world()->r12tech()->ansatz()->diag();
  Ref<R12EnergyIntermediates> r12intermediates=new R12EnergyIntermediates(r12eval(),r12eval()->r12world()->r12tech()->stdapprox());
  // E2 = - Vbar . B^{-1} . Vbar, where Vbar = V + VT
  for(int i=0; i<NSpinCases2; i++) {
    SpinCase2 spincase2i=static_cast<SpinCase2>(i);
    if (r12eval()->dim_oo(spincase2i).n() == 0)
      continue;
    r12intermediates->assign_V(spincase2i,H1_R0[spincase2i]);
  }

  const bool new_energy=(diag==true) ? true : false;
  Ref<MP2R12Energy> r12energy = construct_MP2R12Energy(r12intermediates,debug(),new_energy);

  std::vector<double> E2(NSpinCases2,0.0);
  const int num_unique_spincases2 = (reference_->spin_polarized() ? 3 : 2);
  r12energy->compute();
  for (int s=0; s<num_unique_spincases2; s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    if (r12eval()->dim_oo(spincase2).n() == 0)
      continue;

    E2[s]=r12energy->ef12tot(spincase2);
  }

  double e2;
  ExEnv::out0() << indent << "E2(AB)        = "<< scprintf("%20.15lf",E2[AlphaBeta])
                << std::endl;
  if (num_unique_spincases2 > 2) {
    e2 = E2[AlphaBeta] + E2[AlphaAlpha] + E2[BetaBeta];
    ExEnv::out0() << indent << "E2(BB)        = "<< scprintf("%20.15lf",E2[BetaBeta])
                  << std::endl;
    ExEnv::out0() << indent << "E2(AA)        = "<< scprintf("%20.15lf",E2[AlphaAlpha])
                  << std::endl;
  }
  else {
    e2 = E2[AlphaBeta] + 2.0*E2[AlphaAlpha];
    ExEnv::out0() << indent << "E2(AA)        = "<< scprintf("%20.15lf",2.0*E2[AlphaAlpha])
                  << std::endl;
    ExEnv::out0() << indent << "E2(s)         = "<< scprintf("%20.15lf",E2[AlphaBeta] - E2[AlphaAlpha])
                  << std::endl;
    ExEnv::out0() << indent << "E2(t)         = "<< scprintf("%20.15lf",3.0*E2[AlphaAlpha])
                  << std::endl;
  }
  // e2 will include cabs singles energy, if mbpt2-r12 includes it
  e2 += cabs_singles_energy();

  ExEnv::out0() << indent << "E2            = "<< scprintf("%20.15lf",e2)
                << std::endl;
  ExEnv::out0() << indent << "ECCSD         = "<< scprintf("%20.15lf",eccsd_)
                << std::endl;
  ExEnv::out0() << indent << "ECCSD_PT2R12  = "<< scprintf("%20.15lf",e2 + eccsd_)
                << std::endl;
  ExEnv::out0() << indent << "ECCSD_PT2R12+REF = "<< scprintf("%20.15lf", reference_energy() + e2 + eccsd_)
                << std::endl;

  set_energy(reference_energy() + e2 + eccsd_);
}

double
PsiCCSD_PT2R12::cabs_singles_energy()
{
  return cabs_singles_energy_;
}

void PsiCCSD_PT2R12::print(std::ostream&o) const {
  o << indent << "PsiCCSD_PT2R12:" << std::endl;
  o << incindent;
  o << indent << "Spin-adapted algorithm: " << (spinadapted_ ? "true" : "false") << std::endl;
  o << indent << "Include CABS singles? : " << (cabs_singles_ ? "true" : "false") << std::endl;
  if (cabs_singles_) {
    o << indent << "  E(CABS singles) = " << scprintf("%25.15lf", cabs_singles_energy_)
                                          << std::endl;
  }

  r12world()->print(o);
  PsiCC::print(o);
  o << decindent;
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_PT2R12T_cd(typeid(PsiCCSD_PT2R12T), "PsiCCSD_PT2R12T", 1,
                                    "public PsiCC", 0, create<PsiCCSD_PT2R12T>,
                                    create<PsiCCSD_PT2R12T>);

PsiCCSD_PT2R12T::PsiCCSD_PT2R12T(const Ref<KeyVal>&keyval) :
  PsiCCSD_PT2R12(keyval), e_t_(NAN) {
  if (!replace_Lambda_with_T_)
    throw FeatureNotImplemented("PsiCCSD_PT2R12T::PsiCCSD_PT2R12T() -- cannot properly use Lambdas yet",__FILE__,__LINE__);
}

PsiCCSD_PT2R12T::~PsiCCSD_PT2R12T() {
}

PsiCCSD_PT2R12T::PsiCCSD_PT2R12T(StateIn&s) :
  PsiCCSD_PT2R12(s) {
  s.get(e_t_);
}

int PsiCCSD_PT2R12T::gradient_implemented() const {
  return 0;
}

void PsiCCSD_PT2R12T::save_data_state(StateOut&s) {
  PsiCCSD_PT2R12::save_data_state(s);
  s.put(e_t_);
}

void PsiCCSD_PT2R12T::write_input(int convergence) {
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiCorrWavefunction::write_input(convergence);
  input->write_keyword("psi:wfn", "ccsd_t");
  input->write_keyword("ccenergy:convergence", convergence);
  input->write_keyword("ccenergy:maxiter", maxiter_);
  input->close();
}

void PsiCCSD_PT2R12T::compute() {
  PsiCCSD_PT2R12::compute();
  const double e_ccsd_pt2r12 = value() - reference_energy();
  e_t_ = exenv()->chkpt().rd_e_t();
  const double e_ccsd_pt2r12t = e_ccsd_pt2r12 + e_t_;
  ExEnv::out0() << indent << "E(T)          = " << scprintf("%20.15lf",e_t_) << std::endl;
  ExEnv::out0() << indent << "ECCSD_PT2R12T = " << scprintf("%20.15lf",e_ccsd_pt2r12t)
                << std::endl;
  ExEnv::out0() << indent << "ECCSD_PT2R12T+REF = "<< scprintf("%20.15lf", reference_energy() + e_ccsd_pt2r12t)
                << std::endl;
  set_energy(e_ccsd_pt2r12t + reference_energy());
}

void PsiCCSD_PT2R12T::print(std::ostream&o) const {
  o << indent << "PsiCCSD_PT2R12T:" << std::endl;
  o << incindent;
  PsiCCSD_PT2R12::print(o);
  o << decindent;
}

////////////////////////

namespace {
  template <sc::fastpairiter::PairSymm PSymm_pq, sc::fastpairiter::PairSymm PSymm_ab>
  void xypq_to_xyab(const RefSCMatrix& xypq,
                    RefSCMatrix& xyab,
                    const MOIndexMap& a_to_p,
                    const MOIndexMap& b_to_q,
                    const int np,
                    const int nq)
  {
    const int nxy = xypq.rowdim().n();
    if (nxy != xyab.rowdim().n())
      throw ProgrammingError("xypq_to_xyab() -- row dimensions do not match",__FILE__,__LINE__);
    const int na = a_to_p.size();
    const int nb = b_to_q.size();

    typedef sc::fastpairiter::MOPairIter<PSymm_pq> PQIter;
    typedef sc::fastpairiter::MOPairIter<PSymm_ab> ABIter;
    PQIter pqiter(np,nq);
    ABIter abiter(na,nb);

    for (int xy=0; xy<nxy; ++xy) {
      for(abiter.start(); int(abiter); abiter.next()) {
        const int a = abiter.i();
        const int b = abiter.j();
        const int ab = abiter.ij();
        const int p = a_to_p[a];
        const int q = b_to_q[b];
        const int pq = pqiter.ij(p,q);
        const double elem = xypq.get_element(xy, pq);
        xyab.set_element(xy, ab, elem);
      }
    }
  }

  /// if part2 == false, VT1_ji = T1_ja . V_ai
  /// if part2 == true, VT1_ij = V_ia . (T_ja)^t
  /// note that the result is not symmetric with respect to permutation of particles 1 and 2
  RefSCMatrix contract_Via_T1ja(bool part2, const RefSCMatrix& V, const RefSCMatrix& T1,
                                const Ref<OrbitalSpace>& i, const Ref<OrbitalSpace>& a, const Ref<OrbitalSpace>& j)
  {
    const unsigned int ni = i->rank();
    const unsigned int nj = j->rank();
    const unsigned int na = a->rank();
    // allocate space for the result
    const unsigned int nij = ni*nj;
    RefSCMatrix R(V.rowdim(), new SCDimension(nij), V.kit());

    // Convert T1 to a raw matrix
    double* T1_raw = new double[j->rank() * a->rank()];
    T1.convert(T1_raw);

    const unsigned int nxy = V.rowdim().n();
    const unsigned int nia = V.coldim().n();
    double* raw_ia_blk = new double[V.coldim().n()];
    double* raw_ij_blk = new double[nij];
    RefSCVector Rrow(R.coldim(), R.kit());

    for(unsigned int xy=0; xy<nxy; ++xy) {
      RefSCVector row = V.get_row(xy);
      row.convert(raw_ia_blk);
      if (part2) {
        // V_ia . (T1_ja)^t : raw_ia_blk * T1_raw^t = raw_ij_blk
        C_DGEMM('n','t',ni,nj,na,1.0,raw_ia_blk,na,T1_raw,na,0.0,raw_ij_blk,nj);
      }
      else {
        // T1_ja . Vai : T1_raw * raw_ia_blk = raw_ij_blk
        C_DGEMM('n','n',nj,ni,na,1.0,T1_raw,na,raw_ia_blk,ni,0.0,raw_ij_blk,ni);
      }
      Rrow.assign(raw_ij_blk);
      R.assign_row(Rrow,xy);
    }

    return R;
  }

} // anonymous namespace

