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

#include <ccfiles.h>
#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/utils.h>
#include <chemistry/qc/mbptr12/blas.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/mbptr12/utils.impl.h>
#include <chemistry/qc/mbptr12/pairiter.impl.h>
#include <chemistry/qc/psi/psicc_pt2r12.h>

#define TEST_ViaT1 1

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
                                const Ref<MOIndexSpace>& i, const Ref<MOIndexSpace>& a, const Ref<MOIndexSpace>& j);
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_PT2R12_cd(typeid(PsiCCSD_PT2R12), "PsiCCSD_PT2R12", 1,
                                   "public PsiCC", 0, create<PsiCCSD_PT2R12>,
                                   create<PsiCCSD_PT2R12>);

PsiCCSD_PT2R12::PsiCCSD_PT2R12(const Ref<KeyVal>&keyval) :
  PsiCC(keyval), eccsd_(1.0) {
  if (!replace_Lambda_with_T_)
    throw FeatureNotImplemented("PsiCCSD_PT2R12::PsiCCSD_PT2R12() -- cannot properly use Lambdas yet",__FILE__,__LINE__);
  
  if (mp2_only_)
    PsiCC::do_test_t2_phases();
  
  mbptr12_ = require_dynamic_cast<MBPT2_R12*>(keyval->describedclassvalue("mbpt2r12").pointer(), "PsiCCSD_PT2R12::PsiCCSD_PT2R12\n");
  
  const Ref<R12IntEvalInfo> r12info = mbptr12_->r12evalinfo();
  const Ref<R12Technology> r12tech = r12info->r12tech();
  // cannot do gbc = false yet
  if (!r12tech->gbc())
    throw FeatureNotImplemented("PsiCCSD_PT2R12::PsiCCSD_PT2R12() -- gbc = false is not yet implemented",__FILE__,__LINE__);
  // cannot do ebc=false either
  if (!r12tech->ebc())
    throw FeatureNotImplemented("PsiCCSD_PT2R12::PsiCCSD_PT2R12() -- ebc = false is not yet implemented",__FILE__,__LINE__);
  
}

PsiCCSD_PT2R12::~PsiCCSD_PT2R12() {
}

PsiCCSD_PT2R12::PsiCCSD_PT2R12(StateIn&s) :
  PsiCC(s) {
  mbptr12_ << SavableState::restore_state(s);
  s.get(eccsd_);
}

int PsiCCSD_PT2R12::gradient_implemented() const {
  return 0;
}

void PsiCCSD_PT2R12::save_data_state(StateOut&s) {
  PsiCC::save_data_state(s);
  SavableState::save_state(mbptr12_.pointer(), s);
  s.put(eccsd_);
}

void PsiCCSD_PT2R12::write_input(int convergence) {
  Ref<PsiInput> input = get_psi_input();
  input->open();
  PsiCorrWavefunction::write_input(convergence);
  // if testing T2 transform, obtain MP1 amplitudes
  if (test_t2_phases_)
    input->write_keyword("psi:wfn", "mp2");
  else
    input->write_keyword("psi:wfn", "ccsd");
  input->close();
}

void PsiCCSD_PT2R12::compute() {
  const Ref<R12IntEvalInfo> r12info = mbptr12_->r12evalinfo();
  const Ref<R12Technology> r12tech = r12info->r12tech();
  // params
  const bool ebc = r12tech->ebc();
  
  // compute CCSD wave function
  PsiWavefunction::compute();
  // read CCSD energy
  if (!mp2_only_) {
    psi::PSIO& psio = exenv()->psio();
    psio.open(CC_INFO, PSIO_OPEN_OLD);
    psio.read_entry(CC_INFO, "CCSD Energy", reinterpret_cast<char*>(&eccsd_),
                    sizeof(double));
    psio.close(CC_INFO, 1);
  }
  
  // Compute intermediates
  const double mp2r12_energy = mbptr12_->value();
  Ref<R12IntEval> r12eval = mbptr12_->r12eval();
  RefSCMatrix Vpq[NSpinCases2];
  RefSCMatrix Vab[NSpinCases2];
  RefSCMatrix Via[NSpinCases2];
  RefSCMatrix Vai[NSpinCases2];
  RefSCMatrix Vij[NSpinCases2];
  RefSymmSCMatrix B[NSpinCases2];
  RefSymmSCMatrix X[NSpinCases2];
  RefSCMatrix T2_MP1[NSpinCases2];
  
  const int nspincases2 = r12eval->nspincases2();
  for (int s=0; s<nspincases2; s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    if (r12eval->dim_oo(spincase2).n() == 0)
      continue;

    Ref<R12IntEvalInfo> r12info = r12eval->r12info();
    const Ref<MOIndexSpace>& p1 = r12info->refinfo()->orbs(spin1);
    const Ref<MOIndexSpace>& p2 = r12info->refinfo()->orbs(spin2);
    const unsigned int np1 = p1->rank();
    const unsigned int np2 = p2->rank();
    
    const Ref<MOIndexSpace>& occ1_act = r12eval->occ_act(spin1);
    const Ref<MOIndexSpace>& occ2_act = r12eval->occ_act(spin2);
    const Ref<MOIndexSpace>& vir1_act = r12eval->vir_act(spin1);
    const Ref<MOIndexSpace>& vir2_act = r12eval->vir_act(spin2);
    
    Vpq[s] = r12eval->V(spincase2, p1, p2);
    Vij[s] = r12eval->V(spincase2);
    X[s] = r12eval->X(spincase2);
    B[s] = r12eval->B(spincase2);
    if (test_t2_phases_)
      T2_MP1[s] = r12eval->T2(spincase2);
    
    const Ref<MOIndexSpace>& v1 = r12eval->vir_act(spin1);
    const Ref<MOIndexSpace>& v2 = r12eval->vir_act(spin2);
    const unsigned int nv1 = v1->rank();
    const unsigned int nv2 = v2->rank();
    const Ref<MOIndexSpace>& o1 = r12eval->occ_act(spin1);
    const Ref<MOIndexSpace>& o2 = r12eval->occ_act(spin2);
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
    if (spincase2 != AlphaBeta)
      xypq_to_xyab<AntiSymm,AntiSymm>(Vpq[s],Vab[s],v1_to_p1,v2_to_p2,np1,np2);
    else
      xypq_to_xyab<ASymm,ASymm>(Vpq[s],Vab[s],v1_to_p1,v2_to_p2,np1,np2);
    
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
    if (spincase2 == AlphaBeta && r12info->refinfo()->spin_polarized()) {
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
      if (test_t2_phases_)
        T2_MP1[s].print("MP1 T2 amplitudes");
    }
    if (debug() >= DefaultPrintThresholds::mostO4) {
      Vij[s].print("Vij matrix");
    }
#if TEST_V
    RefSCMatrix Vab_test = r12eval->V(spincase2,vir1_act,vir2_act);
    Vab_test.print("Vab matrix (test)");
    (Vab[s] - Vab_test).print("Vab - Vab (test): should be 0");
#endif
#if TEST_V
    RefSCMatrix Via_test = r12eval->V(spincase2,occ1_act,vir2_act);
    Via_test.print("Via matrix (test)");
    (Via[s] - Via_test).print("Via - Via (test): should be 0");
#endif
  } // end of spincase2 loop
  Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

  //
  // Convert CC amplitudes to MPQC orbitals
  //
  RefSCMatrix T1[NSpinCases1];
  RefSCMatrix T2[NSpinCases2];
  RefSCMatrix Tau2[NSpinCases2];

  // Form transforms and transform T1
  RefSCMatrix MPQC2PSI_tform_oa[NSpinCases1];
  RefSCMatrix MPQC2PSI_tform_va[NSpinCases1];
  const int nspincases1 = r12eval->nspincases1();
  for(int s=0; s<nspincases1; ++s) {
    const SpinCase1 spin = static_cast<SpinCase1>(s);
    // print out MPQC orbitals to compare to Psi orbitals below;
    const Ref<MOIndexSpace>& orbs_sb_mpqc = r12eval->r12info()->refinfo()->orbs_sb(spin);
    if (debug() >= DefaultPrintThresholds::mostN2) {
      orbs_sb_mpqc->coefs().print(prepend_spincase(spin,"MPQC eigenvector").c_str());
      orbs_sb_mpqc->evals().print(prepend_spincase(spin,"MPQC eigenvalues").c_str());
    }
    const Ref<MOIndexSpace>& occ_act = r12eval->occ_act(spin);
    const Ref<MOIndexSpace>& vir_act = r12eval->vir_act(spin);
    
    // Psi stores amplitudes in Pitzer (symmetry-blocked) order. Construct such spaces for active occupied and virtual spaces
    Ref<MOIndexSpace> occ_act_sb_psi = occ_act_sb(spin);
    Ref<MOIndexSpace> vir_act_sb_psi = vir_act_sb(spin);
    
    MPQC2PSI_tform_oa[spin] = transform(*occ_act_sb_psi, *occ_act,
                                        localkit);
    MPQC2PSI_tform_va[spin] = transform(*vir_act_sb_psi, *vir_act,
                                        localkit);
    RefSCMatrix T1_psi = this->T1(spin);
    T1[spin] = transform_T1(MPQC2PSI_tform_oa[spin], MPQC2PSI_tform_va[spin], T1_psi,
                            localkit);
    if (debug() >= DefaultPrintThresholds::mostN2) {
      T1[spin].print(prepend_spincase(spin,"CCSD T1 amplitudes from Psi3 (in MPQC orbitals, obtained by transform):").c_str());
    }
  }
  if (nspincases1 == 1) {
    MPQC2PSI_tform_oa[Beta] = MPQC2PSI_tform_oa[Alpha];
    MPQC2PSI_tform_va[Beta] = MPQC2PSI_tform_va[Alpha];
  }

#define TRY_USING_SPARSEMAP 0
#if TRY_USING_SPARSEMAP
  // map MPQC (energy-ordered) orbitals to Psi (symmetry-ordered) orbitals
  SparseMOIndexMap MPQC2PSI_occ1_act;
  SparseMOIndexMap MPQC2PSI_vir1_act;
  bool use_sparsemap = true;
  {
    bool can_map_occ_act = true;
    bool can_map_vir_act = true;
    try {
      SparseMOIndexMap tmpmap_o(sparsemap(*occ1_act_sb_psi,*occ1_act,1e-9));
      std::swap(MPQC2PSI_occ1_act,tmpmap_o);
    }
    catch (CannotConstructMap& e) {
      can_map_occ_act = false;
    }
    try {
      SparseMOIndexMap tmpmap_v(sparsemap(*vir1_act_sb_psi,*vir1_act,1e-9));
      std::swap(MPQC2PSI_vir1_act,tmpmap_v);
    }
    catch (CannotConstructMap& e) {
      can_map_vir_act = false;
    }
    
    use_sparsemap = use_sparsemap && can_map_occ_act && can_map_vir_act;
    if (use_sparsemap_only_) {
      if (!can_map_occ_act)
        throw ProgrammingError("PsiCCSD_PT2R12::compute() -- cannot map MPQC occ act. orbitals to Psi occ act. orbitals",__FILE__,__LINE__);
      if (!can_map_vir_act)
        throw ProgrammingError("PsiCCSD_PT2R12::compute() -- cannot map MPQC vir act. orbitals to Psi vir act. orbitals",__FILE__,__LINE__);
    }
  }
  {
    bool can_map_occ_act = true;
    bool can_map_vir_act = true;
    try {
      SparseMOIndexMap tmpmap_o(sparsemap(*occ2_act_sb_psi,*occ2_act,1e-9));
      std::swap(MPQC2PSI_occ2_act,tmpmap_o);
    }
    catch (CannotConstructMap& e) {
      can_map_occ_act = false;
    }
    try {
      SparseMOIndexMap tmpmap_v(sparsemap(*vir2_act_sb_psi,*vir2_act,1e-9));
      std::swap(MPQC2PSI_vir2_act,tmpmap_v);
    }
    catch (CannotConstructMap& e) {
      can_map_vir_act = false;
    }

    use_sparsemap = use_sparsemap && can_map_occ_act && can_map_vir_act;
    if (use_sparsemap_only_) {
      if (!can_map_occ_act)
      throw ProgrammingError("PsiCCSD_PT2R12::compute() -- cannot map MPQC occ act. orbitals to Psi occ act. orbitals",__FILE__,__LINE__);
      if (!can_map_vir_act)
      throw ProgrammingError("PsiCCSD_PT2R12::compute() -- cannot map MPQC vir act. orbitals to Psi vir act. orbitals",__FILE__,__LINE__);
    }
  }
  
  if (use_sparsemap) {
    T1[spin1] = transform_T1(MPQC2PSI_occ1_act, MPQC2PSI_vir1_act, T1_psi_1,
                             localkit);
    T1[spin2] = transform_T1(MPQC2PSI_occ2_act, MPQC2PSI_vir2_act, T1_psi_2,
                             localkit);
    T2[spincase2] = transform_T2(MPQC2PSI_occ1_act, MPQC2PSI_occ2_act,
                                 MPQC2PSI_vir1_act, MPQC2PSI_vir2_act, T2_psi,
                                 localkit);
    Tau2[spincase2] = transform_T2(MPQC2PSI_occ1_act, MPQC2PSI_occ2_act,
                                   MPQC2PSI_vir1_act, MPQC2PSI_vir2_act,
                                   Tau2_psi, localkit);
    if (test_t2_phases_) {
      compare_T2(T2[spincase2], T2_MP1[spincase2], spincase2, occ1_act->rank(),
                 occ2_act->rank(), vir1_act->rank(), vir2_act->rank());
    }
    if (debug() >= DefaultPrintThresholds::mostO2N2) {
      T1[spin1].print(prepend_spincase(spin1,"CCSD T1 amplitudes from Psi3 (in MPQC orbitals):").c_str());
      T1[spin2].print(prepend_spincase(spin2,"CCSD T1 amplitudes from Psi3 (in MPQC orbitals):").c_str());
      T2[spincase2].print(prepend_spincase(spincase2,"CCSD T2 amplitudes from Psi3 (in MPQC orbitals):").c_str());
      Tau2[spincase2].print(prepend_spincase(spincase2,"CCSD Tau2 amplitudes from Psi3 (in MPQC orbitals):").c_str());
    }
  }
  else
#endif

  for(int s=0; s<nspincases2; ++s) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    
    if (r12eval->dim_oo(spincase2).n() == 0)
      continue;

    const Ref<MOIndexSpace>& occ1_act = r12eval->occ_act(spin1);
    const Ref<MOIndexSpace>& occ2_act = r12eval->occ_act(spin2);
    const Ref<MOIndexSpace>& vir1_act = r12eval->vir_act(spin1);
    const Ref<MOIndexSpace>& vir2_act = r12eval->vir_act(spin2);
    
    // grab CC amplitudes
    RefSCMatrix T2_psi = this->T2(spincase2);
    T2[spincase2] = transform_T2(MPQC2PSI_tform_oa[spin1], MPQC2PSI_tform_oa[spin2],
                                 MPQC2PSI_tform_va[spin1], MPQC2PSI_tform_va[spin2], T2_psi,
                                 localkit);
    // same-spin T2 from Psi had unrestricted indices to make the tranform easier
    // convert to i<j, a<b by antisymmetrizing and scaling by 0.5
    if (spincase2 != AlphaBeta) {
      RefSCMatrix T2_tmp = T2[spincase2].clone();
      T2_tmp.assign(T2[spincase2]);
      T2[spincase2] = 0;
      T2[spincase2] = localkit->matrix(r12eval->dim_oo(spincase2),r12eval->dim_vv(spincase2));
      antisymmetrize<false>(T2[spincase2],T2_tmp,occ1_act,occ2_act,vir1_act,vir2_act);
      T2[spincase2].scale(0.5);
    }
    
    // Tau = T2 + T1*T1 has same 1st through 3rd order contributions as T2
    if (mp2_only_ || completeness_order_for_intermediates_ < 4) {
      Tau2[spincase2] = T2[spincase2].clone();
      Tau2[spincase2].assign(T2[spincase2]);
    }
    else {
      RefSCMatrix Tau2_psi = this->Tau2(spincase2);
      Tau2[spincase2] = transform_T2(MPQC2PSI_tform_oa[spin1], MPQC2PSI_tform_oa[spin2],
                                     MPQC2PSI_tform_va[spin1], MPQC2PSI_tform_va[spin2], Tau2_psi,
                                     localkit);
      // see comment above for T2
      if (spincase2 != AlphaBeta) {
        RefSCMatrix Tau2_tmp = Tau2[spincase2].clone();
        Tau2_tmp.assign(Tau2[spincase2]);
        Tau2[spincase2] = 0;
        Tau2[spincase2] = localkit->matrix(r12eval->dim_oo(spincase2),r12eval->dim_vv(spincase2));
        antisymmetrize<false>(Tau2[spincase2],Tau2_tmp,occ1_act,occ2_act,vir1_act,vir2_act);
        Tau2[spincase2].scale(0.5);
      }
    }

    if (test_t2_phases_) {
      const Ref<MOIndexSpace>& occ1_act = r12eval->occ_act(spin1);
      const Ref<MOIndexSpace>& vir1_act = r12eval->vir_act(spin1);
      const Ref<MOIndexSpace>& occ2_act = r12eval->occ_act(spin2);
      const Ref<MOIndexSpace>& vir2_act = r12eval->vir_act(spin2);
      compare_T2(T2[spincase2], T2_MP1[spincase2], spincase2, occ1_act->rank(),
                 occ2_act->rank(), vir1_act->rank(), vir2_act->rank());
    }
    if (debug() >= DefaultPrintThresholds::mostO2N2) {
      T2[spincase2].print(prepend_spincase(spincase2,"CCSD T2 amplitudes from Psi3 (in MPQC orbitals, obtained by transform):").c_str());
      Tau2[spincase2].print(prepend_spincase(spincase2,"CCSD Tau2 amplitudes from Psi3 (in MPQC orbitals, obtained by transform):").c_str());
    }
  }
  
  // Compute Hamiltonian matrix elements
  RefSCMatrix H1_0R[NSpinCases2];
  RefSCMatrix H1_R0[NSpinCases2];
  RefSymmSCMatrix H0_RR[NSpinCases2];
  for (int s=0; s<nspincases2; s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);
    if (r12eval->dim_oo(spincase2).n() == 0)
      continue;

    const Ref<MOIndexSpace>& occ1_act = r12eval->occ_act(spin1);
    const Ref<MOIndexSpace>& occ2_act = r12eval->occ_act(spin2);
    const Ref<MOIndexSpace>& vir1_act = r12eval->vir_act(spin1);
    const Ref<MOIndexSpace>& vir2_act = r12eval->vir_act(spin2);
    const bool p1_equiv_p2 = (occ1_act == occ2_act) && (vir1_act == vir2_act);
    
    
    // H1_0R is just Vij
    H1_0R[s] = Vij[s].clone();
    H1_0R[s].assign(Vij[s]);
    if (debug() >= DefaultPrintThresholds::O2N2) {
      H1_0R[s].print(prepend_spincase(spincase2,"<0|Hb|R>").c_str());
    }
    
    // Vij is also the leading term in H1_R0
    H1_R0[s] = Vij[s].clone();
    H1_R0[s].assign(Vij[s]);
    
    // if not testing MP2, compute other terms in <R|Hb|0>
    if (!mp2_only_) {
      
      // the leading term in <R|(HT)|0> is Tau2.Vab
      RefSCMatrix VT2 = Vab[s] * Tau2[s].t();
      // If using the symmetric form of the correction, Lambdas are replaced with T
      if (replace_Lambda_with_T_)
        VT2.scale(2.0);
      RefSCMatrix HT = VT2;  VT2 = 0;
      H1_R0[s].accumulate(HT);

      // the next term is T1.Vai. It's third-order if BC hold, second-order otherwise
      if ( completeness_order_for_intermediates_ >= 3 ||
          (completeness_order_for_intermediates_ >= 2 && !r12eval->r12info()->bc()) ) {

        // Via . T1
        RefSCMatrix VT1_2 = contract_Via_T1ja(true,Via[s],T1[spin2],occ1_act,vir2_act,occ2_act);
        VT1_2.print("VT1_2");
        RefSCMatrix VT1 = VT1_2.clone(); VT1.assign(VT1_2);
        if (p1_equiv_p2) {
          // NOTE: V_{xy}^{ia} T_a^j + V_{xy}^{aj} T_a^i = 2 * symm(V_{xy}^{ia} T_a^j
          VT1.scale(2.0);
          const Ref<MOIndexSpace>& xspace1 = r12eval->xspace(spin1);
          const Ref<MOIndexSpace>& xspace2 = r12eval->xspace(spin2);
          if (spincase2 == AlphaBeta)
            symmetrize12<false,ASymm,ASymm>(VT1, VT1, xspace1, xspace2, occ1_act, occ2_act);
          else
            symmetrize12<false,AntiSymm,ASymm>(VT1, VT1, xspace1, xspace2, occ1_act, occ2_act);
        }
        else {
          // Vai^t . T1
          RefSCMatrix VT1_1 = contract_Via_T1ja(false,Vai[s],T1[spin1],occ2_act,vir1_act,occ1_act);
          VT1_1.print("VT1_1");
          VT1.accumulate(VT1_1);
        }
        // If needed, antisymmetrize VT1
        if (spincase2 != AlphaBeta) {
          RefSCMatrix VT1Anti = HT.clone();
          const Ref<MOIndexSpace>& xspace1 = r12eval->xspace(spin1);
          const Ref<MOIndexSpace>& xspace2 = r12eval->xspace(spin2);
          symmetrize<false,AntiSymm,ASymm,AntiSymm,AntiSymm>(VT1Anti,VT1,xspace1,xspace2,occ1_act,occ2_act);
          VT1 = VT1Anti;
        }
        if (debug() >= DefaultPrintThresholds::mostO2N2) {
          VT1.print("Via.T1 + T1.Vai");
        }
        
        // test against V evaluator
#if TEST_ViaT1
        if (spincase2 == AlphaBeta) {
          RefSCMatrix VIj, ViJ;
          RefSCMatrix vir2_act_coefs = vir2_act->coefs();
          RefSCMatrix to_t1;
          {
            RefSCDimension rowdim = vir2_act->dim();
            RefSCDimension coldim = occ2_act->dim();
            to_t1 = vir2_act_coefs.kit()->matrix(rowdim, coldim);
            double* tmp = new double[rowdim.n() * coldim.n()];
            T1[spin2].t().convert(tmp);
            to_t1.assign(tmp);
            delete[] tmp;
          }
          RefSCMatrix Tocc2_act_coefs = vir2_act_coefs * to_t1;
          Ref<MOIndexSpace> Tocc2_act = new MOIndexSpace("i(t1)", "active occupied orbitals (weighted by T1)",
              occ2_act, Tocc2_act_coefs, occ2_act->basis());
          ViJ = r12eval->V(spincase2, occ1_act, Tocc2_act);
          if (p1_equiv_p2) {
            // NOTE: V_{xy}^{ia} T_a^j + V_{xy}^{aj} T_a^i = 2 * symm(V_{xy}^{ia} T_a^j
            ViJ.scale(2.0);
            symmetrize<false>(ViJ, ViJ, r12eval->xspace(spin1), occ1_act);
          }
          else {
            RefSCMatrix vir1_act_coefs = vir1_act->coefs();
            RefSCMatrix to_t1;
            {
              RefSCDimension rowdim = vir1_act->dim();
              RefSCDimension coldim = occ1_act->dim();
              to_t1 = vir1_act_coefs.kit()->matrix(rowdim, coldim);
              double* tmp = new double[rowdim.n() * coldim.n()];
              T1[spin1].t().convert(tmp);
              to_t1.assign(tmp);
              delete[] tmp;
            }
            RefSCMatrix Tocc1_act_coefs = vir1_act_coefs * to_t1;
            Ref<MOIndexSpace> Tocc1_act = new MOIndexSpace("I(t1)", "active occupied orbitals (weighted by T1)",
                occ1_act, Tocc1_act_coefs, occ1_act->basis());
            VIj = r12eval->V(spincase2, Tocc1_act, occ2_act);
          }
          (VIj+ViJ).print("Via.T1 computed with R12IntEvalInfo::V()");
          VT1.print("Via.T1 + T1.Vai");
          //(VT1 - ViJ - VIj).print("Via.T1 - Via.T1 (test): should be zero");
        }
#endif

        // If using the symmetric form of the correction, Lambdas are replaced with T
        if (replace_Lambda_with_T_)
          VT1.scale(2.0);
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
    
    // H0_RR is just B
    H0_RR[s] = B[s].clone();
    H0_RR[s].assign(B[s]);
  }
  // Make H1_R0[AlphaAlpha] if this is a closed-shell case (it isn't computed then)
  if (nspincases2 == 1) {
    H1_R0[AlphaAlpha] = r12eval->V(AlphaAlpha).clone();
    antisymmetrize(H1_R0[AlphaAlpha], H1_R0[AlphaBeta], r12eval->xspace(Alpha),
                   r12eval->occ_act(Alpha));
  }
  
  // compute the second-order correction: E2 = - H1_0R . H0_RR^{-1} . H1_R0 = C_MP1 . H1_R0
  const int debug = 0;
  Ref<MP2R12Energy> r12energy = new MP2R12Energy(r12eval,r12eval->r12info()->r12tech()->stdapprox(),debug);
  std::vector<double> E2(NSpinCases2,0.0);
  // WARNING only RHF and UHF are considered
  const int num_unique_spincases2 = (reference_->spin_polarized() ? 3 : 2);
  for (int s=0; s<num_unique_spincases2; s++) {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(s);
    if (r12eval->dim_oo(spincase2).n() == 0)
      continue;
    
    RefSCMatrix C_MP1 = r12energy->C(spincase2);
    RefSCMatrix E2_mat = C_MP1.t() * H1_R0[s];
    E2[s] = E2_mat.trace();
  }
  
  double e2;
  ExEnv::out0() << "E2(AB)        = "<< scprintf("%15.10lf",E2[AlphaBeta])
                << endl;
  if (num_unique_spincases2 > 2) {
    e2 = E2[AlphaBeta] + E2[AlphaAlpha] + E2[BetaBeta];
    ExEnv::out0() << "E2(BB)        = "<< scprintf("%15.10lf",E2[BetaBeta])
                  << endl;
    ExEnv::out0() << "E2(AA)        = "<< scprintf("%15.10lf",E2[AlphaAlpha])
                  << endl;
  }
  else {
    e2 = E2[AlphaBeta] + 2.0*E2[AlphaAlpha];
    ExEnv::out0() << "E2(AA)        = "<< scprintf("%15.10lf",2.0*E2[AlphaAlpha])
                  << endl;
    ExEnv::out0() << "E2(s)         = "<< scprintf("%15.10lf",E2[AlphaBeta] - E2[AlphaAlpha])
                  << endl;
    ExEnv::out0() << "E2(t)         = "<< scprintf("%15.10lf",3.0*E2[AlphaAlpha])
                  << endl;
  }
  
  ExEnv::out0() << "E2            = "<< scprintf("%15.10lf",e2)
                << endl;
  ExEnv::out0() << "ECCSD         = "<< scprintf("%15.10lf",eccsd_)
                << endl;
  ExEnv::out0() << "ECCSD_PT2R12  = "<< scprintf("%15.10lf",e2 + eccsd_)
                << endl;
  ExEnv::out0() << "E2(MP2)       = "<< scprintf("%15.10lf",mbptr12_->r12_corr_energy())
                << endl;
  ExEnv::out0() << "EMP2          = "<< scprintf("%15.10lf",mbptr12_->corr_energy() - mbptr12_->r12_corr_energy())
                << endl;
  ExEnv::out0() << "EMP2R12       = "<< scprintf("%15.10lf",mbptr12_->corr_energy())
                << endl;
  
  set_energy(reference_energy() + e2 + eccsd_);
}

//////////////////////////////////////////////////////////////////////////

static ClassDesc PsiCCSD_PT2R12T_cd(typeid(PsiCCSD_PT2R12T), "PsiCCSD_PT2R12T", 1,
                                    "public PsiCC", 0, create<PsiCCSD_PT2R12T>,
                                    create<PsiCCSD_PT2R12T>);

PsiCCSD_PT2R12T::PsiCCSD_PT2R12T(const Ref<KeyVal>&keyval) :
  PsiCCSD_PT2R12(keyval), e_t_(1.0) {
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
  // if testing T2 transform, obtain MP1 amplitudes
  if (test_t2_phases_)
    input->write_keyword("psi:wfn", "mp2");
  else
    input->write_keyword("psi:wfn", "ccsd_t");
  input->close();
}

void PsiCCSD_PT2R12T::compute() {
  PsiCCSD_PT2R12::compute();
  const double e_ccsd_pt2r12 = value() - reference_energy();
  e_t_ = exenv()->chkpt().rd_e_t();
  const double e_ccsd_pt2r12t = e_ccsd_pt2r12 + e_t_;
  ExEnv::out0() << "E(T)          = " << scprintf("%15.10lf",e_t_) << endl;
  ExEnv::out0() << "ECCSD_PT2R12T = " << scprintf("%15.10lf",e_ccsd_pt2r12t)
                << endl;
  set_energy(e_ccsd_pt2r12t);
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
                                const Ref<MOIndexSpace>& i, const Ref<MOIndexSpace>& a, const Ref<MOIndexSpace>& j)
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

