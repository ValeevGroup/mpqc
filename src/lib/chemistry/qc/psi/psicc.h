//
// psicc.h
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

#ifndef _chemistry_qc_psi_psicc_h
#define _chemistry_qc_psi_psicc_h

#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/wfn/spin.h>
#include <math/distarray4/distarray4.h>

namespace sc {
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCC is a Psi coupled cluster wave function

  class PsiCC : public PsiCorrWavefunction {

      Ref<OrbitalSpace> occ_act_sb_[NSpinCases1];
      Ref<OrbitalSpace> vir_act_sb_[NSpinCases1];
      RefSCMatrix T1_[NSpinCases1];
      RefSCMatrix T2_[NSpinCases2];
      Ref<DistArray4> T2_da4_[NSpinCases2];
      RefSCMatrix Tau2_[NSpinCases2];
      RefSCMatrix Lambda1_[NSpinCases1];
      RefSCMatrix Lambda2_[NSpinCases2];
      Ref<DistArray4> Lambda2_da4_[NSpinCases2];

      Ref<OrbitalSpace> occ_sb_[NSpinCases1];

    protected:
      // set to true if want to run only if Psi3 and MPQC orbitals match exactly up to a phase
      static const bool use_sparsemap_only_ = false;

      /// read in T1-like quantity of spincase spin using DPD label L
      RefSCMatrix T1(SpinCase1 spin, const std::string& L);
      /// read in T2-like quantity of spincase spin using DPD label L
      RefSCMatrix T2(SpinCase2 spin, const std::string& L);
      /// read in T2-like quantity of spincase spin using DPD label L
      Ref<DistArray4> T2_distarray4(SpinCase2 spin, const std::string& L);

      // get the T1 amplitudes (frozen core size) from
      // non frozen core PSI CCSD density calculation as compute_1rdm = true
      RefSCMatrix T1_fzc(SpinCase1 spin, const std::string& dpdlabel);
      // get the T2 amplitudes (frozen core size) from
      // non frozen core PSI CCSD density calculation as compute_1rdm = true
      Ref<DistArray4> T2_distarray4_fzc(SpinCase2 spin12,
                                        const std::string& dpdlabel);

      // read PSI3 ccsd one-particle density of spincase spin
      RefSCMatrix Onerdm(SpinCase1 spin);
      // read PSI3 relaxation 1e density of spincase spin
      RefSCMatrix Onerdm_relax_X(SpinCase1 spin);
      // test
      RefSCMatrix Onerdm_relax_D(SpinCase1 spin);
//      void compute_onerdm_relax(RefSCMatrix& Dorbs_alpha,
//                                RefSCMatrix& Dorbs_beta);

#if 0
      /// transform T1 to the new basis using sparse maps
      RefSCMatrix
          transform_T1(
                       const SparseMOIndexMap& occ_act_map,
                       const SparseMOIndexMap& vir_act_map,
                       const RefSCMatrix& T1,
                       const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit()) const;
      /// transform T2 to the new basis using sparse maps
      RefSCMatrix
          transform_T2(
                       const SparseMOIndexMap& occ1_act_map,
                       const SparseMOIndexMap& occ2_act_map,
                       const SparseMOIndexMap& vir1_act_map,
                       const SparseMOIndexMap& vir2_act_map,
                       const RefSCMatrix& T2,
                       const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit()) const;
      /// transform T1 to the new basis using dense transforms
      RefSCMatrix
          transform_T1(
                       const RefSCMatrix& occ_act_tform,
                       const RefSCMatrix& vir_act_tform,
                       const RefSCMatrix& T1,
                       const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit()) const;
      /// transform T2 to the new basis using dense transforms
      RefSCMatrix
          transform_T2(
                       const RefSCMatrix& occ1_act_tform,
                       const RefSCMatrix& occ2_act_tform,
                       const RefSCMatrix& vir1_act_tform,
                       const RefSCMatrix& vir2_act_tform,
                       const RefSCMatrix& T2,
                       const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit()) const;
#endif
      /// compare T2 and T2_ref (check that elements < zero are in the same place and elements > soft_zero have the same sign)
      void compare_T2(const RefSCMatrix& T2, const RefSCMatrix& T2_ref, SpinCase2 spin12,
                      unsigned int no1, unsigned int no2, unsigned int nv1,
                      unsigned int nv2, double zero = 1e-8) const;

      /// initialize Psi3 DPD library
      void dpd_start();
      /// stop Psi3 DPD library
      void dpd_stop();

      /// return active occupied orbital space (symmetry-blocked)
      const Ref<OrbitalSpace>& occ_act_sb(SpinCase1 spin);
      /// return active virtual orbital space (symmetry-blocked)
      const Ref<OrbitalSpace>& vir_act_sb(SpinCase1 spin);

      /// return occupied orbital space (symmetry-blocked)
      const Ref<OrbitalSpace>& occ_sb(SpinCase1 spin);

      int maxiter_;   //< max number of CC iterations
      static const int default_maxiter = 50;
      
      bool diis_;
      int diis_nvector_;

    public:
      PsiCC(const Ref<KeyVal>&);
      PsiCC(StateIn&);
      ~PsiCC();
      void save_data_state(StateOut&);

      /// return T amplitudes of rank 1. The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked).
      virtual const RefSCMatrix& T1(SpinCase1 spin1);
      /// return T amplitudes of rank 2. The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked).
      virtual const RefSCMatrix& T2(SpinCase2 spin2);
      /// return T amplitudes of rank 2. The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked)
      /// @param dpdkey libdpd key specifying the data; the default is "t" (seek T2 amplitudes), replace with "tau" to obtain Tau2 amplitudes (t2 + t1 ^ t1)
      virtual Ref<DistArray4> T2_da4(SpinCase2 spin2, std::string dpdkey = "t");
      /// return Tau2 amplitudes. The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked).
      virtual const RefSCMatrix& Tau2(SpinCase2 spin2);
      /// return Lambda amplitudes of rank 1
      virtual const RefSCMatrix& Lambda1(SpinCase1 spin1);
      /// return Lambda amplitudes of rank 2
      virtual const RefSCMatrix& Lambda2(SpinCase2 spin2);
      /// return T amplitudes of rank 2. The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked)
      virtual Ref<DistArray4> Lambda2_da4(SpinCase2 spin2);

      void obsolete();
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCCSD is a concrete implementation of Psi CCSD wave function

  class PsiCCSD : public PsiCC {
    protected:
      void write_input(int conv);
      double pccsd_alpha_;
      double pccsd_beta_;
      double pccsd_gamma_;
    public:
      PsiCCSD(const Ref<KeyVal>&);
      PsiCCSD(StateIn&);
      ~PsiCCSD();
      void save_data_state(StateOut&);
      bool analytic_gradient_implemented() const;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCCSD_T is a concrete implementation of Psi CCSD(T) wave function

  class PsiCCSD_T : public PsiCC {
    protected:
      void write_input(int conv);
    public:
      PsiCCSD_T(const Ref<KeyVal>&);
      PsiCCSD_T(StateIn&);
      ~PsiCCSD_T();

      void save_data_state(StateOut&);
  };

  ///////////////////////////////////////////////////////////////////
  /// PsiCC2 is a concrete implementation of Psi ground-state CC2 wave function
  class PsiCC2 : public PsiCC {
    protected:
      void write_input(int conv);
    public:
      PsiCC2(const Ref<KeyVal>&);
      PsiCC2(StateIn&);
      ~PsiCC2();

      void save_data_state(StateOut&);
  };
  

  ///////////////////////////////////////////////////////////////////
  /// PsiCC3 is a concrete implementation of Psi ground-state CC3 wave function
  class PsiCC3 : public PsiCC {
    protected:
      void write_input(int conv);
    public:
      PsiCC3(const Ref<KeyVal>&);
      PsiCC3(StateIn&);
      ~PsiCC3();

      void save_data_state(StateOut&);
  };

} // namespace

#endif /* header guard */
