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

#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_psicc_h
#define _chemistry_qc_psi_psicc_h

#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/mbptr12/spin.h>

namespace sc {
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCC is a Psi coupled cluster wave function

  class PsiCC : public PsiCorrWavefunction {
      RefSCMatrix T1_[NSpinCases1];
      RefSCMatrix T2_[NSpinCases2];
      RefSCMatrix Tau2_[NSpinCases2];
      RefSCMatrix Lambda1_[NSpinCases1];
      RefSCMatrix Lambda2_[NSpinCases2];

    protected:
      // set to true if want to run only if Psi3 and MPQC orbitals match exactly up to a phase
      static const bool use_sparsemap_only_ = false;

      /// set to true to test whether T2 transform from Psi3 to MPQC orbitals works correctly (causes Psi3 to do MP1 instead of CCSD)
      static bool test_t2_phases_;
      static void do_test_t2_phases() {
        test_t2_phases_ = true;
      }
      
      /// read in T1-like quantity using DPD label L
      RefSCMatrix T1(const std::string& L);
      /// read in T2-like quantity using DPD label L
      RefSCMatrix T2(const std::string& L);

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
      /// compare T2 and T2_ref (check that elements < zero are in the same place and elements > soft_zero have the same sign)
      void compare_T2(const RefSCMatrix& T2, const RefSCMatrix& T2_ref,
                      unsigned int no1, unsigned int no2, unsigned int nv1,
                      unsigned int nv2, double zero = 1e-8) const;

    public:
      PsiCC(const Ref<KeyVal>&);
      PsiCC(StateIn&);
      ~PsiCC();
      void save_data_state(StateOut&);

      /// return T amplitudes of rank 1. The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked).
      virtual const RefSCMatrix& T1(SpinCase1 spin1);
      /// return T amplitudes of rank 2. The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked).
      virtual const RefSCMatrix& T2(SpinCase2 spin2);
      /// return Tau2 amplitudes. The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked).
      virtual const RefSCMatrix& Tau2(SpinCase2 spin2);
      /// return Lambda amplitudes of rank 1
      virtual const RefSCMatrix& Lambda1(SpinCase1 spin1);
      /// return Lambda amplitudes of rank 2
      virtual const RefSCMatrix& Lambda2(SpinCase2 spin2);
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCCSD is a concrete implementation of Psi CCSD wave function

  class PsiCCSD : public PsiCC {
    protected:
      void write_input(int conv);
    public:
      PsiCCSD(const Ref<KeyVal>&);
      PsiCCSD(StateIn&);
      ~PsiCCSD();
      void save_data_state(StateOut&);
      int gradient_implemented() const;
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
      int gradient_implemented() const;
  };
  
} // namespace

#endif /* header guard */
