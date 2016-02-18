//
// psicc_pt2r12.h
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

#ifndef _chemistry_qc_psi_psiccpt2r12_h
#define _chemistry_qc_psi_psiccpt2r12_h

#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/psi/psicc.h>

namespace sc {

  ///////////////////////////////////////////////////////////////////
  /// PsiCC_PT2R12 is used to implement \f$\mathrm{CC}-(2)_{\overline{R12}}\f$ methods
  class PsiCC_PT2R12 : public PsiCC {
      Ref<R12IntEval> r12eval_;           // the R12 intermediates evaluator
      Ref<R12WavefunctionWorld> r12world_;   // parameters for r12eval_
      Ref<MP2R12Energy> mp2r12_energy_;
      bool spinadapted_;
      bool cabs_singles_;
      double cabs_singles_energy_;
      double pccsd_alpha_;
      double pccsd_beta_;
      double pccsd_gamma_;
    protected:
      /// set to true to use Ts instead of Lambdas
      static const bool replace_Lambda_with_T_ = true;
      /// EXPERTS-ONLY: if you want to enable TA-based evaluation of higher-order terms
      /// turn this on
      static const bool need_lambda_ = false;
      /** default was to include up to 3rd-order terms in the energy (V.T2)
          current default is to include higher-order terms in the energy also -- this means that
          V.T1 terms are included even in closed-shell calculations,
          in contrast to the original formulation ... the max order for the intermediates
          is thus 4 (V is second order, T1 is 1st(ROHF)/2nd(RHF,UHF) */
      static const unsigned int completeness_order_ = 4;

      void write_basic_input(int conv);

      void compute_ept2r12();

    public:
      /** The KeyVal constructor uses keywords of PsiCC, WavefunctionWorld, and R12WavefunctionWorld, and the following keywords
          <dl>

    <dt><tt>spinadapted</tt><dd> This boolean specifies whether to compute spin-adapted
    or spin-orbital pair energies. Default is to compute spin-adapted energies for closed-shell
    systems and spin-orbital energies for open-shell systems. For some references, e.g. UHF, this keyword
    is not used.

      <dt><tt>cabs_singles</tt><dd> Evaluate the second-order energy contribution from
      CABS singles and include it into the CC-R12 energy. The default is true.

      <dt><tt>pccsd_alpha</tt><dd> The default is 1.0 .

      <dt><tt>pccsd_beta</tt><dd> The default is 1.0 .

      <dt><tt>pccsd_gamma</tt><dd> The default is 1.0 .

      </dl> */
      PsiCC_PT2R12(const Ref<KeyVal>&);
      PsiCC_PT2R12(StateIn&);
      ~PsiCC_PT2R12();
      void save_data_state(StateOut&);

      const Ref<R12WavefunctionWorld>& r12world() const { return r12world_; }
      const Ref<R12IntEval>& r12eval() const { return r12eval_; }
      // CABS singles contribution to the total energy
      double cabs_singles_energy();

      // CCSD_F12 orbital relaxation contribution to 1rdm
      void compute_onerdm_relax(const Ref<R12EnergyIntermediates>& r12intermediates,
                                RefSCMatrix& Dorbs_alpha,
                                RefSCMatrix& Dorbs_beta);

      /// print
      void print(std::ostream&o=ExEnv::out0()) const;

      void obsolete();
  };

  ///////////////////////////////////////////////////////////////////
  /// PsiCCSD_PT2R12 is a concrete implementation of the \f$\mathrm{CCSD}-(2)_{\overline{R12}}\f$ method
  class PsiCCSD_PT2R12 : public PsiCC_PT2R12 {
      double eccsd_;
      void write_input(int conv);
    public:
      /** The KeyVal constructor uses keywords of PsiCC_PT2R12.
      */
      PsiCCSD_PT2R12(const Ref<KeyVal>&);
      PsiCCSD_PT2R12(StateIn&);
      ~PsiCCSD_PT2R12();
      void save_data_state(StateOut&);
      void compute();
      /// print
      void print(std::ostream&o=ExEnv::out0()) const;
  };

  ///////////////////////////////////////////////////////////////////
  /// PsiCCSD_PT2R12T is a concrete implementation of the \f$\mathrm{CCSD}(T)_{\overline{R12}}\f$ method
  class PsiCCSD_PT2R12T : public PsiCC_PT2R12 {
      double eccsd_;
      double e_t_;
      void write_input(int conv);
    public:
      /** The KeyVal constructor uses keywords of PsiCC_PT2R12.
      */
      PsiCCSD_PT2R12T(const Ref<KeyVal>&);
      PsiCCSD_PT2R12T(StateIn&);
      ~PsiCCSD_PT2R12T();
      void save_data_state(StateOut&);
      void compute();
      /// print
      void print(std::ostream&o=ExEnv::out0()) const;
  };

  ///////////////////////////////////////////////////////////////////
  /// PsiCC3_PT2R12 is a concrete implementation of the ground-state \f$\mathrm{CC3}-(2)_{\overline{R12}}\f$ method
  class PsiCC3_PT2R12 : public PsiCC_PT2R12 {
      double ecc3_;
      void write_input(int conv);
    public:
      /** The KeyVal constructor uses keywords of PsiCC_PT2R12.
      */
      PsiCC3_PT2R12(const Ref<KeyVal>&);
      PsiCC3_PT2R12(StateIn&);
      ~PsiCC3_PT2R12();
      void save_data_state(StateOut&);
      void compute();
      /// print
      void print(std::ostream&o=ExEnv::out0()) const;
  };

}

#endif /* header guard */
