//
// psiwfn.h
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

#ifndef _chemistry_qc_psi_psiwfn_h
#define _chemistry_qc_psi_psiwfn_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/psi/psiexenv.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/moindexspace.h>

namespace sc {
  
  ///////////////////////////////////////////////////////////////////
  /** PsiWavefunction is an abstract base for all Psi wave functions.
   Its KeyVal constructor is invoked by all KeyVal constructors of
   concrete implementations of PsiWavefunction.
   */

  class PsiWavefunction : public Wavefunction {
      
      Ref<PsiExEnv> exenv_;
      /// All Psi wave functions can at least compute the energy
      int value_implemented() const {
        return 1;
      }
      
    protected:
      int nirrep_;
      char *memory_;
      /// Prepares a complete Psi input file. The input file is assumed to have been opened.
      virtual void write_input(int conv) =0;

      std::vector<int> read_occ(const Ref<KeyVal> &keyval, const char *name,
                                int nirrep);

      /// return the debug level
      int debug() const;

    public:
      /** The KeyVal constructor.

       <dl>

       <dt><tt>psienv</tt><dd> Specifies a PsiExEnv object.  There
       is no default.

       <dt><tt>memory</tt><dd> This integer specifies the amount of memory
       (in bytes) for Psi to use. The default is 2000000.

       <dt><tt>debug</tt><dd> This integer can be used to produce output
       for debugging.  The default is 0.

       </dl> */
      PsiWavefunction(const Ref<KeyVal>&);
      PsiWavefunction(StateIn&);
      ~PsiWavefunction();

      void save_data_state(StateOut&);

      /** Writes out Psi input file entries specific to this PsiWavefunction.
       The input file is assumed to have been opened. */
      virtual void write_basic_input(int conv);
      void compute();
      void print(std::ostream&o=ExEnv::out0()) const;
      RefSymmSCMatrix density();
      int nirrep() const { return nirrep_; }

      /// Return an associated PsiExEnv object
      Ref<PsiExEnv> exenv() const {
        return exenv_;
      }
      ;
      /// Return an associated PsiInput object
      Ref<PsiInput> get_psi_input() const {
        return exenv_->get_psi_input();
      }
      ;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiSCF is an abstract base for all Psi SCF wave functions

  class PsiSCF : public PsiWavefunction {
      RefDiagSCMatrix evals_[NSpinCases1];
      RefSCMatrix coefs_[NSpinCases1];
      std::vector<unsigned int> occpi_[NSpinCases1];
      std::vector<unsigned int> uoccpi_[NSpinCases1];
      std::vector<unsigned int> mopi_;
      
    protected:
      std::vector<int> docc_;
      std::vector<int> socc_;
      int multp_;
      int charge_;

    public:
      PsiSCF(const Ref<KeyVal>&);
      PsiSCF(StateIn&);
      ~PsiSCF();
      void save_data_state(StateOut&);

      // use spin datatypes defined in spin.h
      typedef sc::SpinCase1 SpinCase1;
      enum RefType {rhf, hsoshf, uhf};
      /// Returns the PsiSCF::RefType of this particular Psi SCF wave function
      virtual PsiSCF::RefType reftype() const =0;
      /// Returns the eigenvalues matrix
      virtual const RefDiagSCMatrix& evals(SpinCase1 spin = Alpha);
      /// Returns the coefficient matrix
      virtual const RefSCMatrix& coefs(SpinCase1 spin = Alpha);
      /// Number of occupied orbitals of spin S per irrep
      const std::vector<unsigned int>& occpi(SpinCase1 S);
      /// Number of unoccupied orbitals of spin S per irrep
      const std::vector<unsigned int>& uoccpi(SpinCase1 S);
      /// Number of orbitals per irrep
      const std::vector<unsigned int>& mopi();
      /// Number of electrons
      int nelectron();

      /// number of MOs
      unsigned int nmo();
      /// number of occupied MOs of spin
      unsigned int nocc(SpinCase1 spin);
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCLHF is a concrete implementation of Psi RHF wave function

  class PsiCLHF : public PsiSCF {
    protected:
      void write_input(int conv);
    public:
      PsiCLHF(const Ref<KeyVal>&);
      PsiCLHF(StateIn&);
      ~PsiCLHF();

      void write_basic_input(int conv);
      int spin_polarized() {
        return 0;
      }
      ;
      int gradient_implemented() const {
        return 1;
      }
      ;
      PsiSCF::RefType reftype() const {
        return rhf;
      }
      ;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiHSOSHF is a concrete implementation of Psi ROHF wave function

  class PsiHSOSHF : public PsiSCF {
    protected:
      void write_input(int conv);
    public:
      PsiHSOSHF(const Ref<KeyVal>&);
      PsiHSOSHF(StateIn&);
      ~PsiHSOSHF();

      void write_basic_input(int conv);
      int spin_polarized() {
        return 0;
      }
      ;
      int gradient_implemented() const {
        return 1;
      }
      ;
      PsiSCF::RefType reftype() const {
        return hsoshf;
      }
      ;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiUHF is a concrete implementation of Psi UHF wave function

  class PsiUHF : public PsiSCF {
    protected:
      void write_input(int conv);
    public:
      PsiUHF(const Ref<KeyVal>&);
      PsiUHF(StateIn&);
      ~PsiUHF();

      void write_basic_input(int conv);
      int spin_polarized() {
        return 1;
      }
      ;
      int gradient_implemented() const {
        return 1;
      }
      ;
      PsiSCF::RefType reftype() const {
        return uhf;
      }
      ;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCorrWavefunction is a Psi correlated wave function

  class PsiCorrWavefunction : public PsiWavefunction {
    protected:
      Ref<PsiSCF> reference_;
      Ref<MOIndexSpace> occ_act_sb_[NSpinCases1];
      Ref<MOIndexSpace> vir_act_sb_[NSpinCases1];
      unsigned int nfzc_;
      unsigned int nfzv_;
      std::vector<unsigned int> frozen_docc_;
      std::vector<unsigned int> frozen_uocc_;
      void write_input(int conv);

    public:
      PsiCorrWavefunction(const Ref<KeyVal>&);
      PsiCorrWavefunction(StateIn&);
      ~PsiCorrWavefunction();
      void save_data_state(StateOut&);
      int spin_polarized() {
        return reference_->spin_polarized();
      }
      const Ref<PsiSCF>& reference() const { return reference_; }
      /// Number of electrons
      int nelectron();
      /// symmetry-blocked space of active occupied orbitals from Psi3
      const Ref<MOIndexSpace>& occ_act_sb(SpinCase1);
      /// symmetry-blocked space of active virtual orbitals from Psi3
      const Ref<MOIndexSpace>& vir_act_sb(SpinCase1);
      /// # of frozen doubly-occupied orbitals per irrep
      const std::vector<unsigned int>& frozen_docc();
      /// # of frozen unoccupied orbitals per irrep
      const std::vector<unsigned int>& frozen_uocc();

      /// reference energy
      virtual double reference_energy();
      
  };
  
}
#endif
