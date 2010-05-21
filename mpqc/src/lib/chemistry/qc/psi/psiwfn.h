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
#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/psi/psiexenv.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/orbitalspace.h>
#include <chemistry/qc/mbptr12/mbptr12.h>

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
      size_t memory_;
      char *memory_str_;
      /// Prepares a complete Psi input file. The input file is assumed to have been opened.
      virtual void write_input(int conv) =0;

      std::vector<unsigned int> read_occ(const Ref<KeyVal> &keyval, const char *name,
                                         size_t nirrep);

      /// return the debug level
      int debug() const;

    public:
      /** The KeyVal constructor.

       <dl>

       <dt><tt>psienv</tt><dd> Specifies a PsiExEnv object.  There
       is no default.

       <dt><tt>memory</tt><dd> This integer specifies the amount of memory
       (in bytes) for Psi to use. The default is 32000000.

       <dt><tt>debug</tt><dd> This integer can be used to produce output
       for debugging.  The default is 0.

       </dl> */
      PsiWavefunction(const Ref<KeyVal>&);
      PsiWavefunction(StateIn&);
      ~PsiWavefunction();

      void save_data_state(StateOut&);

      /// returns the Psi3 convention for the ordering of the cartesian functions
      static Integral::CartesianOrdering cartesian_ordering();

      /** Writes out Psi input file entries specific to this PsiWavefunction.
       The input file is assumed to have been opened. */
      virtual void write_basic_input(int conv);
      void compute();
      void print(std::ostream&o=ExEnv::out0()) const;
      int nirrep() const { return nirrep_; }

      RefSymmSCMatrix density();
      /// Returns the MO basis density (blocked by symmetry)
      virtual RefSymmSCMatrix mo_density(SpinCase1 spin = AnySpinCase1) =0;

      /// Return an associated PsiExEnv object
      Ref<PsiExEnv> exenv() const {
        return exenv_;
      }

      /// Return an associated PsiInput object
      Ref<PsiInput> get_psi_input() const {
        return exenv_->get_psi_input();
      }

#if 0
      /// Returns a map from shells in Psi3 basis to std::pair<shell,contraction> in MPQC basis (note that Psi3 does not handle general contractions)
      std::vector< std::pair<unsigned int,unsigned int> > shell_map();
      /// Returns a map from AO in Psi3 basis to AO in MPQC basis
      std::vector<unsigned int> ao_map();
#endif

      /// return Psi3 nuclear repulsion energy
      double nuclear_repulsion_energy() const;
  };

  class PsiSCF;

  ///////////////////////////////////////////////////////////////////
  /// PsiCorrWavefunction is a Psi correlated wave function

  class PsiCorrWavefunction : public PsiWavefunction {
    protected:
      Ref<PsiSCF> reference_;
      unsigned int nfzc_;
      unsigned int nfzv_;
      mutable std::vector<unsigned int> frozen_docc_;
      mutable std::vector<unsigned int> frozen_uocc_;
      RefSymmSCMatrix mo_density_[NSpinCases1];

      void write_input(int conv);
      // set fz2restr to use RESTRICTED instead of FROZEN keywords
      void write_input_frozen2restricted(int conv, bool fz2restr);

      double valacc_to_refacc() const { return 100.0; }

      /// returns the index map that transforms indices in which densities are reported in Psi
      /// to the symmetry-blocked indices. Single-reference correlated methods in Psi
      /// report densities using QT-ordered indices, whereas multireference methods
      /// use RAS ordering. The default implementation of this function assumes the former.
      /// Overload as necessary.
      virtual std::vector<unsigned int> map_density_to_sb();

    public:
      /** A KeyVal constructor is used to generate a PsiCorrWavefunction
          object from the input. The full list of keywords
          that are accepted is below.

          <table border="1">

          <tr><td>%Keyword<td>Type<td>Default<td>Description

          <tr><td><tt>reference</tt><td>PsiSCF<td>none<td>specifies the reference wavefunction.

          <tr><td><tt>frozen_docc</tt><td>vector<int><td>none<td>specifies the number of doubly-occupied
          orbitals in each irreducible representation that will be kept "frozen" in correlated calculations;
          the number of elements in frozen_docc must equal the order of the computational point group.

          <tr><td><tt>frozen_uocc</tt><td>vector<int><td>none<td>specifies the number of non-occupied
          orbitals in each irreducible representation that will be kept "frozen" in correlated calculations;
          the number of elements in frozen_uocc must equal the order of the computational point group.

          <tr><td><tt>nfzc</tt><td>int<td>0<td>specifies the total number of doubly-occupied
          orbitals that will be kept "frozen" in correlated calculations;
          this keyword is ignored if <tt>frozen_docc</tt> is given. The orbitals to be kept frozen will
          be those of the lowest energy in reference object.

          <tr><td><tt>nfzv</tt><td>int<td>0<td>specifies the total number of non-occupied
          orbitals that will be kept "frozen" in correlated calculations;
          this keyword is ignored if <tt>frozen_uocc</tt> is given. The orbitals to be kept frozen will
          be those of the highest energy in reference object.

          </table>
       */
      PsiCorrWavefunction(const Ref<KeyVal>&);
      PsiCorrWavefunction(StateIn&);
      ~PsiCorrWavefunction();
      void save_data_state(StateOut&);

      void print(std::ostream& os) const;
      int spin_polarized();
      /// sets the desired value accuracy
      void set_desired_value_accuracy(double acc);

      void compute(); // compute is overloaded because reference object needs to be marked computed

      const Ref<PsiSCF>& reference() const;
      /// Number of electrons
      int nelectron();

#if 0   // these should be implemented in PsiCC classes
      /// symmetry-blocked space of active occupied orbitals from Psi3
      const Ref<OrbitalSpace>& occ_act_sb(SpinCase1);
      /// symmetry-blocked space of active virtual orbitals from Psi3
      const Ref<OrbitalSpace>& vir_act_sb(SpinCase1);
#endif
      /// total # of frozen doubly-occupied orbitals
      unsigned int nfzc() const;
      /// total # of frozen unoccupied orbitals
      unsigned int nfzv() const;
      /// symmetry-blocked space of MO's from Psi3
      /// the default implementation returns the orbitals from reference()
      /// can be overridden if this wfn changes reference orbitals
      virtual const Ref<OrbitalSpace>&  orbs_sb(SpinCase1 spin);

      /// # of frozen doubly-occupied orbitals per irrep
      const std::vector<unsigned int>& frozen_docc() const;
      /// # of frozen unoccupied orbitals per irrep
      const std::vector<unsigned int>& frozen_uocc() const;
      /// # of occupied active orbitals per irrep
      const std::vector<unsigned int> docc_act();
      const std::vector<unsigned int> socc();
      const std::vector<unsigned int> uocc_act();
      const std::vector<unsigned int> docc();
      const std::vector<unsigned int> uocc();

      /// reference energy
      virtual double reference_energy();

      /// return one-particel density matrix in symmetry-blocked orbitals \sa orbs_sb()
      RefSymmSCMatrix mo_density(SpinCase1 spin);
#if 0
      /// return one-particle density matrix as a symmetric matrix indexed by (moindex1,moindex2).
      RefSymmSCMatrix onepdm(const SpinCase1 &spin);
      RefSymmSCMatrix onepdm();
      /// this twopdm is stored in chemist's (Mulliken) notation order.
      /// Access element (ij|km) by get_element(ordinary_INDEX(i,j), ordinary_INDEX(k,m))
      RefSymmSCMatrix twopdm();
      /// this twopdm is stored in Dirac notation order.
      RefSymmSCMatrix twopdm_dirac();
      RefSymmSCMatrix twopdm_dirac_from_components();
#endif
      RefSymmSCMatrix twopdm_dirac(const SpinCase2 &pairspin);
      void print_onepdm_vec(FILE *output,const RefSCVector &opdm,double TOL);
      void print_onepdm_mat(FILE *output,const RefSymmSCMatrix &opdm,double TOL);
      void print_twopdm_mat(FILE *output,const RefSymmSCMatrix &tpdm, double TOL);
      void print_twopdm_arr(FILE *output,double *tpdm,double TOL);
  };


  ///////////////////////////////////////////////////////////////////
  /// PsiSCF is an abstract base for all Psi SCF wave functions

  class PsiSCF : public PsiWavefunction {
      Ref<OrbitalSpace> orbs_sb_[NSpinCases1];
      RefDiagSCMatrix evals_[NSpinCases1];
      RefSCMatrix coefs_[NSpinCases1];
      RefSymmSCMatrix mo_density_[NSpinCases1];
      std::vector<unsigned int> occpi_[NSpinCases1];
      std::vector<unsigned int> uoccpi_[NSpinCases1];
      std::vector<unsigned int> mopi_;
      std::vector<double> occupation_[NSpinCases1];
      void compute_occupations(SpinCase1 spin);

      // PsiCorrWavefunction needs to be able to set value of PsiSCF
      friend void PsiCorrWavefunction::compute();
    protected:
      std::vector<unsigned int> docc_;
      std::vector<unsigned int> socc_;
      int multp_;
      int charge_;
      static const int maxiter = 200;
      /// how much lower is the desired accuracy of the guess?
      static double guess_acc_ratio() { return 1e4; }

      /// guess wave function is only used to get the occupations
      Ref<OneBodyWavefunction> guess_wfn_;

    public:
      PsiSCF(const Ref<KeyVal>&);
      PsiSCF(StateIn&);
      ~PsiSCF();
      void save_data_state(StateOut&);

      /// imports occupations from obwfn. Will throw if docc_ and socc_ had been initialized
      /// and do not match obwfn.
      void import_occupations(const Ref<OneBodyWavefunction>& obwfn);
      enum RefType {rhf, hsoshf, uhf};
      /// Returns the PsiSCF::RefType of this particular Psi SCF wave function
      virtual PsiSCF::RefType reftype() const =0;
      /// Returns the eigenvalues matrix
      virtual const RefDiagSCMatrix& evals(SpinCase1 spin = AnySpinCase1);
      /// Returns the coefficient matrix in AO basis
      virtual const RefSCMatrix& coefs(SpinCase1 spin = AnySpinCase1);
      /// Number of occupied orbitals of spin S per irrep
      const std::vector<unsigned int>& occpi(SpinCase1 S);
      /// Number of unoccupied orbitals of spin S per irrep
      const std::vector<unsigned int>& uoccpi(SpinCase1 S);
      /// Number of orbitals per irrep
      const std::vector<unsigned int>& mopi();
      /// symmetry-blocked space
      const Ref<OrbitalSpace>&  orbs_sb(SpinCase1 spin);
      /// Number of electrons
      int nelectron();
      /// Returns the total occupation for orbital mo
      double occupation(int mo);
      /// Returns the occupation for alpha orbitals
      double alpha_occupation(int mo);
      /// Returns the occupation for beta orbitals
      double beta_occupation(int mo);
      RefSymmSCMatrix density();
      RefSymmSCMatrix alpha_density();
      RefSymmSCMatrix beta_density();

      RefSymmSCMatrix mo_density(SpinCase1 spin = AnySpinCase1);

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
      void print(std::ostream& os = ExEnv::out0()) const;

      void write_basic_input(int conv);
      int spin_polarized() {
        return 0;
      }
      int gradient_implemented() const {
        return 1;
      }
      PsiSCF::RefType reftype() const {
        return rhf;
      }
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
      void print(std::ostream& os = ExEnv::out0()) const;

      void write_basic_input(int conv);
      int spin_polarized() {
        return 1;
      }
      int gradient_implemented() const {
        return 1;
      }
      PsiSCF::RefType reftype() const {
        return hsoshf;
      }
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
      void print(std::ostream& os = ExEnv::out0()) const;

      void write_basic_input(int conv);
      int spin_polarized() {
        return 1;
      }
      int gradient_implemented() const {
        return 1;
      }
      PsiSCF::RefType reftype() const {
        return uhf;
      }
  };

#if 0
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /// PsiCorrWavefunction_PT2R12: a corrlated wave function with a perturbational explicitly correlated correction.

  class PsiCorrWavefunction_PT2R12 : public PsiCorrWavefunction {
    protected:
      size_t memory_r12_;
      Ref<SCF> reference_mpqc_;
      Ref<R12IntEval> r12eval_;
      Ref<R12WavefunctionWorld> r12world_;
      enum onepdm_type { HF = 0, correlated = 1 };
      onepdm_type opdm_type_;
      bool tpdm_from_opdms_;
      bool exactgamma_phi_twoelec_;  /// use the exact gamma for two-electron systems.
      int debug_;
    protected:
      /// MPQC to Psi transform matrix if Psi uses QT ordering.
      virtual RefSCMatrix MPQC2PSI_transform_matrix(SpinCase1 spin) = 0;
      /// MPQC to Psi transform matrix if Psi uses ras ordering.
      //RefSCMatrix MPQC2PSI_transform_matrix(SpinCase1 spin,
      //                                      const std::vector<unsigned int> &ras1,
      //                                      const std::vector<unsigned int> &ras2,
      //                                      const std::vector<unsigned int> &ras3);
      void write_input(int convergence);
      /// Returns Hcore in MO basis
      RefSymmSCMatrix hcore_mo();
      RefSymmSCMatrix hcore_mo(SpinCase1 spin);
      /// molecular integrals in chemist's notation
      RefSymmSCMatrix moints();  // closed shell case
      RefSCMatrix moints(SpinCase2 pairspin);
      /// This is the g (in Dirac notation) that should be used.
      RefSCMatrix g(SpinCase2 pairspin);
      /**
       * general-reference Fock atrix (cf. eqn. (71b) of W. Kutzelnigg, D. Mukherjee, J. Chem Phys. 107 (1997) p. 432)
       * if opdm_ in R12WavefunctionWorld is set, otherwise it returns the ordinary Fock operator.
       */
      RefSCMatrix f(SpinCase1 spin);
      /// phi of a single-reference Hartree-Fock reference (needed only for debugging)
      RefSCMatrix phi_HF(SpinCase2 pairspin);
      /**
       * phi not truncated in lambda, but without the three-pdm. This should be always used for
       * two-electron systems.
       */
      RefSCMatrix phi_twoelec(SpinCase2 pairspin);
      /*
       * phi truncated in lambda: terms with three-particle lambda's or higher or terms with
       * squares (or higher) of two-particle lambda's are neglected.
       */
      RefSCMatrix phi(SpinCase2 pairspin);
      /*
       * phi truncated in lambda: terms with three-particle lambda's or higher or terms with
       * squares (or higher) of two-particle lambda's are neglected.
       *
       * \sa PsiCorrWavefunction::phi()
       *  this version does not use 2-pdm, but a 2-body cumulant (2-lambda)
       */
      RefSymmSCMatrix phi_cumulant(SpinCase2 pairspin);

      /// Returns the Psi3 computed correlated onepdm's transformed to MPQC MO's. Function called from function onepdm_refmo.
      RefSymmSCMatrix onepdm_transformed(const SpinCase1 &spin);
      /// Returns the closed shell total onepdm_transformed, i.e. Alpha plus Beta, transformed to MPQC MO's.
      RefSymmSCMatrix onepdm_transformed();
      /// Returns the Hartree-Fock onepdm (Just 1's on the occ-occ block diagonal).
      RefSymmSCMatrix onepdm_refmo_HF(const SpinCase1 &spin);
      /// Returns the closed shell total twopdm_transformed, i.e. Alpha plus Beta, in chemist's (Mulliken) notation, transformed to MPQC MO's.
      RefSymmSCMatrix twopdm_transformed();
      /// Returns the closed shell total twopdm_transformed, i.e. Alpha plus Beta, in Dirac notation, transformed to MPQC MO's.
      RefSymmSCMatrix twopdm_transformed_dirac();
      /// Returns the Psi3 computed correlated twopdm's tranformed to MPQC MO's. Function called from function twopdm_refmo.
      RefSymmSCMatrix twopdm_transformed_dirac(const SpinCase2 &pairspin);
      /// Returns the Psi3 computed correlated onepdm's tranformed to MPQC MO's
      RefSymmSCMatrix twopdm_refmo_from_onepdm(const SpinCase2 &pairspin);
    public:
      PsiCorrWavefunction_PT2R12(const Ref<KeyVal> &keyval);
      PsiCorrWavefunction_PT2R12(StateIn &s);
      ~PsiCorrWavefunction_PT2R12();
      void save_data_state(StateOut &s);
      /// This function should be used to get the onepdm's from Psi, transformed to reference MO's.
      RefSymmSCMatrix onepdm_refmo(const SpinCase1 &spin);
      /// This function should be used to get the onepdm's from Psi, transformed to reference MO's.
      RefSymmSCMatrix twopdm_refmo(const SpinCase2 &pairspin);
      /// This function should be used to get the two-particle lambda's from Psi, transformed to reference MO's.
      RefSymmSCMatrix lambda_refmo(const SpinCase2 &pairspin);
      /// Returns the geminal coefficients.
      RefSCMatrix C(SpinCase2 S);
      RefSCMatrix V_genref_projector2(SpinCase2 pairspin);
      RefSCMatrix V_transformed_by_C(SpinCase2 pairspin);
      RefSymmSCMatrix X_transformed_by_C(SpinCase2 pairspin);
      RefSymmSCMatrix B_transformed_by_C(SpinCase2 pairspin);
      /// computes the projected contribution to the energy.
      double compute_DC_energy_GenRefansatz2();
      /// Returns the dipole moments for closed shell wave functions from onepdm's (mainly for testing).
      RefSCVector dipolemoments(const Ref<DipoleData> &dipoledata = new DipoleData());
      double energy_HF();
      double energy_conventional();
      double energy_conventional_so();
      /// This function computes the "old" General_PT2R12 correction, i.e. the one invoking projector 1.
      double energy_PT2R12_projector1(SpinCase2 pairspin);
      double energy_PT2R12_projector2(SpinCase2 pairspin);
      /// Prints the pair energies as a trace of the given matrix Contrib_mat
      void print_pair_energies(const RefSCMatrix &Contrib_mat,SpinCase2 pairspin);

      const Ref<R12WavefunctionWorld>& r12world() const { return r12world_; }
  };
#endif

  namespace detail {
    /// returns 1-RDM in symmetry-blocked form.
    /// Depending on the type of a calculation, on-disk Psi densities will be reported in different formats.
    /// Hence the user must provide a map from the on-disk order (e.g. RAS)
    /// to the target symmetry-blocked order (can be a subset of the full symmetry-blocked order)
    /// @param spin the spincase
    /// @param mopi the number of orbitals in each block of the target dimension. The output matrix will have blocked dimensions.
    /// @param dmap maps on-disk indices (calculation-specific) to the target indices
    /// @param kit optional SCMatrixKit that will be used to create the target 1-RDM matrix. If omitted, the default BlockedSCMatrixKit will be used
    RefSymmSCMatrix rdopdm(SpinCase1 spin,
                           const std::vector<unsigned int>& mopi,
                           const std::vector<unsigned int>& dmap,
                           Ref<SCMatrixKit> kit = 0);
    /// this function reads 2-rdm in physicists format. Depending on the
    /// calculation, on-disk Psi densities will be reported in different formats.
    /// Hence the user must provide the index map (dmap). @sa sc::detail::rdopdm
    RefSymmSCMatrix rdtpdm(SpinCase2 pairspin,
                           const std::vector<unsigned int>& dmap,
                           Ref<SCMatrixKit> kit = 0);
  }
}
#endif
