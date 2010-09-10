//
// psici.h
//
// Copyright (C) 2008 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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

#ifndef _chemistry_qc_psi_psici_h
#define _chemistry_qc_psi_psici_h

#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/mbptr12/spin.h>

namespace sc {

  /** PsiRASCI is a general (RAS) CI PsiWavefunction.
   */
  class PsiRASCI : public PsiCorrWavefunction {
    public:
    /** A KeyVal constructor is used to generate a PsiRASCI
        object from the input. It recognizes all keywords of
        PsiCorrWavefunction class and the following keywords:

        <table border="1">

        <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

        <dt><tt>prerequisite</tt><dd> Specifies an object of class PsiWavefunction that must be "computed" prior to computing this object.
        The default is null object, which means clean compute, i.e. start running Psi by invoking "input" module which will overwrite the
        checkpoint file.

        <dt><tt>root</tt><dd> Specifies which state to solve for. The default is 1, i.e. the ground state. Value of <tt>root</tt>
        cannot be greater than <tt>nroots</tt>.

        <dt><tt>nroots</tt><dd> Specifies the number of CI vectors to seek. The default is the value specified with keyword
        root (see keyword <tt>nroots</tt>).

        <dt><tt>multiplicity</tt><dd> Specifies the multiplicity of the state to solve for. The default is the same as the multiplicity of
        the object specified with the <tt>reference</tt> keyword. \sa PsiCorrWavefunction

        <tr><td><tt>valence_obwfn</tt><td>OneBodyWavefunction<td>null<td>This optional keyword specifies
        an object that will provide the orbital ordering for the initial guess. It is recommended to use
        an SCF object with the minimal basis needed to express the orbitals used in defining the RAS spaces.
        For example, for a valence CASSCF this means that SCF with an STO-3G basis will suffice. For states
        with Rydberg character one may want to choose an appropriate ANO basis set.

        </table>
     */
      PsiRASCI(const Ref<KeyVal> &keyval);
      PsiRASCI(StateIn &s);
      ~PsiRASCI();
      void save_data_state(StateOut &s);
      void compute();
      void print(std::ostream&) const;

      /// returns vector that specifies the number of RAS1 orbitals in each irrep
      const std::vector<unsigned int>& ras1() const { return ras1_; }
      /// returns vector that specifies the number of RAS2 orbitals in each irrep
      const std::vector<unsigned int>& ras2() const { return ras2_; }
      /// returns vector that specifies the number of RAS3 orbitals in each irrep
      const std::vector<unsigned int>& ras3() const { return ras3_; }
      /// returns the maximum number of electrons allowed in RAS3 space
      unsigned int ras3_max() const { return ras3_max_; }

      RefSymmSCMatrix mo_density(SpinCase1 spin); // mo_density is overloaded because detci
                                                  // reports density in active orbitals only
      /// if this is PsiRASSCF this will return RASSCF orbitals
      const Ref<OrbitalSpace>& orbs_sb(SpinCase1 spin);
      /// returns occupied OrbitalSpace. If ras3_max=0 this is a subset
      /// of the space reported by orbs_sp(). This space is symmetry-blocked.
      const Ref<OrbitalSpace>& occ(SpinCase1 spin);
      /// 1-pdm in the space reported by occ()
      RefSymmSCMatrix onepdm_occ(SpinCase1 spin);
      /// 2-pdm in the space reported by occ()
      RefSymmSCMatrix twopdm_occ(SpinCase2 spin);

    protected:

      bool opdm_print_;   /// print the one-particle density matrix
      bool tpdm_print_;   /// print the two-particle density matrix
      int root_;          /// compute a specific root of the wave function
      int multiplicity_;  /// the spin multiplicity of the target state
      int nroots_;        /// number of roots for detci calculations
      int target_sym_;    /// the symmetry (irrep) of the target root
      int h0_blocksize_;  /// block size for the H0 guess
      int ex_lvl_;        /// CI excitation level
      bool repl_otf_;     /// do CI string replacements on the fly. saves memory, but is slower.

      // this data may need to be modified by RASSCF
      int energy_convergence_;
      int convergence_;
      int maxiter_;       /// maxiter for detci

      Ref<OrbitalSpace> orbs_sb_[NSpinCases1];
      Ref<OrbitalSpace> occ_[NSpinCases1];

      RefSymmSCMatrix onepdm_occ_[NSpinCases1];
      RefSymmSCMatrix twopdm_occ_[NSpinCases1];

      // optional RAS info
      // it is initialized automatically
      std::vector<unsigned int> ras1_;
      std::vector<unsigned int> ras2_;
      std::vector<unsigned int> ras3_;
      int ras3_max_;

      double scf_levelshift_;      /// Psi3 cscf levelshift
      int scf_stop_levelshift_;    /// number of iterations, for which the levelshift is applied

      /** This OneBodyWavefunction defines valence orbitals.
       *
       * Purpose:
       * If a basis set with diffuse functions is used for a CASSCF calculation,
       * there may be energetically low lying diffuse obitals entering the active
       * space if the orbitals are not reordered. The reference wave function helps
       * to find the appropriate reordering by a 'black box' procedure which can be
       * explained as follows. If the basis set of the reference wave function is
       * chosen to be a small or a minimal basis set, it's virtual orbitals a aready
       * a qualitatively good description of the virtual orbitals of the system,
       * without having the problem of orbitals that have to be reordered. Thus, by
       * computing the overlap between the reference and the original wave function,
       * the 'correct' active space orbitals of the orginal wave function are those
       * which show the largest overlap with the corresponding active space orbitals
       * of the reference wave function. This procedure is used as soon as the
       * the "valence_obwfn" keyword is specified.
       */
      Ref<OneBodyWavefunction> valence_obwfn_;

      /// orbital reordering
      std::vector<unsigned int> moorder_;

      void write_input(int convergence);
      void write_rasci_input(int convergence, bool rasscf);

      std::vector<unsigned int> map_density_to_sb();
  };

  /// PsiRASSCF is a type of a PsiRASCI wavefunction that implements orbital optimization.
  class PsiRASSCF : public PsiRASCI {
    public:
      /** A KeyVal constructor is used to generate a PsiRASSCF
          object from the input. It recognizes all keywords of
          PsiRASCI class and the following keywords:

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          <tr><td><tt>state_average</tt><td>boolean<td>false<td>whether to do state-averaging. The default is to
          compute optimal orbitals for the average of all states (see keyword <tt>num_states</tt>).

          <tr><td><tt>relax_core</tt><td>boolean<td>false<td>whether to keep the occupied orbitals that are not part of
          RAS I or II spaces fixed in RASSCF or to relax them.

          <tr><td><tt>valence_obwfn</tt><td>OneBodyWavefunction<td>null<td>This optional keyword specifies
          an object that will provide the orbital ordering for the initial guess. It is recommended to use
          an SCF object with the minimal basis needed to express the orbitals used in defining the RAS spaces.
          For example, for a valence CASSCF this means that SCF with an STO-3G basis will suffice. For states
          with Rydberg character one may want to choose an appropriate ANO basis set.

          </table>
       */
      PsiRASSCF(const Ref<KeyVal>& kv);
      PsiRASSCF(StateIn&);
      ~PsiRASSCF();
      void save_data_state(StateOut&);
      void compute();
      void print(std::ostream&) const;

    private:
      static ClassDesc class_desc_;

      int rasscf_energy_convergence_;
      int rasscf_convergence_;
      int rasscf_maxiter_;    /// max number of iterations in rasscf
      int diis_start_; /// after X cycles, diis starts

      bool state_average_;   //< state average?
      bool relax_core_;

      bool run_detci_only_;  // hack to allow running state-averaged RASSCF

      void write_input(int convergence);

  };


}

#endif /*_chemistry_qc_psi_psici_h*/
