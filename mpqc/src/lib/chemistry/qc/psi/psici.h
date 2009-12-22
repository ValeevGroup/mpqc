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

  /** PsiCI is a general (RAS) CI PsiWavefunction.
   */
  class PsiCI : public PsiCorrWavefunction {
    public:
      /**
       * KeyVal constructor uses the following keywords
        <dl>

        <dt><tt>valence_obwfn</tt><dd> Specifies the OneBodyWavefunction object used to determine the valence orbitals.
        Recommended to use the minimal-basis HF wavefunction.

        <dt><tt>root</tt><dd> Specifies which CI vector to pick in DETCI. The default is 1, i.e. the lowest root.

        <dt><tt>detci_num_roots</tt><dd> Specifies the number of CI vectors to seek in DETCI. The default is the value specified with keyword
        root. \sa keyword "detcas_detci_num_roots"

       */
      PsiCI(const Ref<KeyVal> &keyval);
      PsiCI(StateIn &s);
      ~PsiCI();
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
      /// if rasscf = true then return MCSCF orbitals
      const Ref<OrbitalSpace>& orbs_sb(SpinCase1 spin);
      /// returns occupied OrbitalSpace. For CAS methods this is a subset
      /// of the space reported by orbs_sp(). This space is symmetry-blocked.
      const Ref<OrbitalSpace>& occ(SpinCase1 spin);
      /// 1-pdm in the space reported by occ()
      RefSymmSCMatrix onepdm_occ(SpinCase1 spin);
      /// 2-pdm in the space reported by occ()
      RefSymmSCMatrix twopdm_occ(SpinCase2 spin);

    private:
      bool opdm_print_;   /// print the one-particle density matrix
      bool tpdm_print_;   /// print the two-particle density matrix
      int root_;          /// compute a specific root of the wave function
      int detci_num_roots_;           /// number of roots for detci calculations
      int detcas_detci_num_roots_;    /// number of roots for detci in detcas calculations
      int h0_blocksize_;  /// block size for the H0 guess
      int ex_lvl_;        /// CI excitation level
      bool repl_otf_;     /// do CI string replacements on the fly. saves memory, but is slower.
      int detci_energy_convergence_;
      int detcas_energy_convergence_;
      int detcas_detci_energy_convergence_; /// energy convergence of the detci energy in each step of detcas calculation
      int detci_convergence_;
      int detcas_convergence_;
      int detcas_detci_convergence_;  /// convergence of the detci wave function in each step of detcas calculation
      bool detcas_diis_;   /// use DIIS for detcas calculations
      int detcas_detci_maxiter_;    /// maxiter for detci in detcas calculations
      int detci_maxiter_;     /// maxiter for detci (not in detcas calculations)
      std::vector<unsigned int> detcas_detci_average_states_;   /// vector of states over which averaging is performed in a detci of a detcas calculation

      // do orbital optimization first?
      bool rasscf_;
      std::string wfn_type_;  /// wfn keyword is set to this. can be detci or detcas

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
      std::string reorder_;
      std::vector<unsigned int> moorder_;

      void write_input(int convergence);

      std::vector<unsigned int> map_density_to_sb();
  };

}

#endif /*_chemistry_qc_psi_psici_h*/
