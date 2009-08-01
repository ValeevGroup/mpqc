//
// psici_pt2r12.h
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

#ifndef PSICI_PT2R12_H_
#define PSICI_PT2R12_H_

#include <map>
#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/mbptr12/spin.h>

namespace sc {

  /** PsiCI_PT2R12 is a Psi3 CI computation with a following computation of a perturbative R12 correction.
   *
   */
  class PsiCI_PT2R12 : public PsiCorrWavefunction_PT2R12
  {
    private:
      double eci_;    /// the conventional CI energy

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
      std::vector<int> detcas_detci_average_states_;   /// vector of states over which averaging is performed in a detci of a detcas calculation

      // do orbital optimization first?
      bool rasscf_;
      std::string wfn_type_;  /// wfn keyword is set to this. can be detci or detcas

      // optional RAS info
      std::vector<int> ras1_;
      std::vector<int> ras2_;
      std::vector<int> ras3_;
      int ras3_max_;

      double scf_levelshift_;      /// Psi3 cscf levelshift
      int scf_stop_levelshift_;    /// number of iterations, for which the levelshift is applied

      Ref<OneBodyWavefunction> valence_obwfn_;   /** This OneBodyWavefunction defines valence orbitals.
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
                                   * the "reference_wfn" keyword is specified.
                                   */

      /// orbital reordering
      std::string reorder_;
      std::vector<int> moorder_;

      void write_input(int convergence);

      /** Computes the overlap matrix between the valence_obwfn_ (rowdim) and the
       *  wave function. Both dimensions refer to MOs in symmetry blocked order. */
      RefSCMatrix overlap_with_valence_obwfn();
      /** Returns a map between the MOs of the reference wave function (first index) and
       *  the MOs of the wave function between which the overlap is the largest. */
      std::map<int,int> reference_index_map(const RefSCMatrix &overlap);
      /** Returns a map between the valence_obwfn_->nmo MOs and the lowest MO's of the
       *  corresponding symmetries of the wave function in symmetry blocked MO order. */
      std::map<int,int> standard_index_map(const RefSCMatrix &overlap);
      /** Index map from standard symmetry blocked order (as it is used in Psi3) to
       *  a symmetry blocked order where the lowest orbitals correspond to the valence orbitals
       *  from valence_obwfn_. */
      std::map<int,int> psi_index_map(const RefSCMatrix &overlap);
      /** Returns a vector of orbital indices in symmetry blocked order where the lowest
       *  MOs correspond to the respective MOs in the valence wfn. This vector
       *  can be used as it is in the Psi3 cscf keyword "reorder" to avoid the inclusion
       *  of unphysical diffuse function into the active space of CASSCF calculations.
       *
       *  \param overlap is the overlap between a reference wfn (e.g. minimal basis HF) and a given wfn
       */
      std::vector<int> map_to_valence_order(const RefSCMatrix &overlap);
      /** Returns a vector of the corresponding symmetry (block) indices for each MO
       *  in energy order. */
      std::vector<int> mo_symmetries_in_energetic_order();
      virtual RefSCMatrix MPQC2PSI_transform_matrix(SpinCase1 spin);
      std::vector<unsigned int> mo_blocks();
      std::vector<unsigned int> frzc_blocks();
      std::vector<unsigned int> docc_blocks();
      std::vector<unsigned int> docc_act_blocks();
      std::vector<unsigned int> socc_blocks();
      std::vector<unsigned int> uocc_act_blocks();
      std::vector<unsigned int> uocc_blocks();
      std::vector<unsigned int> frzv_blocks();
      void print_all_blocks(std::ostream &o=ExEnv::out0());

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
      PsiCI_PT2R12(const Ref<KeyVal> &keyval);
      PsiCI_PT2R12(StateIn &s);
      ~PsiCI_PT2R12();
      void save_data_state(StateOut &s);
      void compute();
  };

}

#endif /*PSICI_PT2R12_H_*/
