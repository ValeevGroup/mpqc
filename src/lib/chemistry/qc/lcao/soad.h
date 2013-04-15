//
// soad.h
//
// Copyright (C) 2013 Edward Valeev
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _mpqc_src_lib_chemistry_qc_lcao_soad_h
#define _mpqc_src_lib_chemistry_qc_lcao_soad_h

#include <chemistry/qc/wfn/obwfn.h>
#include <chemistry/qc/lcao/wfnworld.h>

namespace sc {

  /** SuperpositionOfAtomicDensities is a OneBodyWavefunction
   *  useful as a guess for other OneBodyWavefunction objects.
   *  It generates the guess by combining approximate spin-averaged atomic densities
   *  to produce an approximate density of the neutral molecule.
   *  That can be followed up by an additional SCF-type density relaxation step.
   *  */
  class SuperpositionOfAtomicDensities : public OneBodyWavefunction {
    public:

      /** A KeyVal constructor is used to generate a SuperpositionOfAtomicDensities
          object from the input. In addition to all keywords of OneBodyWavefunction,
          the following list of keywords
          is accepted:

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

            <tr><td><tt>relax</tt><td>boolean<td>true<td>if true, one iteration of SCF will be performed and the resulting
            density will be idempotent. Otherwise, the density is projected from the minimal basis and
            is (in general) not idempotent.

          </table>

          Since this meant primarily to provide orbitals, not just
          densities, a guess density is used to produce a Fock matrix that is diagonalized;
          the resulting orbitals are populated and this produces the (idempotent) density reported by density() and ao_density().
          To compute a guess density (non-idempotent)
          use guess_minimal_density() and guess_density() .

          Since (nonrelativistic) minimal bases in MPQC are only available for atoms up to
          Rn (Z=86), this constructor will throw if the Molecule object contains atoms with Z>86.
       */
      SuperpositionOfAtomicDensities(const Ref<KeyVal>& kv);
      SuperpositionOfAtomicDensities(StateIn&);
      ~SuperpositionOfAtomicDensities();
      void save_data_state(StateOut&);

      /**
       * This is provided for greater composability. Generates a minimal basis set that can be used to generate a guess density using
       * guess_minimal_density(). STO-6G (for H-Kr) and WTBS (Rb-Rn) bases are used to specify
       * the atomic orbitals.
       * @param mol the Molecule object
       * @return a minimal basis set
       */
      static Ref<GaussianBasisSet> minimal_basis_set(const Ref<Molecule>& mol);

      /**
       *  This is provided for greater composability. For a given minimal basis generates a total
       *  (spin-free) non-idempotent density by superposition of average atomic densities.
       * @param minimal_basis_set a minimal basis set with basis functions of a given l on each atom in aufbau order (\sa guess_minimal_basis() )
       * @param intf Integral factory
       * @return the AO basis density
       */
      static RefSymmSCMatrix guess_minimal_density(const Ref<GaussianBasisSet>& minimal_basis_set,
                                                   const Ref<Integral>& intf);
      /**
       *  This is provided for greater composability. Similar to guess_minimal_density(), but takes any basis.
       *  The density is obtained by projection from a minimal basis.
       * @param basis_set a basis set
       * @param intf Integral factory
       * @return the AO basis density
       */
      static RefSymmSCMatrix guess_density(const Ref<GaussianBasisSet>& basis_set,
                                           const Ref<Integral>& intf);

      // implement pure virtual methods of Wavefunction and below
      /// reports the (idempotent) density corresponding to the orbitals and occupancies reported by oso_eigenvectors() and occupation()
      RefSymmSCMatrix density();
      int spin_polarized();
      int nelectron();
      int value_implemented() const;
      void compute();
      void obsolete();

      // implement pure virtual methods of OneBodyWavefunction
      /// not spin-unrestricted
      int spin_unrestricted();
      /// eigenvectors (and eigenvalues) obtained by diagonalizing the Fock matrix constructed from the minimal basis
      RefSCMatrix oso_eigenvectors();
      RefDiagSCMatrix eigenvalues();
      double occupation(int irrep, int vectornum);

    private:
      static ClassDesc class_desc_;

      bool relax_;
      RefDiagSCMatrix occs_;

      //Ref<GaussianBasisSet> minbasis_;   //< STO-6G and/or WTBS basis (or other minimal basis)
      //RefSymmSCMatrix minbasis_density_; //< in AO basis

  };


} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
