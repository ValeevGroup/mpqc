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

#ifndef _mpqc_src_lib_chemistry_qc_wfn_soad_h
#define _mpqc_src_lib_chemistry_qc_wfn_soad_h

#include <chemistry/qc/wfn/obwfn.h>

namespace sc {

  /** SuperpositionOfAtomicDensities is a OneBodyWavefunction
   *  useful as a guess for other OneBodyWavefunction objects.
   *  It uses WTBS and/or STO-6G orbitals as atomic orbitals.
   *
   *  */
  class SuperpositionOfAtomicDensities : public OneBodyWavefunction {
    public:

    /* these keywords may not even be necessary!
    <tr><td><tt>total_charge</tt><td>integer<td>0<td>the total charge of the molecule
    <tr><td><tt>spin_unrestricted</tt><td>boolean<td>false<td>produce same spin-up and and spin-down orbitals?
     */

      /** A KeyVal constructor is used to generate a SuperpositionOfAtomicDensities
          object from the input. In addition to all keywords of OneBodyWavefunction,
          the following list of keywords
          is accepted:

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

          </table>

          Implementation notes. WTBS basis (and STO-6G for H) is used to specify the atomic orbitals and construct the
          guess density; this guess density is then projected unto the basis given by the
          the keyword <tt>basis</tt>. Since WTBS basis is only available for atoms up to
          Rn, constructor will throw if the Molecule object contains heavier atoms.
       */
      SuperpositionOfAtomicDensities(const Ref<KeyVal>& kv);
      SuperpositionOfAtomicDensities(StateIn&);
      ~SuperpositionOfAtomicDensities();
      void save_data_state(StateOut&);

      // implement pure virtual methods of Wavefunction and below
      /// reports the density projected onto the given basis set
      RefSymmSCMatrix density();
      int spin_polarized();
      int nelectron();
      int value_implemented() const;
      void compute();
      void obsolete();

      // implement pure virtual methods of OneBodyWavefunction
      /// not spin-unrestricted
      int spin_unrestricted();
      RefSCMatrix oso_eigenvectors();
      RefDiagSCMatrix eigenvalues();
      double occupation(int irrep, int vectornum);

    private:
      static ClassDesc class_desc_;

      int total_charge_;
      bool spin_unrestricted_;

      Ref<GaussianBasisSet> minbasis_;   //< STO-6G and/or WTBS basis (or other minimal basis)
      RefSymmSCMatrix minbasis_density_; //< in AO basis

  };


} // end of namespace sc

#endif // end of header guard


// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
