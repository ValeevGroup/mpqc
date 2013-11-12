//
// eht.h
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
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

#ifndef _chemistry_qc_wfn_eht_h
#define _chemistry_qc_wfn_eht_h

#include <chemistry/qc/wfn/obwfn.h>

namespace sc {

  /// @addtogroup ChemistryElectronicStructureOneBody
  /// @{

/** This computes the extended Huckel energy and wavefunction.  It is useful
   as a quick initial guess for other one body wavefunctions.  */
class ExtendedHuckelWfn: public OneBodyWavefunction {
  private:
    int nirrep_;
    int *docc_;
    int *socc_;
    int total_charge_;
    int user_occ_;

    /// guesses occupations by minimizing total energy of free-electron model and maximizing multiplicity
    void fill_occ(const RefDiagSCMatrix &evals,
                  int nelectron, int *docc, int *socc);

    void compute();

    RefSymmSCMatrix h_eht_oso();

  public:
    ExtendedHuckelWfn(StateIn&);
    /**
     *  The KeyVal constructor accepts all keywords of OneBodyWavefunction class, plus the following additional keywords:
        <dl>

        <dt><tt>total_charge</tt><dd> Specifies the total charge
           of the system. This charge is defined without taking into account custom nuclear
           charges or classical charges, i.e. total charge = sum of atomic numbers of nuclei
           - number of electrons.  The default is 0.

        <dt><tt>socc</tt><dd> This vector of integers gives the total
        number of singly occupied orbitals of each irreducible
        representation. All electrons are assumed to have m_s = +1/2.
        If socc is given, then docc must be given.

        <dt><tt>docc</tt><dd> This vector of integers gives the total
        number of doubly occupied orbitals of each irreducible
        representation.  If docc is given, then socc must be given.

        </dl>
     */
    ExtendedHuckelWfn(const Ref<KeyVal>&);
    ~ExtendedHuckelWfn();

    void save_data_state(StateOut&);

    double occupation(int irrep, int vectornum);

    RefSCMatrix oso_eigenvectors();
    RefDiagSCMatrix eigenvalues();
    RefSymmSCMatrix density();
    double magnetic_moment() const;
    int spin_unrestricted();

    int value_implemented() const;
};

/// @}
// end of addtogroup ChemistryElectronicStructureOneBody

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
