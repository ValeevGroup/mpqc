//
// moindexspace.h
//
// Copyright (C) 2004 Edward Valeev
//
// Author: Edward Valeev <edward.valeev@chemistry.gatech.edu>
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

#ifndef _chemistry_qc_mbptr12_moindexspace_h
#define _chemistry_qc_mbptr12_moindexspace_h

#include <vector>
#include <util/ref/ref.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <math/scmat/abstract.h>
#include <util/state/statein.h>
#include <chemistry/qc/basis/basis.h>

using namespace std;

namespace sc {

  /** Class MOIndexSpace describes a range of molecular orbitals or
      similar objects that are linear combinations of basis functions
      (e.g. atomic orbitals). Examples of sets that can be described
      using MOIndexSpace are occupied MOs and virtual MOs. */

class MOIndexSpace : virtual public SavableState {

public:

  /// Describes the ordering of indices
  enum IndexOrder { symmetry = 0, energy = 1, other = 2 };

private:

  Ref<GaussianBasisSet> basis_;    // The AO basis
  RefSCMatrix coefs_;              // AO-> MO transformation coefficients (nao by nmo matrix)
  vector<int> mosym_;              // irrep of each orbital

  int ng_;            // Order of the point group
  int rank_;          // The rank of this space
  int full_rank_;     // Rank of the full space, i.e. number of MOs
  int nblocks_;       // Number of blocks
  vector<int> first_mo_;   // Index of the first MO in each block
  vector<int> nmo_;        // Number of MOs in each block
  vector<int> offsets_;    // Full-space index of the first MO in each block

  IndexOrder moorder_;

  // checks mosym_ for irrep indices outside of the allowed range
  void check_mosym() const;

  // initialize the object
  void init();

public:
  MOIndexSpace(StateIn&);
  /** Constructs an MOIndexSpace object using a symmetry-blocked transformation coefficient matrix,
      basis set, and block offsets */
  MOIndexSpace(const RefSCMatrix& coefs, const Ref<GaussianBasisSet> basis, const vector<int>& offsets);
  /** Constructs an MOIndexSpace object using a non-blocked transformation coefficient matrix,
      basis set, block offsets, and MO irreps */
  MOIndexSpace(const RefSCMatrix& coefs, const Ref<GaussianBasisSet> basis, const vector<int>& offsets,
               const vector<int>& mosym, IndexOrder moorder = other);
  ~MOIndexSpace();

  void save_data_state(StateOut&);

  /// Returns the AO basis set
  Ref<GaussianBasisSet> basis() const;  
  /// Returns the coefficient matrix
  RefSCMatrix coefs() const;
  /// Returns the orbital symmetry array
  vector<int> mosym() const;
  /// Returns the order of the orbitals
  IndexOrder moorder() const;
  /// Returns the rank of the space
  int rank() const;
  /// Returns the rank of the full space
  int full_rank() const;
  /// Returns the number of blocks
  int nblocks() const;
  /// Returns the first orbital index in each block
  vector<int> first_mo() const;
  /// Returns the number of orbitals in each block
  vector<int> nmo() const;
  /// Returns the full-space index of the first orbital in each block
  vector<int> offsets() const;

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


