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
      (e.g. atomic orbitals). In general, such sets are subspaces of
      a full space of orbitals supported by the given basis. Orbitals
      can be symmetry-blocked, ordered by energy, etc.
      Examples of sets that can be described
      using MOIndexSpace are occupied MOs and virtual MOs. */

class MOIndexSpace : virtual public SavableState {

public:

  /// Describes the ordering of indices
  enum IndexOrder { symmetry = 0, energy = 1, undefined = 2 };

private:

  std::string name_;                      // String identifier for the orbital space
  
  Ref<GaussianBasisSet> basis_;    // The AO basis
  RefSCMatrix coefs_;              // AO-> MO transformation coefficients (nao by nmo matrix)
  RefDiagSCMatrix evals_;          // "eigenvalues" associated with the MOs
  RefSCDimension modim_;           // The MO dimension
  vector<int> mosym_;              // irrep of each orbital

  int rank_;          // The rank of this space
  int full_rank_;     // Rank of the full space, i.e. number of MOs
  int nblocks_;       // Number of blocks
  vector<int> nmo_;        // Number of MOs in each block
  vector<int> offsets_;    // Index of the first MO in each block relative to the first MO of that block in full space
  vector<int> map_to_full_space_;  // Full-space index

  IndexOrder moorder_;

  // checks mosym_ for irrep indices outside of the allowed range
  void check_mosym() const;

  // determines offsets_ and nmo_ from nfzc, nfzv, and evals
  void frozen_to_blockinfo(const int nfzc, const int nfzv, const RefDiagSCMatrix& evals);

  // computes coefficient matrix from the full coefficient matrix. If moorder_ == energy
  // then the MO vectors will be sorted by their eigenvalues
  void full_coefs_to_coefs(const RefSCMatrix& full_coefs, const RefDiagSCMatrix& evals);

  // initialize the object
  void init();
  
  // sorting functions borrowed from mbpt.cc
  static void dquicksort(double *item,int *index,int n);
  static void dqs(double *item,int *index,int left,int right);
  

public:
  MOIndexSpace(StateIn&);
  /** This function constructs an MOIndexSpace from (blocked) space full_coefs.
      Block i will contain vectors [ offsets[i], offsets[i]+nmopi[i]-1 ] . By default,
      the space maintains the same blocked structure and the same order within blocks
      as the original space (moorder=symmetry). If moorder=energy and eigenvalues
      evals are provided, then all vectors will be put in one block and
      sorted according to ascending evals.

      \param full_coefs -- symmetry-blocked transformation coefficient matrix
      (AO by MO) for the full space
      \param basis -- basis set
      \param offsets -- block offsets
      \param nmopi -- new block sizes
      \param moorder -- specifies new ordering of vectors
      \param evals -- used to sort the vectors
      */
  MOIndexSpace(std::string name, const RefSCMatrix& full_coefs, const Ref<GaussianBasisSet> basis,
               const vector<int>& offsets, const vector<int>& nmopi, IndexOrder moorder = symmetry,
               const RefDiagSCMatrix& evals = 0);
  /** This constructor should be used when the MOIndexSpace object is a subspace of a full orbital space.
      Similarly to the previous constructor, it constructs an MOIndexSpace object using a symmetry-blocked
      transformation coefficient matrix (AO by MO) for the full space,
      basis set, "eigenvalues" and the number of orbitals with lowest (nfzc) and highest (nfzv) eigenvalues
      to be dropped. The orbitals in the constructed space are ordered by energy. */
  MOIndexSpace(std::string name, const RefSCMatrix& full_coefs, const Ref<GaussianBasisSet> basis,
               const RefDiagSCMatrix& evals, int nfzc, int nfzv, IndexOrder moorder = energy);
  /** This constructor should be used when the MOIndexSpace object is the full orbital space.
      The orbitals will be symmetry-blocked. */
  MOIndexSpace(std::string name, const RefSCMatrix& full_coefs, const Ref<GaussianBasisSet> basis);

  /* This constructor should be used when the MOIndexSpace object is the full orbital space.
     Constructs an MOIndexSpace object using a non-blocked transformation coefficient matrix
     (AO by MO) for this space, basis set, and optional order info and MO irreps. */
  /*MOIndexSpace(std::string name, const RefSCMatrix& coefs, const Ref<GaussianBasisSet> basis,
               IndexOrder moorder = undefined, const vector<int>& mosym = 0,
               const RefDiagSCMatrix& evals = 0);*/
  /** This constructor is a true hack introduced because I have no idea how to construct what I need.
      It will copy orig_space but replace it's coefs with new_coefs, and its basis with new_basis. */
  MOIndexSpace(std::string name, const Ref<MOIndexSpace>& orig_space, const RefSCMatrix& new_coefs,
               const Ref<GaussianBasisSet>& new_basis);
  ~MOIndexSpace();

  void save_data_state(StateOut&);

  /// Returns the AO basis set
  const std::string& name() const;  
  /// Returns the AO basis set
  const Ref<GaussianBasisSet>& basis() const;  
  /// Returns the coefficient matrix
  const RefSCMatrix& coefs() const;
  /// Returns the "eigenvalues" matrix
  const RefDiagSCMatrix& evals() const;
  /// Returns the orbital symmetry array
  const vector<int>& mosym() const;
  /// Returns the order of the orbitals
  IndexOrder moorder() const;
  /// Returns the rank of the space
  int rank() const;
  /// Returns the rank of the full space
  int full_rank() const;
  /// Returns the number of blocks
  int nblocks() const;
  /// Returns the number of orbitals in each block
  const vector<int>& nmo() const;
  /// Returns the full-space index of the first orbital in each block
  const vector<int>& offsets() const;
  /// Returns the full-space index
  int to_full_space(const int i) const;

  /// Returns how much "significant" (i.e. O^2) memory this object uses
  size_t memory_in_use() const;

  /// Prints out this
  void print(std::ostream&o=ExEnv::out0()) const;
  /// Produces a short summary
  void print_summary(std::ostream& os) const;

};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:


