//
// moindexspace.h
//
// Copyright (C) 2004 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_moindexspace_h
#define _chemistry_qc_mbptr12_moindexspace_h

#include <vector>
#include <stdexcept>
#include <util/ref/ref.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/group/thread.h>
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

  /// max length of string returned by id()
  static const unsigned int max_id_length = 10;
  std::string id_;                        // see documentation for id()
  std::string name_;                      // String identifier for the orbital space

  Ref<GaussianBasisSet> basis_;    // The AO basis
  Ref<Integral> integral_;
  RefSCMatrix coefs_;              // AO-> MO transformation coefficients (nao by nmo matrix)
  RefDiagSCMatrix evals_;          // "eigenvalues" associated with the MOs
  RefSCDimension modim_;           // The MO dimension
  std::vector<unsigned int> mosym_;              // irrep of each orbital

  unsigned int rank_;          // The rank of this space
  unsigned int nblocks_;       // Number of blocks
  vector<unsigned int> nmo_;        // Number of MOs in each block

  IndexOrder moorder_;

  // checks mosym_ for irrep indices outside of the allowed range
  void check_mosym() const;

  // determines nmo_ from nfzc, nfzv, and evals and returns offsets for each block
  std::vector<unsigned int> frozen_to_blockinfo(unsigned int nfzc, unsigned int nfzv,
                                                const RefDiagSCMatrix& evals);

  // computes coefficient matrix from the full coefficient matrix. If moorder_ == energy
  // then the MO vectors will be sorted by their eigenvalues
  void full_coefs_to_coefs(const RefSCMatrix& full_coefs, const RefDiagSCMatrix& evals,
                           const std::vector<unsigned int>& offsets);

  // initialize the object
  void init();

  // sorting functions borrowed from mbpt.cc
  static void dquicksort(double *item,unsigned int *index,unsigned int n);


public:
  MOIndexSpace(StateIn&);
  /// Copy constructor
  MOIndexSpace(const MOIndexSpace&);
  /** This function constructs an MOIndexSpace from (blocked) space full_coefs.
      Block i will contain vectors [ offsets[i], offsets[i]+nmopi[i]-1 ] . By default,
      the space maintains the same blocked structure and the same order within blocks
      as the original space (moorder=symmetry). If moorder=energy and eigenvalues
      evals are provided, then all vectors will be put in one block and
      sorted according to ascending evals.

      \param full_coefs -- symmetry-blocked transformation coefficient matrix
      (AO by MO) for the full space
      \param basis -- basis set
      \param integral -- integral factory
      \param offsets -- block offsets
      \param nmopi -- new block sizes
      \param moorder -- specifies new ordering of vectors
      \param evals -- used to sort the vectors
      */
  MOIndexSpace(const std::string& id, const std::string& name, const RefSCMatrix& full_coefs,
               const Ref<GaussianBasisSet>& basis, const Ref<Integral>& integral,
               const std::vector<unsigned int>& offsets, const std::vector<unsigned int>& nmopi,
               const IndexOrder& moorder = symmetry,
               const RefDiagSCMatrix& evals = 0);
  /** This constructor should be used when the MOIndexSpace object is a subspace of a full orbital space.
      Similarly to the previous constructor, it constructs an MOIndexSpace object using a symmetry-blocked
      transformation coefficient matrix (AO by MO) for the full space,
      basis set, "eigenvalues" and the number of orbitals with lowest (nfzc) and highest (nfzv) eigenvalues
      to be dropped. The orbitals in the constructed space are ordered by energy. */
  MOIndexSpace(const std::string& id, const std::string& name, const RefSCMatrix& full_coefs,
               const Ref<GaussianBasisSet>& basis, const Ref<Integral>& integral,
               const RefDiagSCMatrix& evals, unsigned int nfzc, unsigned int nfzv,
               const IndexOrder& moorder = energy);
  /** This constructor should be used when the MOIndexSpace object is the full orbital space.
      The orbitals will be symmetry-blocked. */
  MOIndexSpace(const std::string& id, const std::string& name, const RefSCMatrix& full_coefs,
               const Ref<GaussianBasisSet>& basis, const Ref<Integral>& integral);

  /** This constructor is a true hack introduced because I have no idea how to construct what I need.
      It will copy orig_space but replace its coefs with new_coefs, and its basis with new_basis. */
  MOIndexSpace(const std::string& id, const std::string& name,
               const Ref<MOIndexSpace>& orig_space, const RefSCMatrix& new_coefs,
               const Ref<GaussianBasisSet>& new_basis);
  ~MOIndexSpace();

  void save_data_state(StateOut&);

  /// Returns a self-contained expressive label
  const std::string& name() const;
  /// Returns a short (preferably, one, max 10 character) identifier for the space
  /** Suggested convention:
      i -- active occupied orbitals
      a -- active unoccupied orbitals
      m -- all occupied orbitals
      e -- all virtual orbitals
      p -- all Hartree-Fock orbitals
      p' -- all RI functions
  */
  const std::string& id() const;
  /// returns the dimension correspondign to this space
  const RefSCDimension& dim() const;
  /// Returns the AO basis set
  const Ref<GaussianBasisSet>& basis() const;
  /// Returns the integral factory used to instantiate the coefficient matrix
  const Ref<Integral>& integral() const;
  /// Returns the coefficient matrix
  const RefSCMatrix& coefs() const;
  /// Returns the "eigenvalues" matrix
  const RefDiagSCMatrix& evals() const;
  /// Returns the orbital symmetry array
  const vector<unsigned int>& mosym() const;
  /// Returns the order of the orbitals
  IndexOrder moorder() const;
  /// Returns the rank of the space
  unsigned int rank() const;
  /// Returns the number of blocks
  unsigned int nblocks() const;
  /// Returns the number of orbitals in each block
  const vector<unsigned int>& nmo() const;

  /// Returns how much "significant" (i.e. O^2) memory this object uses
  size_t memory_in_use() const;

  /// Prints out this
  void print(std::ostream&o=ExEnv::out0()) const;
  /// Produces a short summary
  void print_summary(std::ostream& os) const;
  /// Prints out this in details (coefficients, etc.)
  void print_detail(std::ostream&o=ExEnv::out0()) const;

};

class CannotConstructMap: public std::logic_error {
  public:
  CannotConstructMap() : std::logic_error("Cannot map given MOIndexSpaces") {}
};

/** s2<<s1 returns map from s1 to s2. Throws CannotConstructMap if map cannot be constructed.
    Map can be constructed if and only if:
    1) s1.basis() == s2.basis()
    2) s1.integral() == s2.integral()
    2) s1.rank() <= s2.rank()
    3) for every MO in s1 there is an identical (for now, including phase) MO in s2
  */
typedef std::vector<unsigned int> MOIndexMap;
MOIndexMap
operator<<(const MOIndexSpace& space2, const MOIndexSpace& space1);

/** sparsemap(s2,s1) returns a sparse one-to-one map from s1 to s2. Throws CannotConstructMap if map cannot be constructed.
    Map can be constructed if and only if:
    1) s1.basis() == s2.basis()
    2) s1.integral() == s2.integral()
    2) s1.rank() <= s2.rank()
    3) for every MO in s1 there is an MO in s2 that differs by at most the sign.
  */
typedef std::vector< std::pair<unsigned int,double> > SparseMOIndexMap;
SparseMOIndexMap
sparsemap(const MOIndexSpace& space2, const MOIndexSpace& space1, double hardzero=1e-12);

/** transform(s2,s1) returns matrix that transforms s1 to s2.
    Throws CannotConstructMap if the transform cannot be constructed.
    The transform can only be constructed if
    1) s1.integral() == s2.integral()
    2) s1.rank() >= s2.rank()
    3) overlap of s1.basis() and s2.basis() is not zero

    the returned matrix has dimensions of s2.coefs().coldim() and s1.coefs().coldim() and is allocated using kit.
  */
RefSCMatrix
transform(const MOIndexSpace& space2, const MOIndexSpace& space1, const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit());

/** overlap(s2,s1) returns the overlap matrix between s2 and s1. It can be computed if
    1) s1.integral() is compatible with s2.integral()
    The matrix has dimensions of s2.coefs().coldim() and s1.coefs().coldim() and is allocated using kit.
    Throws if the overlap cannot be computed.
  */
RefSCMatrix
overlap(const MOIndexSpace& space2, const MOIndexSpace& space1, const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit());

/** in(s1,s2) returns true if s1 is in s2 */
bool
in(const MOIndexSpace& s1, const MOIndexSpace& s2);


/** Registry of globally-known MOIndexSpaces associates MOIndexSpace S with a string key given by S->id().
    It is a singleton. All operations on it are thread-safe.
 */
class MOIndexSpaceRegistry : virtual public SavableState {
  public:
    static const Ref<MOIndexSpaceRegistry>& instance();

    /// returns MOIndexSpace that corresponds to this key. If key is not known, returns null pointer
    Ref<MOIndexSpace> find(const std::string& key) const;
    /// registers this MOIndexSpace
    void add(const Ref<MOIndexSpace>& space);

  private:
    // this is a Singleton: access only through instance()
    MOIndexSpaceRegistry();

    static Ref<MOIndexSpaceRegistry> instance_;

    typedef std::map< std::string, Ref<MOIndexSpace> > MOIndexSpaceMap;
    MOIndexSpaceMap space_map_;
    // std::map's operations are not reentrant, hence lock the map every time
    Ref<ThreadLock> lock_;
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:


