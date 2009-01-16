//
// orbitalspace.h
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

#ifndef _chemistry_qc_mbptr12_orbitalspace_h
#define _chemistry_qc_mbptr12_orbitalspace_h

#include <vector>
#include <stdexcept>
#include <util/ref/ref.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/group/thread.h>
#include <math/scmat/abstract.h>
#include <util/state/statein.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/registry.h>

using namespace std;

namespace sc {

  /** Class OrbitalSpace describes a range of orbitals
   that are linear combinations of Gaussian basis functions
   (e.g. atomic orbitals). In general, such sets are subspaces of
   a full space of orbitals supported by the given basis. Orbitals
   can be symmetry-blocked, ordered by energy, etc.
   Examples of sets that can be described
   using OrbitalSpace are occupied MOs and virtual MOs. */

  class OrbitalSpace: virtual public SavableState {

    public:

      ///Describes the ordering of indices
      enum IndexOrder {
        symmetry = 0, ///< symmetry corresponds to orbitals ordered by symmetry, then (i.e. within each symmetry block) by energy.
        energy = 1, ///< energy corresponds to orbitals ordered by energy.
        correlated = 2
      /**< correlated corresponds to orbitals ordered as: frozen occupied, active occupied, active virtual, frozen virtual;
       within each block orbitals are ordered by symmetry, and within each symmetry block orbitals are
       ordered by energy.*/
      };

    private:

      static const unsigned int max_id_length = 30; // max length of string returned by id()
      std::string id_; // see documentation for id()
      std::string name_; // String identifier for the orbital space

      Ref<GaussianBasisSet> basis_; // The AO basis
      Ref<Integral> integral_;
      RefSCMatrix coefs_; // AO->MO transformation coefficients (nao by nmo matrix)
      RefDiagSCMatrix evals_; // "eigenvalues" associated with the MOs
      RefSCDimension dim_; // only here to allow dim() return const &
      std::vector<unsigned int> orbsym_; // irrep of each orbital
      std::vector<unsigned int> block_sizes_; // Number of orbitals in each block

      // checks orbsym_ for irrep indices outside of the allowed range
      void check_orbsym() const;

      // determines block_sizes_ from nfzc, nfzv, and evals and returns offsets for each block
      std::vector<unsigned int>
          frozen_to_blockinfo(unsigned int nfzc, unsigned int nfzv,
                              const RefDiagSCMatrix& evals);

      // computes coefficient matrix from the full coefficient matrix. If moorder_ == energy
      // then the MO vectors will be sorted by their eigenvalues
      void full_coefs_to_coefs(const RefSCMatrix& full_coefs,
                               const RefDiagSCMatrix& evals, const std::vector<
                                   unsigned int>& offsets, IndexOrder moorder);

      // initialize the object
      void init();

      // sorting functions borrowed from mbpt.cc
      static void dquicksort(double *item, unsigned int *index, unsigned int n);

    public:
      OrbitalSpace(StateIn&);
      /// Copy constructor
      OrbitalSpace(const OrbitalSpace&);
      /** This function constructs an OrbitalSpace from a set of vectors whose coefficients are given by
       full_coefs. Block i will contain vectors [ block_offsets[i], block_offsets[i]+block_sizes[i] ).
       The number of blocks is given by block_offsets.size(); it must be identical to
       block_sizes.size().

       \param full_coefs -- symmetry-blocked transformation coefficient matrix
       (AO by MO) for the full space
       \param evals -- the orbital energies
       \param basis -- basis set
       \param integral -- integral factory
       \param block_offsets -- block offsets
       \param block_sizes -- new block sizes
       \param moorder -- describes the ordering of vectors
       */
      OrbitalSpace(const std::string& id, const std::string& name,
                   const RefSCMatrix& full_coefs, const RefDiagSCMatrix& evals,
                   const Ref<GaussianBasisSet>& basis,
                   const Ref<Integral>& integral,
                   const std::vector<unsigned int>& block_offsets,
                   const std::vector<unsigned int>& block_sizes,
                   const IndexOrder& moorder);

      /** This function constructs an OrbitalSpace from (blocked) space full_coefs.
       Block i will contain vectors [ offsets[i], offsets[i]+nmopi[i]-1 ] . By default,
       the space maintains the same blocked structure and the same order within blocks
       as the original space (moorder=symmetry). If moorder=energy and eigenvalues
       evals are provided, then all vectors will be put in one block and
       sorted according to ascending evals.

       \param full_coefs -- symmetry-blocked transformation coefficient matrix
       (AO by MO) for the full space
       \param basis -- basis set
       \param integral -- integral factory
       \param block_offsets -- block offsets
       \param block_sizes -- new block sizes
       \param moorder -- specifies new ordering of vectors
       \param evals -- used to sort the vectors
       */
      OrbitalSpace(const std::string& id, const std::string& name,
                   const RefSCMatrix& full_coefs,
                   const Ref<GaussianBasisSet>& basis,
                   const Ref<Integral>& integral,
                   const std::vector<unsigned int>& block_offsets,
                   const std::vector<unsigned int>& block_sizes,
                   const IndexOrder& moorder = symmetry,
                   const RefDiagSCMatrix& evals = 0);
      /** This constructor should be used when the OrbitalSpace object is a subspace of a full orbital space.
       Similarly to the previous constructor, it constructs an OrbitalSpace object using a symmetry-blocked
       transformation coefficient matrix (AO by MO) for the full space,
       basis set, "eigenvalues" and the number of orbitals with lowest (nfzc) and highest (nfzv) eigenvalues
       to be dropped. The orbitals in the constructed space are ordered by energy. */
      OrbitalSpace(const std::string& id, const std::string& name,
                   const RefSCMatrix& full_coefs,
                   const Ref<GaussianBasisSet>& basis,
                   const Ref<Integral>& integral, const RefDiagSCMatrix& evals,
                   unsigned int nfzc, unsigned int nfzv,
                   const IndexOrder& moorder = energy);
      /** This constructor should be used when the OrbitalSpace object is the full orbital space.
       The orbitals will be symmetry-blocked. */
      OrbitalSpace(const std::string& id, const std::string& name,
                   const RefSCMatrix& full_coefs,
                   const Ref<GaussianBasisSet>& basis,
                   const Ref<Integral>& integral);

      /** This constructor is a true hack introduced because I have no idea how to construct what I need.
       It will copy orig_space but replace its coefs with new_coefs, and its basis with new_basis. */
      OrbitalSpace(const std::string& id, const std::string& name, const Ref<
          OrbitalSpace>& orig_space, const RefSCMatrix& new_coefs, const Ref<
          GaussianBasisSet>& new_basis);
      ~OrbitalSpace();

      void save_data_state(StateOut&);

      OrbitalSpace& operator=(const OrbitalSpace& other);

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
      /// returns the dimension corresponding to this space
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
      const vector<unsigned int>& orbsym() const;
      /// Returns the rank of the space
      unsigned int rank() const;
      /// Returns the number of blocks
      unsigned int nblocks() const;
      /// Returns the number of orbitals in each block
      const vector<unsigned int>& block_sizes() const;

      /// Returns how much "significant" (i.e. O^2) memory this object uses
      size_t memory_in_use() const;

      /// Prints out this
      void print(std::ostream&o = ExEnv::out0()) const;
      /// Produces a short summary
      void print_summary(std::ostream& os) const;
      /// Prints out this in details (coefficients, etc.)
      void print_detail(std::ostream&o = ExEnv::out0()) const;

  };

  bool operator==(const OrbitalSpace& space1, const OrbitalSpace& space2);

  class CannotConstructMap: public std::logic_error {
    public:
      CannotConstructMap() :
        std::logic_error("Cannot map given OrbitalSpaces") {
      }
  };

  /** s2<<s1 returns map from s1 to s2. Throws CannotConstructMap if map cannot be constructed.
   Map can be constructed if and only if:
   1) s1.basis() == s2.basis()
   2) s1.integral() == s2.integral()
   2) s1.rank() <= s2.rank()
   3) for every MO in s1 there is an identical (for now, including phase) MO in s2
   */
  typedef std::vector<unsigned int> MOIndexMap;
  MOIndexMap operator<<(const OrbitalSpace& space2, const OrbitalSpace& space1);

  /** sparsemap(s2,s1) returns a sparse one-to-one map from s1 to s2. Throws CannotConstructMap if map cannot be constructed.
   Map can be constructed if and only if:
   1) s1.basis() == s2.basis()
   2) s1.integral() == s2.integral()
   2) s1.rank() <= s2.rank()
   3) for every MO in s1 there is an MO in s2 that differs by at most the sign.
   */
  typedef std::vector<std::pair<unsigned int, double> > SparseMOIndexMap;
  SparseMOIndexMap sparsemap(const OrbitalSpace& space2,
                             const OrbitalSpace& space1, double hardzero =
                                 1e-12);

  /** transform(s2,s1) returns matrix that transforms s1 to s2.
   Throws CannotConstructMap if the transform cannot be constructed.
   The transform can only be constructed if
   1) s1.integral() == s2.integral()
   2) s1.rank() >= s2.rank()
   3) overlap of s1.basis() and s2.basis() is not zero

   the returned matrix has dimensions of s2.coefs().coldim() and s1.coefs().coldim() and is allocated using kit.
   */
  RefSCMatrix transform(const OrbitalSpace& space2, const OrbitalSpace& space1,
                        const Ref<SCMatrixKit>& kit =
                            SCMatrixKit::default_matrixkit());

  /** overlap(s2,s1) returns the overlap matrix between s2 and s1. It can be computed if
   1) s1.integral() is compatible with s2.integral()
   The matrix has dimensions of s2.coefs().coldim() and s1.coefs().coldim() and is allocated using kit.
   Throws if the overlap cannot be computed.
   */
  RefSCMatrix overlap(const OrbitalSpace& space2, const OrbitalSpace& space1,
                      const Ref<SCMatrixKit>& kit =
                          SCMatrixKit::default_matrixkit());

  /** in(s1,s2) returns true if s1 is in s2 */
  bool in(const OrbitalSpace& s1, const OrbitalSpace& s2);

  /**
   Parses keys of OrbitalSpace. Although OrbitalSpace object with arbitrary keys can be created,
   other components of MPQC assume certain rules for the keys. The simplest rule is that a key for
   a spin-polarized OrbitalSpace (i.e., that has spin-alpha and spin-beta variants) must encode the spin information.
   This class can encode and decode this information.

   Another rule is how transformed OrbitalSpace objects are labeled. \sa ParsedTranformedOrbitalSpaceKey
   */
  class ParsedOrbitalSpaceKey {
    public:
      /// throws if key is not properly formatted
      ParsedOrbitalSpaceKey(const std::string& key);

      const std::string& key() const {
        return key_;
      }
      const std::string& label() const {
        return label_;
      }
      SpinCase1 spin() const {
        return spin_;
      }

      static std::string key(const std::string& label, SpinCase1 spin);

    private:
      std::string key_;
      std::string label_;
      SpinCase1 spin_;
  };

  /**
   Parses keys of a "transformed" OrbitalSpace. Although OrbitalSpace object with arbitrary keys can be created,
   other components of MPQC assume certain rules for the keys. "Transformed" OrbitalSpace
   represents space U obtained by transformation of space V as
   U_i^j = C_i^k V_k^j, where V is the coefficient matrix of the original space, C is the transformation
   matrix, and U is the transformed space. Transformation matrices typically represent common
   one-electron operators, e.g., Fock, Coulomb, exchange, kinetic energy, etc.
   */
  class ParsedTransformedOrbitalSpaceKey {
    public:
      typedef enum {
        Fock = 0,
        Coulomb = 1,
        Exchange = 2,
        KineticEnergy = 3,
        CoreHamiltonian = 4,
        CorePlusCoulomb = 5,
        InvalidOneBodyOperator = 6
      } OneBodyOperator;
#define ParsedTransformedOrbitalSpaceKey_OneBodyOperatorLabelInitializer {"F", "J", "K", "T", "h", "h+J"}
      static const char* OneBodyOperatorLabel[]; // = ParsedTransformedOrbitalSpaceKey_OneBodyOperatorLabelInitializer

      /// throws if key is not properly formatted
      ParsedTransformedOrbitalSpaceKey(const std::string& key);

      const std::string& key() const {
        return key_;
      }
      const std::string& label() const {
        return label_;
      }
      SpinCase1 spin() const {
        return spin_;
      }
      const std::string& original_label() const {
        return original_label_;
      }
      SpinCase1 original_spin() const {
        return spin_;
      }
      OneBodyOperator transform_operator() const {
        return transform_operator_;
      }

      static std::string key(const std::string& label, SpinCase1 spin,
                             const std::string& original_label,
                             SpinCase1 original_spin, OneBodyOperator oper);

    private:
      std::string key_;
      std::string label_;
      SpinCase1 spin_;
      std::string original_label_;
      SpinCase1 original_spin_;
      OneBodyOperator transform_operator_;
  };

  /// registry of globally-known OrbitalSpace objects
  typedef Registry<std::string, Ref<OrbitalSpace> ,
      detail::SingletonCreationPolicy> OrbitalSpaceRegistry;
  /// helper function to form a key/space pair from a OrbitalSpace
  std::pair<std::string, Ref<OrbitalSpace> > make_keyspace_pair(const Ref<
      OrbitalSpace>& space, SpinCase1 spin = AnySpinCase1);

  /// registry of globally-known OrbitalSpace objects that describe AO basis spaces
  typedef Registry<Ref<GaussianBasisSet> , Ref<OrbitalSpace> ,
      detail::SingletonCreationPolicy> AOSpaceRegistry;

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:


