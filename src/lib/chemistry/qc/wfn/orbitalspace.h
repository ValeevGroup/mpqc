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

#ifndef _chemistry_qc_mbptr12_orbitalspace_h
#define _chemistry_qc_mbptr12_orbitalspace_h

#include <bitset>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <util/ref/ref.h>
#include <util/state/statein.h>
#include <util/state/stateout.h>
#include <util/group/thread.h>
#include <math/scmat/abstract.h>
#include <util/state/statein.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/basis/operator.h>
#include <chemistry/qc/wfn/spin.h>
#include <util/misc/registry.h>

namespace sc {

  /// @addtogroup ChemistryElectronicStructureOneBody
  /// @{

  /// Orbital = index + attributes
  template <typename Attributes>
  class DecoratedOrbital {
    public:
      DecoratedOrbital() {}
      DecoratedOrbital(size_t index, const Attributes& attr) : index_(index), attr_(attr) {}

      size_t index() const { return index_; }
      const Attributes& attr() const { return attr_; }

    private:
      size_t index_;
      Attributes attr_;
  };

  /// Orbital in a blocked space
  typedef DecoratedOrbital< unsigned int > BlockedOrbital;

    /// MO is irrep, energy, occupation number
    struct MolecularOrbitalAttributes {
      public:
        MolecularOrbitalAttributes(unsigned int irrep, double energy,
                                   double occnum) :
          irrep_(irrep), energy_(energy), occnum_(occnum) {
        }

        unsigned int irrep() const {
          return irrep_;
        }
        double energy() const {
          return energy_;
        }
        double occnum() const {
          return occnum_;
        }

      private:
        unsigned int irrep_;
        double energy_;
        double occnum_;
    };
    /// Same as MolecularOrbitalAttributes, plus spin
    struct MolecularSpinOrbitalAttributes : public MolecularOrbitalAttributes {
      public:
        MolecularSpinOrbitalAttributes(unsigned int irrep,
                                       double energy,
                                       double occnum,
                                       SpinCase1 spin) :
          MolecularOrbitalAttributes(irrep, energy, occnum), spin_(spin)
          {}

        SpinCase1 spin() const { return spin_; }

      private:
        SpinCase1 spin_;
    };


  /**
   * Describes particle-hole attributes of orbitals
   */
  struct ParticleHoleOrbitalAttributes : public std::bitset<2> {
    ParticleHoleOrbitalAttributes() {}
    ParticleHoleOrbitalAttributes(unsigned long val) : std::bitset<2>(val) {}
    ParticleHoleOrbitalAttributes(const std::bitset<2>& bs) : std::bitset<2>(bs) {}

    static ParticleHoleOrbitalAttributes Hole;          //!< only holes can be created
    static ParticleHoleOrbitalAttributes Particle;      //!< only particles can be created
    static ParticleHoleOrbitalAttributes Any;           //!< holes and particles can be created
    static ParticleHoleOrbitalAttributes None;          //!< neither holes nor particles can be created
  };

  typedef DecoratedOrbital< MolecularOrbitalAttributes > MolecularOrbital;
  typedef DecoratedOrbital< MolecularSpinOrbitalAttributes > MolecularSpinOrbital;

  /// order by symmetry first, then by energy, then by occ num
  struct SymmetryMOOrder {
    public:
      SymmetryMOOrder(unsigned int nirreps) : nirreps_(nirreps) {}

      bool operator()(const MolecularOrbital& o1, const MolecularOrbital& o2) const {
        if (o1.attr().irrep() < o2.attr().irrep())
          return true;
        else if (o1.attr().irrep() == o2.attr().irrep()) {
          if (o1.attr().energy() < o2.attr().energy())
            return true;
          else if (o1.attr().energy() == o2.attr().energy()) {
            if (o1.attr().occnum() < o2.attr().occnum())
              return true;
          }
        }
        return false;
      }
      unsigned int block(const MolecularOrbital& o) const {
        return o.attr().irrep();
      }
      unsigned int nblocks() const {
        return nirreps_;
      }

    private:
      unsigned int nirreps_;
  };

  /// order by energy first, then by symmetry. EnergyCompare specifies the weak strict ordering of orbitals wrt energy
  template <typename EnergyCompare = std::less<double> > struct EnergyMOOrder {
    public:
      bool operator()(const MolecularOrbital& o1, const MolecularOrbital& o2) const {
        if ( ecompare(o1.attr().energy(), o2.attr().energy()) )
          return true;
        else if (o1.attr().energy() == o2.attr().energy()) {
          if (o1.attr().irrep() < o2.attr().irrep())
            return true;
        }
        return false;
      }
      unsigned int block(const MolecularOrbital& o) const {
        return 0;
      }
      unsigned int nblocks() const {
        return 1;
      }
    private:
      EnergyCompare ecompare;
  };

  /// order by occupation first, then by symmetry, then by energy
  struct CorrelatedMOOrder {
    public:
      CorrelatedMOOrder(unsigned int nirreps) : nirreps_(nirreps) {}

      bool operator()(const MolecularOrbital& o1, const MolecularOrbital& o2) const {
        // occupieds come before virtuals
        if (o1.attr().occnum() > o2.attr().occnum())
          return true;
        else if (o1.attr().occnum() == o2.attr().occnum()) {
          if (o1.attr().irrep() < o2.attr().irrep())
            return true;
          else if (o1.attr().irrep() == o2.attr().irrep()) {
            if (o1.attr().energy() < o2.attr().energy())
              return true;
          }
        }
        return false;
      }

      unsigned int block(const MolecularOrbital& o) const {
        const unsigned int irrep = o.attr().irrep();
        const double occnum = o.attr().occnum();
        // occupieds come before virtuals
        const int occblock = (occnum == 1.0) ? 0 : 1;
        return occblock * nirreps_ + irrep;
      }
      unsigned int nblocks() const {
        return nirreps_ * 2;
      }

    private:
      unsigned int nirreps_;
  };

  /// order by occupation first, then by spin, then by symmetry, then by energy
  struct CorrelatedSpinMOOrder {
    public:
      CorrelatedSpinMOOrder(unsigned int nirreps) : nirreps_(nirreps) {}

      bool operator()(const MolecularSpinOrbital& o1, const MolecularSpinOrbital& o2) const {
        // occupieds come before virtuals
        if (o1.attr().occnum() > o2.attr().occnum())
          return true;
        else if (o1.attr().occnum() == o2.attr().occnum()) {
          if (o1.attr().spin() < o2.attr().spin())
            return true;
          else if (o1.attr().spin() == o2.attr().spin()) {
            if (o1.attr().irrep() < o2.attr().irrep())
              return true;
            else if (o1.attr().irrep() == o2.attr().irrep()) {
              if (o1.attr().energy() < o2.attr().energy())
                return true;
            }
          }
        }
        return false;
      }

      unsigned int block(const MolecularSpinOrbital& o) const {
        const unsigned int irrep = o.attr().irrep();
        const SpinCase1 spin = o.attr().spin();
        const unsigned int spincase = (spin == Alpha) ? 0 : 1;
        const double occnum = o.attr().occnum();
        const unsigned int occblock = (occnum == 1.0) ? 0 : 1;
        // occupieds come before virtuals, Alpha before Beta
        const unsigned int result = (occblock * 2 + spincase) * nirreps_ + irrep;
        return result;
      }
      unsigned int nblocks() const {
        return nirreps_ * 4;
      }

    private:
      unsigned int nirreps_;
  };

  namespace detail {

    template <typename Container> struct ContainerAdaptor {
      public:
        typedef typename Container::value_type value_type;
        ContainerAdaptor(const Container& cont) : cont_(cont) {}
        size_t size() const { return cont_.size(); }
        value_type elem(size_t i) const { return cont_[i]; }

      private:
        const Container& cont_;
    };

    template<> struct ContainerAdaptor<RefDiagSCMatrix> {
      public:
        typedef RefDiagSCMatrix Container;
        typedef double value_type;
        ContainerAdaptor(const Container& cont) : cont_(cont) {}
        size_t size() const { return cont_.dim().n(); }
        value_type elem(size_t i) const { return cont_(i); }

      private:
        const Container& cont_;
    };

  } // namespace detail

  /// mask out first n MOs in the order defined by Compare. By default mask the n lowest-energy MOs
  template <typename Attribute,
            typename AttributeContainer,
            typename Compare = std::less<Attribute> >
  struct MolecularOrbitalMask {
    private:
      typedef DecoratedOrbital<Attribute> MO;
      struct _compare {
        bool operator()(const MO& mo1,
            const MO& mo2) const {
          Compare comp;
          return comp(mo1.attr(), mo2.attr());
        }
      };

    public:
      MolecularOrbitalMask(unsigned int n,
                           const AttributeContainer& attributes) :
        mask_(detail::ContainerAdaptor<AttributeContainer>(attributes).size(),true)
        {
          // validate input
          if (n == 0) return;
          const size_t nmos = mask_.size();
          MPQC_ASSERT(n < nmos);

          // copy attributes to vector of MOs
          std::vector<MO> mos;
          detail::ContainerAdaptor<AttributeContainer> contadaptor(attributes);
          for(size_t mo=0; mo<nmos; ++mo) {
            mos.push_back(MO(mo,contadaptor.elem(mo)));
          }

          // sort
          _compare comp;
          std::stable_sort(mos.begin(), mos.end(), comp);

          // copy to mask
          for(unsigned int t=0; t<n; ++t) {
            mask_[ mos[t].index() ] = false;
          }
        }

      const std::vector<bool>& mask() const { return mask_; }

      bool operator[](size_t o) const { return mask_[o]; }

    private:
      std::vector<bool> mask_;
  };


  ///////////////////////////////////////////////////////////////////////////

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
        correlated = 2, /**< correlated corresponds to orbitals ordered as: frozen occupied, active occupied, active virtual, frozen virtual;
       within each block orbitals are ordered by symmetry, and within each symmetry block orbitals are
       ordered by energy.*/
        general = 3 ///< any other
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
      std::vector<unsigned int> block_offsets_; // Start of each block
      std::vector<unsigned int> block_sizes_; // Number of orbitals in each block

      // checks orbsym_ for irrep indices outside of the allowed range
      void check_orbsym() const;

      // determines block_sizes_ from nfzc, nfzv, and evals and returns offsets for each block
      std::vector<unsigned int>
          frozen_to_blockinfo(unsigned int nfzc, unsigned int nfzv,
                              const RefDiagSCMatrix& evals);

      // computes coefficient matrix from the full coefficient matrix. If moorder == energy
      // then the MO vectors will be sorted by their eigenvalues
      void full_coefs_to_coefs(const RefSCMatrix& full_coefs,
                               const RefDiagSCMatrix& evals,
                               const std::vector<unsigned int>& offsets,
                               IndexOrder moorder);
      /// initialize the object
      void init();

      // sorting functions borrowed from mbpt.cc
      static void dquicksort(double *item, unsigned int *index, unsigned int n);

    protected:
      /// Empty constructor only useful for derived classes -- don't forget to call init()
      OrbitalSpace();
      /// initialize the object by mapping the original space to a space with indexmap
      void init(const std::string& id, const std::string& name,
                const Ref<GaussianBasisSet>& basis, const Ref<Integral>& integral,
                const RefSCMatrix& coefs,
                const RefDiagSCMatrix& evals,
                const std::vector<unsigned int>& orbsyms,
                unsigned int nblocks,
                const std::vector<BlockedOrbital>& indexmap);

    public:
      OrbitalSpace(StateIn&);
      /// Copy constructor
      OrbitalSpace(const OrbitalSpace&);
      /** This function constructs an OrbitalSpace from a set of vectors whose coefficients are given by
       full_coefs. Block i will contain vectors [ block_offsets[i], block_offsets[i]+block_sizes[i] ).
       The number of blocks is given by block_offsets.size(); it must be identical to
       block_sizes.size().

       \param id -- identifier for this orbital space
       \param name -- descriptive name for this orbital space
       \param full_coefs -- symmetry-blocked transformation coefficient matrix
       (AO by MO) for the full space
       \param evals -- the orbital energies
       \param basis -- basis set
       \param integral -- integral factory
       \param block_offsets -- block offsets
       \param block_sizes -- new block sizes
       */
      OrbitalSpace(const std::string& id, const std::string& name,
                   const RefSCMatrix& full_coefs, const RefDiagSCMatrix& evals,
                   const Ref<GaussianBasisSet>& basis,
                   const Ref<Integral>& integral,
                   const std::vector<unsigned int>& block_offsets,
                   const std::vector<unsigned int>& block_sizes);

      /** This function constructs an OrbitalSpace from (blocked) space full_coefs.
       Block i will contain vectors [ offsets[i], offsets[i]+nmopi[i]-1 ] . By default,
       the space maintains the same blocked structure and the same order within blocks
       as the original space (moorder=symmetry). If moorder=energy and eigenvalues
       evals are provided, then all vectors will be put in one block and
       sorted according to ascending evals.

       \param id -- identifier for this orbital space
       \param name -- descriptive name for this orbital space
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
      OrbitalSpace(const std::string& id, const std::string& name,
                   const Ref<OrbitalSpace>& orig_space,
                   const RefSCMatrix& new_coefs,
                   const Ref<GaussianBasisSet>& new_basis);
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
      /// Returns the coefficient matrix built with a non-blocked kit
      RefSCMatrix coefs_nb() const;
      /// Returns the "eigenvalues" matrix
      const RefDiagSCMatrix& evals() const;
      /// Returns the orbital symmetry array
      const std::vector<unsigned int>& orbsym() const;
      /// Returns the rank of the space
      unsigned int rank() const;
      /// Returns the number of blocks
      unsigned int nblocks() const;
      /// Returns the number of orbitals in each block
      const std::vector<unsigned int>& block_sizes() const;

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

  typedef std::vector<unsigned int> MOIndexMap;
  /** s2<<s1 returns map from s1 to s2. Throws CannotConstructMap if map cannot be constructed.
   Map can be constructed if and only if:
   1) s1.basis() == s2.basis()
   2) s1.integral() == s2.integral()
   2) s1.rank() <= s2.rank()
   3) for every MO in s1 there is an identical (for now, including phase) MO in s2
   */
  MOIndexMap operator<<(const OrbitalSpace& space2, const OrbitalSpace& space1);

  /** same as operator<<(), except if some orbital in space1 is not contained in space2, map it to -1.
      \param expect_same_bases setting to true will enforce space1.basis() == space2.basis()
    */
  std::vector<int> map(const OrbitalSpace& space2, const OrbitalSpace& space1, bool expect_same_bases = true);

  typedef std::vector<std::pair<unsigned int, double> > SparseMOIndexMap;
  /** sparsemap(s2,s1) returns a sparse one-to-one map from s1 to s2. Throws CannotConstructMap if map cannot be constructed.
   Map can be constructed if and only if:
   1) s1.basis() == s2.basis()
   2) s1.integral() == s2.integral()
   2) s1.rank() <= s2.rank()
   3) for every MO in s1 there is an MO in s2 that differs by at most the sign.
   */
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
   Parses keys of OrbitalSpace. Although OrbitalSpace objects with arbitrary keys can be created,
   other components of MPQC assume certain rules for the keys. The simplest rule is that a key for
   a spin-polarized OrbitalSpace (i.e., that has spin-alpha and spin-beta variants)
   must encode the spin information.
   This class can encode and decode this information.

   Another rule is how transformed OrbitalSpace objects are labeled. \sa ParsedTranformedOrbitalSpaceKey
   */
  class ParsedOrbitalSpaceKey {
    public:
      struct exception {};

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
   Parses keys of a "transformed" OrbitalSpace. Although OrbitalSpace objects with arbitrary keys can be created,
   other components of MPQC assume certain rules for the keys. "Transformed" OrbitalSpace
   represents space U obtained by transformation of space V as
   U_i^j = C_i^k V_k^j, where V is the coefficient matrix of the original space, C is the transformation
   matrix, and U is the transformed space. Transformation matrices typically represent common
   one-electron operators, e.g., Fock, Coulomb, exchange, kinetic energy, etc. (\sa OneBodyOper::type)

   The key format is "idU_operC(idV)".
   */
  class ParsedTransformedOrbitalSpaceKey {
    public:
      struct exception { };

      /// throws ParsedTransformedOrbitalSpaceKey::exception if the key is not properly formatted
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
      const std::string& support_label() const {
        return support_label_;
      }
      SpinCase1 support_spin() const {
        return spin_;
      }
      OneBodyOper::type transform_operator() const {
        return transform_operator_;
      }

      static std::string key(const std::string& label, SpinCase1 spin,
                             const std::string& original_label,
                             SpinCase1 original_spin, OneBodyOper::type oper);

      static bool valid_key(const std::string& key);

    private:
      std::string key_;
      std::string label_;
      SpinCase1 spin_;
      std::string support_label_;
      SpinCase1 support_spin_;
      OneBodyOper::type transform_operator_;
  };

  /// registry of globally-known OrbitalSpace objects
  typedef Registry<std::string, Ref<OrbitalSpace> ,
                   detail::NonsingletonCreationPolicy,
                   std::equal_to<std::string>,
                   RefObjectEqual<OrbitalSpace> > OrbitalSpaceRegistry;
  /// helper function to form a key/space pair from a OrbitalSpace
  std::pair<std::string, Ref<OrbitalSpace> > make_keyspace_pair(const Ref<
      OrbitalSpace>& space, SpinCase1 spin = AnySpinCase1);
  /// helper function to create a key basename (i.e. without the spin label)
  /// that is guaranteed for any spin to not be found in this OrbitalSpaceRegistry
  std::string new_unique_key(const Ref<OrbitalSpaceRegistry>& oreg);

  /// registry of globally-known OrbitalSpace objects that describe AO basis spaces
  typedef Registry<Ref<GaussianBasisSet>, Ref<OrbitalSpace>,
                   detail::NonsingletonCreationPolicy,
                   std::equal_to< Ref<GaussianBasisSet> >,
                   RefObjectEqual<OrbitalSpace> > AOSpaceRegistry;

  ////////////////////////////////

  /** This is an OrbitalSpace ordered according to the Order type.
      Order defines strict weak ordering for Orbital objects
      in this space and the blocking scheme.
   */
  template <typename Order>
  class OrderedOrbitalSpace : public OrbitalSpace {
    public:
      OrderedOrbitalSpace(const std::string& id, const std::string& name,
                          const Ref<GaussianBasisSet>& basis,
                          const Ref<Integral>& integral,
                          const RefSCMatrix& coefs, const RefDiagSCMatrix& evals,
                          const RefDiagSCMatrix& occnums,
                          const std::vector<unsigned int>& orbsyms,
                          const Order& order);

      OrderedOrbitalSpace(StateIn&);
      void save_data_state(StateOut&);
      ~OrderedOrbitalSpace();

    private:

      typedef OrderedOrbitalSpace this_type;
      // ClassDesc object
      static ClassDesc class_desc_;
  };

  /** Same as OrderedOrbitalSpace, except for spin-orbitals
   */
  template <typename Order>
  class OrderedSpinOrbitalSpace : public OrbitalSpace {
    public:
      OrderedSpinOrbitalSpace(const std::string& id, const std::string& name,
                          const Ref<GaussianBasisSet>& basis,
                          const Ref<Integral>& integral,
                          const RefSCMatrix& coefs_alpha,
                          const RefSCMatrix& coefs_beta,
                          const RefDiagSCMatrix& evals_alpha,
                          const RefDiagSCMatrix& evals_beta,
                          const RefDiagSCMatrix& occnums_alpha,
                          const RefDiagSCMatrix& occnums_beta,
                          const std::vector<unsigned int>& orbsyms_alpha,
                          const std::vector<unsigned int>& orbsyms_beta,
                          const Order& order);

      OrderedSpinOrbitalSpace(StateIn&);
      void save_data_state(StateOut&);
      ~OrderedSpinOrbitalSpace();

    private:

      typedef OrderedSpinOrbitalSpace this_type;
      // ClassDesc object
      static ClassDesc class_desc_;
  };

  ////////////////////////////////

  /** This is an OrbitalSpace produced from an existing one by masking out some Orbitals
   */
  class MaskedOrbitalSpace : public OrbitalSpace {
    public:
      MaskedOrbitalSpace(const std::string& id, const std::string& name,
                         const Ref<OrbitalSpace>& orig_space,
                         const std::vector<bool>& mask);

      MaskedOrbitalSpace(StateIn&);
      void save_data_state(StateOut&);

  };

  ////////////////////////////////

  /** This is an OrbitalSpace produced from an existing one by getting rid of the blocking.
   * The number of blocks in the result space is always one.
   */
  class NonblockedOrbitalSpace : public OrbitalSpace {
    public:
      NonblockedOrbitalSpace(const std::string& id, const std::string& name,
                       const Ref<OrbitalSpace>& orig_space);

      NonblockedOrbitalSpace(StateIn&);
      void save_data_state(StateOut&);

  };

  ////////////////////////////////

  /** This is an OrbitalSpace describing a set of atomic orbitals.
      The resulting dimension has 1 block, with the subdimension provided by GaussianBasisSet::basisdim()
   */
  class AtomicOrbitalSpace : public OrbitalSpace {
    public:
      AtomicOrbitalSpace(const std::string& id, const std::string& name,
                         const Ref<GaussianBasisSet>& basis,
                         const Ref<Integral>& integral);

      AtomicOrbitalSpace(StateIn&);
      void save_data_state(StateOut&);

  };

  ////////////////////////////////

  /** This is a union of two OrbitalSpaces s1 and s2. s1 and s2 may be supported by different basis sets.
   */
  class OrbitalSpaceUnion : public OrbitalSpace {
    public:
      /// \param merge_blocks if set to true, will make sure that s1 and s2 have the same blocking structure and the result
      /// will have the same number of blocks. If set to false, the total number of blocks in the result will be
      /// sum of the number of blocks in each.
      OrbitalSpaceUnion(const std::string& id, const std::string& name,
                        const OrbitalSpace& s1, const OrbitalSpace& s2,
                        bool merge_blocks = true);

      OrbitalSpaceUnion(StateIn&);
      void save_data_state(StateOut&);
  };

  ////////////////////////////////

  /** This is an empty OrbitalSpace.
   */
  class EmptyOrbitalSpace : public OrbitalSpace {
    public:
      EmptyOrbitalSpace(const std::string& id, const std::string& name,
                        const Ref<GaussianBasisSet>& basis,
                        const Ref<Integral>& integral,
                        const IndexOrder& moorder = symmetry
                        );

      EmptyOrbitalSpace(StateIn&);
      void save_data_state(StateOut&);
    private:
      static ClassDesc class_desc_;
  };

  /// @}
  // end of addtogroup ChemistryElectronicStructureOneBody

}

#include <chemistry/qc/wfn/orbitalspace.timpl.h>

#endif

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:


