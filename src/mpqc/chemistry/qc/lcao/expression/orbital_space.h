//
// Created by Chong Peng on 2/16/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_SPACE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_SPACE_H_

#include <memory>

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/expression/operator.h"
#include "mpqc/chemistry/qc/lcao/expression/orbital_index.h"
#include "mpqc/math/groups/group.h"
#include "mpqc/util/misc/exception.h"

namespace mpqc {
namespace lcao {

/// @addtogroup ChemistryESLCAO
/// @{

/**
 *  \brief OrbitalSpace represents a set of LCAO
 *
 *  \tparam Array the type that represents the expansion coefficients
 *
 */
template <typename Array>
class OrbitalSpace {
 public:
  using Group = ::mpqc::math::Group;

  /**
    *  every class that can provide OrbitalSpace (e.g. RHF) will publicly
    *  inherit from OrbitalSpace::Provider
    *
    *  @sa Provides
    */
  class Provider {
   public:
    static constexpr const std::size_t default_lcao_blocksize = 20;

    /// @return true if \c OrbitalSpace can be computed.
    virtual bool can_evaluate(OrbitalSpace* ospace = nullptr) = 0;

    /// @return true if \c OrbitalSpace is available without additional
    /// computation.
    virtual bool is_available(OrbitalSpace* ospace = nullptr) = 0;

    /// @brief computes an OrbitalSpace and assigns to \c *ospace

    /// @param ospace pointer to the OrbitalSpace object that will contain the
    /// result
    /// @param target_precision optional precision parameter, its use
    ///        depends on the type of Provider (some Providers will ignore it,
    ///        some will use as the guidance on the precision of the energy,
    ///        etc.)
    /// @param target_lcao_blocksize optional parameter to determine the average
    ///        blocksize of the LCAO dimension (blocking of the AO dimension is
    ///        determined by the clustering of the Molecule). May be ignored,
    ///        hence it is recommended to reblock the coefficients as needed.
    virtual void evaluate(
        OrbitalSpace* ospace, double target_precision = 0,
        std::size_t target_lcao_blocksize = default_lcao_blocksize) = 0;
  };

  /// makes a default-initialized array
  OrbitalSpace() = default;

  /**
   * Constructs an OrbitalSpace for a system without symmetry.
   *
   *  @param idx     an OrbitalIndex that represents this space; it is converted
   *                 to the base index
   *  @param ao_idx  an OrbitalIndex that represents the AO space supporting
   *                 this space; it is converted to the base index
   *  @param tarray  a TiledArray::DistArray type
   */
  OrbitalSpace(const OrbitalIndex& idx, const OrbitalIndex& ao_idx,
               const Array& tarray)
      : index_(make_base_index(idx)),
        ao_index_(make_base_index(ao_idx)),
        coefs_(tarray),
        group_(std::make_shared<math::groups::Z1>()),
        irrep_indices_(coefs_.trange().tiles_range().extent_data()[1], 0) {
    if (coefs_.elements_range().rank() != 2)
      throw ProgrammingError(
          "OrbitalSpace ctor expects a matrix of coefficients", __FILE__,
          __LINE__);
  }

  /**
   * Constructs an OrbitalSpace for a system with symmetry. LCAOs are
   * grouped according to the irreducible irreps of the symmetry group.
   *
   *  @param idx     an OrbitalIndex that represents this space; it is converted
   *                 to the base index
   *  @param ao_idx  an OrbitalIndex that represents the AO space supporting
   *                 this space; it is converted to the base index
   *  @param tarray  a TiledArray::DistArray type
   *  @param group   the symmetry group object
   *  @param irrep_indices the irrep (ordinal) indices for each block of the
   *                 column dimension of tarray; if omitted, assume that the
   *                 blocks are spread uniformly among the irreps,
   *                 in the order of increasing irrep indices.
   */
  OrbitalSpace(const OrbitalIndex& idx, const OrbitalIndex& ao_idx,
               const Array& tarray, std::shared_ptr<const Group> group,
               std::vector<Group::ordinal_type> irrep_indices =
                   std::vector<Group::ordinal_type>())
      : index_(make_base_index(idx)),
        ao_index_(make_base_index(ao_idx)),
        coefs_(tarray),
        group_(group),
        irrep_indices_(std::move(irrep_indices)) {
    // expect a matrix of coefficients
    if (coefs_.elements_range().rank() != 2)
      throw ProgrammingError(
          "OrbitalSpace ctor expects a matrix of coefficients", __FILE__,
          __LINE__);

    const auto nirreps = group_->order();
    const auto ntiles_col = coefs_.trange().tiles_range().extent_data()[1];

    // if not given irrep indices, assume the blocks are evenly divisible among
    // irreps
    if (irrep_indices_.empty()) {
      if (ntiles_col % nirreps != 0)
        throw ProgrammingError("OrbitalSpace: cannot deduce irrep indices",
                               __FILE__, __LINE__);

      irrep_indices_.resize(ntiles_col);
      auto ntiles_per_irrep = ntiles_col / nirreps;
      for (auto irrep = 0, tile = 0; irrep != nirreps; ++irrep)
        for (auto t = 0; t != ntiles_per_irrep; ++t, ++tile)
          irrep_indices_[tile] = irrep;
    } else if (irrep_indices_.size() != ntiles_col)
      throw ProgrammingError("OrbitalSpace: # of irrep indices != # of tiles",
                             __FILE__, __LINE__);
  }

  /// constructs this OrbitalSpace using a visiting Provider
  /// @tparam Visitor a class derived from OrbitalSpace::Provider
  template <typename Visitor>
  void evaluate(std::shared_ptr<Visitor> visitor) {
    if (auto provider = std::dynamic_pointer_cast<Provider>(visitor)) {
      provider->evaluate(*this);
    } else
      throw ProgrammingError(
          "OrbitalSpace::evaluate: visitor does not provide OrbitalSpace "
          "objects",
          __FILE__, __LINE__);
  }

  virtual ~OrbitalSpace() = default;

  /// @return the base OrbitalIndex object for this space
  const OrbitalIndex& index() const { return index_; }

  /// @return the base OrbitalIndex object for the AO space supporting this
  /// space
  const OrbitalIndex& ao_index() const { return ao_index_; }

  /// @return a const reference to the coefficient matrix (an \c Array object,
  ///         whose rows are AOs, and columns are LCAOs).
  const Array& coefs() const { return coefs_; }

  /// @return rank of this space
  size_t rank() const { return coefs_.elements_range().extent_data()[1]; }

  /// @return rank of the AO space that supports this
  size_t ao_rank() const { return coefs_.elements_range().extent_data()[0]; }

  /// @return the TiledRange1 object for this space
  const TA::TiledRange1& trange() const { return coefs_.trange().data()[1]; }

  /// @return the tiled range object for the AO space that supports this
  const TA::TiledRange1& ao_trange() const { return coefs_.trange().data()[0]; }

  /// @return the \c std::string object that contains a brief description of
  /// this space
  const std::string& descriptor() const { return descriptor_; }

  /// @return the symmetry group
  std::shared_ptr<const Group> group() const { return group_; }

  /// @return the vector of irrep (ordinal) indices for each block of the column
  /// dimension of coefs
  const std::vector<Group::ordinal_type>& irrep_indices() const {
    return irrep_indices_;
  }

  /// interface to TA::Array () function
  TA::expressions::TsrExpr<Array, true> operator()(const std::string& vars) {
    return coefs_(vars);
  };

  /// interface to TA::Array () function
  TA::expressions::TsrExpr<const Array, true> operator()(
      const std::string& vars) const {
    return coefs_(vars);
  };

 private:
  OrbitalIndex index_;
  OrbitalIndex ao_index_;
  std::string descriptor_;
  Array coefs_;
  std::shared_ptr<const Group> group_;
  std::vector<Group::ordinal_type> irrep_indices_;
};  // class OrbitalSpace

/// @brief an OrbitalSpace where each orbital in addition to irrep
///        has additional attributes

/// @tparam Array the type that represents the expansion coefficients
/// @tparam Attribute the type of decoration
/// @tparam AttributeTag useful to distinguish instances of
/// DecoratedOrbitalSpace with built-in
///         Attributes
template <typename Array, typename OrbitalAttribute,
          typename AttributeTag = std::true_type>
class DecoratedOrbitalSpace : virtual public OrbitalSpace<Array> {
 public:
  /**
    *  every class that can provide DecoratedOrbitalSpace will publicly
    *  inherit from DecoratedOrbitalSpace::Provider
    *
    *  @sa Provides
    */
  class Provider {
   public:
    static constexpr const std::size_t default_lcao_blocksize =
        OrbitalSpace<Array>::Provider::default_lcao_blocksize;

    /// @return true if \c DecoratedOrbitalSpace can be computed.
    virtual bool can_evaluate(DecoratedOrbitalSpace* space = nullptr) = 0;

    /// @return true if \c DecoratedOrbitalSpace is available without additional
    /// computation.
    virtual bool is_available(DecoratedOrbitalSpace* ospace = nullptr) = 0;

    /// @brief computes an DecoratedOrbitalSpace and assigns to \c *space

    /// @param space pointer to the DecoratedOrbitalSpace object that will
    ///        contain the result
    /// @param target_precision optional precision parameter, its use
    ///        depends on the type of Provider (some Providers will ignore it,
    ///        some will use as the guidance on the precision of the energy,
    ///        etc.)
    /// @param target_lcao_blocksize optional parameter to determine the average
    ///        blocksize of the LCAO dimension (blocking of the AO dimension is
    ///        determined by the clustering of the Molecule). May be ignored,
    ///        hence it is recommended to reblock the coefficients as needed.
    virtual void evaluate(
        DecoratedOrbitalSpace* space, double target_precision = 0,
        std::size_t target_lcao_blocksize = default_lcao_blocksize) = 0;
  };

  DecoratedOrbitalSpace() = default;

  /**
   * Constructor
   *
   *  @param idx     an OrbitalIndex that represents this space; it is converted
   *                 to the base index
   *  @param ao_idx  an OrbitalIndex that represents the AO space supporting
   *                 this space; it is converted to the base index
   *  @param tarray  a TiledArray::DistArray type
   *  @param attributes a vector of OrbitalAttribute objects
   *  @param particle_type_tag an optional particle type tag (default is 0)
   */
  DecoratedOrbitalSpace(const OrbitalIndex& idx, const OrbitalIndex& ao_idx,
                        const Array& tarray,
                        const std::vector<OrbitalAttribute>& attributes)
      : OrbitalSpace<Array>(idx, ao_idx, tarray), attributes_(attributes) {}

  ~DecoratedOrbitalSpace() = default;

  /// constructs this DecoratedOrbitalSpace using a visiting Provider
  /// @tparam Visitor a class derived from OrbitalSpace::Provider
  template <typename Visitor>
  void evaluate(std::shared_ptr<Visitor> visitor) {
    if (auto provider = std::dynamic_pointer_cast<Provider>(visitor)) {
      provider->evaluate(*this);
    } else
      throw ProgrammingError(
          "DecoratedOrbitalSpace::evaluate: visitor does not provide "
          "DecoratedOrbitalSpace "
          "objects",
          __FILE__, __LINE__);
  }

  const std::vector<OrbitalAttribute>& attributes() const {
    return attributes_;
  }

 private:
  std::vector<OrbitalAttribute> attributes_;

};  // class DecoratedOrbitalSpace

namespace detail {
template <std::size_t tag>
struct AttributeTag {};
using CanonicalAttributeTag = AttributeTag<0>;
using PopulatedAttributeTag = AttributeTag<1>;
}  // namespace detail

/// Canonical orbitals are decorated by energies, in each irrep block appear in
/// the order of increasing enegry
template <typename Array>
using CanonicalOrbitalSpace =
    DecoratedOrbitalSpace<Array, double, detail::CanonicalAttributeTag>;

/// Populated orbitals are decorated by occupancies, no order assumed
template <typename Array>
using PopulatedOrbitalSpace =
    DecoratedOrbitalSpace<Array, double, detail::PopulatedAttributeTag>;

/// @}

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_SPACE_H_
