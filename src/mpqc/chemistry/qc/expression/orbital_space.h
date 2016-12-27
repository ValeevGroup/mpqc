//
// Created by Chong Peng on 2/16/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_SPACE_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_SPACE_H_

#include <memory>

#include "mpqc/chemistry/qc/basis/basis.h"
#include "mpqc/chemistry/qc/expression/orbital_index.h"
#include "operator.h"

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
  OrbitalSpace() = default;

  /**
   * Constructor
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
        coefs_(tarray) {}

  ~OrbitalSpace() = default;

  /// @return the base OrbitalIndex object for this space
  const OrbitalIndex& index() const { return index_; }

  /// @return the base OrbitalIndex object for the AO space supporting this space
  const OrbitalIndex& ao_index() const { return ao_index_; }

  /// @return a const reference to the coefficient matrix (an \c Array object,
  ///         whose rows are AOs, and columns are LCAOs).
  const Array& coefs() const { return coefs_; }

  /// @return rank of this space
  size_t rank() const {
    return coefs_.trange().elements_range().extent_data()[1];
  }

  /// @return rank of the AO space
  size_t ao_rank() const {
    return coefs_.trange().elements_range().extent_data()[0];
  }

  /// @return the \c std::string object that contains a brief description of this space
  const std::string& descriptor() const {
    return descriptor_;
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
};  // class mpqc::lcao::OrbitalSpace

/// @}

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_SPACE_H_
