//
// Created by Chong Peng on 2/16/16.
//

#ifndef MPQC_ORBITAL_SPACE_H
#define MPQC_ORBITAL_SPACE_H

#include <memory>

#include "operator.h"
#include <mpqc/chemistry/qc/basis/basis.h>
#include <mpqc/chemistry/qc/expression/orbital_index.h>

namespace mpqc {

/**
 *  \brief Class that represent a set of LCAO
 */

template <typename Array>
class OrbitalSpace {
 public:
  OrbitalSpace() = default;

  /**
   * Constructor
   *
   *  @param mo_index  OrbitalIndex that represent mo space
   *  @param ao_index  OrbitalIndex that represent ao space
   *  @param tarray    a TiledArray::DistArray type
   */

  OrbitalSpace(const OrbitalIndex& mo_index, const OrbitalIndex& ao_index,
               const Array& tarray)
      : mo_index_(mo_index), ao_index_(ao_index), coefs_(tarray) {}

  ~OrbitalSpace() = default;

  /// @return OrbitalIndex object for this space
  const OrbitalIndex& mo_key() const { return mo_index_; }

  /// @return OrbitalIndex object for the AO space supporting this space
  const OrbitalIndex& ao_key() const { return ao_index_; }

  /// @return the coefficient matrix (rows = AOs, columns = LCAOs)
  const Array& array() const { return coefs_; }

  /// @return rank of this space
  size_t rank() const {
    return coefs_.trange().elements_range().extent_data()[1];
  }

  /// @return rank of the AO space
  size_t ao_rank() const {
    return coefs_.trange().elements_range().extent_data()[0];
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
  OrbitalIndex mo_index_;
  OrbitalIndex ao_index_;
  Array coefs_;
};
}

#endif  // MPQC_ORBITAL_SPACE_H
