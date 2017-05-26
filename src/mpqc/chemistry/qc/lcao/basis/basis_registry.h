//
// Created by Chong Peng on 9/9/16.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_LCAO_BASIS_BASIS_REGISTRY_H_
#define SRC_MPQC_CHEMISTRY_QC_LCAO_BASIS_BASIS_REGISTRY_H_

#include "mpqc/chemistry/qc/lcao/basis/basis.h"
#include "mpqc/chemistry/qc/lcao/expression/orbital_registry.h"

namespace mpqc {
namespace lcao {

namespace gaussian {

/**
 * Typedef of OrbitalBasisRegistry
 * A Registry that maps OrbitalIndex to a Gaussian Basis
 */
using OrbitalBasisRegistry = OrbitalRegistry<std::shared_ptr<Basis>>;

}  // namespace gaussian

// clang-format off
/**
 * \brief The KeyVal constructor
 *
 * \param kv the KeyVal object, it will be queried for the following keywords:
 *
 * | Keyword | Type | Default| Description |
 * |---------|------|--------|-------------|
 * | basis | AtomicBasis | none | the orbital basis set|
 * | df_basis | AtomicBasis | none | an optional density-fitting basis set |
 * | aux_basis | AtomicBasis | none | an optional auxiliary basis set for the F12 theories |
 * | vir_basis | AtomicBasis | none | an optional basis set for supporting unoccupied orbitals in dual-basis methods |
 **/
// clang-format on
template <>
gaussian::OrbitalBasisRegistry::OrbitalRegistry(const KeyVal& kv);

/**
 * specialization of OrbitalRegistry::clear() only removes elements that do not
 * update their state automatically. Only gaussian::AtomicBasis objects update with
 * atomic positions.
 */
template <>
void gaussian::OrbitalBasisRegistry::clear();

}  // namespace lcao
}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_LCAO_BASIS_BASIS_REGISTRY_H_
