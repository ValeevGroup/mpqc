//
// Created by Chong Peng on 9/9/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_REGISTRY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_REGISTRY_H_

#include "mpqc/chemistry/qc/basis/basis.h"
#include "mpqc/chemistry/qc/expression/orbital_registry.h"

namespace mpqc {
namespace lcao {

namespace gaussian {

/**
 * Typedef of OrbitalBasisRegistry
 * A Registry that maps OrbitalIndex to a Gaussian Basis
 */
using OrbitalBasisRegistry = OrbitalRegistry<std::shared_ptr<Basis>>;

}  // namespace gaussian

/**
 * KeyVal Constructor for OrbitalBasisRegistry
 */
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

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_REGISTRY_H_
