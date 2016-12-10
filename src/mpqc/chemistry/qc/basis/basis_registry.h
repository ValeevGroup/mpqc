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
using OrbitalBasisRegistry = OrbitalRegistry<gaussian::Basis>;

}  // namespace gaussian

/**
 * KeyVal Constructor for OrbitalBasisRegistry
 */

template <>
OrbitalRegistry<gaussian::Basis>::OrbitalRegistry(const KeyVal& kv);

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_BASIS_BASIS_REGISTRY_H_
