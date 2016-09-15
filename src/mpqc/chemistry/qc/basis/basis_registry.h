//
// Created by Chong Peng on 9/9/16.
//

#ifndef MPQC_BASIS_REGISTRY_H
#define MPQC_BASIS_REGISTRY_H

#include <mpqc/chemistry/qc/basis/basis.h>
#include <mpqc/chemistry/qc/expression/orbital_registry.h>

namespace mpqc {
namespace basis {

/**
 * Typedef of OrbitalBasisRegistry
 * A Registry that map OrbitalIndex to Basis
 */
using OrbitalBasisRegistry = OrbitalRegistry<mpqc::basis::Basis>;

}  // end of namespace basis

/**
 * KeyVal Constructor for OrbitalBasisRegistry
 */

template <>
mpqc::OrbitalRegistry<mpqc::basis::Basis>::OrbitalRegistry(const KeyVal& kv);

}  // end of namespace mpqc

#endif  // MPQC_BASIS_REGISTRY_H
