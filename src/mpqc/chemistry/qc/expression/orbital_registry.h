//
// Created by Chong Peng on 2/16/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_REGISTRY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_REGISTRY_H_

#include "mpqc/chemistry/qc/expression/formula_registry.h"
#include "mpqc/chemistry/qc/expression/orbital_space.h"

namespace mpqc {
namespace lcao {
/**
 *  \brief map OrbitalIndex to Value object
 */
template <typename Value>
class OrbitalRegistry : public Registry<OrbitalIndex, Value> {
 public:
  using Key = OrbitalIndex;
  using value_type = typename Registry<Key, Value>::value_type;
  using element_type = typename Registry<Key, Value>::element_type;
  using iterator = typename Registry<Key, Value>::iterator;
  using const_iterator = typename Registry<Key, Value>::const_iterator;

  OrbitalRegistry() = default;
  OrbitalRegistry(const element_type& map) : Registry<Key, Value>(map) {}

  // for interface in OrbitalBasisRegistry
  OrbitalRegistry(const KeyVal& kv)  {}

  /// add Value that has index() function as key type
  void add(const Value& val) { this->insert(val.index(), val); }

  /// add by Key and Value
  void add(const Key& key, const Value& val) { this->insert(key, val); }
};

/**
 * Typedef of OrbitalSpaceRegistry
 * A Registry that map OrbitalIndex to OrbitalSpace class
 */
template <typename Array>
using OrbitalSpaceRegistry = OrbitalRegistry<OrbitalSpace<Array>>;

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_REGISTRY_H_
