//
// Created by Chong Peng on 2/16/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_REGISTRY_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_REGISTRY_H_

#include <memory>

#include "mpqc/chemistry/qc/lcao/expression/formula_registry.h"
#include "mpqc/chemistry/qc/lcao/expression/orbital_space.h"
#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"

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
  using container_type = typename Registry<Key, Value>::container_type;
  using iterator = typename Registry<Key, Value>::iterator;
  using const_iterator = typename Registry<Key, Value>::const_iterator;

  OrbitalRegistry() = default;
  OrbitalRegistry(const container_type& map) : Registry<Key, Value>(map) {}

  /// The KeyVal ctor must be provided as a specialization
  OrbitalRegistry(const KeyVal& kv) {}

  /// add Value that has index() function as key type
  void add(const Value& val) { this->insert(val.index(), val); }

  /// add by Key and Value
  /// @note \c val is copied
  void add(const Key& key, const Value& val) { this->insert(key, val); }

  /// specialize, if needed
  void clear() override { Registry<OrbitalIndex, Value>::clear(); }
};

/**
 * @brief OrbitalSpaceRegistry is an OrbitalRegistry that maps OrbitalIndex to
 * OrbitalSpace class.
 *
 * OrbitalSpaceRegistry contains additional functionality (e.g.
 * trange1_engine())
 * specific to this type of registry.
 *
 */
template <typename Array>
class OrbitalSpaceRegistry : public OrbitalRegistry<OrbitalSpace<Array>> {
 public:
  using TRange1Engine = ::mpqc::utility::TRange1Engine;

  /// accessor to the TRange1Engine object that provides shortcuts to the
  /// OrbitalSpace info already in the registry
  /// @warning this method will likely go away, the users are encouraged
  ///          to access the OrbitalSpace info directly, i.e. instead of
  ///          ```trange1_engine()->get_occ_tr()``` do
  ///          ```find("i").trange()```
  const std::shared_ptr<const TRange1Engine>& trange1_engine() const {
    if (!trange1_engine_) {
      try {
        auto nobs = this->retrieve("p").rank();
        const auto& occ_space = this->retrieve("m");
        auto ndocc = occ_space.rank();
        auto ndocc_act = this->retrieve("i").rank();
        auto n_frozen_core = ndocc - ndocc_act;
        const auto& uocc_space = this->retrieve("a");

        // can't always deduce user's occ and unocc blk sizes, just grab dimensions of the first tile in each
        const auto& occ_tile0 = occ_space.trange().tile(0);
        auto occ_blksize = occ_tile0.second - occ_tile0.first;
        const auto& uocc_tile0 = uocc_space.trange().tile(0);
        auto uocc_blksize = uocc_tile0.second - uocc_tile0.first;

        trange1_engine_ = std::make_shared<const TRange1Engine>(
            ndocc, nobs, occ_blksize, uocc_blksize, n_frozen_core);
      } catch (ProgrammingError&) {
        throw ProgrammingError(
            "OrbitalSpaceRegistry::trange1_engine failed -- likely some "
            "OrbitalSpace objects are missing",
            __FILE__, __LINE__);
      }
    }
    return trange1_engine_;
  }

 private:
  mutable std::shared_ptr<const TRange1Engine> trange1_engine_ = nullptr;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_EXPRESSION_ORBITAL_REGISTRY_H_
