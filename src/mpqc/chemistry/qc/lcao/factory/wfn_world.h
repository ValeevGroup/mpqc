//
// Created by Eduard Valeyev on 8/5/18.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_WFN_WORLD_H
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_WFN_WORLD_H

#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/chemistry/qc/lcao/basis/basis_registry.h"

namespace mpqc {
namespace lcao {

/// specialization of ::mpqc::WavefunctionWorld to the LCAO case
class WavefunctionWorld : public ::mpqc::WavefunctionWorld {
 public:
  /**
   * \brief The KeyVal constructor
   *
   * \param kv The KeyVal object; it will be queried for all keywords of OrbitalBasisRegistry, and the following keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c atoms | Molecule or UnitCell | none | |
   * | \c molecule | Molecule | none | This will be queried only if \c atoms is not given. This keyword is deprecated and may be removed in the future |
   **/
  explicit WavefunctionWorld(KeyVal const &kv);
  ~WavefunctionWorld() override = default;

  /// Return a reference to the basis registry
  const std::shared_ptr<gaussian::OrbitalBasisRegistry> &basis_registry() { return basis_registry_; }

 private:
  std::shared_ptr<gaussian::OrbitalBasisRegistry> basis_registry_;
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_LCAO_FACTORY_WFN_WORLD_H
