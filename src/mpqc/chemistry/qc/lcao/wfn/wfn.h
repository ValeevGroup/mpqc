/*
 * wfn.h
 *
 *  Created on: Apr 27, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_

#include <functional>
#include <memory>

#include "mpqc/chemistry/qc/wfn/wfn.h"
#include "mpqc/chemistry/qc/lcao/basis/basis_registry.h"
#include "mpqc/util/misc/observer.h"

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
  const std::shared_ptr<gaussian::OrbitalBasisRegistry>& basis_registry() { return basis_registry_; }

 private:
  std::shared_ptr<gaussian::OrbitalBasisRegistry> basis_registry_;
};


/// Wavefunction computes a wave function (or a wave function-like quantity,
/// like Green's function or reduced density matrix) in a Gaussian basis.

/// \todo elaborate Wavefunction documentation
class Wavefunction : public ::mpqc::Wavefunction {
 public:
  // clang-format off
  /**
   *  \brief The KeyVal constructor
   *
   * @param[in] kv the KeyVal object will be queried for all keywords of ::mpqc::Wavefunction
   */
  // clang-format on
  Wavefunction(const KeyVal& kv);
  virtual ~Wavefunction();

  /// @return shared_ptr to the WavefunctionWorld object that this Wavefunction belongs to
  std::shared_ptr<WavefunctionWorld> wfn_world() const {
    return std::static_pointer_cast<WavefunctionWorld>(::mpqc::Wavefunction::wfn_world());
  }
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_
