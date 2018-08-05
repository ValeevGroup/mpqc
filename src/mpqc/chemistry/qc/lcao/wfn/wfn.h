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
#include "mpqc/chemistry/qc/lcao/factory/wfn_world.h"
#include "mpqc/util/misc/observer.h"

namespace mpqc {
namespace lcao {

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
  virtual ~Wavefunction() = default;

  /// @return shared_ptr to the WavefunctionWorld object that this Wavefunction belongs to
  std::shared_ptr<WavefunctionWorld> wfn_world() const {
    return std::static_pointer_cast<WavefunctionWorld>(::mpqc::Wavefunction::wfn_world());
  }
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_
