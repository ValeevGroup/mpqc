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
#include "mpqc/chemistry/qc/lcao/wfn/wfn_world.h"
#include "mpqc/util/misc/observer.h"

namespace mpqc {
namespace lcao {

/// Wavefunction computes a wave function (or a wave function-like quantity,
/// like
/// Green's function or reduced density matrix) in a Gaussian basis.

/// \todo elaborate Wavefunction documentation
class Wavefunction : public ::mpqc::Wavefunction {
 private:
  /** Pointer to the WfnWorld */
  std::shared_ptr<WavefunctionWorld> wfn_world_;

 public:
  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for the following keywords:
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "wfn_world" OR \c "..:wfn_world" OR \c "$:wfn_world" |
   * WavefunctionWorld | none | This specifies the WavefunctionWorld object in
   * which this object will live initially. If not found, the contents of this
   * KeyVal object will be used to construct a new WavefunctionWorld object |
   *
   */
  Wavefunction(const KeyVal& kv);
  virtual ~Wavefunction();

  virtual void obsolete() { ::mpqc::Wavefunction::obsolete(); };

  const std::shared_ptr<WavefunctionWorld>& wfn_world() const {
    return wfn_world_;
  }
};

}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_
