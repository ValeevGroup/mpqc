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

#include "mpqc/chemistry/qc/wfn/wfn_world.h"
#include "mpqc/util/misc/observer.h"

namespace mpqc {

/// Wavefunction = opaque function of atoms, only has 2 states: computed and not
/// computed.
/// TODO It needs some sort of precision tracking to facilitate reuse.
class Wavefunction : virtual public DescribedClass, public utility::Observer {
 public:
  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | atoms | Molecule or UnitCell | none | the collection of Atoms |
   *
   * For compatibility, if keyword \c atoms not found keyword \c molecule will
   * be queried.
   *
   */
  Wavefunction(const KeyVal& kv);
  virtual ~Wavefunction();

  virtual void obsolete() { computed_ = false; }
  bool computed() const { return computed_; }

  std::shared_ptr<Molecule> atoms() { return atoms_; }
  std::shared_ptr<const Molecule> atoms() const { return atoms_; }

  virtual void print(std::ostream& os = ExEnv::out0()) const {
    os << indent << "Wavefunction (type = " << this->class_key() << "):\n" << incindent;
    os << *atoms_;
    os << decindent;
    os << std::endl;
  }

 protected:
  bool computed_ = false;

 private:
  std::shared_ptr<Molecule> atoms_;
};  // class Wavefunction

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
   * | KeyWord | Type | Default| Description |
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
