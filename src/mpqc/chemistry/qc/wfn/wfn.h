/*
 * wfn.h
 *
 *  Created on: Jan 28, 2017
 *      Author: evaleev
 */

#ifndef SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_
#define SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/util/misc/exenv.h"
#include "mpqc/util/keyval/keyval.h"

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

}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_ */
