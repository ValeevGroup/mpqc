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
<<<<<<< HEAD
namespace lcao {

class PropertyBase;

/// Wavefunction computes a wave function (or a wave function-like quantity, like
/// Green's function or reduced density matrix) in a Gaussian basis.

/// \todo elaborate Wavefunction documentation
class Wavefunction : public mpqc::Energy {
 private:
  /** Pointer to the WfnWorld
   *
   * \note No need to make this shared Wfn is just a member of the world it
   *lives in so no ownership here.
   *
   * \warning Wfn should never delete or allocate this pointer.
   *
   * \note by chong I changed this to shared pointer, for example, MP2 and HF
   *          will share the same wfn_world
   */
  std::shared_ptr<WavefunctionWorld> wfn_world_;
=======
>>>>>>> d616d50c1009a3c300a6a51d0494c9d9d8b1713f

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

<<<<<<< HEAD
  virtual void compute(PropertyBase* pb) = 0;
  virtual double value() = 0;
  virtual void obsolete() {
    mpqc::Energy::obsolete();
  };
=======
  virtual void obsolete() { computed_ = false; }
  bool computed() const { return computed_; }

  std::shared_ptr<Molecule> atoms() { return atoms_; }
  std::shared_ptr<const Molecule> atoms() const { return atoms_; }
>>>>>>> d616d50c1009a3c300a6a51d0494c9d9d8b1713f

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
