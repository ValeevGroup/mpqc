/*
 * wfn.h
 *
 *  Created on: Jan 28, 2017
 *      Author: evaleev
 */

#ifndef SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_
#define SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_

#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/util/core/exenv.h"
#include "mpqc/util/keyval/keyval.h"

namespace mpqc {

/// WavefunctionWorld is an environment for one or more collaborating Wavefunction objects.

/// It provides an execution context (madness::World) and a set of atoms
class WavefunctionWorld : virtual public DescribedClass {
 public:

  /**
   * \brief The KeyVal constructor
   *
   * \param kv The KeyVal object; it will be queried for the following keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c atoms | Molecule or UnitCell | none | |
   * | \c molecule | Molecule | none | This will be queried only if \c atoms is not given. This keyword is deprecated and may be removed in the future |
   **/
  explicit WavefunctionWorld(KeyVal const &kv);
  ~WavefunctionWorld() override = default;

  /// Return a reference to the MADNESS world object
  madness::World &world() { return world_; }

  /// Return a reference to the molecule in the world
  const std::shared_ptr<Molecule>& atoms() { return atoms_; }
  std::shared_ptr<const Molecule> atoms() const { return atoms_; }

 private:
  madness::World &world_;
  std::shared_ptr<Molecule> atoms_;
};

/// Wavefunction = opaque function of atoms, only has 2 states: computed and not
/// computed.
class Wavefunction : virtual public DescribedClass, public utility::Observer {
 public:
  /**
   *  \brief The KeyVal constructor
   *
   * @param[in] kv The KeyVal object; it will be queried for the following keywords:
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c "wfn_world" OR \c "..:wfn_world" OR \c "$:wfn_world" | WavefunctionWorld | none | This specifies the WavefunctionWorld object in which this object will live initially. If not found, the contents of this KeyVal object will be used to construct a new WavefunctionWorld object |
   */
  Wavefunction(const KeyVal& kv);
  virtual ~Wavefunction();

  /// @return shared_ptr to the WavefunctionWorld object that this Wavefunction belongs to
  const std::shared_ptr<WavefunctionWorld>& wfn_world() const {
    return wfn_world_;
  }

  virtual void obsolete() { computed_ = false; }
  bool computed() const { return computed_; }

  std::shared_ptr<Molecule> atoms() { return wfn_world_->atoms(); }
  std::shared_ptr<const Molecule> atoms() const { return wfn_world_->atoms(); }

  virtual void print(std::ostream& os = ExEnv::out0()) const {
    os << indent << "Wavefunction (type = " << this->class_key() << "):\n" << incindent;
    os << *atoms();
    os << decindent;
    os << std::endl;
  }

 protected:
  bool computed_ = false;

  /// sets wfn_world_ to @c wfn_world
  /// @param[in] wfn_world pointer to the new WavefunctionWorld object
  void reset_wfn_world(const std::shared_ptr<WavefunctionWorld>& wfn_world) {
    wfn_world_ = wfn_world;
  }

 private:
  std::shared_ptr<WavefunctionWorld> wfn_world_;
};  // class Wavefunction

}  // namespace mpqc

#endif /* SRC_MPQC_CHEMISTRY_QC_WFN_WFN_H_ */
