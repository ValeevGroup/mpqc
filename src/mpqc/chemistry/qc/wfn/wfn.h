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

namespace mpqc {

/// Wavefunction = opaque function of atoms, only has 2 states: computed and not
/// computed.
/// TODO It needs some sort of precision tracking to facilitate reuse.
class Wavefunction : public DescribedClass {
 public:
  /**
   *  \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for the following keywords:
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | atoms | Molecule or UnitCell | none | the collection of Atoms |
   *
   */
  Wavefunction(const KeyVal& kv);
  virtual ~Wavefunction();

  virtual void obsolete() { computed_ = false; }
  bool computed() const { return computed_; }

  const std::shared_ptr<Molecule>& atoms() const { return atoms_; }

 protected:
  bool computed_ = false;

  void set_atoms(std::shared_ptr<Molecule> atoms) { atoms_ = atoms; }

 private:
  std::shared_ptr<Molecule> atoms_;
};  // class Wavefunction

class Property;

namespace lcao {

/// Wavefunction computes a wave function (or a wave function-like quantity,
/// like
/// Green's function or reduced density matrix) in a Gaussian basis.

/// \todo elaborate Wavefunction documentation
class Wavefunction : public ::mpqc::Wavefunction {
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
   * KeyVal object will be used to construct WavefunctionWorld object |
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
