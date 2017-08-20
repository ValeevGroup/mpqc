/*
 * wfn_world.h
 *
 *  Created on: Aug 18, 2016
 *      Author: Drew Lewis
 */
#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_

#include <memory>

#include <tiledarray.h>

#include "mpqc/util/keyval/keyval.h"
#include "mpqc/chemistry/molecule/molecule.h"
#include "mpqc/chemistry/qc/lcao/basis/basis_registry.h"

namespace mpqc {
namespace lcao {

/// WavefunctionWorld is an environment for one or more collaborating Wavefunction objects.

/// It provides an execution context (madness::World), a molecule, and
/// basis and operator registries.
class WavefunctionWorld : virtual public DescribedClass {
 public:

 private:
  madness::World &world_;
  std::shared_ptr<Molecule> atoms_;
  std::shared_ptr<gaussian::OrbitalBasisRegistry> basis_registry_;

 public:
  /**
   * \brief The KeyVal constructor
   *
   * \param kv The KeyVal object will be queried for all keywords of OrbitalBasisRegistry, and the following keywords:
   *
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | \c atoms | Molecule or UnitCell | none | |
   * | \c molecule | Molecule | none | This will be queried only if \c atoms is not given. This keyword is deprecated and may be removed in the future |
   **/
  explicit WavefunctionWorld(KeyVal const &kv);
  ~WavefunctionWorld() override = default;

  /// Return a reference to the madness world
  madness::World &world() { return world_; }

  /// Return a reference to the molecule in the world
  const std::shared_ptr<Molecule>& atoms() const { return atoms_; }

  /// Return a reference to Basis Registry
  const std::shared_ptr<gaussian::OrbitalBasisRegistry>& basis_registry() { return basis_registry_; }

};


}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
