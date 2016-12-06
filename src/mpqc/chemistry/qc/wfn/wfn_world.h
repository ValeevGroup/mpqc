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
#include "mpqc/chemistry/qc/basis/basis_registry.h"

namespace mpqc {
namespace qc {

/// WavefunctionWorld is an environment for one or more collaborating Wavefunction objects.

/// It provides an execution context (madness::World), a molecule, and
/// basis and operator registries.
class WavefunctionWorld : public DescribedClass {
 public:

 private:
  madness::World &world_;
  std::shared_ptr<Molecule> mol_;
  std::shared_ptr<basis::OrbitalBasisRegistry> basis_registry_;

 public:
  /**
   * \brief The KeyVal constructor
   *
   * The KeyVal object will be queried for all keywords of OrbitalBasisRegistry, and the following keywords:
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | molecule | Molecule | none | |
   **/
  WavefunctionWorld(KeyVal const &kv);
  ~WavefunctionWorld();

  /// Return a reference to the madness world
  madness::World &world() { return world_; }

  /// Return a reference to the molecule in the world
  Molecule const &molecule() const { return *mol_; }

  /// Return Basis Registry
  std::shared_ptr<basis::OrbitalBasisRegistry> const basis_registry() {return basis_registry_;}

};


}  // namespace qc
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
