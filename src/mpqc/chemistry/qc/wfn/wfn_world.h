/*
 * wfn_world.h
 *
 *  Created on: Aug 18, 2016
 *      Author: Drew Lewis
 */
#ifndef MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
#define MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_

#include <memory>

#include <tiledarray.h>

#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/molecule/molecule.h>
#include <mpqc/chemistry/qc/basis/basis_registry.h>

namespace mpqc {
namespace qc {

/**
 * WfnWorld is a work environment for one or more collaborating Wfn.
 * It provides an execution context (madness::World), a molecule, and
 * basis and operator registries.
 */
class WavefunctionWorld : public DescribedClass {
 public:

 private:
  madness::World &world_;
  std::shared_ptr<molecule::Molecule> mol_;
  std::shared_ptr<basis::OrbitalBasisRegistry> basis_registry_;

 public:
  /**
   * \brief KeyVal constructor
   *
   * it take all keys to construct OrbitalBasisRegistry
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
  molecule::Molecule const &molecule() const { return *mol_; }

  /// Return Basis Registry
  std::shared_ptr<basis::OrbitalBasisRegistry> const basis_registry() {return basis_registry_;}

};


}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
