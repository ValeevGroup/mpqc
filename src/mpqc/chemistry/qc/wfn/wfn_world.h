/*
 * wfn_world.h
 *
 *  Created on: Aug 18, 2016
 *      Author: Drew Lewis
 */
#ifndef MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
#define MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_

#include "../../../../../include/tiledarray.h"
#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/molecule/molecule.h>
#include <mpqc/chemistry/qc/integrals/atomic_integral.h>

#include <memory>

namespace mpqc {
namespace qc {

/**
 * WfnWorld is a work environment for one or more collaborating Wfn.
 * It provides an execution context (madness::World), a molecule, and
 * basis and operator registries.
 */
class WfnWorld : public DescribedClass {
 public:
  using AOIntegral = integrals::AtomicIntegral<TA::TensorD, TA::SparsePolicy>;
  using ArrayType = AOIntegral::TArray;

 private:
  madness::World &world_;
  AOIntegral ao_ints_;
  std::shared_ptr<molecule::Molecule> mol_;
  std::shared_ptr<basis::OrbitalBasisRegistry> basis_registry_;

 public:
  WfnWorld(KeyVal const &kv);
  ~WfnWorld();

  /// Return a reference to the madness world
  madness::World &world() { return world_; }

  /// Return a reference to the molecule in the world
  molecule::Molecule const &molecule() const { return *mol_; }

  /*! Return a reference to the AtomicIntegral Library
   *
   * \note This reference can't be made const without modifying the
   * AtomicIntegral library so that certain members are mutable.
   */
  AOIntegral &ao_integrals() { return ao_ints_; }
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
