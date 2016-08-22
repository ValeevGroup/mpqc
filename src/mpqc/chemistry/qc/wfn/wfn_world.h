/*
 * wfn_world.h
 *
 *  Created on: Aug 18, 2016
 *      Author: Drew Lewis
 */
#ifndef MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
#define MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_

// #include "../../../../../include/tiledarray.h"
#include <mpqc/util/keyval/keyval.hpp>
#include <mpqc/chemistry/molecule/molecule.h>
#include <mpqc/chemistry/qc/integrals/atomic_integral.h>

#include <memory>

namespace mpqc {
namespace qc {

class WfnWorld : public DescribedClass {
 private:
 // madness::World &world_;
  std::shared_ptr<molecule::Molecule> mol_;
 //  std::shared_ptr<integrals::AtomicIntegral<TA::TensorD, TA::SparsePolicy>>
 //      ao_ints_;

 public:
  WfnWorld(KeyVal const &kv);

  /// Return a reference to the madness world
  // madness::World &world() { return world_ }

  /// Return a reference to the molecule in the world
  molecule::Molecule const &molecule() const { return *mol_; }

  /*! Return a reference to the AtomicIntegral Library
   *
   * \note This reference can't be made const without modifying the
   * AtomicIntegral library so that certain members are mutable.
   */
  // AtomicIntegral &ao_integrals() { return *ao_ints_; }
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_WFN_WORLD_H_
