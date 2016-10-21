/*
 * ao_wfn.h
 *
 *  Created on: Aug 17, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
#define MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_

#include <mpqc/chemistry/qc/wfn/wfn.h>
#include <mpqc/util/keyval/keyval.hpp>

namespace mpqc {
namespace qc {

template<typename Tile>
class AOWavefunction : public Wavefunction {
 public:
  using AOIntegral = integrals::AtomicIntegral<Tile, TA::SparsePolicy>;
  using ArrayType = typename AOIntegral::TArray;

  AOWavefunction(const KeyVal &kv) : Wavefunction(kv)
  {
    ao_ints_ = integrals::detail::construct_atomic_integral<TA::TensorD, TA::SparsePolicy>(kv);
    ao_ints_->set_orbital_basis_registry(this->wfn_world()->basis_registry());

  }
  ~AOWavefunction() = default;

  void compute(PropertyBase *pb) override {
    throw std::logic_error("Not Implemented!");
  }

  void obsolete() override {
    ao_integrals().registry().purge(wfn_world()->world());
  }

  /*! Return a reference to the AtomicIntegral Library
   *
   * \note This reference can't be made const without modifying the
   * AtomicIntegral library so that certain members are mutable.
   */
  AOIntegral &ao_integrals() { return *ao_ints_; }
private:
  std::shared_ptr<AOIntegral> ao_ints_;
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
