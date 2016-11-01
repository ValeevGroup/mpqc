/*
 * ao_wfn.h
 *
 *  Created on: Aug 17, 2016
 *      Author: Drew Lewis
 */

#ifndef MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
#define MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_

#include <mpqc/chemistry/qc/wfn/wfn.h>
#include <mpqc/chemistry/qc/integrals/atomic_integral.h>
#include <mpqc/chemistry/qc/integrals/direct_atomic_integral.h>

namespace mpqc {
namespace qc {

template<typename Tile, typename Policy>
class AOWavefunction : public Wavefunction {
 public:
  using AOIntegral = integrals::AtomicIntegral<Tile, Policy>;
  using DirectAOIntegral = integrals::DirectAtomicIntegral<Tile, Policy>;
  using ArrayType = typename AOIntegral::TArray;

  AOWavefunction(const KeyVal &kv) : Wavefunction(kv)
  {
    ao_ints_ = integrals::detail::construct_atomic_integral<Tile, Policy>(kv);
    ao_ints_->set_orbital_basis_registry(this->wfn_world()->basis_registry());

    direct_ao_ints_ = integrals::detail::construct_direct_atomic_integral<Tile,Policy>(kv);
    direct_ao_ints_->set_orbital_basis_registry(this->wfn_world()->basis_registry());

  }
  virtual ~AOWavefunction() = default;

  void compute(PropertyBase *pb) override {
    throw std::logic_error("Not Implemented!");
  }

  /// obsolete, purge the registry in AOIntegral and DirectAOInetgral
  void obsolete() override {
    ao_integrals().registry().purge(wfn_world()->world());
    direct_ao_integrals().registry().purge(wfn_world()->world());
    Wavefunction::obsolete();
  }

  double value() override {
    return 0.0;
  };

  /*! Return a reference to the AtomicIntegral Library
   *
   * \note This reference can't be made const without modifying the
   * AtomicIntegral library so that certain members are mutable.
   */
  AOIntegral &ao_integrals() { return *ao_ints_; }

  /// return a reference to DirectAtomicIntegral
  DirectAOIntegral& direct_ao_integrals() { return *direct_ao_ints_;}
private:
  std::shared_ptr<AOIntegral> ao_ints_;
  std::shared_ptr<DirectAOIntegral> direct_ao_ints_;
};

}  // namespace qc
}  // namespace mpqc

#endif  // MPQC_CHEMISTRY_QC_WFN_AO_WFN_H_
