//
// Created by Chong Peng on 3/2/17.
//

#ifndef SRC_MPQC_CHEMISTRY_QC_PROPERTIES_EXCITATION_ENERGY_H_
#define SRC_MPQC_CHEMISTRY_QC_PROPERTIES_EXCITATION_ENERGY_H_

#include "mpqc/chemistry/qc/properties/property.h"

namespace mpqc {

/**
 * Property that computes Excitation Energy
 */

class ExcitationEnergy : public WavefunctionProperty<std::vector<double>> {
 public:
  using typename WavefunctionProperty<std::vector<double>>::function_base_type;

  class Provider : public math::FunctionVisitorBase<function_base_type> {
   public:
    virtual bool can_evaluate(ExcitationEnergy* ex_energy) = 0;
    /// evaluate Excitation Energy and use set_value to assign the values to \c
    /// ex_energy
    virtual void evaluate(ExcitationEnergy* ex_energy) = 0;
  };

  // clang-format off
  /**
   * @brief The KeyVal constructor
   * @param kv the KeyVal object, it will be queried for all
   *        keywords of the WavefunctionProperty class, as well as the following
   * keywords:
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | n_roots | unsigned int | 3 | number of states to compute |
   * | n_guess | unsigned int | n_roots | number of guess states used to compute |
   * | singlets | bool | true | compute singlet excitation energy, only apply to closed-shell system|
   * | triplets | bool | false | compute triplet excitation energy, only apply to closed-shell system|
   *
   * @note This constructor overrides the default target precision to 1e-5 \f$ \sim 0.3~{\rm meV} \f$.
   */
  // clang-format on
  ExcitationEnergy(const KeyVal& kv);

  ~ExcitationEnergy() = default;

  /// @return number of roots
  std::size_t n_roots() const;

  /// @return number of guess
  std::size_t n_guess() const;

  /// set the number of roots
  void set_n_roots(unsigned int n_roots_);

  /// @return if do singlets
  bool singlets() const;

  /// @return if do triplets
  bool triplets() const;

 private:
  void do_evaluate() override;

  std::size_t n_roots_;
  std::size_t n_guess_;
  bool singlets_;
  bool triplets_;
};

}  // namespace mpqc

#endif  // SRC_MPQC_CHEMISTRY_QC_PROPERTIES_EXCITATION_ENERGY_H_
