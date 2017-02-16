#ifndef SRC_MPQC_CHEMISTRY_QC_PROPERTIES_GFPOLE_H_
#define SRC_MPQC_CHEMISTRY_QC_PROPERTIES_GFPOLE_H_

#include "mpqc/chemistry/qc/properties/property.h"

namespace mpqc {

/** Computes the real poles of a Green's function
 */
class GFRealPole : public WavefunctionProperty<double> {
public:
  using typename WavefunctionProperty<double>::function_base_type;

  class Provider : public math::FunctionVisitorBase<function_base_type> {
  public:
    virtual bool can_evaluate(GFRealPole* pole) = 0;
    /// Provider::evaluate computes the taylor expansion of the real pole
    /// with respect to the molecular coordinates
    /// and uses set_value to assign the values to \c energy
    virtual void evaluate(GFRealPole* pole) = 0;
  };

  // clang-format off
  /**
   * @brief The KeyVal constructor
   * @param kv the KeyVal object, it will be queried for all
   *        keywords of the WavefunctionProperty class, as well as the following keywords:
   * | Keyword | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | target | int | -1 | target pole index, defined relative to the Fermi level; e.g., -1 denotes the first pole below the Fermi level, +1 denotes the first one above. In practice the Wavefunction class may use this as a initial guess rather than try to guarantee the index of the computed pole; refer to the documentation of the appropriate Wavefunction class|
   *
   * @note This constructor overrides the default target precision to 1e-4 \f$ \sim 3~{\rm meV} \f$.
   */
  // clang-format on
  explicit GFRealPole(const KeyVal& kv);

  /// @return the target pole index (see the KeyVal ctor of GFRealPole for the explanation).
  int target() const;

private:
  int target_;
  void do_evaluate() override;
};

} // namespace mpqc

#endif //  MPQC4_SRC_MPQC_CHEMISTRY_QC_PROPERTIES_GFPOLE_H_
