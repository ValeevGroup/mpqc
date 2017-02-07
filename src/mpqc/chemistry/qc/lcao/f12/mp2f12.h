//
// Created by Chong Peng on 3/31/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_MP2F12_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_MP2F12_H_

#include <string>

#include "mpqc/chemistry/qc/lcao/f12/cabs_singles.h"
#include "mpqc/chemistry/qc/lcao/f12/f12_intermediates.h"
#include "mpqc/chemistry/qc/lcao/integrals/f12_utility.h"
#include "mpqc/chemistry/qc/lcao/scf/mo_build.h"
#include "mpqc/chemistry/qc/lcao/wfn/lcao_wfn.h"
#include "mpqc/chemistry/qc/properties/energy.h"

namespace mpqc {
namespace lcao {

/**
 *  \brief MP2F12 method for closed shell
 */

template <typename Tile>
class RMP2F12 : public LCAOWavefunction<Tile, TA::SparsePolicy>,
                public Provides<Energy> {
 public:
  using TArray = TA::DistArray<Tile, TA::SparsePolicy>;

  // clang-format off
  /**
   * KeyVal constructor
   * @param kv
   * keywords: takes all keywords from LCAOWavefunction
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | ref | Wavefunction | none | reference Wavefunction, RHF for example |
   * | approaximation | char | C | approaximation to use (C or D) |
   * | cabs_singles | bool | true | if do CABSSingles calculation |
   *
   */
  // clang-format on
  RMP2F12(const KeyVal& kv);
  ~RMP2F12() = default;

  void obsolete() override;

  const std::shared_ptr<Eigen::VectorXd> orbital_energy() const {
    return orbital_energy_;
  }
  const double mp2_corr_energy() const { return mp2_corr_energy_; }
  const double f12_energy() const { return mp2_f12_energy_; }
  const double cabs_singles_energy() const { return singles_energy_; }

 protected:
  bool can_evaluate(Energy* energy) override;

  void evaluate(Energy* result) override;

  /// compute mp2f12 energy
  /// @return tuple of E_MP2 matrix and E_F12 matrix
  std::tuple<RowMatrix<double>, RowMatrix<double>> compute();

 private:
  /// initialize mp2f12
  virtual void init(double ref_precision);

  /// initialize orbital energy
  virtual void init_orbital_energy();

  /// function to compute B intermediate
  virtual TArray compute_B();

  /// function to compute V intermediate
  virtual TArray compute_V();

  /// function to compute X intermediate
  virtual TArray compute_X();

  /// function to compute C intermediate
  virtual TArray compute_C();

  /// function to compute T2 amplitudes
  virtual std::tuple<TArray, TArray> compute_T();

  /// function to compute CABS singles correction
  virtual double compute_cabs_singles();

 protected:
  char approximation_;
  TA::SparseShape<float> ijij_ijji_shape_;
  bool cabs_singles_;
  std::shared_ptr<Wavefunction> ref_wfn_;
  std::shared_ptr<Eigen::VectorXd> orbital_energy_;

 private:
  /// MP2 correlation energy
  double mp2_corr_energy_ = 0.0;
  /// F12 contribution to MP2F12 energy
  double mp2_f12_energy_ = 0.0;
  /// CABS singles energy
  double singles_energy_ = 0.0;
  /// computed precision
  double computed_precision_ = std::numeric_limits<double>::max();
};

/**
 *  \brief MP2F12 method for closed shell with RI
 */

template <typename Tile>
class RIRMP2F12 : public RMP2F12<Tile> {
 public:
  using TArray = TA::DistArray<Tile, TA::SparsePolicy>;
  /**
 * KeyVal constructor
 * @param kv
 * keywords: takes all keywords from LCAOWavefunction
 *
 * | KeyWord | Type | Default| Description |
 * |---------|------|--------|-------------|
 * | ref | Wavefunction | none | reference Wavefunction, RHF for example |
 * | approaximation | char | C | approaximation to use (C or D) |
 * | cabs_singles | bool | true | if do CABSSingles calculation |
 *
 */
  RIRMP2F12(const KeyVal& kv);
  ~RIRMP2F12() = default;

 private:
  void init_orbital_energy() override;
  TArray compute_B() override;
  TArray compute_V() override;
  TArray compute_X() override;
  TArray compute_C() override;
  std::tuple<TArray, TArray> compute_T() override;
  double compute_cabs_singles() override;
};

extern template class RMP2F12<TA::TensorD>;
extern template class RIRMP2F12<TA::TensorD>;

}  // namespace lcao
}  // mpqc

#include "mp2f12_impl.h"

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_MP2F12_H_
