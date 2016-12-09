//
// Created by Chong Peng on 3/31/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_MP2F12_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_MP2F12_H_

#include <string>

#include "mpqc/chemistry/qc/f12/cabs_singles.h"
#include "mpqc/chemistry/qc/f12/f12_intermediates.h"
#include "mpqc/chemistry/qc/integrals/f12_utility.h"
#include "mpqc/chemistry/qc/mbpt/mp2.h"
#include "mpqc/chemistry/qc/wfn/trange1_engine.h"
#include "mpqc/math/external/eigen/eigen.h"

namespace mpqc {
namespace f12 {

/**
 *  \brief MP2F12 method for closed shell
 */

template <typename Tile>
class RMP2F12 : public qc::LCAOWavefunction<Tile, TA::SparsePolicy> {
 public:
  using TArray = TA::DistArray<Tile,TA::SparsePolicy>;

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
  RMP2F12(const KeyVal& kv);
  ~RMP2F12() = default;

  double value() override;
  std::tuple<RowMatrix<double>, RowMatrix<double>> compute();
  void compute(qc::PropertyBase* pb) override;
  void obsolete() override;

 private:
  virtual TArray compute_B();
  virtual TArray compute_V();
  virtual TArray compute_X();
  virtual TArray compute_C();
  virtual std::tuple<TArray, TArray> compute_T();
  virtual double compute_cabs_singles();

 protected:
  char approximation_;
  TA::SparseShape<float> ijij_ijji_shape_;
  bool cabs_singles_;
  std::shared_ptr<qc::Wavefunction> ref_wfn_;
};


/**
 *  \brief MP2F12 method for closed shell with RI
 */

template <typename Tile>
class RIRMP2F12 : public RMP2F12<Tile> {
 public:
  using TArray = TA::DistArray<Tile,TA::SparsePolicy>;
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
  virtual TArray compute_B() override;
  virtual TArray compute_V() override;
  virtual TArray compute_X() override;
  virtual TArray compute_C() override;
  virtual std::tuple<TArray, TArray> compute_T() override;
  virtual double compute_cabs_singles() override;
};

extern template class RMP2F12<TA::TensorD>;
extern template class RIRMP2F12<TA::TensorD>;

}  // end of namespace f12
}  // mpqc

#include "mp2f12_impl.h"

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_MP2F12_H_
