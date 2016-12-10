//
// Created by Chong Peng on 7/12/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_DBCCSD_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_DBCCSD_H_

#include "mpqc/chemistry/qc/cc/ccsd.h"

namespace mpqc {
namespace lcao {

/**
 *  \breif Dual basis CCSD method
 */
template <typename Tile, typename Policy>
class DBCCSD : public CCSD<Tile, Policy> {
  using TArray = TA::DistArray<Tile, Policy>;

 public:
  DBCCSD() = default;
  virtual ~DBCCSD() {}

  /**
   *  KeyVal constructor
   *  keywords: all keywords of CCSD class
   *
   * | KeyWord | Type | Default| Description |
   * |---------|------|--------|-------------|
   * | method | string | standard | method to compute ccsd (standard or df) |
   *
   */

  DBCCSD(const KeyVal& kv) : CCSD<Tile, Policy>(kv) {
    std::string method = kv.value<std::string>("method", "standard");
    if (method == "direct") {
      throw std::invalid_argument(
          "Integral Direct Dual Basis CCSD is not Implemented!!\n");
    }
  };

 private:
  void init() override {
    if (this->orbital_energy() == nullptr ||
        this->trange1_engine() == nullptr) {
      auto& lcao_factory = this->lcao_factory();
      auto mol = lcao_factory.ao_factory().molecule();
      Eigen::VectorXd orbital_energy;
      this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(
          lcao_factory, orbital_energy, mol, this->is_frozen_core(),
          this->occ_block(), this->unocc_block());
      this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
    }
  }
};

extern template
class DBCCSD<TA::TensorD, TA::SparsePolicy>;

}  // namespace  lcao
}  // namespace  mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_CC_DBCCSD_H_
