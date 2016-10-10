//
// Created by Chong Peng on 7/12/16.
//

#ifndef MPQC_DBCCSD_H
#define MPQC_DBCCSD_H

#include <mpqc/chemistry/qc/cc/ccsd.h>

namespace mpqc {
namespace cc {

/**
 * Takes all the argument of CCSD class
 */
template <typename Tile, typename Policy>
class DBCCSD : public CCSD<Tile, Policy> {
  using TArray = TA::DistArray<Tile, Policy>;

 public:
  DBCCSD() = default;
  DBCCSD(integrals::LCAOFactory<Tile, Policy> &lcao_factory,
         rapidjson::Document &options)
      : CCSD<Tile, Policy>(lcao_factory, options) {
    auto direct = this->options_.HasMember("Direct")
                      ? this->options_["Direct"].GetBool()
                      : true;
    if (direct == true) {
      throw std::runtime_error(
          "Integral Direct Dual Basis CCSD is not Implemented!!\n");
    }
  }
  /// compute function
  virtual double compute() {
    // initialize
    init(this->options_);
    TArray t1;
    TArray t2;

    double ccsd_corr = 0.0;
    ccsd_corr = this->compute_ccsd_conventional(t1, t2);

    this->T1_ = t1;
    this->T2_ = t2;

    //                ccsd_intermediate_->clean_two_electron();

    return ccsd_corr;
  }

 private:
  void init(const rapidjson::Document &in) {
    if (this->orbital_energy_ == nullptr || this->trange1_engine_ == nullptr) {
      auto &lcao_factory = this->ccsd_intermediate_->lcao_factory();
      auto mol = lcao_factory.atomic_integral().molecule();
      Eigen::VectorXd orbital_energy;
      this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(
          lcao_factory, orbital_energy, in, mol);
      this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
    }
  }
};
}
}

#endif  // MPQC_DBCCSD_H
