//
// Created by Chong Peng on 7/12/16.
//

#ifndef MPQC_DBCCSD_H
#define MPQC_DBCCSD_H

#include <mpqc/chemistry/qc/cc/ccsd.h>

namespace mpqc{
namespace cc{

/**
 * Takes all the argument of CCSD class
 */
template <typename Tile, typename Policy>
class DBCCSD : public CCSD<Tile,Policy>{

using TArray = TA::DistArray<Tile,Policy>;
public:
  DBCCSD() = default;
  DBCCSD(integrals::MolecularIntegral<Tile,Policy>& mo_int, rapidjson::Document &options)
          : CCSD<Tile,Policy>(mo_int,options)
  { }

  /// compute function
  virtual double compute(){
    // initialize
    init(this->options_);
    TArray t1;
    TArray t2;

    auto direct = this->options_.HasMember("Direct") ? this->options_["Direct"].GetBool(): true;
    double ccsd_corr = 0.0;
    if(direct){
//                    double ccsd_corr = compute_ccsd_direct_ao(t1, t2);
      ccsd_corr = this->compute_ccsd_direct(t1, t2);
    }
    else {
      ccsd_corr = this->compute_ccsd_nondirect(t1,t2);
    }

    this->T1_ = t1;
    this->T2_ = t2;

//                ccsd_intermediate_->clean_two_electron();

    return ccsd_corr;
  }
private:

  void init(const rapidjson::Document &in) {
    if(this->orbital_energy_== nullptr || this->trange1_engine_ == nullptr) {
      auto &mo_int = this->ccsd_intermediate_->mo_integral();
      auto mol = mo_int.atomic_integral().molecule();
      int occ = mol.occupation(0) / 2;
      Eigen::VectorXd orbital_energy;
      this->trange1_engine_ = closed_shell_dualbasis_mo_build_eigen_solve_svd(mo_int, orbital_energy, in, mol, occ);
      this->orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
    }
  }

};

}
}

#endif //MPQC_DBCCSD_H
