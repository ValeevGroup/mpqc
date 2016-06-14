//
// Created by Chong Peng on 6/13/16.
//

#ifndef MPQC_DBMP2_H
#define MPQC_DBMP2_H

#include <mpqc/chemistry/qc/mbpt/mp2.h>

namespace mpqc{
namespace mbpt{


template <typename Tile, typename Policy>
class DBMP2{
public:
  using TArray = TA::DistArray<Tile,Policy>;
  using MolecularIntegral = integrals::MolecularIntegral<Tile,Policy>;

  DBMP2(MolecularIntegral& mo_int, const rapidjson::Document& in){

    auto& ao_int = mo_int.atomic_integral();
    auto orbital_registry = mo_int.orbital_space();
    auto mol = mo_int.atomic_integral().molecule();
    int occ = mol->occupation(0)/2;
    Eigen::VectorXd orbital_energy;
    auto trange1_engine = closed_shell_dual_basis_mo_build_eigen_solve(ao_int, *orbital_registry, orbital_energy, in, *mol, occ);

      mp2_ = std::make_shared<MP2<Tile,Policy>>(MP2<Tile,Policy>(mo_int,orbital_energy,trange1_engine));

  }

  double compute(const rapidjson::Document& in){
    return mp2_->compute(in);
  }

private:
  std::shared_ptr<MP2<Tile,Policy>> mp2_;

};

}
}

#endif //MPQC_DBMP2_H
