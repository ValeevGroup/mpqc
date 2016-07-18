//
// Created by Chong Peng on 7/12/16.
//

#ifndef MPQC_DB_CCSDF12_H
#define MPQC_DB_CCSDF12_H

#include <mpqc/chemistry/qc/cc/dbccsd.h>
#include <mpqc/chemistry/qc/f12/ccsdf12.h>
#include <mpqc/chemistry/qc/f12/db_f12_intermediates.h>

namespace mpqc{
namespace f12{

template <typename Tile>
class DBCCSDF12 : public CCSDF12<Tile>{

public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using MolecularIntegralClass = integrals::MolecularIntegral<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  DBCCSDF12(MolecularIntegralClass & mo_int, rapidjson::Document& options) : CCSDF12<Tile>(construct_db_ccsd(mo_int,options))
  {
  }

  virtual real_t compute(){

    // compute ccsd
    real_t ccsd = this->ccsd_->compute();

    auto& option = this->ccsd_->options();
    // initialize CABS orbitals
    closed_shell_dualbasis_cabs_mo_build_svd(this->mo_int_, option,
                                   this->ccsd_->trange1_engine());

    Matrix Eij_F12;

    std::string method = option.HasMember("Method") ? option["Method"].GetString() : "df";

    if (method == "df") {
      Eij_F12 = compute_c_df();
    } else {
      throw std::runtime_error("Wrong CCSDF12 Method");
    }

    return ccsd + Eij_F12.sum();

    return 0;
  }

  using CCSDF12<Tile>::debug;

private:

  Matrix compute_c_df();

  std::shared_ptr<cc::CCSD<Tile,Policy>> construct_db_ccsd(MolecularIntegralClass & mo_int, rapidjson::Document& options){
    std::shared_ptr<cc::CCSD<Tile,Policy>> dbccsd = std::make_shared<cc::DBCCSD<Tile,Policy>>(mo_int,options);
    return dbccsd;
  }

};

template <typename Tile>
typename DBCCSDF12<Tile>::Matrix DBCCSDF12<Tile>::compute_c_df()
{
  auto& mo_integral = this->mo_int_;
  auto& world = mo_integral.get_world();
  Matrix Eij_F12;

  // clean MO integrals
  mo_integral.registry().clear();

  auto nocc = this->ccsd_->trange1_engine()->get_active_occ();

  // create shape
  auto occ_tr1 = this->ccsd_->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // compute V_ijij_ijji
  TArray V_ijij_ijji = compute_V_ijij_ijji_db_df(mo_integral, ijij_ijji_shape);

  // VT2 contribution
  TArray tmp = compute_VT2_ijij_ijji_db_df(mo_integral, this->ccsd_->t2(), ijij_ijji_shape);
  V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");

  // VT1 contribution
  {
    TArray tmp =
            compute_VT1_ijij_ijji_db_df(mo_integral, this->ccsd_->t1(), ijij_ijji_shape);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // V contribution to energy
  Eij_F12 = V_ijij_ijji("i1,j1,i2,j2")
          .reduce(f12::F12PairEnergyReductor<Tile>(2 * C_ijij_bar,
                                                   2 * C_ijji_bar, nocc));
  if (debug()) utility::print_par(world, "E_V: ", Eij_F12.sum(), "\n");

  // compute X term
  TArray X_ijij_ijji = compute_X_ijij_ijji_db_df(mo_integral, ijij_ijji_shape);

  auto Fij = this->mo_int_.compute(L"(i|F|j)[df]");
  auto Fij_eigen = array_ops::array_to_eigen(Fij);
  f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);
  {
    Matrix eij = X_ijij_ijji("i1,j1,i2,j2")
            .reduce(f12::F12PairEnergyReductor<Tile>(
                    CC_ijij_bar, CC_ijji_bar, nocc));
    eij *= -1;
    if (debug()) utility::print_par(world, "E_X: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  // compute B term
  {
    TArray B_ijij_ijji = compute_B_ijij_ijji_db_df(mo_integral, ijij_ijji_shape);
    Matrix eij = B_ijij_ijji("i1,j1,i2,j2")
            .reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,
                                                CC_ijji_bar, nocc));
    if (debug()) utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  if (debug()) utility::print_par(world, "E_F12: ", Eij_F12.sum(), "\n");

  return Eij_F12;

}

}// end of namespace f12
}// end of namespace mpqc

#endif //MPQC_DB_CCSDF12_H
