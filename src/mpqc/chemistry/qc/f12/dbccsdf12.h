//
// Created by Chong Peng on 7/12/16.
//

#ifndef MPQC_DB_CCSDF12_H
#define MPQC_DB_CCSDF12_H

#include <mpqc/chemistry/qc/cc/dbccsd.h>
#include <mpqc/chemistry/qc/f12/ccsdf12.h>
#include <mpqc/chemistry/qc/f12/db_f12_intermediates.h>
#include <mpqc/chemistry/qc/f12/cabs_singles.h>

namespace mpqc {
namespace f12 {

/***
 * \class DBCCSDF12 Dual basis CCSD F12 class
 *
 * takes all options from CCSDF12 and the flowing options
 *
 * @param RIMethod: OBS or VBS, if OBS, RI Basis is union of OBS+CABS; if VBS,
 * RI Basis is union of VBS+CABS, default VBS
 *
 */

template <typename Tile>
class DBCCSDF12 : public CCSDF12<Tile> {
 public:
  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using LCAOFactoryType = integrals::LCAOFactory<Tile, Policy>;

  using real_t = typename Tile::scalar_type;
  using Matrix = RowMatrix<real_t>;

  DBCCSDF12(LCAOFactoryType& lcao_factory, rapidjson::Document& options)
      : CCSDF12<Tile>(construct_db_ccsd(lcao_factory, options)) {}

  virtual real_t compute() {
    auto& world = this->lcao_factory_.get_world();
    // compute ccsd
    real_t ccsd = this->ccsd_->compute();

    auto& option = this->ccsd_->options();
    // initialize CABS orbitals
    closed_shell_dualbasis_cabs_mo_build_svd(this->lcao_factory_, option,
                                             this->ccsd_->trange1_engine());

    Matrix Eij_F12;

    std::string method =
        option.HasMember("Method") ? option["Method"].GetString() : "df";

    if (method == "df") {
      Eij_F12 = compute_c_df();
    } else {
      throw std::runtime_error("Wrong DBCCSDF12 Method");
    }

    real_t e_f12 = Eij_F12.sum();
    if (debug()) utility::print_par(world, "E_F12: ", e_f12, "\n");

    // compute cabs singles
    real_t e_s = 0.0;
    bool singles =
        option.HasMember("Singles") ? option["Singles"].GetBool() : true;
    if (singles) {
      auto single_time0 = mpqc_time::fenced_now(world);

      // non-canonical, don't include F_m^a
      CABSSingles<Tile> cabs_singles(this->lcao_factory_, false);
      e_s = cabs_singles.compute();
      if (debug()) {
        utility::print_par(world, "E_S: ", e_s, "\n");
      }
      auto single_time1 = mpqc_time::fenced_now(world);
      auto single_time = mpqc_time::duration_in_s(single_time0, single_time1);
      mpqc::utility::print_par(world, "Total CABS Singles Time:  ", single_time,
                               "\n");
    }

    return ccsd + e_f12 + e_s;
  }

  using CCSDF12<Tile>::debug;

 private:
  Matrix compute_c_df();

  std::shared_ptr<cc::CCSD<Tile, Policy>> construct_db_ccsd(
      LCAOFactoryType& lcao_factory, rapidjson::Document& options) {
    std::shared_ptr<cc::CCSD<Tile, Policy>> dbccsd =
        std::make_shared<cc::DBCCSD<Tile, Policy>>(lcao_factory, options);
    return dbccsd;
  }
};

template <typename Tile>
typename DBCCSDF12<Tile>::Matrix DBCCSDF12<Tile>::compute_c_df() {
  auto& lcao_factory = this->lcao_factory_;
  auto& world = lcao_factory.get_world();
  Matrix Eij_F12;

  // clean LCAO Integrals
  lcao_factory.registry().clear();

  auto nocc = this->ccsd_->trange1_engine()->get_active_occ();

  // create shape
  auto occ_tr1 = this->ccsd_->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1, occ_tr1, occ_tr1, occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // compute V_ijij_ijji
  TArray V_ijij_ijji = compute_V_ijij_ijji_db_df(lcao_factory, ijij_ijji_shape);

  // VT2 contribution
  TArray tmp = compute_VT2_ijij_ijji_db_df(lcao_factory, this->ccsd_->t2(),
                                           ijij_ijji_shape);
  V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");

  // VT1 contribution
  {
    TArray tmp = compute_VT1_ijij_ijji_db_df(lcao_factory, this->ccsd_->t1(),
                                             ijij_ijji_shape);
    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
  }

  // V contribution to energy
  Eij_F12 = V_ijij_ijji("i1,j1,i2,j2")
                .reduce(f12::F12PairEnergyReductor<Tile>(2 * C_ijij_bar,
                                                         2 * C_ijji_bar, nocc));
  if (debug()) utility::print_par(world, "E_V: ", Eij_F12.sum(), "\n");

  // compute X term
  TArray X_ijij_ijji = compute_X_ijij_ijji_db_df(lcao_factory, ijij_ijji_shape);

  auto Fij = this->lcao_factory_.compute(L"(i|F|j)[df]");
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
    TArray B_ijij_ijji =
        compute_B_ijij_ijji_db_df(lcao_factory, ijij_ijji_shape);
    Matrix eij = B_ijij_ijji("i1,j1,i2,j2")
                     .reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,
                                                         CC_ijji_bar, nocc));
    if (debug()) utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  return Eij_F12;
}

}  // end of namespace f12
}  // end of namespace mpqc

#endif  // MPQC_DB_CCSDF12_H
