//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC_DBMP2F12_H
#define MPQC_DBMP2F12_H

#include <mpqc/chemistry/qc/f12/mp2f12.h>
#include <mpqc/chemistry/qc/mbpt/dbmp2.h>
#include <mpqc/chemistry/qc/f12/db_f12_intermediates.h>


namespace mpqc{
namespace f12{


template <typename Tile>
class DBMP2F12 : protected MP2F12<Tile>{
public:

  using Policy = TA::SparsePolicy;
  using TArray = TA::DistArray<Tile, Policy>;
  using MolecularIntegralClass = integrals::MolecularIntegral<Tile,Policy>;

  using real_t = typename MP2F12<Tile>::real_t;
  using Matrix = typename MP2F12<Tile>::Matrix;
  using MP2F12<Tile>::debug;
  using MP2F12<Tile>::mo_integral;
  using MP2F12<Tile>::trange1_engine;
  using MP2F12<Tile>::orbital_energy;

  DBMP2F12() = default;

  DBMP2F12(MolecularIntegralClass& mo_int) : MP2F12<Tile>(construct_db_mp2(mo_int)) {}

  virtual real_t compute(const rapidjson::Document& in){
    auto& world = this->mo_integral().get_world();
    // dbmp2 time
    auto mp2_time0 = mpqc_time::fenced_now(world);

    double dbmp2 = this->mp2_->compute(in);

    auto mp2_time1 = mpqc_time::fenced_now(world);
    auto mp2_time = mpqc_time::duration_in_s(mp2_time0, mp2_time1);
    mpqc::utility::print_par(world, "Total DBMP2 Time:  ", mp2_time, "\n");

    // start f12
    auto f12_time0 = mpqc_time::fenced_now(world);

    // solve cabs orbitals
    auto orbital_registry = this->mo_integral().orbital_space();
    closed_shell_dualbasis_cabs_mo_build_svd(this->mo_integral(),in,this->mp2_->trange1_engine());

    std::string method = in.HasMember("Method") ? in["Method"].GetString() : "df";

    Matrix mp2_eij, f12_eij;

    if(method == "df"){
      std::tie(mp2_eij, f12_eij) = compute_db_mp2_f12_c_df();
    }
    else{
      throw std::runtime_error("Wrong DBMP2F12 Method");
    }

    if (world.rank() == 0) {
      auto nocc = this->mp2_->trange1_engine()->get_active_occ();
      printf(
          "  i0     i1       eij(mp2)        eij(f12)      eij(mp2-f12) \n"
          "====== ====== =============== =============== ===============\n");
      for (int i = 0; i != nocc; ++i)
        for (int j = i; j != nocc; ++j)
          printf("%4d   %4d   %15.12lf %15.12lf %15.12lf\n", i, j, mp2_eij(i, j),
                 f12_eij(i, j), mp2_eij(i, j) + f12_eij(i, j));
    }

    auto ef12 = f12_eij.sum();
    if (debug()) {
      utility::print_par(this->mo_integral().get_world(), "E_DBMP2: ", dbmp2, "\n");
      utility::print_par(this->mo_integral().get_world(), "E_DBMP2F12: ", ef12, "\n");
    }

    auto f12_time1 = mpqc_time::fenced_now(world);
    auto f12_time = mpqc_time::duration_in_s(f12_time0, f12_time1);

    mpqc::utility::print_par(world, "Total DBF12 Time:  ", f12_time, "\n");

    return dbmp2 + ef12;
  }

private:

  std::shared_ptr<mbpt::MP2<Tile,Policy>> construct_db_mp2(MolecularIntegralClass& mo_int){
    std::shared_ptr<mbpt::MP2<Tile,Policy>> dbmp2 = std::make_shared<mbpt::DBMP2<Tile,Policy>>(mo_int);
    return dbmp2;
  }

  /// DBMP2-F12 C approach with density fitting
  std::tuple<Matrix,Matrix> compute_db_mp2_f12_c_df();

};


template <typename Tile>
std::tuple<typename DBMP2F12<Tile>::Matrix, typename DBMP2F12<Tile>::Matrix>
DBMP2F12<Tile>::compute_db_mp2_f12_c_df() {

  auto& world = this->mo_integral().get_world();

  Matrix Eij_MP2, Eij_F12;

  auto n_active_occ = this->mp2_->trange1_engine()->get_active_occ();
  auto n_occ = this->mp2_->trange1_engine()->get_occ();
  auto n_frozen = this->mp2_->trange1_engine()->get_nfrozen();

  // create shape
  auto occ_tr1 = this->mp2_->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  // T Term
  TArray t2;
  {
    utility::print_par(world, "Compute T_abij With DF \n" );

    TArray g_abij;
    g_abij("a,b,i,j") = this->mo_integral().compute(L"<i j|G|a b>[df]")("i,j,a,b");
    t2 = mpqc::cc::d_abij(g_abij, *(this->mp2_->orbital_energy()), n_occ, n_frozen);

    // compute MP2 energy and pair energies
    TArray TG_ijij_ijji;
    TG_ijij_ijji("i1,j1,i2,j2") =
        (t2("a,b,i1,j1") * g_abij("a,b,i2,j2"))
            .set_shape(ijij_ijji_shape);
    Eij_MP2 = TG_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(2, -1, n_active_occ));
  }

  //compute V term
  TArray V_ijij_ijji = compute_V_ijij_ijji_db_df(mo_integral(), ijij_ijji_shape);
  {
    // G integral in MO not needed, still need G integral in AO to compute F, K, hJ
    this->mo_integral().registry().remove_operation(world, L"G");

    //contribution from V_ijij_ijji
    Matrix eij = V_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(2 * C_ijij_bar,2 * C_ijji_bar,n_active_occ));
    if (debug()) utility::print_par(world, "E_V: ", eij.sum(), "\n");
    Eij_F12 = eij;
  }

  // compute C term
  TArray C_ijab = compute_C_ijab_df(mo_integral());

  {
    utility::print_par(world, "Compute CT With DF \n" );
    V_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

    Matrix eij = V_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(2 * C_ijij_bar,2 * C_ijji_bar,n_active_occ));
    if (debug()) utility::print_par(world, "E_CT: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  // compute X term
  TArray X_ijij_ijji = compute_X_ijij_ijji_db_df(mo_integral(), ijij_ijji_shape);
  {

    auto Fij = this->mo_integral().compute(L"<i|F|j>[df]");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

    Matrix eij = X_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar,n_active_occ));
    eij *= -1;
    if (debug()) utility::print_par(world, "E_X: ", eij.sum(), "\n");
    Eij_F12 += eij;

  }

  // compute B term
  TArray B_ijij_ijji = compute_B_ijij_ijji_db_df(mo_integral(), ijij_ijji_shape);
  {
    Matrix eij = B_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar,n_active_occ));
    if (debug()) utility::print_par(world, "E_B: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  {
    utility::print_par(world, "Compute CC Term With DF \n");
    auto C_bar_ijab = f12::convert_C_ijab(C_ijab, n_occ, n_frozen, *(this->mp2_->orbital_energy()));
    B_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

    Matrix eij = B_ijij_ijji("i1,j1,i2,j2").reduce(F12PairEnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar,n_active_occ));
    if (debug()) utility::print_par(world, "E_CC: ", eij.sum(), "\n");
    Eij_F12 += eij;
  }

  return std::make_tuple(Eij_MP2,Eij_F12);
}

}  // namespace f12
}  // namespace mpqc

#endif //MPQC_DBMP2F12_H
