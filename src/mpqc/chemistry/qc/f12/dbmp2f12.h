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

  DBMP2F12() = default;

  DBMP2F12(MolecularIntegralClass& mo_int) : MP2F12<Tile>(construct_db_mp2(mo_int)) {}

  virtual double compute(const rapidjson::Document& in){
    double result = 0.0;
    auto& world = this->mo_int_.get_world();
    // dbmp2 time
    auto mp2_time0 = mpqc_time::fenced_now(world);

    result += this->mp2_->compute(in);

    auto mp2_time1 = mpqc_time::fenced_now(world);
    auto mp2_time = mpqc_time::duration_in_s(mp2_time0, mp2_time1);
    mpqc::utility::print_par(world, "Total DBMP2 Time:  ", mp2_time, "\n");

    // start f12
    auto f12_time0 = mpqc_time::fenced_now(world);

    // solve cabs orbitals
    auto& ao_int = this->mo_int_.atomic_integral();
    auto orbital_registry = this->mo_int_.orbital_space();
    closed_shell_dualbasis_cabs_mo_build_svd(this->mo_int_,in,this->mp2_->trange1_engine());

    std::string method = in.HasMember("Method") ? in["Method"].GetString() : "df";

    if(method == "df"){
      result += compute_db_mp2_f12_c_df();
    }
    else{
      throw std::runtime_error("Wrong DBMP2F12 Method");
    }
    auto f12_time1 = mpqc_time::fenced_now(world);
    auto f12_time = mpqc_time::duration_in_s(f12_time0, f12_time1);

    mpqc::utility::print_par(world, "Total DBF12 Time:  ", f12_time, "\n");


    return result;
  }

private:

  std::shared_ptr<mbpt::MP2<Tile,Policy>> construct_db_mp2(MolecularIntegralClass& mo_int){
    std::shared_ptr<mbpt::MP2<Tile,Policy>> dbmp2 = std::make_shared<mbpt::DBMP2<Tile,Policy>>(mo_int);
    return dbmp2;
  }

  /// DBMP2-F12 C approach with density fitting
  double compute_db_mp2_f12_c_df();

};


template <typename Tile>
double DBMP2F12<Tile>::compute_db_mp2_f12_c_df() {

  auto& world = this->mo_int_.get_world();

  double E = 0.0;

  auto& mo_integral = this->mo_int_;

  auto occ = this->mp2_->trange1_engine()->get_active_occ();

  // create shape
  auto occ_tr1 = this->mp2_->trange1_engine()->get_occ_tr1();
  TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
  auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

  //compute V term
  TArray V_ijij_ijji = compute_V_ijij_ijji_db_df(mo_integral, ijij_ijji_shape);
  {
    // G integral in MO not needed, still need G integral in AO to compute F, K, hJ
    this->mo_int_.registry().remove_operation(world, L"G");

    //contribution from V_ijij_ijji
    double E_v = V_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(2 * C_ijij_bar,2 * C_ijji_bar));
    utility::print_par(world, "E_V: ", E_v, "\n");
    E += E_v;
  }

  // T Term
  TArray t2;
  {
    utility::print_par(world, "Compute T_abij With DF \n" );

    TArray g_iajb;
    g_iajb = this->mo_int_.compute(L"<i j|G|a b>[df]");
    g_iajb("a,b,i,j") = g_iajb("i,j,a,b");
    t2 = mpqc::cc::d_abij(g_iajb,*(this->mp2_->orbital_energy()),occ);

  }

  // compute C term
  TArray C_ijab = compute_C_ijab_df(mo_integral);

  {
    utility::print_par(world, "Compute CT With DF \n" );
    V_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

    double E_ct = V_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(2 * C_ijij_bar,2 * C_ijji_bar));
    utility::print_par(world, "E_CT: ", E_ct, "\n");
    E += E_ct;
  }

  // compute X term
  TArray X_ijij_ijji = compute_X_ijij_ijji_db_df(mo_integral, ijij_ijji_shape);
  {

    auto Fij = this->mo_int_.compute(L"<i|F|j>[df]");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

    double E_x = -X_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
    utility::print_par(world, "E_X: ", E_x, "\n");
    E += E_x;

  }

  // compute B term
  TArray B_ijij_ijji = compute_B_ijij_ijji_db_df(mo_integral, ijij_ijji_shape);
  {
    double E_b = B_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
    utility::print_par(world, "E_B: ", E_b, "\n");
    E += E_b;
  }

  {
    utility::print_par(world, "Compute CC Term With DF \n");
    auto C_bar_ijab = f12::convert_C_ijab(C_ijab, occ, *(this->mp2_->orbital_energy()));
    B_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

    double E_cc = B_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
    utility::print_par(world, "E_CC: ", E_cc, "\n");
    E += E_cc;
  }

  utility::print_par(world, "E_F12: ", E, "\n");
  return E;
}

} // end of namespace f12
} // end of namespace mpqc



#endif //MPQC_DBMP2F12_H
