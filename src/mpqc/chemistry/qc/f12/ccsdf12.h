//
// Created by Chong Peng on 4/12/16.
//

#ifndef MPQC_CCSDF12_H
#define MPQC_CCSDF12_H


#include "../../../../../include/tiledarray.h"
#include "../../../../../common/namespaces.h"
#include <mpqc/chemistry/qc/integrals/molecular_integral.h>
#include "../../../../../utility/trange1_engine.h"
#include <mpqc/chemistry/qc/f12/f12_intermediates.h>
#include <mpqc/chemistry/qc/cc/ccsd.h>

namespace mpqc{
namespace f12{

template <typename Tile>
class CCSDF12 {
public:

    using Policy = TA::SparsePolicy;
    using TArray = TA::DistArray<Tile,Policy>;
    using MolecularIntegralClass = integrals::MolecularIntegral<Tile,Policy>;


    CCSDF12(integrals::MolecularIntegral<Tile, Policy> &mo_int, rapidjson::Document& options) : mo_int_(mo_int)
    {
        ccsd_ = std::make_shared<cc::CCSD<Tile, Policy>>(mo_int,options);
    }

    double compute(){
        // compute ccsd
        double ccsd = ccsd_->compute();

        auto& option = ccsd_->options();
        // initialize CABS orbitals
        auto& ao_int = this->mo_int_.atomic_integral();
        auto orbital_registry = this->mo_int_.orbital_space();
        closed_shell_cabs_mo_build_svd(this->mo_int_,option,this->ccsd_->trange1_engine());


        auto lazy_two_electron_int = ccsd_->intermediate()->direct_ao();
        double f12;

        bool df;
        std::string method = option.HasMember("Method") ? option["Method"].GetString() : "df";
        if(method == "four center"){
                f12 = compute_c(lazy_two_electron_int);
        }
        else if(method == "df"){
                f12 = compute_c_df(lazy_two_electron_int);
        }
        else{
            throw std::runtime_error("Wrong CCSDF12 Method");
        }

        return ccsd + f12;
    }

private:
    /// standard approach
    template<typename DirectArray>
    double compute_c_df(const DirectArray &darray);

    template<typename DirectArray>
    double compute_c(const DirectArray &darray);

private:

    MolecularIntegralClass& mo_int_;
    std::shared_ptr<cc::CCSD<Tile,Policy>> ccsd_;

//    TArray ccsd_->t1(); /// t1 amplitude
//    TArray ccsd_->t2(); /// t2 amplitude

};


template <typename Tile>
template <typename DirectArray>
double CCSDF12<Tile>::compute_c_df(const DirectArray &darray) {

    auto& mo_integral = mo_int_;
    auto& world = mo_integral.get_world();
    double E = 0.0;

    // clean MO integrals
    mo_integral.registry().clear();

    // create shape
    auto occ_tr1 = ccsd_->trange1_engine()->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    // compute V_ijij_ijji
    TArray V_ijij_ijji = compute_V_ijij_ijji_df(mo_integral, ijij_ijji_shape);

    // VT2 contribution
    if(darray.is_initialized()){
        TArray tmp = compute_VT2_ijij_ijji_df_direct(mo_integral, ccsd_->t2(), ijij_ijji_shape, darray);
        V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    }else{
        TArray tmp = compute_VT2_ijij_ijji_df(mo_integral,ccsd_->t2(),ijij_ijji_shape);
        V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    }

    // VT1 contribution
    {
        TArray tmp = compute_VT1_ijij_ijji_df(mo_integral,ccsd_->t1(),ijij_ijji_shape);
        V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    }

    // V contribution to energy
    double E_v = V_ijij_ijji("i1,j1,i2,j2").reduce(f12::F12EnergyReductor<Tile>(2 * C_ijij_bar, 2 * C_ijji_bar));
    utility::print_par(world, "E_V: ", E_v, "\n");
    E += E_v;

    // compute X term
    TArray X_ijij_ijji = compute_X_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    // R_ipjq not needed
    mo_int_.registry().remove_formula(world, L"(i1 p|R|j1 q)[df]");

    auto Fij = mo_int_.compute(L"(i|F|j)[df]");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

    double E_x = -X_ijij_ijji("i1,j1,i2,j2").reduce(f12::F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));

    utility::print_par(world, "E_X: ", E_x, "\n");
    E += E_x;

    // compute B term
    TArray B_ijij_ijji = compute_B_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    double E_b = B_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
    utility::print_par(world, "E_B: ", E_b, "\n");
    E += E_b;

    utility::print_par(world, "E_F12: ", E, "\n");

    return E;
}


template <typename Tile>
template <typename DirectArray>
double CCSDF12<Tile>::compute_c(const DirectArray &darray) {

    auto& mo_integral = mo_int_;
    auto& world = mo_integral.get_world();
    double E = 0.0;

    // clean MO integrals
    mo_integral.registry().clear();

    // create shape
    auto occ_tr1 = ccsd_->trange1_engine()->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    // compute V_ijij_ijji
    TArray V_ijij_ijji = compute_V_ijij_ijji(mo_integral, ijij_ijji_shape);

//    std::cout << "V_ijij_ijji" << std::endl;
//    std::cout << V_ijij_ijji << std::endl;

    // VT2 contribution
//    if(darray.is_initialized()){
//        TArray tmp = compute_VT2_ijij_ijji_df_direct(mo_integral, ccsd_->t2(), ijij_ijji_shape, darray);
//        V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
//    }else
    {
        TArray tmp = compute_VT2_ijij_ijji(mo_integral,ccsd_->t2(),ijij_ijji_shape);
        V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");

    }

    // VT1 contribution
    {
        TArray tmp = compute_VT1_ijij_ijji(mo_integral,ccsd_->t1(),ijij_ijji_shape);
        V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    }

    // V contribution to energy
    double E_v = V_ijij_ijji("i1,j1,i2,j2").reduce(f12::F12EnergyReductor<Tile>(2 * C_ijij_bar, 2 * C_ijji_bar));
    utility::print_par(world, "E_V: ", E_v, "\n");
    E += E_v;

//    {
//        utility::print_par(world, "Compute CC Term Without DF \n");
//        auto C_ijab = compute_C_ijab(mo_integral);
//        auto C_bar_ijab = f12::convert_C_ijab(C_ijab, occ, *orbital_energy_);
//        V_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);
//
//        double E_cc = V_ijij_ijji("i1,j1,i2,j2").reduce(f12::CLF12Energy<Tile>(CC_ijij_bar,CC_ijji_bar));
//        utility::print_par(world, "E_CC: ", E_cc, "\n");
//        E += E_cc;
//    }
    // compute X term
    TArray X_ijij_ijji = compute_X_ijij_ijji(mo_integral, ijij_ijji_shape);
//    std::cout << "X_ijij_ijji" << std::endl;
//    std::cout << X_ijij_ijji << std::endl;


    // R_ipjq not needed
    mo_int_.registry().remove_formula(world, L"(i1 p|R|j1 q)");

    auto Fij = mo_int_.compute(L"(i|F|j)");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

    double E_x = -X_ijij_ijji("i1,j1,i2,j2").reduce(f12::F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));

    utility::print_par(world, "E_X: ", E_x, "\n");
    E += E_x;

    // compute B term
    TArray B_ijij_ijji = compute_B_ijij_ijji(mo_integral, ijij_ijji_shape);
//    std::cout << "B_ijij_ijji" << std::endl;
//    std::cout << B_ijij_ijji << std::endl;

    double E_b = B_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
    utility::print_par(world, "E_B: ", E_b, "\n");
    E += E_b;

    utility::print_par(world, "E_F12: ", E, "\n");

    return E;
}

}//end of namespace f12
}//end of namespace mpqc


#endif //MPQC_CCSDF12_H
