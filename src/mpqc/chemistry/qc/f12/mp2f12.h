//
// Created by Chong Peng on 3/31/16.
//

#ifndef MPQC_MP2F12_H
#define MPQC_MP2F12_H

#include <string>

#include <mpqc/chemistry/qc/f12/f12_utility.h>
#include <mpqc/chemistry/qc/f12/f12_intermediates.h>
#include "../../../../../utility/cc_utility.h"
#include "../../../../../utility/trange1_engine.h"

namespace mpqc{
namespace f12{

template <typename Tile>
class MP2F12 {

public:
    using Policy = TA::SparsePolicy;
    using TArray = TA::DistArray<Tile, Policy>;
    using MolecularIntegralClass = integrals::MolecularIntegral<Tile,Policy>;

    MP2F12() = default;

    MP2F12(MolecularIntegralClass& mo_int, std::shared_ptr<TRange1Engine> tre, const Eigen::VectorXd& ens)
            : mo_int_(mo_int), tre_(tre), orbital_energy_(std::make_shared<Eigen::VectorXd>(ens)) {}

    double compute_mp2_f12_c_df();

    double compute_mp2_f12_c();

    double compute(const rapidjson::Document& in){

        std::string method = in.HasMember("Method") ? in["Method"].GetString() : "df";

        double mp2_f12_energy = 0.0;

        if(method == "four center"){
            mp2_f12_energy = compute_mp2_f12_c();
        }
        else if(method == "df"){
            mp2_f12_energy = compute_mp2_f12_c_df();
        }
        else{
            throw std::runtime_error("Wrong MP2F12 Method");
        }

        return mp2_f12_energy;
    }

private:

private:

    MolecularIntegralClass& mo_int_;
    std::shared_ptr<TRange1Engine> tre_;
    std::shared_ptr<Eigen::VectorXd> orbital_energy_;
};

template <typename Tile>
double MP2F12<Tile>::compute_mp2_f12_c_df() {

    auto& world = mo_int_.get_world();

    double E = 0.0;

    auto& mo_integral = mo_int_;

    auto occ = tre_->get_actual_occ();

    TArray t2;
    {
        utility::print_par(world, "Compute T_abij With DF \n" );

        TArray g_iajb;
        g_iajb = mo_int_.compute(L"<i j|G|a b>[df]");
        g_iajb("a,b,i,j") = g_iajb("i,j,a,b");
        t2 = mpqc::cc::d_abij(g_iajb,*orbital_energy_,occ);

    }

    // create shape
    auto occ_tr1 = tre_->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    //compute V term
    TArray V_ijij_ijji = compute_V_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    {

        // G integral in MO not needed, still need G integral in AO to compute F, K, hJ
        mo_int_.registry().remove_operation(world, L"G");

        auto V_map = V_ijij_ijji.get_pmap();
        auto local = V_map->local_size();
        std::cout << "V PMap Local Size, Rank " << world.rank() << " Size " << local << std::endl;

        //contribution from V_ijij_ijji
        // NB factor of 2 from the Hylleraas functional
        double E_v = V_ijij_ijji("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(2 * C_ijij_bar,2 * C_ijji_bar));
        utility::print_par(world, "E_V: ", E_v, "\n");
        E += E_v;
    }

    // compute C term
    TArray C_ijab = compute_C_ijab_df(mo_integral);

    {
        utility::print_par(world, "Compute CT With DF \n" );
        V_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

        // NB factor of 2 from the Hylleraas functional
        double E_ct = V_ijij_ijji("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(2 * C_ijij_bar,2 * C_ijji_bar));
        utility::print_par(world, "E_CT: ", E_ct, "\n");
        E += E_ct;
    }


    // compute X term
    TArray X_ijij_ijji = compute_X_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    {

        // R_ipjq not needed
        mo_int_.registry().remove_formula(world, L"<i1 j1|R|p q>[df]");

        auto Fij = mo_int_.compute(L"<i|F|j>[df]");
        auto Fij_eigen = array_ops::array_to_eigen(Fij);
        f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

        double E_x = -X_ijij_ijji("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_X: ", E_x, "\n");
        E += E_x;

    }

    // compute B term
    TArray B_ijij_ijji = compute_B_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    {
        double E_b = B_ijij_ijji("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_B: ", E_b, "\n");
        E += E_b;
    }

    {

        utility::print_par(world, "Compute CC Term With DF \n");
        auto C_bar_ijab = f12::convert_C_ijab(C_ijab, occ, *orbital_energy_);
        B_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

        double E_cc = B_ijij_ijji("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_CC: ", E_cc, "\n");
        E += E_cc;
    }

    utility::print_par(world, "E_F12: ", E, "\n");
    return E;
}


template <typename Tile>
double MP2F12<Tile>::compute_mp2_f12_c(){

    auto& world = mo_int_.get_world();

    double E = 0.0;

    auto& mo_integral = mo_int_;

    auto occ = tre_->get_actual_occ();

    TArray t2_nodf;
    {
        utility::print_par(world, "Compute T_abij Without DF \n" );
        TArray g_iajb;
        g_iajb = mo_integral.compute(L"<i j|G|a b>");
        g_iajb("a,b,i,j") = g_iajb("i,j,a,b");
        t2_nodf = mpqc::cc::d_abij(g_iajb,*orbital_energy_,occ);
    }

    // create shape
    auto occ_tr1 = tre_->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    TArray V_ijij_ijji_nodf = compute_V_ijij_ijji(mo_integral, ijij_ijji_shape);
    {
        double E_v = V_ijij_ijji_nodf("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(2 * C_ijij_bar, 2 * C_ijji_bar));
        utility::print_par(world, "E_V: ", E_v, "\n");
        E += E_v;
    }

    TArray C_ijab_nodf = compute_C_ijab(mo_integral);

    {
        utility::print_par(world, "Compute CT Without DF \n" );
        V_ijij_ijji_nodf("i1,j1,i2,j2") = (C_ijab_nodf("i1,j1,a,b")*t2_nodf("a,b,i2,j2")).set_shape(ijij_ijji_shape);

        double E_ct = V_ijij_ijji_nodf("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(2 * C_ijij_bar, 2 * C_ijji_bar));
        utility::print_par(world, "E_CT: ", E_ct, "\n");
        E += E_ct;
    }

    TArray X_ijij_ijji_nodf = compute_X_ijij_ijji(mo_integral, ijij_ijji_shape);
    {

        // compute energy contribution
        auto Fij = mo_integral.compute(L"<i|F|j>");
        auto Fij_eigen = array_ops::array_to_eigen(Fij);
        f12::convert_X_ijkl(X_ijij_ijji_nodf, Fij_eigen);

        double E_x = -X_ijij_ijji_nodf("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_X: ", E_x, "\n");
        E += E_x;
    }

    TArray B_ijij_ijji_nodf = compute_B_ijij_ijji(mo_integral,ijij_ijji_shape);
    {

        double E_b = B_ijij_ijji_nodf("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_B: ", E_b, "\n");
        E += E_b;
    }

    {
        utility::print_par(world, "Compute CC Term Without DF \n");
        auto C_bar_ijab = f12::convert_C_ijab(C_ijab_nodf, occ, *orbital_energy_);
        B_ijij_ijji_nodf("i1,j1,i2,j2") = (C_ijab_nodf("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

        double E_cc = B_ijij_ijji_nodf("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_CC: ", E_cc, "\n");
        E += E_cc;
    }

    utility::print_par(world, "E_F12: ", E, "\n");
    return E;
}

}// end of namespce f12
} // end of namespace mpqc



#endif //MPQC_MP2F12_H
