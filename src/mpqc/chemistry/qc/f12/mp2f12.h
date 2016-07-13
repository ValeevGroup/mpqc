//
// Created by Chong Peng on 3/31/16.
//

#ifndef MPQC_MP2F12_H
#define MPQC_MP2F12_H

#include <string>

#include <mpqc/chemistry/qc/f12/f12_utility.h>
#include <mpqc/chemistry/qc/mbpt/mp2.h>
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

    MP2F12(MolecularIntegralClass& mo_int) : mo_int_(mo_int)
    {
        mp2_ = std::make_shared<mbpt::MP2<Tile,Policy>>(mo_int);
    }

    MP2F12(std::shared_ptr<mbpt::MP2<Tile,Policy>> mp2) : mo_int_(mp2->mo_integral()), mp2_(mp2) {}

    std::tuple<double, double> compute_mp2_f12_c_df();

    std::tuple<double, double> compute_mp2_f12_c();

    double compute(const rapidjson::Document& in){

        auto& world = this->mo_int_.get_world();
        auto f12_time0 = mpqc_time::fenced_now(world);
        // solve mo orbitals
        mp2_->init(in);

        // solve cabs orbitals
        auto& ao_int = this->mo_int_.atomic_integral();
        auto orbital_registry = this->mo_int_.orbital_space();
        closed_shell_cabs_mo_build_svd(this->mo_int_,in,this->mp2_->trange1_engine());

        std::string method = in.HasMember("Method") ? in["Method"].GetString() : "df";

        double mp2_energy, f12_energy;

        if(method == "four center"){
          std::tie(mp2_energy, f12_energy) = compute_mp2_f12_c();
        }
        else if(method == "df"){
          std::tie(mp2_energy, f12_energy) = compute_mp2_f12_c_df();
        }
        else{
            throw std::runtime_error("Wrong MP2F12 Method");
        }

        utility::print_par(mo_int_.get_world(), "E_MP2: ", mp2_energy, "\n");
        utility::print_par(mo_int_.get_world(), "E_F12: ", f12_energy, "\n");

        auto f12_time1 = mpqc_time::fenced_now(world);
        auto f12_time = mpqc_time::duration_in_s(f12_time0, f12_time1);
        mpqc::utility::print_par(world, "Total MP2F12 Time:  ", f12_time, "\n");

        return mp2_energy + f12_energy;
    }

protected:

    MolecularIntegralClass& mo_int_;
    std::shared_ptr<mbpt::MP2<Tile,Policy>> mp2_;
};

template <typename Tile>
std::tuple<double,double> MP2F12<Tile>::compute_mp2_f12_c_df() {

    auto& world = mo_int_.get_world();

    double E_MP2, E_F12 = 0.0;

    auto& mo_integral = mo_int_;

    auto occ = mp2_->trange1_engine()->get_active_occ();

    // create shape
    auto occ_tr1 = mp2_->trange1_engine()->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    TArray t2;
    {
        utility::print_par(world, "Compute T_abij With DF \n" );

        TArray g_abij;
        g_abij("a,b,i,j") = mo_int_.compute(L"<i j|G|a b>[df]")("i,j,a,b");
        t2 = mpqc::cc::d_abij(g_abij,*(mp2_->orbital_energy()),occ);

        // compute MP2 energy and pair energies
        TArray TG_ijij_ijji;
        TG_ijij_ijji("i1,j1,i2,j2") =
            (t2("a,b,i1,j1") * g_abij("a,b,i2,j2"))
                .set_shape(ijij_ijji_shape);
        E_MP2 = TG_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(2, -1));
    }

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
        double E_v = V_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(2 * C_ijij_bar,2 * C_ijji_bar));
        utility::print_par(world, "E_V: ", E_v, "\n");
        E_F12 += E_v;
    }

    // compute C term
    TArray C_ijab = compute_C_ijab_df(mo_integral);

    {
        utility::print_par(world, "Compute CT With DF \n" );
        V_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

        // NB factor of 2 from the Hylleraas functional
        double E_ct = V_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(2 * C_ijij_bar,2 * C_ijji_bar));
        utility::print_par(world, "E_CT: ", E_ct, "\n");
        E_F12 += E_ct;
    }


    // compute X term
    TArray X_ijij_ijji = compute_X_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    {

        // R_ipjq not needed
        mo_int_.registry().remove_formula(world, L"<i1 j1|R|p q>[df]");

        auto Fij = mo_int_.compute(L"<i|F|j>[df]");
        auto Fij_eigen = array_ops::array_to_eigen(Fij);
        f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

        double E_x = -X_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_X: ", E_x, "\n");
        E_F12 += E_x;

    }

    // compute B term
    TArray B_ijij_ijji = compute_B_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    {
        double E_b = B_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_B: ", E_b, "\n");
        E_F12 += E_b;
    }

    {

        utility::print_par(world, "Compute CC Term With DF \n");
        auto C_bar_ijab = f12::convert_C_ijab(C_ijab, occ, *(mp2_->orbital_energy()));
        B_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

        double E_cc = B_ijij_ijji("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_CC: ", E_cc, "\n");
        E_F12 += E_cc;
    }

    return std::make_tuple(E_MP2, E_F12);
}


template <typename Tile>
std::tuple<double,double> MP2F12<Tile>::compute_mp2_f12_c(){

    auto& world = mo_int_.get_world();

    double E_MP2 = 0.0, E_F12 = 0.0;

    auto& mo_integral = mo_int_;

    auto occ = mp2_->trange1_engine()->get_active_occ();

    // create shape
    auto occ_tr1 = mp2_->trange1_engine()->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    TArray t2_nodf;  // t2_abij
    {
        utility::print_par(world, "Compute T_abij Without DF \n" );
        TArray g_abij;
        g_abij("a,b,i,j") = mo_integral.compute(L"<i j|G|a b>")("i,j,a,b");
        t2_nodf = mpqc::cc::d_abij(g_abij,*(mp2_->orbital_energy()),occ);
        TArray TG_ijij_ijji_nodf;
        TG_ijij_ijji_nodf("i1,j1,i2,j2") =
            (t2_nodf("a,b,i1,j1") * g_abij("a,b,i2,j2"))
                .set_shape(ijij_ijji_shape);
        E_MP2 = TG_ijij_ijji_nodf("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(2, -1));
    }

    TArray V_ijij_ijji_nodf = compute_V_ijij_ijji(mo_integral, ijij_ijji_shape);
    {
        double E_v = V_ijij_ijji_nodf("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(2 * C_ijij_bar, 2 * C_ijji_bar));
        utility::print_par(world, "E_V: ", E_v, "\n");
        E_F12 += E_v;
    }

    TArray C_ijab_nodf = compute_C_ijab(mo_integral);

    {
        utility::print_par(world, "Compute CT Without DF \n" );
        V_ijij_ijji_nodf("i1,j1,i2,j2") = (C_ijab_nodf("i1,j1,a,b")*t2_nodf("a,b,i2,j2")).set_shape(ijij_ijji_shape);

        double E_ct = V_ijij_ijji_nodf("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(2 * C_ijij_bar, 2 * C_ijji_bar));
        utility::print_par(world, "E_CT: ", E_ct, "\n");
        E_F12 += E_ct;
    }

    TArray X_ijij_ijji_nodf = compute_X_ijij_ijji(mo_integral, ijij_ijji_shape);
    {

        // compute energy contribution
        auto Fij = mo_integral.compute(L"<i|F|j>");
        auto Fij_eigen = array_ops::array_to_eigen(Fij);
        f12::convert_X_ijkl(X_ijij_ijji_nodf, Fij_eigen);

        double E_x = -X_ijij_ijji_nodf("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_X: ", E_x, "\n");
        E_F12 += E_x;
    }

    TArray B_ijij_ijji_nodf = compute_B_ijij_ijji(mo_integral,ijij_ijji_shape);
    {

        double E_b = B_ijij_ijji_nodf("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_B: ", E_b, "\n");
        E_F12 += E_b;
    }

    {
        utility::print_par(world, "Compute CC Term Without DF \n");
        auto C_bar_ijab = f12::convert_C_ijab(C_ijab_nodf, occ, *(mp2_->orbital_energy()));
        B_ijij_ijji_nodf("i1,j1,i2,j2") = (C_ijab_nodf("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

        double E_cc = B_ijij_ijji_nodf("i1,j1,i2,j2").reduce(F12EnergyReductor<Tile>(CC_ijij_bar,CC_ijji_bar));
        utility::print_par(world, "E_CC: ", E_cc, "\n");
        E_F12 += E_cc;
    }

    return std::make_tuple(E_MP2, E_F12);
}

}// end of namespace f12
} // end of namespace mpqc



#endif //MPQC_MP2F12_H
