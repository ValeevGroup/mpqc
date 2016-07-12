//
// Created by Chong Peng on 3/31/16.
//

#ifndef MPQC_MP2F12_H
#define MPQC_MP2F12_H

#include <string>

#include <mpqc/chemistry/qc/mbpt/mp2.h>
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

    /// constructor using MO Integral with orbitals computed
    MP2F12(std::shared_ptr<mbpt::MP2<Tile,Policy>>& mp2) : mo_int_(mp2->mo_integral()), mp2_(mp2)
    { }

    MP2F12(std::shared_ptr<mbpt::MP2<Tile,Policy>>&& mp2) : mo_int_(mp2->mo_integral()), mp2_(mp2)
    { }

    /// constructfor using MO Integral without orbitals computed
    MP2F12(MolecularIntegralClass& mo_int) : mo_int_(mo_int){

        mp2_ = std::make_shared<mbpt::MP2<Tile,Policy>>(mo_int_);

    }

    virtual double compute(const rapidjson::Document& in){

        auto& world = mo_int_.get_world();
        // mp2 time
        auto mp2_time0 = mpqc_time::fenced_now(world);

        double mp2_e = mp2_->compute(in);

        auto mp2_time1 = mpqc_time::fenced_now(world);
        auto mp2_time = mpqc_time::duration_in_s(mp2_time0, mp2_time1);
        mpqc::utility::print_par(world, "Total MP2 Time:  ", mp2_time, "\n");


        auto f12_time0 = mpqc_time::fenced_now(world);
        // solve cabs orbitals
        auto& ao_int = mo_int_.atomic_integral();
        auto orbital_registry = mo_int_.orbital_space();
        closed_shell_cabs_mo_build_svd(ao_int, *orbital_registry, in, mp2_->trange1_engine());

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
        auto f12_time1 = mpqc_time::fenced_now(world);
        auto f12_time = mpqc_time::duration_in_s(f12_time0, f12_time1);
        mpqc::utility::print_par(world, "Total F12 Time:  ", f12_time, "\n");

        return mp2_e + mp2_f12_energy;
    }

private:

  /// MP2-F12 C approach with density fitting
  double compute_mp2_f12_c_df();

  /// MP2-F12 C approach
  double compute_mp2_f12_c();

protected:
    struct MP2F12Energy {
        using result_type =  double;
        using argument_type =  Tile;

        double iiii;
        double ijij;
        double ijji;

        MP2F12Energy() = default;
        MP2F12Energy(MP2F12Energy const &) = default;
        MP2F12Energy(double c1, double c2, double c3) : iiii(c1), ijij(c2), ijji(c3) {}

        result_type operator()() const { return 0.0; }

        result_type operator()(result_type const &t) const { return t; }

        void operator()(result_type &me, result_type const &other) const {
            me += other;
        }

        void operator()(result_type &me, argument_type const &tile) const {
            auto const &range = tile.range();

            TA_ASSERT(range.rank() == 4);

            auto const st = range.lobound_data();
            auto const fn = range.upbound_data();
            auto tile_idx = 0;

            auto sti = st[0];
            auto fni = fn[0];
            auto stj = st[1];
            auto fnj = fn[1];
            auto stk = st[2];
            auto fnk = fn[2];
            auto stl = st[3];
            auto fnl = fn[3];

            for (auto i = sti; i < fni; ++i) {
                for (auto j = stj; j < fnj; ++j) {
                    for (auto k = stk; k < fnk; ++k) {
                        for (auto l = stl; l < fnl; ++l, ++tile_idx) {
                            // iiii
                            if( (i==j) && (i==k) && (k==l)){
                                me += iiii*tile.data()[tile_idx];
                            }
                            // ijij
                            else if( j > i && k==i && l==j){
                                me += ijij*tile.data()[tile_idx];
                            }
                            // ijji
                            else if( j > i && l==i && k==j){
                                me += ijji*tile.data()[tile_idx];
                            }
                        }
                    }
                }
            }
        }
    };

protected:

    MolecularIntegralClass& mo_int_;
    std::shared_ptr<mbpt::MP2<Tile,Policy>> mp2_;
};



template <typename Tile>
double MP2F12<Tile>::compute_mp2_f12_c_df() {

    auto& world = mo_int_.get_world();

    double E = 0.0;

    auto& mo_integral = mo_int_;

    auto occ = mp2_->trange1_engine()->get_actual_occ();

    TArray t2;
    {
        utility::print_par(world, "Compute T_abij With DF \n" );

        TArray g_iajb;
        g_iajb = mo_int_.compute(L"<i j|G|a b>[df]");
        g_iajb("a,b,i,j") = g_iajb("i,j,a,b");
        t2 = mpqc::cc::d_abij(g_iajb,*(mp2_->orbital_energy()),occ);

    }

    // create shape
    auto occ_tr1 = mp2_->trange1_engine()->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    //compute V term
    TArray V_ijij_ijji = compute_V_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    {

        // G integral in MO not needed, still need G integral in AO to compute F, K, hJ
        mo_int_.registry().remove_operation(world, L"G");

        //contribution from V_ijij_ijji
        double E_v = V_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(1.0,2.5,-0.5));
        utility::print_par(world, "E_V: ", E_v, "\n");
        E += E_v;
    }

    // compute C term
    TArray C_ijab = compute_C_ijab_df(mo_integral);

    {
        utility::print_par(world, "Compute CT With DF \n" );
        V_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

        double E_ct = V_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(1.0,2.5,-0.5));
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

        double E_x = -X_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_X: ", E_x, "\n");
        E += E_x;

    }

    // compute B term
    TArray B_ijij_ijji = compute_B_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    {
        double E_b = B_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_B: ", E_b, "\n");
        E += E_b;
    }

    {

        utility::print_par(world, "Compute CC Term With DF \n");
        auto C_bar_ijab = f12::convert_C_ijab(C_ijab, occ, *(mp2_->orbital_energy()));
        B_ijij_ijji("i1,j1,i2,j2") = (C_ijab("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

        double E_cc = B_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
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

    auto occ = mp2_->trange1_engine()->get_actual_occ();

    TArray t2_nodf;
    {
        utility::print_par(world, "Compute T_abij Without DF \n" );
        TArray g_iajb;
        g_iajb = mo_integral.compute(L"<i j|G|a b>");
        g_iajb("a,b,i,j") = g_iajb("i,j,a,b");
        t2_nodf = mpqc::cc::d_abij(g_iajb,*(mp2_->orbital_energy()),occ);
    }

    // create shape
    auto occ_tr1 = mp2_->trange1_engine()->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    TArray V_ijij_ijji_nodf = compute_V_ijij_ijji(mo_integral, ijij_ijji_shape);
    {
        double E_v = V_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(1.0,2.5,-0.5));
        utility::print_par(world, "E_V: ", E_v, "\n");
        E += E_v;
    }

    TArray C_ijab_nodf = compute_C_ijab(mo_integral);

    {
        utility::print_par(world, "Compute CT Without DF \n" );
        V_ijij_ijji_nodf("i1,j1,i2,j2") = (C_ijab_nodf("i1,j1,a,b")*t2_nodf("a,b,i2,j2")).set_shape(ijij_ijji_shape);

        double E_ct = V_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(1.0,2.5,-0.5));
        utility::print_par(world, "E_CT: ", E_ct, "\n");
        E += E_ct;
    }

    TArray X_ijij_ijji_nodf = compute_X_ijij_ijji(mo_integral, ijij_ijji_shape);
    {

        // compute energy contribution
        auto Fij = mo_integral.compute(L"<i|F|j>");
        auto Fij_eigen = array_ops::array_to_eigen(Fij);
        f12::convert_X_ijkl(X_ijij_ijji_nodf, Fij_eigen);

        double E_x = -X_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_X: ", E_x, "\n");
        E += E_x;
    }

    TArray B_ijij_ijji_nodf = compute_B_ijij_ijji(mo_integral,ijij_ijji_shape);
    {

        double E_b = B_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_B: ", E_b, "\n");
        E += E_b;
    }

    {
        utility::print_par(world, "Compute CC Term Without DF \n");
        auto C_bar_ijab = f12::convert_C_ijab(C_ijab_nodf, occ, *(mp2_->orbital_energy()));
        B_ijij_ijji_nodf("i1,j1,i2,j2") = (C_ijab_nodf("i1,j1,a,b")*C_bar_ijab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

        double E_cc = B_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_CC: ", E_cc, "\n");
        E += E_cc;
    }

    utility::print_par(world, "E_F12: ", E, "\n");
    return E;
}

}// end of namespce f12
} // end of namespace mpqc



#endif //MPQC_MP2F12_H
