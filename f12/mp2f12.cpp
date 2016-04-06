//
// Created by Chong Peng on 3/31/16.
//

#include "mp2f12.h"
#include "f12_utility.h"
#include "../utility/cc_utility.h"

void mpqc::f12::MP2F12::compute_mp2_f12_c_df() {

    auto& world = mo_int_.get_world();

    double E = 0.0;

    auto& mo_integral = mo_int_;
    auto& ao_integral = mo_int_.atomic_integral();

    auto occ = tre_->get_occ();

    TArray t2;
    {
        utility::print_par(world, "Compute T_abij With DF \n" );

        TArray g_iajb;
        g_iajb = mo_int_.compute(L"(i a|G|j b)[df]");
        g_iajb("a,b,i,j") = g_iajb("i,a,j,b");
        t2 = mpqc::cc::d_abij(g_iajb,orbital_energy_,occ);

    }

    // create shape
    auto occ_tr1 = tre_->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    //compute V term
    TArray V_ijij_ijji;
    {
        utility::print_par(world, "Compute V_ijij_ijji With DF \n" );
        V_ijij_ijji("i1,j1,i2,j2") = (mo_integral(L"(Κ |GR|i2 i1)")*mo_integral(L"(Κ|GR|Λ)[inv]")*mo_integral(L"(Λ |GR|j1 j2)")).set_shape(ijij_ijji_shape);

        // all types of GR integral not needed
        mo_int_.remove_operation_all(world, L"GR");

//        std::cout << V_ijij_ijji << std::endl;
        V_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|G|j1 q)[df]")*mo_integral(L"(i2 p|R|j2 q)[df]")).set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_ijji << std::endl;
        V_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 m|G|j1 a')[df]")*mo_integral(L"(i2 m|R|j2 a')[df]")).set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_ijji << std::endl;
        V_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|G|i1 a')[df]")*mo_integral(L"(j2 m|R|i2 a')[df]")).set_shape(ijij_ijji_shape);
//        std::cout << V_ijij_ijji << std::endl;

        // G integral in MO not needed, still need G integral in AO to compute F, K, hJ
        mo_int_.registry().remove_operation(world, L"G");


        //contribution from V_ijij_ijji
        double E_v = V_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(1.0,2.5,-0.5));
        utility::print_par(world, "E_V: ", E_v, "\n");
        E += E_v;
    }

    // compute C term
    TArray C_iajb;
    {
        utility::print_par(world, "Compute C_iajb With DF \n" );
        C_iajb("i,a,j,b") = mo_integral(L"(i a|R|j a')[df]")*mo_integral(L"(b|F|a')[df]");
        C_iajb("i,a,j,b") += mo_integral(L"(j b|R|i a')[df]")*mo_integral(L"(a|F|a')[df]");
    }


    {
        utility::print_par(world, "Compute CT With DF \n" );
        V_ijij_ijji("i1,j1,i2,j2") = (C_iajb("i1,a,j1,b")*t2("a,b,i2,j2")).set_shape(ijij_ijji_shape);

        double E_ct = V_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(1.0,2.5,-0.5));
        utility::print_par(world, "E_CT: ", E_ct, "\n");
//        E += E_ct;
    }


    // compute X term
    TArray X_ijij_ijji;
    {
        utility::print_par(world, "Compute X_ijij_ijji With DF \n" );
        X_ijij_ijji("i1,j1,i2,j2") = (mo_integral(L"(Κ |R2|i1 i2)")*ao_integral(L"(Κ|R2|Λ)[inv]")*mo_integral(L"(Λ |R2|j1 j2)")).set_shape(ijij_ijji_shape);


//        std::cout << X_ijij_ijji << std::endl;
        X_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 q)[df]")*mo_integral(L"(i2 p|R|j2 q)[df]")).set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_ijji << std::endl;
        X_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 m|R|j1 a')[df]")*mo_integral(L"(i2 m|R|j2 a')[df]")).set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_ijji << std::endl;
        X_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|R|i1 a')[df]")*mo_integral(L"(j2 m|R|i2 a')[df]")).set_shape(ijij_ijji_shape);
//        std::cout << X_ijij_ijji << std::endl;

        // R_ipjq not needed
        mo_int_.registry().remove_formula(world, L"(i1 p|R|j1 q)[df]");

        auto Fij = mo_int_.compute(L"(i|F|j)[df]");
        auto Fij_eigen = array_ops::array_to_eigen(Fij);
        f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

        double E_x = -X_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_X: ", E_x, "\n");
        E += E_x;

    }

    // compute B term
    TArray B_ijij_ijji;
    {

        utility::print_par(world, "Compute B_ijij_ijji With DF \n");

        B_ijij_ijji("i1,j1,i2,j2") = (mo_integral(L"(Κ |dR2|i1 i2)")*ao_integral(L"(Κ|dR2|Λ)[inv]")*mo_integral(L"(Λ |dR2|j1 j2)")).set_shape(ijij_ijji_shape);

        mo_int_.remove_operation_all(world, L"dR2");
//        std::cout << B_ijij_ijji << std::endl;
        auto hJ = mo_int_.compute(L"(P' | hJ | i)[df]");
        B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(i1 P'|R2|j1 j2)[df]")*hJ("P',i2")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(j1 P'|R2|i1 i2)[df]")*hJ("P',j2")).set_shape(ijij_ijji_shape);

        mo_int_.remove_operation_all(world, L"R2");
        mo_int_.remove_operation_all(world, L"hJ");
//        std::cout << B_ijij_ijji << std::endl;

        B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(i2 R'|R|j2 Q')[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(j2 R'|R|i2 Q')[df]")).set_shape(ijij_ijji_shape);

        // AO R integral not needed
        mo_int_.atomic_integral().registry().remove_operation(world, L"R");

//        std::cout << B_ijij_ijji << std::endl;
        B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(i2 R'|R|j2 m)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(j2 R'|R|i2 m)[df]")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_ijji << std::endl;

        B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(i2 P'|R|j2 b')[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(j2 P'|R|i2 b')[df]")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_ijji << std::endl;
        // P' doesn't appear later
        mo_int_.registry().remove_orbital(world, L"P'");

        B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(i2 r|R|j2 a)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(j2 r|R|i2 a)[df]")).set_shape(ijij_ijji_shape);


//        std::cout << B_ijij_ijji << std::endl;
        B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(i2 n|R|j2 b')[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(j2 n|R|i2 b')[df]")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_ijji << std::endl;


        B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(j2 a|R|i2 a')[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(i2 a|R|j2 a')[df]")).set_shape(ijij_ijji_shape);

//        std::cout << B_ijij_ijji << std::endl;

        double E_b = B_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_B: ", E_b, "\n");
        E += E_b;
    }

    {

        utility::print_par(world, "Compute CC Term With DF \n");
        auto C_bar_iajb = f12::convert_C_iajb(C_iajb, occ, orbital_energy_);
        B_ijij_ijji("i1,j1,i2,j2") = (C_iajb("i1,a,j1,b")*C_bar_iajb("i2,a,j2,b")).set_shape(ijij_ijji_shape);

        double E_cc = B_ijij_ijji("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_CC: ", E_cc, "\n");
//        E += E_cc;
    }

    utility::print_par(world, "E_F12: ", E, "\n");
}


void mpqc::f12::MP2F12::compute_mp2_f12_c(){

    auto& world = mo_int_.get_world();

    double E = 0.0;

    auto& mo_integral = mo_int_;
    auto& ao_integral = mo_int_.atomic_integral();

    auto occ = tre_->get_occ();

    TArray t2_nodf;
    {
        utility::print_par(world, "Compute T_abij Without DF \n" );
        TArray g_iajb;
        g_iajb = mo_integral.compute(L"(i a|G|j b)");
        g_iajb("a,b,i,j") = g_iajb("i,a,j,b");
        t2_nodf = mpqc::cc::d_abij(g_iajb,orbital_energy_,occ);
    }

    // create shape
    auto occ_tr1 = tre_->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    TArray V_ijij_ijji_nodf;
    {

        utility::print_par(world, "Compute V_ijij_ijji Without DF \n" );
        // all types of GR integral not needed
        V_ijij_ijji_nodf("i1,j1,i2,j2") = mo_integral(L"(i1 i2|GR|j1 j2)").set_shape(ijij_ijji_shape);
        mo_int_.remove_operation_all(world, L"GR");
        V_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|G|j1 q)")*mo_integral(L"(i2 p|R|j2 q)")).set_shape(ijij_ijji_shape);
        V_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 m|G|j1 a')")*mo_integral(L"(i2 m|R|j2 a')")).set_shape(ijij_ijji_shape);
        V_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|G|i1 a')")*mo_integral(L"(j2 m|R|i2 a')")).set_shape(ijij_ijji_shape);

        // G integral in MO not needed, still need G integral in AO to compute F, K, hJ
        mo_int_.registry().remove_operation(world, L"G");
        double E_v = V_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(1.0,2.5,-0.5));
        utility::print_par(world, "E_V: ", E_v, "\n");
        E += E_v;
    }

    TArray C_iajb_nodf;
    {
        utility::print_par(world, "Compute C_iajb Without DF \n" );
        C_iajb_nodf("i,a,j,b") = mo_integral(L"(i a|R|j a')")*mo_integral(L"(b|F|a')");
        C_iajb_nodf("i,a,j,b") += mo_integral(L"(j b|R|i a')")*mo_integral(L"(a|F|a')");
    }

    {
        utility::print_par(world, "Compute CT Without DF \n" );
        V_ijij_ijji_nodf("i1,j1,i2,j2") = (C_iajb_nodf("i1,a,j1,b")*t2_nodf("a,b,i2,j2")).set_shape(ijij_ijji_shape);

        double E_ct = V_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(1.0,2.5,-0.5));
        utility::print_par(world, "E_CT: ", E_ct, "\n");
    }

    TArray X_ijij_ijji_nodf;
    {
        utility::print_par(world, "Compute X_ijij_ijji Without DF \n" );
        X_ijij_ijji_nodf("i1,j1,i2,j2") = mo_integral(L"(i1 i2 |R2|j1 j2)").set_shape(ijij_ijji_shape);
        X_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 q)")*mo_integral(L"(i2 p|R|j2 q)")).set_shape(ijij_ijji_shape);
        X_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 m|R|j1 a')")*mo_integral(L"(i2 m|R|j2 a')")).set_shape(ijij_ijji_shape);
        X_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|R|i1 a')")*mo_integral(L"(j2 m|R|i2 a')")).set_shape(ijij_ijji_shape);

        // R_ipjq not needed
        mo_int_.registry().remove_formula(world, L"(i1 p|R|j1 q)[df]");

        // compute energy contribution
        auto Fij = mo_integral.compute(L"(i|F|j)");
        auto Fij_eigen = array_ops::array_to_eigen(Fij);
        f12::convert_X_ijkl(X_ijij_ijji_nodf, Fij_eigen);

        double E_x = -X_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_X: ", E_x, "\n");
        E += E_x;
    }

    TArray B_ijij_ijji_nodf;
    {
        utility::print_par(world, "Compute B_ijij_ijji Without DF \n");

        B_ijij_ijji_nodf("i1,j1,i2,j2") = (mo_integral(L"(i1 i2 |dR2|j1 j2)")).set_shape(ijij_ijji_shape);
        mo_int_.remove_operation_all(world, L"dR2");

        auto hJ = mo_integral.compute(L"(P' | hJ | i)");
        B_ijij_ijji_nodf("i1,j1,i2,j2") += (mo_integral(L"(i1 P'|R2|j1 j2)")*hJ("P',i2")).set_shape(ijij_ijji_shape);
        B_ijij_ijji_nodf("i1,j1,i2,j2") += (mo_integral(L"(j1 P'|R2|i1 i2)")*hJ("P',j2")).set_shape(ijij_ijji_shape);

        mo_int_.remove_operation_all(world, L"R2");
        mo_int_.remove_operation_all(world, L"hJ");

        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 Q')")*mo_integral(L"(P'|K|R')")*mo_integral(L"(i2 R'|R|j2 Q')")).set_shape(ijij_ijji_shape);
        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 Q')")*mo_integral(L"(P'|K|R')")*mo_integral(L"(j2 R'|R|i2 Q')")).set_shape(ijij_ijji_shape);

        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 m)")*mo_integral(L"(P'|F|R')")*mo_integral(L"(i2 R'|R|j2 m)")).set_shape(ijij_ijji_shape);
        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 m)")*mo_integral(L"(P'|F|R')")*mo_integral(L"(j2 R'|R|i2 m)")).set_shape(ijij_ijji_shape);

        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 m|R|j1 b')")*mo_integral(L"(m|F|P')")*mo_integral(L"(i2 P'|R|j2 b')")).set_shape(ijij_ijji_shape);
        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 m|R|i1 b')")*mo_integral(L"(m|F|P')")*mo_integral(L"(j2 P'|R|i2 b')")).set_shape(ijij_ijji_shape);
        // P' doesn't appear later
        mo_int_.registry().remove_orbital(world, L"P'");

        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 a)")*mo_integral(L"(p|F|r)")*mo_integral(L"(i2 r|R|j2 a)")).set_shape(ijij_ijji_shape);
        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (mo_integral(L"(j1 p|R|i1 a)")*mo_integral(L"(p|F|r)")*mo_integral(L"(j2 r|R|i2 a)")).set_shape(ijij_ijji_shape);

        B_ijij_ijji_nodf("i1,j1,i2,j2") += (mo_integral(L"(i1 m|R|j1 b')")*mo_integral(L"(m|F|n)")*mo_integral(L"(i2 n|R|j2 b')")).set_shape(ijij_ijji_shape);
        B_ijij_ijji_nodf("i1,j1,i2,j2") += (mo_integral(L"(j1 m|R|i1 b')")*mo_integral(L"(m|F|n)")*mo_integral(L"(j2 n|R|i2 b')")).set_shape(ijij_ijji_shape);

        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 p|R|j1 a)")*mo_integral(L"(p|F|a')")*mo_integral(L"(j2 a|R|i2 a')")).set_shape(ijij_ijji_shape);
        B_ijij_ijji_nodf("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 p|R|i1 a)")*mo_integral(L"(p|F|a')")*mo_integral(L"(i2 a|R|j2 a')")).set_shape(ijij_ijji_shape);


        double E_b = B_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_B: ", E_b, "\n");
        E += E_b;
    }

    {
        utility::print_par(world, "Compute CC Term Without DF \n");
        auto C_bar_iajb = f12::convert_C_iajb(C_iajb_nodf, occ, orbital_energy_);
        B_ijij_ijji_nodf("i1,j1,i2,j2") = (C_iajb_nodf("i1,a,j1,b")*C_bar_iajb("i2,a,j2,b")).set_shape(ijij_ijji_shape);

        double E_cc = B_ijij_ijji_nodf("i1,j1,i2,j2").reduce(MP2F12Energy(0.25,0.4375,0.0625));
        utility::print_par(world, "E_CC: ", E_cc, "\n");
//        E += E_cc;
    }

    utility::print_par(world, "E_F12: ", E, "\n");

}
