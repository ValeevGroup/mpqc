//
// Created by Chong Peng on 7/1/15.
//

#ifndef TILECLUSTERCHEM_CCSD_H
#define TILECLUSTERCHEM_CCSD_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"

#include "trange1_engine.h"
#include "mo_block.h"
#include "../ta_routines/tarray_block.h"

#include "ccsd_intermediates.h"
#include "./diis_ccsd.h"
#include "integral_generator.h"
#include "lazy_integral.h"


namespace tcc {
    namespace cc {

        // CCSD class that computed CCSD energy

        template<typename Tile, typename Policy>
        class CCSD {

        public:
            typedef TA::Array <double, 2, Tile, Policy> TArray2;
            typedef TA::Array <double, 3, Tile, Policy> TArray3;
            typedef TA::Array <double, 4, Tile, Policy> TArray4;

            typedef tcc::TArrayBlock<double, 2, Tile, Policy, tcc::MOBlock> TArrayBlock2;
            typedef tcc::TArrayBlock<double, 4, Tile, Policy, tcc::MOBlock> TArrayBlock4;


            CCSD(const TArray2 &fock, const Eigen::VectorXd &ens,
                 const std::shared_ptr<TRange1Engine> &tre,
                 const std::shared_ptr<CCSDIntermediate<Tile, Policy>> &inter) :
                    ens_(ens), tre_(tre), intermediate_(inter)
            {
                auto mo_block = std::make_shared<tcc::MOBlock>(*tre_);
                fock_ = TArrayBlock2(fock, mo_block);
            }

            // compute function
            virtual void compute(){

                TArray2 t1;
                TArray4 t2;

                double ccsd_corr = compute_ccsd(t1, t2);

                intermediate_->clean_two_electron();

            }


            // dummy way of doing CCSD
            // store all the integrals in memory
            // used as reference for development
            double compute_ccsd_dummy(TArray2& t1, TArray4& t2) {

                auto n_occ = tre_->get_actual_occ();

                TArray4 g_abij = intermediate_->get_abij();

                TArray2 f_ai;
                f_ai("a,i") = fock_("a,i");

                g_abij.get_world().gop.fence();

//      std::cout << g_abij << std::endl;

                TArray2 d1(f_ai.get_world(), f_ai.trange(), f_ai.get_shape(),
                           f_ai.get_pmap());
                // store d1 to local
                d_ai(d1, ens_, n_occ);

                TArray4 d2(g_abij.get_world(), g_abij.trange(),
                           g_abij.get_shape(), g_abij.get_pmap());
                // store d2 distributed
                d_abij(d2, ens_, n_occ);

                t1("a,i") = f_ai("a,i") * d1("a,i");
                t2("a,b,i,j") = g_abij("a,b,i,j") * d2("a,b,i,j");

//      std::cout << t1 << std::endl;
//      std::cout << t2 << std::endl;
                TArray4 tau;
                tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

                double E0 = 0.0;
                double E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
                            TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"),
                                    tau("a,b,i,j"));
                double dE = std::abs(E1 - E0);
//      std::cout << E1 << std::endl;

                // get all two electron integrals
                TArray4 g_ijkl = intermediate_->get_ijkl();
                TArray4 g_abcd = intermediate_->get_abcd();
                TArray4 g_iajb = intermediate_->get_iajb();
                TArray4 g_iabc = intermediate_->get_iabc();
                TArray4 g_aibc = intermediate_->get_aibc();
                TArray4 g_ijak = intermediate_->get_ijak();
                TArray4 g_ijka = intermediate_->get_ijka();

//                intermediate_->clean();

                if (g_abij.get_world().rank() == 0) {
                    std::cout << "start iteration" << std::endl;
                }
                //optimize t1 and t2
                std::size_t iter = 0ul;
                double error = 1.0;
                TArray2 r1;
                TArray4 r2;
                TA::DIIS <tcc::cc::T1T2<double, Tile, Policy>> diis(1);
                while ((dE >= 1.0e-7 || error >= 1e-7)) {

                    //start timer
                    auto t0 = tcc::tcc_time::now();

                    // intermediates for t1
                    // external index i and a
                    // vir index a b c d
                    // occ index i j k l
                    TArray2 h_ac, h_ki, h_kc;
                    {
                        h_ac("a,c") = -(2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * tau("a,d,k,l");

                        h_ki("k,i") = (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * tau("c,d,i,l");

                        h_kc("k,c") = f_ai("c,k") + (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * t1("d,l");
                    }

                    // compute residual r1(n) = t1(n+1) - t1(n)
                    // external index i and a
                    {
                        r1("a,i") = -t1("a,i") + d1("a,i") * (
                                //
                                f_ai("a,i") - 2.0 * f_ai("c,k") * t1("a,k") * t1("c,i")
                                //
                                + h_ac("a,c") * t1("c,i") - t1("a,k") * h_ki("k,i")
                                //
                                + h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k") + t1("c,i") * t1("a,k"))
                                //
                                + (2.0 * g_abij("c,a,k,i") - g_iajb("k,a,i,c")) * t1("c,k")
                                //
                                + (2.0 * g_iabc("k,a,c,d") - g_iabc("k,a,d,c")) * tau("c,d,k,i")
                                //
                                - (2.0 * g_ijak("k,l,c,i") - g_ijak("l,k,c,i")) * tau("c,a,k,l")
                        );
                    }

//        g_abij.get_world().gop.fence();
                    // intermediates for t2
                    // external index i j a b

                    TArray4 a_klij, b_abij, j_akic, k_kaic, T;
                    TArray2 g_ki, g_ac;
                    {

                        T("d,b,i,l") =
                                0.5 * t2("d,b,i,l") + t1("d,i") * t1("b,l");

                        a_klij("k,l,i,j") = g_ijkl("k,l,i,j")
                                            + g_ijka("k,l,i,c") * t1("c,j") +
                                            g_ijak("k,l,c,j") * t1("c,i")
                                            +
                                            g_abij("c,d,k,l") * tau("c,d,i,j");

                        b_abij("a,b,i,j") = g_abcd("a,b,c,d") * tau("c,d,i,j")
                                            -
                                            g_aibc("a,k,c,d") * tau("c,d,i,j") *
                                            t1("b,k")
                                            -
                                            g_iabc("k,b,c,d") * tau("c,d,i,j") *
                                            t1("a,k");

                        g_ki("k,i") = h_ki("k,i") + f_ai("c,k") * t1("c,i")
                                      + (2.0 * g_ijka("k,l,i,c") -
                                         g_ijka("l,k,i,c")) * t1("c,l");

                        g_ac("a,c") = h_ac("a,c") - f_ai("c,k") * t1("a,k")
                                      + (2.0 * g_aibc("a,k,c,d") -
                                         g_aibc("a,k,d,c")) * t1("d,k");

                        j_akic("a,k,i,c") = g_abij("a,c,i,k") -
                                            g_ijka("l,k,i,c") * t1("a,l")
                                            + g_aibc("a,k,d,c") * t1("d,i") -
                                            g_abij("c,d,k,l") * T("d,a,i,l")
                                            + 0.5 * (2.0 * g_abij("c,d,k,l") -
                                                     g_abij("d,c,k,l")) *
                                              t2("a,d,i,l");

                        k_kaic("k,a,i,c") = g_iajb("k,a,i,c") -
                                            g_ijka("k,l,i,c") * t1("a,l") +
                                            g_iabc("k,a,d,c") * t1("d,i") -
                                            g_abij("d,c,k,l") * T("d,a,i,l");

                    }

//        g_abij.get_world().gop.fence();
                    // compute residual r2(n) = t2(n+1) - t2(n)
                    {
                        r2("a,b,i,j") = -t2("a,b,i,j") + d2("a,b,i,j") * (
                                //
                                g_abij("a,b,i,j")
                                //
                                + a_klij("k,l,i,j") * tau("a,b,k,l")
                                //
                                + b_abij("a,b,i,j")

                                // permutation part
                                //
                                + (g_ac("a,c") * t2("c,b,i,j") -
                                   g_ki("k,i") * t2("a,b,k,j"))
                                + (g_ac("b,c") * t2("c,a,j,i") -
                                   g_ki("k,j") * t2("b,a,k,i"))

                                + (g_iabc("i,c,a,b") -
                                   g_iajb("k,b,i,c") * t1("a,k")) *
                                  t1("c,j")
                                + (g_iabc("j,c,b,a") -
                                   g_iajb("k,a,j,c") * t1("b,k")) *
                                  t1("c,i")
                                //
                                - (g_ijak("i,j,a,k") +
                                   g_abij("a,c,i,k") * t1("c,j")) *
                                  t1("b,k")
                                - (g_ijak("j,i,b,k") +
                                   g_abij("b,c,j,k") * t1("c,i")) *
                                  t1("a,k")

                                + 0.5 * (2.0 * j_akic("a,k,i,c") -
                                         k_kaic("k,a,i,c")) *
                                  (2.0 * t2("c,b,k,j") - t2("b,c,k,j"))
                                + 0.5 * (2.0 * j_akic("b,k,j,c") -
                                         k_kaic("k,b,j,c")) *
                                  (2.0 * t2("c,a,k,i") - t2("a,c,k,i"))

                                - 0.5 * k_kaic("k,a,i,c") * t2("b,c,k,j") -
                                k_kaic("k,b,i,c") * t2("a,c,k,j")
                                - 0.5 * k_kaic("k,b,j,c") * t2("a,c,k,i") -
                                k_kaic("k,a,j,c") * t2("b,c,k,i")
                        );
                    }
//        g_abij.get_world().gop.fence();

                    t1("a,i") = t1("a,i") + r1("a,i");
                    t2("a,b,i,j") = t2("a,b,i,j") + r2("a,b,i,j");

                    tcc::cc::T1T2<double, Tile, Policy> t(t1, t2);
                    tcc::cc::T1T2<double, Tile, Policy> r(r1, r2);
                    error = r.norm() / size(t);
                    diis.extrapolate(t, r);

                    //update t1 and t2
                    t1("a,i") = t.first("a,i");
                    t2("a,b,i,j") = t.second("a,b,i,j");

                    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

                    // recompute energy
                    E0 = E1;
                    E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i"))
                         +
                         TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")),
                                 tau("a,b,i,j"));
                    dE = std::abs(E0 - E1);
                    iter += 1ul;

                    auto t1 = tcc_time::now();
                    auto duration = tcc_time::duration_in_s(t0, t1);

                    if (g_abij.get_world().rank() == 0) {
                        std::cout << iter << "  " << dE << "  " << error <<
                        "  " << E1 << "  " << duration << std::endl;
                    }

                    g_abij.get_world().gop.fence();
//        std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-20.10f", Delta_E)
//        << scprintf("%-15.10f", E_1) << std::endl;

                }
                if (g_abij.get_world().rank() == 0) {
                    std::cout << "CCSD Energy  " << E1 << std::endl;
                }
                return E1;
            }


            // ccsd energy for performance calculation
            double compute_ccsd(TArray2& t1, TArray4& t2) {

                auto n_occ = tre_->get_actual_occ();

                TArray4 g_abij = intermediate_->get_abij();

                TArray2 f_ai;
                f_ai("a,i") = fock_("a,i");

                g_abij.get_world().gop.fence();

//      std::cout << g_abij << std::endl;

                TArray2 d1(f_ai.get_world(), f_ai.trange(), f_ai.get_shape(),
                           f_ai.get_pmap());
                d_ai(d1, ens_, n_occ);

                TArray4 d2(g_abij.get_world(), g_abij.trange(),
                           g_abij.get_shape(), g_abij.get_pmap());
                d_abij(d2, ens_, n_occ);

                t1("a,i") = f_ai("a,i") * d1("a,i");
                t2("a,b,i,j") = g_abij("a,b,i,j") * d2("a,b,i,j");

//      std::cout << t1 << std::endl;
//      std::cout << t2 << std::endl;
                TArray4 tau;
                tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

                double E0 = 0.0;
                double E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i")) +
                            TA::dot(2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j"),
                                    tau("a,b,i,j"));
                double dE = std::abs(E1 - E0);
//      std::cout << E1 << std::endl;

                // get all two electron integrals
                TArray4 g_ijkl = intermediate_->get_ijkl();
                TArray4 g_iajb = intermediate_->get_iajb();
                TArray4 g_ijak = intermediate_->get_ijak();
                TArray4 g_ijka = intermediate_->get_ijka();

                // get three center integrals
                TArray3 Xab = intermediate_->get_Xab();
                TArray3 Xij = intermediate_->get_Xij();
                TArray3 Xai = intermediate_->get_Xai();

                // get mo coefficient
                TArray2 ca = intermediate_->get_Ca();
                TArray2 ci = intermediate_->get_Ci();

                if (g_abij.get_world().rank() == 0) {
                    std::cout << "start iteration" << std::endl;
                    std::cout << "iter " << "    deltaE    " << "            residual       "
                    << "      energy     " << "    U/second  " << " total/second "<<std::endl;
                }
                //optimize t1 and t2
                std::size_t iter = 0ul;
                double error = 1.0;
                TArray2 r1;
                TArray4 r2;
                TA::DIIS <tcc::cc::T1T2<double, Tile, Policy>> diis(1);
                while ((dE >= 1.0e-7 || error >= 1e-7)) {

                    //start timer
                    auto t0 = tcc::tcc_time::now();

                    TArray4 u2_u11;
                    // compute half transformed intermediates
                    auto tu0 = tcc::tcc_time::now();
                    {
                        u2_u11 = intermediate_->compute_u2_u11(t2, t1);
                    }
                    auto tu1 = tcc_time::now();
                    auto duration_u = tcc_time::duration_in_s(tu0, tu1);

//                    if (g_abij.get_world().rank() == 0) {
//                        std::cout << "Time to compute U intermediates   " << duration << std::endl;
//                    }


                    // intermediates for t1
                    // external index i and a
                    TArray2 h_ac, h_ki, h_kc;
                    {
                        h_ac("a,c") = -(2.0 * g_abij("c,d,k,l") - g_abij("c,d,l,k")) * tau("a,d,k,l");

                        h_ki("k,i") = (2.0 * g_abij("c,d,k,l") - g_abij("d,c,k,l")) * tau("c,d,i,l");

                        h_kc("k,c") = f_ai("c,k")
                                      + (-g_abij("d,c,k,l")+ 2.0 * g_abij("c,d,k,l")) * t1("d,l");
                    }
//        g_abij.get_world().gop.fence();

                    // compute residual r1(n) = t1(n+1) - t1(n)
                    // external index i and a
                    {
                        r1("a,i") = -t1("a,i") + d1("a,i") * (
                                //
                                f_ai("a,i")
                                - 2.0 * f_ai("c,k") * t1("c,i") * t1("a,k")
                                //
                                + h_ac("a,c") * t1("c,i") - t1("a,k") * h_ki("k,i")
                                //
                                + h_kc("k,c") * (2.0 * t2("c,a,k,i") - t2("c,a,i,k")
                                                 + t1("a,k") * t1("c,i") )
                                //
                                + (2.0 * g_abij("c,a,k,i") - g_iajb("k,a,i,c")) * t1("c,k")
                                //
                                + (2.0 * u2_u11("p,r,k,i")- u2_u11("p,r,i,k")) * ci("p,k") * ca("r,a")
                                //
                                - (2.0 * g_ijak("k,l,c,i") - g_ijak("l,k,c,i")) * tau("c,a,k,l")
                        );
                    }

//        g_abij.get_world().gop.fence();
                    // intermediates for t2
                    // external index i j a b

                    TArray4 a_klij, b_abij, j_akic, k_kaic, T;
                    TArray2 g_ki, g_ac;
                    {



                        T("d,b,i,l") =
                                0.5 * t2("d,b,i,l") + t1("d,i") * t1("b,l");

                        a_klij("k,l,i,j") = g_ijkl("k,l,i,j")

                                            + g_ijka("k,l,i,c") * t1("c,j")

                                            + g_ijak("k,l,c,j") * t1("c,i")

                                            + g_abij("c,d,k,l") * tau("c,d,i,j");

                        b_abij("a,b,i,j") = (u2_u11("p,r,i,j")*ca("r,b") - ci("r,k")*t1("b,k")*u2_u11("p,r,i,j")) * ca("p,a")

                                            - u2_u11("p,r,i,j") * ci("p,k") * ca("r,b") * t1("a,k");

                        g_ki("k,i") = h_ki("k,i") + f_ai("c,k") * t1("c,i")

                                      + (2.0 * g_ijka("k,l,i,c")

                                      - g_ijka("l,k,i,c")) * t1("c,l");

                        g_ac("a,c") = h_ac("a,c") - f_ai("c,k") * t1("a,k")

                                      + (2.0*Xai("X,d,k")*t1("d,k"))*Xab("X,a,c")

                                      - (Xab("X,a,d")*t1("d,k"))*Xai("X,c,k");

                        j_akic("a,k,i,c") = g_abij("a,c,i,k")

                                            - g_ijka("l,k,i,c") * t1("a,l")

                                            + (Xab("X,a,d")*t1("d,i"))*Xai("X,c,k")

                                            - g_abij("c,d,k,l") * T("d,a,i,l")

                                            + (g_abij("c,d,k,l") - 0.5*g_abij("d,c,k,l")) * t2("a,d,i,l");

                        k_kaic("k,a,i,c") = g_iajb("k,a,i,c")

                                            - g_ijka("k,l,i,c") * t1("a,l")

                                            + (Xai("X,d,k")*t1("d,i"))*Xab("X,a,c")

                                            - g_abij("d,c,k,l") * T("d,a,i,l");

                    }

//        g_abij.get_world().gop.fence();
                    // compute residual r2(n) = t2(n+1) - t2(n)
                    {
                        r2("a,b,i,j") = -t2("a,b,i,j") + d2("a,b,i,j") * (
                                //
                                g_abij("a,b,i,j")
                                //
                                + a_klij("k,l,i,j") * tau("a,b,k,l")

                                //
                                + b_abij("a,b,i,j")

                                // permutation part
                                //
                                + (g_ac("a,c") * t2("c,b,i,j")
                                   - g_ki("k,i") * t2("a,b,k,j"))

                                + (g_ac("b,c") * t2("c,a,j,i")
                                   - g_ki("k,j") * t2("b,a,k,i"))

                                + Xab("X,b,c")*t1("c,j")*Xai("X,a,i")

                                + Xab("X,a,c")*t1("c,i")*Xai("X,b,j")

                                - g_iajb("k,b,i,c") * t1("c,j") * t1("a,k")

                                - g_iajb("k,a,j,c") * t1("c,i") * t1("b,k")
                                //
                                - (g_ijak("i,j,a,k") + g_abij("a,c,i,k") * t1("c,j")) * t1("b,k")

                                - (g_ijak("j,i,b,k") + g_abij("b,c,j,k") * t1("c,i")) * t1("a,k")

                                + (j_akic("a,k,i,c") - 0.5*k_kaic("k,a,i,c")) *
                                  (2.0 * t2("c,b,k,j") - t2("b,c,k,j"))

                                + (j_akic("b,k,j,c") - 0.5*k_kaic("k,b,j,c")) *
                                  (2.0 * t2("c,a,k,i") - t2("a,c,k,i"))

                                - 0.5 * k_kaic("k,a,i,c") * t2("b,c,k,j")
                                    - k_kaic("k,b,i,c") * t2("a,c,k,j")

                                - 0.5 * k_kaic("k,b,j,c") * t2("a,c,k,i")
                                    - k_kaic("k,a,j,c") * t2("b,c,k,i")
                        );
                    }
//        g_abij.get_world().gop.fence();

                    t1("a,i") = t1("a,i") + r1("a,i");
                    t2("a,b,i,j") = t2("a,b,i,j") + r2("a,b,i,j");

                    tcc::cc::T1T2<double, Tile, Policy> t(t1, t2);
                    tcc::cc::T1T2<double, Tile, Policy> r(r1, r2);
                    error = r.norm() / size(t);
                    diis.extrapolate(t, r);

                    //update t1 and t2
                    t1("a,i") = t.first("a,i");
                    t2("a,b,i,j") = t.second("a,b,i,j");

                    tau("a,b,i,j") = t2("a,b,i,j") + t1("a,i") * t1("b,j");

                    // recompute energy
                    E0 = E1;
                    E1 = 2.0 * TA::dot(f_ai("a,i"), t1("a,i"))
                         +
                         TA::dot((2.0 * g_abij("a,b,i,j") - g_abij("b,a,i,j")),
                                 tau("a,b,i,j"));
                    dE = std::abs(E0 - E1);
                    iter += 1ul;

                    auto t1 = tcc_time::now();
                    auto duration_t = tcc_time::duration_in_s(t0, t1);

                    if (g_abij.get_world().rank() == 0) {
                        std::cout.precision(15);
                        std::cout<< iter << "  " << dE << "  " << error <<
                        "  " << E1 << "  " << duration_u << " " << duration_t
                        <<std::endl;
                    }

                    g_abij.get_world().gop.fence();
//        std::cout << indent << scprintf("%-5.0f", iter) << scprintf("%-20.10f", Delta_E)
//        << scprintf("%-15.10f", E_1) << std::endl;

                }
                if (g_abij.get_world().rank() == 0) {
                    std::cout << "CCSD Energy  " << E1 << std::endl;
                }
                return E1;
            }

        private:

            void d_abij(TArray4 &abij,
                        const Eigen::VectorXd &ens, std::size_t n_occ) {
                typedef typename TArray4::range_type range_type;
                typedef typename TArray2::iterator iterator;

                auto make_tile = [&ens, n_occ](range_type &range) {

                    auto result_tile = Tile(range);

                    // compute index
                    const auto a0 = result_tile.range().lobound()[0];
                    const auto an = result_tile.range().upbound()[0];
                    const auto b0 = result_tile.range().lobound()[1];
                    const auto bn = result_tile.range().upbound()[1];
                    const auto i0 = result_tile.range().lobound()[2];
                    const auto in = result_tile.range().upbound()[2];
                    const auto j0 = result_tile.range().lobound()[3];
                    const auto jn = result_tile.range().upbound()[3];

                    auto tile_idx = 0;
                    typename Tile::value_type tmp = 1.0;
                    for (auto a = a0; a < an; ++a) {
                        const auto e_a = ens[a + n_occ];
                        for (auto b = b0; b < bn; ++b) {
                            const auto e_b = ens[b + n_occ];
                            for (auto i = i0; i < in; ++i) {
                                const auto e_i = ens[i];
                                for (auto j = j0; j < jn; ++j, ++tile_idx) {
                                    const auto e_j = ens[j];
                                    const auto e_iajb = e_i + e_j - e_a - e_b;
                                    const auto result_abij = tmp / (e_iajb);
                                    result_tile[tile_idx] = result_abij;
                                }
                            }
                        }
                    }
                    return result_tile;
                };

                for (iterator it = abij.begin(); it != abij.end(); ++it) {

                    madness::Future<Tile> tile = abij.get_world().taskq.add(
                            make_tile,
                            abij.trange().make_tile_range(it.ordinal()));

                    *it = tile;
                }

            }

            void d_ai(TArray2 &f_ai, const Eigen::VectorXd &ens, int n_occ) {
                typedef typename TArray2::range_type range_type;
                typedef typename TArray2::iterator iterator;

                auto make_tile = [&ens, n_occ](range_type &range) {

                    auto result_tile = Tile(range);
                    const auto a0 = result_tile.range().lobound()[0];
                    const auto an = result_tile.range().upbound()[0];
                    const auto i0 = result_tile.range().lobound()[1];
                    const auto in = result_tile.range().upbound()[1];

                    auto ai = 0;
                    typename Tile::value_type tmp = 1.0;
                    for (auto a = a0; a < an; ++a) {
                        const auto e_a = ens[a + n_occ];
                        for (auto i = i0; i < in; ++i, ++ai) {
                            const auto e_i = ens[i];
                            const auto e_ia = e_i - e_a;
                            const auto result_ai = tmp / (e_ia);
                            result_tile[ai] = result_ai;
                        }
                    }
                    return result_tile;
                };

                typename TArray2::pmap_interface::const_iterator it = f_ai.get_pmap()->begin();
                typename TArray2::pmap_interface::const_iterator end = f_ai.get_pmap()->end();
                for (; it != end; ++it) {

                    madness::Future<Tile> tile = f_ai.get_world().taskq.add(
                            make_tile,
                            f_ai.trange().make_tile_range(*it));

                    f_ai.set(*it, tile);
                }

            }

        private:
            Eigen::VectorXd ens_;
            std::shared_ptr<tcc::TRange1Engine> tre_;
            std::shared_ptr<tcc::cc::CCSDIntermediate<Tile, Policy>> intermediate_;
            TArrayBlock2 fock_;
        };

    } //namespace cc
} //namespace tcc


#endif //TILECLUSTERCHEM_CCSD_H
