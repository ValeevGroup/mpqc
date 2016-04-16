//
// Created by Chong Peng on 8/21/15.
//

#ifndef TILECLUSTERCHEM_CCSD_T_H_H
#define TILECLUSTERCHEM_CCSD_T_H_H

#include "ccsd.h"

namespace mpqc{
    namespace cc{

        // CCSD_T class that compute CCSD(T) triple calculation

        //Options:
        // All options in CCSD class
        // DFExpr = bool, control if use df in compute ccsd, default is True
        // Increase = int, control the increasement in outer loop, default is 2

        template<typename Tile, typename Policy>
        class CCSD_T : public CCSD<Tile,Policy> {

        public:

            using TArray = TA::DistArray<Tile,Policy>;

            CCSD_T(const Eigen::VectorXd &ens,
                 const std::shared_ptr<TRange1Engine> &tre,
                 const std::shared_ptr<CCSDIntermediate<Tile, Policy>> &inter,
                   rapidjson::Document &options):
                    CCSD<Tile,Policy>(ens,tre,inter,options)
            {}

            void compute(){

                TArray t1;
                TArray t2;

                double ccsd_corr = 0.0;
                // compute CCSD first
                auto direct = this->options_.HasMember("Direct") ? this->options_["Direct"].GetBool(): false;
                if(direct){
                    ccsd_corr = CCSD<Tile, Policy>::compute_ccsd_direct(t1, t2);
                    // TODO smarter way to clean integrals not needed
                    // clean integrals not needed
                    this->ccsd_intermediate_->clean_two_electron();

                }
                else {
                    ccsd_corr = CCSD<Tile,Policy>::compute_ccsd_straight(t1, t2);
                    this->ccsd_intermediate_->clean_two_electron();
                }

                // start CCSD(T)
                if(t1.get_world().rank() == 0){
                    std::cout << "\nBegining CCSD(T) " << std::endl;
                }
                auto time0 = mpqc_time::now();
                double ccsd_t = compute_ccsd_t(t1, t2);
                auto time1 = mpqc_time::now();
                auto duration1 = mpqc_time::duration_in_s(time0, time1);

                if (t1.get_world().rank() == 0) {
                    std::cout << std::setprecision(15);
                    std::cout << "(T) Energy      " << ccsd_t << " Time " << duration1 << std::endl;
//                    std::cout << "(T) Energy      " << ccsd_t_d << " Time " << duration2 << std::endl;
                    std::cout << "CCSD(T) Energy  " << ccsd_t + ccsd_corr << std::endl;
//                    std::cout << "CCSD(T) Energy  " << ccsd_t_d + ccsd_corr << std::endl;
                }

            }

            double compute_ccsd_t(TArray& t1, TArray& t2){
                bool df = this->options_.HasMember("DFExpr") ? this->options_["DFExpr"].GetBool() : true;
                if(df && t1.get_world().rank() == 0){
                    std::cout << "Use Density Fitting Expression to avoid storing G_vovv" << std::endl;
                }
                // get integral
                TArray g_jklc = this->ccsd_intermediate_->get_ijka();
                TArray g_abij = this->ccsd_intermediate_->get_abij();

                TArray g_diba;
                TArray Xdb;
                TArray Xai;

                if (df){
                    Xdb = this->ccsd_intermediate_->get_Xab();
                    Xai = this->ccsd_intermediate_->get_Xai();
                }else{
                    g_diba= this->ccsd_intermediate_->get_aibc();
                }
                // get trange1
                auto tr_occ = this->trange1_engine_->get_occ_tr1();
                auto tr_vir = this->trange1_engine_->get_vir_tr1();

                auto n_tr_occ = this->trange1_engine_->get_occ_blocks();
                auto n_tr_vir = this->trange1_engine_->get_vir_blocks();
                std::size_t n_tr_x = 0;
                if (df){
                    n_tr_x = Xdb.trange().data().front().tiles().second;
                }


                double triple_energy = 0.0;


                std::size_t increase = this->options_.HasMember("Increase") ? this->options_["Increase"].GetInt() : 2;
                if (increase > n_tr_vir){
                    increase = n_tr_vir;
                }
                std::size_t a_increase = increase;
                std::size_t b_increase = increase;
                std::size_t c_increase = increase;

                std::size_t block_size = this->trange1_engine_->get_occ_block_size();
                std::size_t n_blocks = increase*increase*increase*n_tr_occ*n_tr_occ*n_tr_occ;
                double mem = (n_blocks*std::pow(block_size,6)*8)/(std::pow(1024.0,3));

                if(t1.get_world().rank() == 0){
                    std::cout << "Increase in the loop " << increase << std::endl;
                    std::cout << "Number of blocks at each iteration " << n_blocks << std::endl;
                    std::cout << std::setprecision(5);
                    std::cout << "Size of T3 or V3 at each iteration " << mem << " GB" << std::endl;
                }

                // index in virtual blocks
                std::size_t a = 0;
                std::size_t b = 0;
                std::size_t c = 0;

                // number of blocks computed
                std::size_t n_blocks_computed = 0;

                // loop over virtual blocks
                while (a < n_tr_vir){
                    b = 0;
                    if( a + increase >= n_tr_vir){
                        a_increase = (n_tr_vir - a);
                    }else{
                        a_increase = increase;
                    }


                    std::size_t a_end = a + a_increase - 1;
                    while(b <= a_end){
                        c = 0;

                        if(b + increase - 1 > a_end){
                            if (b == a_end){
                                b_increase = 1;
                            }else{
                                b_increase = (a_end - b);
                            }
                        }else{
                            b_increase = increase;
                        }

                        std::size_t b_end = b + b_increase - 1;

                        while(c <= b_end){

                            if(c + increase - 1> b_end){
                                if ( c == b_end){
                                    c_increase = 1;
                                }else{
                                    c_increase = (b_end - c);
                                }
                            }else{
                                c_increase = increase;
                            }

                            std::size_t c_end = c + c_increase - 1;

//                            std::cout << a << " " << b << " " << c << std::endl;
                            std::size_t a_low = a;
                            std::size_t a_up = a + a_increase;
                            std::size_t b_low = b;
                            std::size_t b_up = b + b_increase;
                            std::size_t c_low = c;
                            std::size_t c_up = c + c_increase;

                            std::size_t blocks = (a_up-a_low)*(b_up-b_low)*(c_up-c_low)*n_tr_occ*n_tr_occ*n_tr_occ;
//                            if (t1.get_world().rank() == 0){
//                                std::cout << "{" << a_low << " " << b_low << " " << c_low << "}" << " ";
//                                std::cout << "{" << a_up << " " << b_up << " " << c_up << "} " << blocks << std::endl;
//                            }
                            n_blocks_computed += blocks;

                            typedef std::vector<std::size_t> block;

                            // compute t3
                            TArray t3;
                            // abcijk contribution
                            // g^{di}_{ba}*t^{cd}_{kj} - g^{jk}_{lc}*t^{ab}_{il}
                            {

                                TArray block_g_diba, block_t2_cdkj, block_g_jklc, block_t2_abil;
                                if(df){
                                    // block for Xdb
                                    block Xdb_low {0, 0, b_low};
                                    block Xdb_up {n_tr_x, n_tr_vir, b_up};

                                    // block for Xai
                                    block Xai_low {0, a_low, 0};
                                    block Xai_up {n_tr_x, a_up, n_tr_occ};

                                    block_g_diba("d,i,b,a") = Xai("X,a,i").block(Xai_low,Xai_up)*Xdb("X,d,b").block(Xdb_low,Xdb_up);
                                }else{
                                    // block for g_diba
                                    block g_diba_low {0,0,b_low,a_low};
                                    block g_diba_up {n_tr_vir,n_tr_occ,b_up,a_up};

                                    block_g_diba("d,i,b,a") = g_diba("d,i,b,a").block(g_diba_low,g_diba_up);
                                }

                                // block for t2_cdk
                                block t2_cdkj_low {c_low,0,0,0};
                                block t2_cdkj_up {c_up,n_tr_vir,n_tr_occ,n_tr_occ};
                                block_t2_cdkj("c,d,k,j") = t2("c,d,k,j").block(t2_cdkj_low,t2_cdkj_up);

                                // block for g_jklc
                                block g_jklc_low {0,0,0,c_low};
                                block g_jklc_up {n_tr_occ,n_tr_occ,n_tr_occ,c_up};

                                block_g_jklc("j,k,l,c") = g_jklc("j,k,l,c").block(g_jklc_low,g_jklc_up);

                                // block for t2_abil
                                block t2_abil_low {a_low,b_low,0,0};
                                block t2_abil_up {a_up,b_up,n_tr_occ,n_tr_occ};

                                block_t2_abil("a,b,i,l") = t2("a,b,i,l").block(t2_abil_low,t2_abil_up);


                                t3("a,b,c,i,j,k") = block_g_diba("d,i,b,a")*block_t2_cdkj("c,d,k,j")
                                                    - block_g_jklc("l,k,j,c")*block_t2_abil("a,b,i,l");
                            }

                            // bcajki contribution
                            // g^{dj}_{cb}*t^{ad}_{ik} - g^{ki}_{la}*t^{bc}_{jl}
                            {
                                TArray block_g_djcb, block_g_kila, block_t2_adik, block_t2_bcjl;

                                if(df){
                                   // block for Xdc
                                    block Xdc_low {0,0,c_low};
                                    block Xdc_up {n_tr_x,n_tr_vir,c_up};

                                    // block for Xbj
                                    block Xbj_low {0,b_low,0};
                                    block Xbj_up {n_tr_x,b_up,n_tr_occ};

                                    block_g_djcb("d,j,c,b") = Xai("X,b,j").block(Xbj_low,Xbj_up)*Xdb("X,d,c").block(Xdc_low,Xdc_up);
                                }else{
                                    // block for g_djcb
                                    block g_djcb_low {0,0,c_low,b_low};
                                    block g_djcb_up {n_tr_vir, n_tr_occ, c_up, b_up};

                                    block_g_djcb("d,j,c,b") = g_diba("d,j,c,b").block(g_djcb_low,g_djcb_up);
                                }

                                // block for t2_adik
                                block t2_adik_low {a_low,0,0,0};
                                block t2_adik_up {a_up, n_tr_vir, n_tr_occ, n_tr_occ};
                                block_t2_adik("a,d,i,k") = t2("a,d,i,k").block(t2_adik_low,t2_adik_up);

                                // block for g_kila
                                block g_kila_low {0, 0, 0, a_low};
                                block g_kila_up {n_tr_occ, n_tr_occ, n_tr_occ, a_up};

                                block_g_kila("k,i,l,a") = g_jklc("k,i,l,a").block(g_kila_low,g_kila_up);

                                // block for t2_bcjl
                                block t2_bcjl_low {b_low,c_low,0,0};
                                block t2_bcjl_up {b_up, c_up, n_tr_occ, n_tr_occ};
                                block_t2_bcjl("b,c,j,l") = t2("b,c,j,l").block(t2_bcjl_low,t2_bcjl_up);

                                t3("a,b,c,i,j,k") += block_g_djcb("d,j,c,b")*block_t2_adik("a,d,i,k")
                                                     - block_g_kila("k,i,l,a")*block_t2_bcjl("b,c,j,l");
                            }

                            // cabkij contribution
                            // g^{dk}_{ac}*t^{bd}_{ji} - g^{ij}_{lb}*t^{ca}_{kl}
                            {
                                TArray block_g_dkac, block_g_ijlb, block_t2_bdji, block_t2_cakl;

                                if(df){
                                    // block for Xda
                                    block Xda_low {0,0,a_low};
                                    block Xda_up {n_tr_x,n_tr_vir,a_up};

                                    // block for Xck
                                    block Xck_low {0,c_low,0};
                                    block Xck_up {n_tr_x,c_up,n_tr_occ};

                                    block_g_dkac("d,k,a,c") = Xai("X,c,k").block(Xck_low,Xck_up)*Xdb("X,d,a").block(Xda_low,Xda_up);

                                }else{
                                    // block for g_dkac
                                    block g_dkac_low {0,0,a_low,c_low};
                                    block g_dkac_up {n_tr_vir, n_tr_occ, a_up, c_up};

                                    block_g_dkac("d,k,a,c") = g_diba("d,k,a,c").block(g_dkac_low, g_dkac_up);
                                }

                                // block for t2_bdji
                                block t2_bdji_low {b_low,0,0,0};
                                block t2_bdji_up {b_up, n_tr_vir, n_tr_occ, n_tr_occ};
                                block_t2_bdji("b,d,j,i") = t2("b,d,j,i").block(t2_bdji_low,t2_bdji_up);

                                // block for g_ijlb
                                block g_ijlb_low {0,0,0, b_low};
                                block g_ijlb_up {n_tr_occ, n_tr_occ, n_tr_occ, b_up};

                                block_g_ijlb("i,j,l,b") = g_jklc("i,j,l,b").block(g_ijlb_low,g_ijlb_up);

                                // block for t2_cakl
                                block t2_cakl_low {c_low,a_low,0,0};
                                block t2_cakl_up {c_up, a_up, n_tr_occ, n_tr_occ};

                                block_t2_cakl("c,a,k,l") = t2("c,a,k,l").block(t2_cakl_low,t2_cakl_up);

                                t3("a,b,c,i,j,k") += block_g_dkac("d,k,a,c")*block_t2_bdji("b,d,j,i")
                                                     - block_g_ijlb("i,j,l,b")*block_t2_cakl("c,a,k,l");
                            }

                            // bacjik contribution
                            // g^{dj}_{ab}*t^{cd}_{ki} - g^{ik}_{lc}*t^{ba}_{jl}
                            {

                                TArray block_g_djab, block_t2_cdki, block_g_iklc, block_t2_bajl;

                                if (df){

                                    // block for Xda
                                    block Xda_low{0,0,a_low};
                                    block Xda_up {n_tr_x,n_tr_vir,a_up};

                                    // block for Xbj
                                    block Xbj_low {0,b_low,0};
                                    block Xbj_up {n_tr_x,b_up,n_tr_occ};

                                    block_g_djab("d,j,a,b") = Xai("X,b,j").block(Xbj_low,Xbj_up)*Xdb("X,d,a").block(Xda_low,Xda_up);

                                }else{
                                    // block for g_djab
                                    block g_djab_low {0,0,a_low,b_low};
                                    block g_djab_up {n_tr_vir,n_tr_occ,a_up,b_up};

                                    block_g_djab("d,j,a,b") = g_diba("d,j,a,b").block(g_djab_low,g_djab_up);
                                }

                                // block for t2_cdki
                                block t2_cdki_low {c_low,0,0,0};
                                block t2_cdki_up {c_up,n_tr_vir,n_tr_occ,n_tr_occ};

                                block_t2_cdki("c,d,k,i") = t2("c,d,k,i").block(t2_cdki_low,t2_cdki_up);

                                // block for g_iklc
                                block g_iklc_low {0,0,0,c_low};
                                block g_iklc_up {n_tr_occ,n_tr_occ,n_tr_occ,c_up};

                                block_g_iklc("i,k,l,c") = g_jklc("i,k,l,c").block(g_iklc_low,g_iklc_up);

                                // block for t2_bajl
                                block t2_bajl_low {b_low,a_low,0,0};
                                block t2_bajl_up {b_up,a_up,n_tr_occ,n_tr_occ};

                                block_t2_bajl("b,a,j,l") = t2("b,a,j,l").block(t2_bajl_low,t2_bajl_up);


                                t3("a,b,c,i,j,k") += block_g_djab("d,j,a,b")*block_t2_cdki("c,d,k,i")
                                                     - block_g_iklc("i,k,l,c")*block_t2_bajl("b,a,j,l");
                            }

                            // acbikj contribution
                            // g^{di}_{ca}*t^{bd}_{jk} - g^{kj}_{lb}*t^{ac}_{il}
                            {

                                TArray block_g_dica, block_t2_bdjk, block_g_kjlb, block_t2_acil;


                                if(df){
                                    // block for Xdc
                                    block Xdc_low {0,0,c_low};
                                    block Xdc_up {n_tr_x,n_tr_vir,c_up};

                                    // block for Xai
                                    block Xai_low {0,a_low,0};
                                    block Xai_up {n_tr_x,a_up,n_tr_occ};

                                    block_g_dica("d,i,c,a") = Xai("X,a,i").block(Xai_low,Xai_up)*Xdb("X,d,c").block(Xdc_low,Xdc_up);

                                }else{
                                    // block for g_dica
                                    block g_dica_low {0,0,c_low,a_low};
                                    block g_dica_up {n_tr_vir,n_tr_occ,c_up,a_up};

                                    block_g_dica("d,i,c,a") = g_diba("d,i,c,a").block(g_dica_low,g_dica_up);
                                }

                                // block for t2_bdjk
                                block t2_bdjk_low {b_low,0,0,0};
                                block t2_bdjk_up {b_up,n_tr_vir,n_tr_occ,n_tr_occ};

                                block_t2_bdjk("b,d,j,k") = t2("b,d,j,k").block(t2_bdjk_low,t2_bdjk_up);

                                // block for g_kjlb
                                block g_kjlb_low {0,0,0,b_low};
                                block g_kjlb_up {n_tr_occ,n_tr_occ,n_tr_occ,b_up};

                                block_g_kjlb("k,j,l,b") = g_jklc("k,j,l,b").block(g_kjlb_low,g_kjlb_up);

                                // block for t2_acil
                                block t2_acil_low {a_low,c_low,0,0};
                                block t2_acil_up {a_up,c_up,n_tr_occ,n_tr_occ};

                                block_t2_acil("a,c,i,l") = t2("a,c,i,l").block(t2_acil_low,t2_acil_up);


                                t3("a,b,c,i,j,k") += block_g_dica("d,i,c,a")*block_t2_bdjk("b,d,j,k")
                                                     - block_g_kjlb("k,j,l,b")*block_t2_acil("a,c,i,l");
                            }

                            // cbakji contribution
                            // g^{dk}_{bc}*t^{ad}_{ij} - g^{ji}_{la}*t^{cb}_{kl}
                            {

                                TArray block_g_dkbc, block_g_jila, block_t2_adij, block_t2_cbkl;

                                if(df){
                                    // block for Xdb
                                    block Xdb_low {0,0,b_low};
                                    block Xdb_up {n_tr_x,n_tr_vir,b_up};

                                    // block for Xck
                                    block Xck_low {0,c_low,0};
                                    block Xck_up {n_tr_x,c_up,n_tr_occ};

                                    block_g_dkbc("d,k,b,c") = Xai("X,c,k").block(Xck_low,Xck_up)*Xdb("X,d,b").block(Xdb_low,Xdb_up);

                                }else{

                                    // block for g_dkbc
                                    block g_dkbc_low {0,0,b_low,c_low};
                                    block g_dkbc_up {n_tr_vir, n_tr_occ, b_up, c_up};

                                    block_g_dkbc("d,k,b,c") = g_diba("d,k,b,c").block(g_dkbc_low,g_dkbc_up);
                                }
                                // block for t2_adi
                                block t2_adij_low {a_low,0,0,0};
                                block t2_adij_up {a_up, n_tr_vir, n_tr_occ, n_tr_occ};

                                block_t2_adij("a,d,i,j") = t2("a,d,i,j").block(t2_adij_low,t2_adij_up);

                                // block for g_jila
                                block g_jila_low {0, 0, 0, a_low};
                                block g_jila_up {n_tr_occ, n_tr_occ, n_tr_occ, a_up};

                                block_g_jila("j,i,l,a") = g_jklc("j,i,l,a").block(g_jila_low,g_jila_up);

                                // block for t2_cbkl
                                block t2_cbkl_low {c_low,b_low,0,0};
                                block t2_cbkl_up {c_up, b_up, n_tr_occ, n_tr_occ};

                                block_t2_cbkl("c,b,k,l") = t2("c,b,k,l").block(t2_cbkl_low,t2_cbkl_up);

                                t3("a,b,c,i,j,k") += block_g_dkbc("d,k,b,c")*block_t2_adij("a,d,i,j")
                                                     - block_g_jila("j,i,l,a")*block_t2_cbkl("c,b,k,l");
                            }

                            // compute v3
                            TArray v3;
                            // abcijk contribution
                            // g^{ab}_{ij}*t^{c}_{k}
                            {
                                // block for g_abij
                                block g_abij_low {a_low,b_low,0,0};
                                block g_abij_up {a_up,b_up,n_tr_occ,n_tr_occ};

                                // block for t1_ck
                                block t1_ck_low {c_low,0};
                                block t1_ck_up {c_up, n_tr_occ};
                                v3("a,b,c,i,j,k") = g_abij("a,b,i,j").block(g_abij_low,g_abij_up)*t1("c,k").block(t1_ck_low,t1_ck_up);
                            }

                            // acbikj contribution
                            // g^{ac}_{ik}*t^{b}_{j}
                            {
                                // block for g_acik
                                block g_acik_low {a_low,c_low,0,0};
                                block g_acik_up {a_up,c_up,n_tr_occ,n_tr_occ};

                                // block for t1_bj
                                block t1_bj_low {b_low,0};
                                block t1_bj_up {b_up, n_tr_occ};
                                v3("a,b,c,i,j,k") += g_abij("a,c,i,k").block(g_acik_low,g_acik_up)*t1("b,j").block(t1_bj_low,t1_bj_up);
                            }

                            // bcajki contribution
                            // g^{bc}_{jk}*t^{a}_{i}
                            {
                                // block for g_bcjk
                                block g_bcjk_low {b_low, c_low, 0,0};
                                block g_bcjk_up {b_up, c_up, n_tr_occ, n_tr_occ};

                                // block for t1_ai
                                block t1_ai_low {a_low, 0};
                                block t1_ai_up {a_up, n_tr_occ};

                                v3("a,b,c,i,j,k") += g_abij("b,c,j,k").block(g_bcjk_low,g_bcjk_up)*t1("a,i").block(t1_ai_low,t1_ai_up);
                            }

                            // compute offset
                            std::size_t a_offset = tr_vir.tile(a).first;
                            std::size_t b_offset = tr_vir.tile(b).first;
                            std::size_t c_offset = tr_vir.tile(c).first;
//                            std::cout << a_offset << " " << b_offset << " " << c_offset << std::endl;
                            std::array<std::size_t,6> offset{a_offset,b_offset,c_offset,0,0,0};

                            double tmp_energy = 0.0;
                            if (b_end < a && c_end < b) {

                                auto ccsd_t_reduce = CCSD_T_Reduce(
                                        std::make_shared<Eigen::VectorXd>(this->orbital_energy_),
                                        this->trange1_engine_->get_actual_occ(), offset);
                                tmp_energy =
                                        (
                                                (t3("a,b,c,i,j,k")
                                                 + v3("a,b,c,i,j,k")
                                                )
                                                * (4.0 * t3("a,b,c,i,j,k")
                                                   + t3("a,b,c,k,i,j")
                                                   + t3("a,b,c,j,k,i")
                                                   - 2 * (t3("a,b,c,k,j,i") + t3("a,b,c,i,k,j") + t3("a,b,c,j,i,k"))
                                                )
                                        ).reduce(ccsd_t_reduce);

                                tmp_energy *= 2;
                            } else {

                                auto ccsd_t_reduce = CCSD_T_ReduceSymm(
                                        std::make_shared<Eigen::VectorXd>(this->orbital_energy_),
                                        this->trange1_engine_->get_actual_occ(), offset);
                                tmp_energy =
                                        (
                                                (t3("a,b,c,i,j,k")
                                                 + v3("a,b,c,i,j,k")
                                                )
                                                * (4.0 * t3("a,b,c,i,j,k")
                                                   + t3("a,b,c,k,i,j")
                                                   + t3("a,b,c,j,k,i")
                                                   - 2 * (t3("a,b,c,k,j,i")
                                                          + t3("a,b,c,i,k,j")
                                                          + t3("a,b,c,j,i,k"))
                                                )
                                        ).reduce(ccsd_t_reduce);
                            }

                            triple_energy += tmp_energy;

                            c += c_increase;
                        }
                        b += b_increase;
                    }

                    if(t1.get_world().rank() == 0){
                        print_progress(a, a+increase, n_tr_vir);
                    }
                    a += a_increase;
                }

                if (t1.get_world().rank() == 0){
                    std::cout << "Total Blocks Computed  " << n_blocks_computed;
                    std::cout << " from " << std::pow(n_tr_occ,3)*std::pow(n_tr_vir,3) << std::endl;
                }
                return  triple_energy;
            }

            // compute and store all t3 amplitudes, not recommanded for performance computing
            double compute_ccsd_t_straight(const TArray& t1, const TArray& t2){

                // get integral
                TArray g_jklc = this->ccsd_intermediate_->get_ijka();
                TArray g_diba = this->ccsd_intermediate_->get_aibc();
                TArray g_abij = this->ccsd_intermediate_->get_abij();

                // compute t3
                TArray t3;
                t3("a,b,c,i,j,k") = g_diba("d,i,b,a")*t2("c,d,k,j") - g_jklc("l,k,j,c")*t2("a,b,i,l");
                t3("a,b,c,i,j,k") = t3("a,b,c,i,j,k") + t3("a,c,b,i,k,j") + t3("c,a,b,k,i,j") + t3("c,b,a,k,j,i")
                        + t3("b,c,a,j,k,i") + t3("b,a,c,j,i,k");

                // compute v3
                TArray v3;
                v3("a,b,c,i,j,k") = g_abij("a,b,i,j")*t1("c,k");
                v3("a,b,c,i,j,k") = v3("a,b,c,i,j,k") + v3("b,c,a,j,k,i") + v3("a,c,b,i,k,j");

                std::array<std::size_t,6> offset {0,0,0,0,0,0};

                auto ccsd_t_reduce = CCSD_T_Reduce(
                        std::make_shared<Eigen::VectorXd>(this->orbital_energy_),
                        this->trange1_engine_->get_actual_occ(),
                        offset);

                double triple_energy =
                        (
                                (t3("a,b,c,i,j,k")
                                 + v3("a,b,c,i,j,k")
                                )
                                * (4.0 * t3("a,b,c,i,j,k")
                                   + t3("a,b,c,k,i,j")
                                   + t3("a,b,c,j,k,i")
                                   -2*(t3("a,b,c,k,j,i")
                                       +t3("a,b,c,i,k,j")
                                       +t3("a,b,c,j,i,k")))
                        ).reduce(ccsd_t_reduce);
                triple_energy = triple_energy/3.0;
                return triple_energy;
            }

            double compute_ccsd_t_direct(TArray& t1,TArray& t2){

                // get integral
                TArray g_jklc = this->ccsd_intermediate_->get_ijka();
                TArray g_diba = this->ccsd_intermediate_->get_aibc();
                TArray g_abij = this->ccsd_intermediate_->get_abij();

                // get trange1
                auto tr_occ = this->trange1_engine_->get_occ_tr1();
                auto tr_vir = this->trange1_engine_->get_vir_tr1();

                auto n_tr_occ = this->trange1_engine_->get_occ_blocks();
                auto n_tr_vir = this->trange1_engine_->get_vir_blocks();

                // number of blocks computed
                std::size_t n_blocks_computed = 0;

                double triple_energy = 0.0;
                // loop over virtual blocks
                for (std::size_t a = 0; a < n_tr_vir; ++a){
                    for(std::size_t b = 0; b <= a; ++b){
                        for(std::size_t c = 0; c <= b; ++c){

//                            std::cout << a << " " << b << " " << c << std::endl;
                            std::size_t a_low = a;
                            std::size_t a_up = a + 1;
                            std::size_t b_low = b;
                            std::size_t b_up = b + 1;
                            std::size_t c_low = c;
                            std::size_t c_up = c + 1;
                            std::size_t blocks = (a_up-a_low)*(b_up-b_low)*(c_up-c_low);
//                            std::cout << a_up << " " << b_up << " " << c_up << std::endl;
//                            std::cout << a << " " << b << " " << c << std::endl;
                            n_blocks_computed += blocks;

                            typedef std::vector<std::size_t> block;

                            // compute t3
                            TArray t3;
                            // abcijk contribution
                            // g^{di}_{ba}*t^{cd}_{kj} - g^{jk}_{lc}*t^{ab}_{il}
                            {
                                // block for t2_cdkj
                                block t2_cdkj_low {c_low,0,0,0};
                                block t2_cdkj_up {c_up,n_tr_vir,n_tr_occ,n_tr_occ};

                                // block for t2_abil
                                block t2_abil_low {a_low,b_low,0,0};
                                block t2_abil_up {a_up,b_up,n_tr_occ,n_tr_occ};

                                // block for g_diba
                                block g_diba_low {0,0,b_low,a_low};
                                block g_diba_up {n_tr_vir,n_tr_occ,b_up,a_up};

                                // block for g_jklc
                                block g_jklc_low {0,0,0,c_low};
                                block g_jklc_up {n_tr_occ,n_tr_occ,n_tr_occ,c_up};


                                TArray block_g_diba, block_t2_cdkj, block_g_jklc, block_t2_abil;
                                block_g_diba("d,i,b,a") = g_diba("d,i,b,a").block(g_diba_low,g_diba_up);
                                block_t2_cdkj("c,d,k,j") = t2("c,d,k,j").block(t2_cdkj_low,t2_cdkj_up);

                                block_g_jklc("j,k,l,c") = g_jklc("j,k,l,c").block(g_jklc_low,g_jklc_up);
                                block_t2_abil("a,b,i,l") = t2("a,b,i,l").block(t2_abil_low,t2_abil_up);


                                t3("a,b,c,i,j,k") = block_g_diba("d,i,b,a")*block_t2_cdkj("c,d,k,j")
                                                    - block_g_jklc("l,k,j,c")*block_t2_abil("a,b,i,l");
                            }

                            // bcajki contribution
                            // g^{dj}_{cb}*t^{ad}_{ik} - g^{ki}_{la}*t^{bc}_{jl}
                            {
                                // block for t2_adik
                                block t2_adik_low {a_low,0,0,0};
                                block t2_adik_up {a_up, n_tr_vir, n_tr_occ, n_tr_occ};

                                // block for t2_bcjl
                                block t2_bcjl_low {b_low,c_low,0,0};
                                block t2_bcjl_up {b_up, c_up, n_tr_occ, n_tr_occ};

                                // block for g_djcb
                                block g_djcb_low {0,0,c_low,b_low};
                                block g_djcb_up {n_tr_vir, n_tr_occ, c_up, b_up};

                                // block for g_kila
                                block g_kila_low {0, 0, 0, a_low};
                                block g_kila_up {n_tr_occ, n_tr_occ, n_tr_occ, a_up};

                                TArray block_g_djcb, block_g_kila, block_t2_adik, block_t2_bcjl;

                                block_g_djcb("d,j,c,b") = g_diba("d,j,c,b").block(g_djcb_low,g_djcb_up);
                                block_t2_adik("a,d,i,k") = t2("a,d,i,k").block(t2_adik_low,t2_adik_up);

                                block_g_kila("k,i,l,a") = g_jklc("k,i,l,a").block(g_kila_low,g_kila_up);
                                block_t2_bcjl("b,c,j,l") = t2("b,c,j,l").block(t2_bcjl_low,t2_bcjl_up);

                                t3("a,b,c,i,j,k") += block_g_djcb("d,j,c,b")*block_t2_adik("a,d,i,k")
                                                     - block_g_kila("k,i,l,a")*block_t2_bcjl("b,c,j,l");
                            }

                            // cabkij contribution
                            // g^{dk}_{ac}*t^{bd}_{ji} - g^{ij}_{lb}*t^{ca}_{kl}
                            {
                                // block for t2_bdji
                                block t2_bdji_low {b_low,0,0,0};
                                block t2_bdji_up {b_up, n_tr_vir, n_tr_occ, n_tr_occ};

                                // block for t2_cakl
                                block t2_cakl_low {c_low,a_low,0,0};
                                block t2_cakl_up {c_up, a_up, n_tr_occ, n_tr_occ};

                                // block for g_dkac
                                block g_dkac_low {0,0,a_low,c_low};
                                block g_dkac_up {n_tr_vir, n_tr_occ, a_up, c_up};

                                // block for g_ijlb
                                block g_ijlb_low {0,0,0, b_low};
                                block g_ijlb_up {n_tr_occ, n_tr_occ, n_tr_occ, b_up};

                                TArray block_g_dkac, block_g_ijlb, block_t2_bdji, block_t2_cakl;

                                block_g_dkac("d,k,a,c") = g_diba("d,k,a,c").block(g_dkac_low, g_dkac_up);
                                block_t2_bdji("b,d,j,i") = t2("b,d,j,i").block(t2_bdji_low,t2_bdji_up);

                                block_g_ijlb("i,j,l,b") = g_jklc("i,j,l,b").block(g_ijlb_low,g_ijlb_up);
                                block_t2_cakl("c,a,k,l") = t2("c,a,k,l").block(t2_cakl_low,t2_cakl_up);

                                t3("a,b,c,i,j,k") += block_g_dkac("d,k,a,c")*block_t2_bdji("b,d,j,i")
                                                     - block_g_ijlb("i,j,l,b")*block_t2_cakl("c,a,k,l");
                            }

                            // bacjik contribution
                            // g^{dj}_{ab}*t^{cd}_{ki} - g^{ik}_{lc}*t^{ba}_{jl}
                            {
                                // block for t2_cdki
                                block t2_cdki_low {c_low,0,0,0};
                                block t2_cdki_up {c_up,n_tr_vir,n_tr_occ,n_tr_occ};

                                // block for t2_bajl
                                block t2_bajl_low {b_low,a_low,0,0};
                                block t2_bajl_up {b_up,a_up,n_tr_occ,n_tr_occ};

                                // block for g_djab
                                block g_djab_low {0,0,a_low,b_low};
                                block g_djab_up {n_tr_vir,n_tr_occ,a_up,b_up};

                                // block for g_iklc
                                block g_iklc_low {0,0,0,c_low};
                                block g_iklc_up {n_tr_occ,n_tr_occ,n_tr_occ,c_up};


                                TArray block_g_djab, block_t2_cdki, block_g_iklc, block_t2_bajl;
                                block_g_djab("d,j,a,b") = g_diba("d,j,a,b").block(g_djab_low,g_djab_up);
                                block_t2_cdki("c,d,k,i") = t2("c,d,k,i").block(t2_cdki_low,t2_cdki_up);

                                block_g_iklc("i,k,l,c") = g_jklc("i,k,l,c").block(g_iklc_low,g_iklc_up);
                                block_t2_bajl("b,a,j,l") = t2("b,a,j,l").block(t2_bajl_low,t2_bajl_up);


                                t3("a,b,c,i,j,k") += block_g_djab("d,j,a,b")*block_t2_cdki("c,d,k,i")
                                                     - block_g_iklc("i,k,l,c")*block_t2_bajl("b,a,j,l");
                            }

                            // acbikj contribution
                            // g^{di}_{ca}*t^{bd}_{jk} - g^{kj}_{lb}*t^{ac}_{il}
                            {
                                // block for t2_bdjk
                                block t2_bdjk_low {b_low,0,0,0};
                                block t2_bdjk_up {b_up,n_tr_vir,n_tr_occ,n_tr_occ};

                                // block for t2_acil
                                block t2_acil_low {a_low,c_low,0,0};
                                block t2_acil_up {a_up,c_up,n_tr_occ,n_tr_occ};

                                // block for g_dica
                                block g_dica_low {0,0,c_low,a_low};
                                block g_dica_up {n_tr_vir,n_tr_occ,c_up,a_up};

                                // block for g_kjlb
                                block g_kjlb_low {0,0,0,b_low};
                                block g_kjlb_up {n_tr_occ,n_tr_occ,n_tr_occ,b_up};


                                TArray block_g_dica, block_t2_bdjk, block_g_kjlb, block_t2_acil;
                                block_g_dica("d,i,c,a") = g_diba("d,i,c,a").block(g_dica_low,g_dica_up);
                                block_t2_bdjk("b,d,j,k") = t2("b,d,j,k").block(t2_bdjk_low,t2_bdjk_up);

                                block_g_kjlb("k,j,l,b") = g_jklc("k,j,l,b").block(g_kjlb_low,g_kjlb_up);
                                block_t2_acil("a,c,i,l") = t2("a,c,i,l").block(t2_acil_low,t2_acil_up);


                                t3("a,b,c,i,j,k") += block_g_dica("d,i,c,a")*block_t2_bdjk("b,d,j,k")
                                                     - block_g_kjlb("k,j,l,b")*block_t2_acil("a,c,i,l");
                            }

                            // cbakji contribution
                            // g^{dk}_{bc}*t^{ad}_{ij} - g^{ji}_{la}*t^{cb}_{kl}
                            {
                                // block for t2_adij
                                block t2_adij_low {a_low,0,0,0};
                                block t2_adij_up {a_up, n_tr_vir, n_tr_occ, n_tr_occ};

                                // block for t2_cbkl
                                block t2_cbkl_low {c_low,b_low,0,0};
                                block t2_cbkl_up {c_up, b_up, n_tr_occ, n_tr_occ};

                                // block for g_dkbc
                                block g_dkbc_low {0,0,b_low,c_low};
                                block g_dkbc_up {n_tr_vir, n_tr_occ, b_up, c_up};

                                // block for g_jila
                                block g_jila_low {0, 0, 0, a_low};
                                block g_jila_up {n_tr_occ, n_tr_occ, n_tr_occ, a_up};

                                TArray block_g_dkbc, block_g_jila, block_t2_adij, block_t2_cbkl;

                                block_g_dkbc("d,k,b,c") = g_diba("d,k,b,c").block(g_dkbc_low,g_dkbc_up);
                                block_t2_adij("a,d,i,j") = t2("a,d,i,j").block(t2_adij_low,t2_adij_up);

                                block_g_jila("j,i,l,a") = g_jklc("j,i,l,a").block(g_jila_low,g_jila_up);
                                block_t2_cbkl("c,b,k,l") = t2("c,b,k,l").block(t2_cbkl_low,t2_cbkl_up);

                                t3("a,b,c,i,j,k") += block_g_dkbc("d,k,b,c")*block_t2_adij("a,d,i,j")
                                                     - block_g_jila("j,i,l,a")*block_t2_cbkl("c,b,k,l");
                            }

                            // compute v3
                            TArray v3;
                            // abcijk contribution
                            // g^{ab}_{ij}*t^{c}_{k}
                            {
                                // block for g_abij
                                block g_abij_low {a_low,b_low,0,0};
                                block g_abij_up {a_up,b_up,n_tr_occ,n_tr_occ};

                                // block for t1_ck
                                block t1_ck_low {c_low,0};
                                block t1_ck_up {c_up, n_tr_occ};
                                v3("a,b,c,i,j,k") = g_abij("a,b,i,j").block(g_abij_low,g_abij_up)*t1("c,k").block(t1_ck_low,t1_ck_up);
                            }

                            // acbikj contribution
                            // g^{ac}_{ik}*t^{b}_{j}
                            {
                                // block for g_acik
                                block g_acik_low {a_low,c_low,0,0};
                                block g_acik_up {a_up,c_up,n_tr_occ,n_tr_occ};

                                // block for t1_bj
                                block t1_bj_low {b_low,0};
                                block t1_bj_up {b_up, n_tr_occ};
                                v3("a,b,c,i,j,k") += g_abij("a,c,i,k").block(g_acik_low,g_acik_up)*t1("b,j").block(t1_bj_low,t1_bj_up);
                            }

                            // bcajki contribution
                            // g^{bc}_{jk}*t^{a}_{i}
                            {
                                // block for g_bcjk
                                block g_bcjk_low {b_low, c_low, 0,0};
                                block g_bcjk_up {b_up, c_up, n_tr_occ, n_tr_occ};

                                // block for t1_ai
                                block t1_ai_low {a_low, 0};
                                block t1_ai_up {a_up, n_tr_occ};

                                v3("a,b,c,i,j,k") += g_abij("b,c,j,k").block(g_bcjk_low,g_bcjk_up)*t1("a,i").block(t1_ai_low,t1_ai_up);
                            }

                            // compute offset
                            std::size_t a_offset = tr_vir.tile(a).first;
                            std::size_t b_offset = tr_vir.tile(b).first;
                            std::size_t c_offset = tr_vir.tile(c).first;
//                            std::cout << a_offset << " " << b_offset << " " << c_offset << std::endl;
                            std::array<std::size_t,6> offset{a_offset,b_offset,c_offset,0,0,0};

                            double tmp_energy = 0.0;
                            if ( b < a && c < b){

                                auto ccsd_t_reduce = CCSD_T_Reduce(
                                        std::make_shared<Eigen::VectorXd>(this->orbital_energy_),
                                        this->trange1_engine_->get_actual_occ(), offset);
                                tmp_energy =
                                        (
                                                (t3("a,b,c,i,j,k")
                                                 + v3("a,b,c,i,j,k")
                                                )
                                                * (4.0 * t3("a,b,c,i,j,k")
                                                   + t3("a,b,c,k,i,j")
                                                   + t3("a,b,c,j,k,i")
                                                   -2*(t3("a,b,c,k,j,i")+t3("a,b,c,i,k,j")+t3("a,b,c,j,i,k"))
                                                )
                                        ).reduce(ccsd_t_reduce);

                                tmp_energy *= 2;
                            }else{

                                auto ccsd_t_reduce = CCSD_T_ReduceSymm(
                                        std::make_shared<Eigen::VectorXd>(this->orbital_energy_),
                                        this->trange1_engine_->get_actual_occ(), offset);

                                tmp_energy =
                                        (
                                                (t3("a,b,c,i,j,k")
                                                 + v3("a,b,c,i,j,k")
                                                )
                                                * (4.0 * t3("a,b,c,i,j,k")
                                                   + t3("a,b,c,k,i,j")
                                                   + t3("a,b,c,j,k,i")
                                                   -2*(t3("a,b,c,k,j,i")
                                                       +t3("a,b,c,i,k,j")
                                                       +t3("a,b,c,j,i,k"))
                                                )
                                        ).reduce(ccsd_t_reduce);
                            }

                            triple_energy += tmp_energy;
                        }
                    }
                    if(t1.get_world().rank() == 0){
                        print_progress(a, a+1, n_tr_vir);
                    }
                }

                if(t1.get_world().rank() == 0){
                    std::cout << "Total Blocks Computed  " << n_blocks_computed;
                    std::cout << " from " << std::pow(n_tr_occ,3)*std::pow(n_tr_vir,3) << std::endl;
                }
                return  triple_energy;

            }


        private:

            struct ReduceBase{
                typedef double result_type;
                typedef Tile argument_type;

                std::shared_ptr<Eig::VectorXd> vec_;
                unsigned int n_occ_;
                std::array<std::size_t,6> offset_;

                ReduceBase(std::shared_ptr<Eig::VectorXd> vec, int n_occ, std::array<std::size_t,6> offset)
                : vec_(std::move(vec)), n_occ_(n_occ) , offset_(offset){ }

                ReduceBase(ReduceBase const &) = default;

                result_type operator()() const { return 0.0; }

                result_type operator()(result_type const &t) const { return t; }

                void operator()(result_type &me, result_type const &other) const {
                    me += other;
                }
            };

            struct CCSD_T_Reduce : public ReduceBase{
                typedef typename ReduceBase::result_type result_type;
                typedef typename ReduceBase::argument_type argument_type ;

                CCSD_T_Reduce(std::shared_ptr<Eig::VectorXd> vec, int n_occ, std::array<std::size_t,6> offset)
                : ReduceBase(vec,n_occ,offset){ }

                CCSD_T_Reduce(CCSD_T_Reduce const &) = default;

                using ReduceBase::operator();

                void operator()(result_type &me, argument_type const &tile) const {

                    auto const &ens = *this->vec_;
                    std::size_t n_occ = this->n_occ_;
                    auto offset_ = this->offset_;

                    const auto a0 = tile.range().lobound()[0];
                    const auto an = tile.range().upbound()[0];
                    const auto b0 = tile.range().lobound()[1];
                    const auto bn = tile.range().upbound()[1];
                    const auto c0 = tile.range().lobound()[2];
                    const auto cn = tile.range().upbound()[2];
                    const auto i0 = tile.range().lobound()[3];
                    const auto in = tile.range().upbound()[3];
                    const auto j0 = tile.range().lobound()[4];
                    const auto jn = tile.range().upbound()[4];
                    const auto k0 = tile.range().lobound()[5];
                    const auto kn = tile.range().upbound()[5];

                    // get the offset
                    const auto a_offset = offset_[0];
                    const auto b_offset = offset_[1];
                    const auto c_offset = offset_[2];
                    const auto i_offset = offset_[3];
                    const auto j_offset = offset_[4];
                    const auto k_offset = offset_[5];

                    auto tile_idx = 0;
                    for (auto a = a0; a < an; ++a) {
                        const auto e_a = ens[a + n_occ + a_offset];
                        for (auto b = b0; b < bn; ++b) {
                            const auto e_b = ens[b + n_occ + b_offset];
                            for(auto c = c0; c < cn; ++c){
                                const auto e_c = ens[c + n_occ + c_offset];
                                for (auto i = i0; i < in; ++i) {
                                    const auto e_i = ens[i + i_offset];
                                    for (auto j = j0; j < jn; ++j) {
                                        const auto e_j = ens[j + j_offset];
                                        for (auto k = k0; k < kn; ++k, ++tile_idx){
                                            const auto e_k = ens[k + k_offset];

                                            const auto e_abcijk = e_i + e_j + e_k - e_a - e_b - e_c;

                                            me += (1.0/e_abcijk) * tile[tile_idx];

                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }; // structure CCSD_T_Reduce


            struct CCSD_T_ReduceSymm : public ReduceBase {
                typedef typename ReduceBase::result_type result_type;
                typedef typename ReduceBase::argument_type argument_type ;

                CCSD_T_ReduceSymm(std::shared_ptr<Eig::VectorXd> vec, int n_occ, std::array<std::size_t,6> offset)
                : ReduceBase(vec,n_occ,offset){ }

                CCSD_T_ReduceSymm(CCSD_T_ReduceSymm const &) = default;

                using ReduceBase::operator();

                void operator()(result_type &me, argument_type const &tile) const {

                    auto const &ens = *this->vec_;
                    std::size_t n_occ = this->n_occ_;

                    // get the offset
                    const auto a_offset = this->offset_[0];
                    const auto b_offset = this->offset_[1];
                    const auto c_offset = this->offset_[2];
                    const auto i_offset = this->offset_[3];
                    const auto j_offset = this->offset_[4];
                    const auto k_offset = this->offset_[5];

                    // compute index in the whole array
                    const auto a0 = tile.range().lobound()[0] + a_offset;
                    const auto an = tile.range().upbound()[0] + a_offset;
                    const auto b0 = tile.range().lobound()[1] + b_offset;
                    const auto bn = tile.range().upbound()[1] + b_offset;
                    const auto nb = bn - b0;
                    const auto c0 = tile.range().lobound()[2] + c_offset;
                    const auto cn = tile.range().upbound()[2] + c_offset;
                    const auto nc = cn - c0;
                    const auto i0 = tile.range().lobound()[3] + i_offset;
                    const auto in = tile.range().upbound()[3] + i_offset;
                    const auto ni = in - i0;
                    const auto j0 = tile.range().lobound()[4] + j_offset;
                    const auto jn = tile.range().upbound()[4] + j_offset;
                    const auto nj = jn - j0;
                    const auto k0 = tile.range().lobound()[5] + k_offset;
                    const auto kn = tile.range().upbound()[5] + k_offset;
                    const auto nk = kn - k0;

                    const auto njk = nj*nk;
                    const auto nijk = ni*njk;
                    const auto ncijk = nc*nijk;
                    const auto nbcijk = nb*ncijk;

                    typename Tile::value_type tmp = 0.0;

                    // use symmetry in loop, only sum result over c <= b <= a
                    for (auto a = a0; a < an; ++a) {
                        const auto e_a = ens[a + n_occ];
                        for (auto b = b0; b < bn && b <= a; ++b) {
                            const auto e_b = ens[b + n_occ];
                            for(auto c = c0; c < cn && c <= b; ++c){
                                const auto e_c = ens[c + n_occ];
                                for (auto i = i0; i < in; ++i) {
                                    const auto e_i = ens[i];
                                    for (auto j = j0; j < jn; ++j) {
                                        const auto e_j = ens[j];
                                        for (auto k = k0; k < kn; ++k){
                                            const auto e_k = ens[k];

                                            const auto tile_idx = (a-a0)*nbcijk + (b-b0)*ncijk + (c-c0)*nijk + (i-i0)*njk + (j-j0)*nk + (k-k0);

                                            const auto e_abcijk = e_i + e_j + e_k - e_a - e_b - e_c;

                                            tmp = (1.0/e_abcijk) * tile[tile_idx];
                                            // 6 fold symmetry if none in a,b,c equal
                                            if (a!=b && a!=c && b!=c) {
                                                tmp = 2.0 * tmp;
                                            }
                                                // if diagonal, a==b==c no symmetry
                                            else if(a==b && b==c){
                                                tmp = 0;
                                            }
                                                // three fold symmetry if two in a,b,c equal
                                            else{
                                                tmp = tmp;
                                            }
                                            me += tmp;

                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }; // structure CCSD_T_ReduceSymm

        }; // class CCSD_T

    } // namespace cc
}  // namespace mpqc

#endif //TILECLUSTERCHEM_CCSD_T_H_H
