//
// Created by Chong Peng on 8/21/15.
//

#ifndef TILECLUSTERCHEM_CCSD_T_H_H
#define TILECLUSTERCHEM_CCSD_T_H_H

#include "ccsd.h"

namespace tcc{
    namespace cc{

        // CCSD_T class that compute CCSD(T) triple calculation
        template<typename Tile, typename Policy>
        class CCSD_T : public CCSD<Tile,Policy> {

        public:

            typedef TA::Array <double, 2, Tile, Policy> TArray2;
            typedef TA::Array <double, 4, Tile, Policy> TArray4;
            typedef TA::Array <double, 6, Tile, Policy> TArray6;

            CCSD_T(const TArray2 &fock, const Eigen::VectorXd &ens,
                 const std::shared_ptr<TRange1Engine> &tre,
                 const std::shared_ptr<CCSDIntermediate<Tile, Policy>> &inter):
                    CCSD<Tile,Policy>(fock,ens,tre,inter)
            {}

            void compute(){

                TArray2 t1;
                TArray4 t2;

                // compute ccsd first
                double ccsd_corr = CCSD<Tile,Policy>::compute_ccsd(t1,t2);

                double ccsd_t = compute_ccsd_t_straight(t1, t2);

                if (t1.get_world().rank() == 0) {
                    std::cout << "(T) Energy      " << ccsd_t << std::endl;
                    std::cout << "CCSD(T) Energy  " << ccsd_t + ccsd_corr << std::endl;
                }

            }

            double compute_ccsd_t(const TArray2& t1, const TArray4& t2){
                return 0;
            }

            double compute_ccsd_t_straight(const TArray2& t1, const TArray4& t2){

                // get integral
                TArray4 g_jklc = this->ccsd_intermediate_->get_ijka();
                TArray4 g_diba = this->ccsd_intermediate_->get_aibc();
                TArray4 g_abij = this->ccsd_intermediate_->get_abij();

                // compute t3
                TArray6 t3;
                t3("a,b,c,i,j,k") = g_diba("d,i,b,a")*t2("c,d,k,j") - g_jklc("l,k,j,c")*t2("a,b,i,l");
                t3("a,b,c,i,j,k") = t3("a,b,c,i,j,k") + t3("a,c,b,i,k,j") + t3("c,a,b,k,i,j") + t3("c,b,a,k,j,i")
                        + t3("b,c,a,j,k,i") + t3("b,a,c,j,i,k");

                // compute v3
                TArray6 v3;
                v3("a,b,c,i,j,k") = g_abij("a,b,i,j")*t1("c,k");
                v3("a,b,c,i,j,k") = v3("a,b,c,i,j,k") + v3("b,c,a,j,k,i") + v3("a,c,b,i,k,j");

                double triple_energy =
                        (
                                (t3("a,b,c,i,j,k")
                                 + v3("a,b,c,i,j,k"))
                                * (4 * t3("a,b,c,i,j,k")
                                   + t3("a,b,c,k,i,j")
                                  + t3("a,b,c,j,k,i")
                                   - 4 * t3("a,b,c,k,j,i")
                                  - t3("a,b,c,i,k,j")
                                   - t3("a,b,c,j,i,k"))
                        ).reduce(CCSD_TRed(
                                std::make_shared<Eigen::VectorXd>(this->orbital_energy_),
                                this->trange1_engine_->get_actual_occ()));
                triple_energy = triple_energy/3.0;
                return triple_energy;
            }

            double compute_ccsd_t_direct(const TArray2& t1, const TArray4& t2){

                // get integral
                TArray4 g_jklc = this->intermediate_->get_ijka();
                TArray4 g_diba = this->intermediate_->get_aibc();
                TArray4 g_abij = this->intermediate_->get_abij();

                // get trange1
                auto tr_occ = this->tre_->get_occ_tr1();
                auto tr_vir = this->tre_->get_vir_tr1();

                auto n_tr_occ = this->tre_->get_occ_blocks();
                auto n_tr_vir = this->tre_->get_vir_blocks();


            }


        private:

            struct CCSD_TRed {
                using result_type = double;
                using argument_type = Tile;

                std::shared_ptr<Eig::VectorXd> vec_;
                unsigned int n_occ_;

                CCSD_TRed(std::shared_ptr<Eig::VectorXd> vec, int n_occ)
                        : vec_(std::move(vec)), n_occ_(n_occ) { }

                CCSD_TRed(CCSD_TRed const &) = default;

                result_type operator()() const { return 0.0; }

                result_type operator()(result_type const &t) const { return t; }

                void operator()(result_type &me, result_type const &other) const {
                    me += other;
                }

                void operator()(result_type &me, argument_type const &tile) const {

                    auto const &ens = *vec_;
                    std::size_t n_occ = n_occ_;

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

                    auto tile_idx = 0;
                    typename Tile::value_type tmp = 1.0;
                    for (auto a = a0; a < an; ++a) {
                        const auto e_a = ens[a + n_occ];
                        for (auto b = b0; b < bn; ++b) {
                            const auto e_b = ens[b + n_occ];
                            for(auto c = c0; c < cn; ++c){
                                const auto e_c = ens[c + n_occ];
                                for (auto i = i0; i < in; ++i) {
                                    const auto e_i = ens[i];
                                    for (auto j = j0; j < jn; ++j) {
                                        const auto e_j = ens[j];
                                        for (auto k = k0; k < kn; ++k, ++tile_idx){
                                            const auto e_k = ens[k];

                                            const auto e_abcijk = e_i + e_j + e_k - e_a - e_b - e_c;

                                            me += 1 / (e_abcijk)* tile[tile_idx];

                                        }
                                    }
                                }
                            }

                        }
                    }
                }
            }; // structure CCSD_TRed


        }; // class CCSD_T

    } // namespace cc
}  // namespace tcc

#endif //TILECLUSTERCHEM_CCSD_T_H_H
