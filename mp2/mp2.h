//
// Created by Chong Peng on 6/24/15.
//

#ifndef TILECLUSTERCHEM_MP2_H
#define TILECLUSTERCHEM_MP2_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../integrals/molecular_integral.h"
#include "../utility/trange1_engine.h"
#include "../utility/parallel_print.h"


using namespace mpqc;

namespace mpqc{


    template<typename Tile, typename Policy>
    class MP2 {
    public:
        using TArray = TA::DistArray<Tile,Policy>;
        using MolecularIntegral = integrals::MolecularIntegral<Tile,Policy>;


        MP2(MolecularIntegral& mo_int, const Eigen::VectorXd& orbital_energy, const std::shared_ptr<TRange1Engine> tre)
                : mo_int_(mo_int), orbital_energy_(std::make_shared<Eigen::VectorXd>(orbital_energy)), trange1_engine_(tre){}

        MP2() = default;

        void compute_df() {

            auto g_iajb = mo_int_.compute(L"(i a|G|j b)[df]");
            // compute mp2 energy
            double energy_mp2 = (g_iajb("i,a,j,b") * (2 * g_iajb("i,a,j,b") - g_iajb("i,b,j,a"))).reduce(Mp2Red(orbital_energy_, trange1_engine_->get_actual_occ()));

            if (g_iajb.get_world().rank() == 0) {
                std::cout << "MP2 Energy With DF: " << energy_mp2 << std::endl;
            }
        }

        void compute() {

            auto g_iajb = mo_int_.compute(L"(i a|G|j b)");
            // compute mp2 energy
            double energy_mp2 = (g_iajb("i,a,j,b") * (2 * g_iajb("i,a,j,b") - g_iajb("i,b,j,a"))).reduce(Mp2Red(orbital_energy_, trange1_engine_->get_actual_occ()));

            if (g_iajb.get_world().rank() == 0) {
                std::cout << "MP2 Energy  " << energy_mp2 << std::endl;
            }
        }

        const std::shared_ptr<Eigen::VectorXd> get_en() const {
            return orbital_energy_;
        }

    private:

        struct Mp2Red {
            using result_type = double;
            using argument_type = Tile;

            std::shared_ptr<Eig::VectorXd> vec_;
            unsigned int n_occ_;

            Mp2Red(std::shared_ptr<Eig::VectorXd> vec, int n_occ)
                    : vec_(std::move(vec)), n_occ_(n_occ) { }

            Mp2Red(Mp2Red const &) = default;

            result_type operator()() const { return 0.0; }

            result_type operator()(result_type const &t) const { return t; }

            void operator()(result_type &me, result_type const &other) const {
                me += other;
            }

            void operator()(result_type &me, argument_type const &tile) const {
                auto const &range = tile.range();
                auto const &vec = *vec_;
                auto const st = range.lobound_data();
                auto const fn = range.upbound_data();
                auto tile_idx = 0;

                auto sti = st[0];
                auto fni = fn[0];
                auto sta = st[1];
                auto fna = fn[1];
                auto stj = st[2];
                auto fnj = fn[2];
                auto stb = st[3];
                auto fnb = fn[3];

                for (auto i = sti; i < fni; ++i) {
                    const auto e_i = vec[i];
                    for (auto a = sta; a < fna; ++a) {
                        const auto e_ia = e_i - vec[a + n_occ_];
                        for (auto j = stj; j < fnj; ++j) {
                            const auto e_iaj = e_ia + vec[j];
                            for (auto b = stb; b < fnb; ++b, ++tile_idx) {
                                const auto e_iajb = e_iaj - vec[b + n_occ_];
                                me += 1 / (e_iajb) * tile.data()[tile_idx];
                            }
                        }
                    }
                }
            }
        };

    private:
        MolecularIntegral& mo_int_;
        std::shared_ptr<mpqc::TRange1Engine> trange1_engine_;
        std::shared_ptr<Eigen::VectorXd> orbital_energy_;
    };

}
#endif //TILECLUSTERCHEM_MP2_H
