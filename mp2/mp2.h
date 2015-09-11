//
// Created by Chong Peng on 6/24/15.
//

#ifndef TILECLUSTERCHEM_MP2_H
#define TILECLUSTERCHEM_MP2_H

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../ta_routines/deep_filter.h"
#include "../ta_routines/array_to_eigen.h"
#include "../ta_routines/tarray_block.h"
#include "../cc/trange1_engine.h"
#include "../cc/mo_block.h"


using namespace tcc;

namespace tcc {


    template<typename Tile, typename Policy>
    class MP2 {

    public:

        typedef TA::Array<double, 2, Tile, Policy> TArray2;
        typedef TA::Array<double, 3, Tile, Policy> TArray3;
        typedef TA::Array<double, 4, Tile, Policy> TArray4;

        MP2() { };

        MP2(const TArray2 &fock, const TArray2 &s_ab, const TArray3 &Xab,
            const std::shared_ptr<TRange1Engine> tre) : trange1_engine_(tre) {

            // initialize intergral g
            init(fock, s_ab, Xab);
        };

        void compute() {

            // compute mp2 energy
            double energy_mp2 = (g_("i,a,j,b") * (2 * g_("i,a,j,b") - g_("i,b,j,a"))).reduce(Mp2Red(orbital_energy_, trange1_engine_->get_actual_occ()));

            if (g_.get_world().rank() == 0) {
                std::cout << "MP2 Energy  " << energy_mp2 << std::endl;
            }
        }

        const TArray4 &get_g() const {
            return g_;
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
                for (auto i = st[0]; i < fn[0]; ++i) {
                    const auto e_i = vec[i];
                    for (auto a = st[1]; a < fn[1]; ++a) {
                        const auto e_ia = e_i - vec[a + n_occ_];
                        for (auto j = st[2]; j < fn[2]; ++j) {
                            const auto e_iaj = e_ia + vec[j];
                            for (auto b = st[3]; b < fn[3]; ++b, ++tile_idx) {
                                const auto e_iajb = e_iaj - vec[b + n_occ_];
                                me += 1 / (e_iajb) * tile.data()[tile_idx];
                            }
                        }
                    }
                }
            }
        };

        //template <typename Tile, typename Policy>
        //void init(const TArray2& fock, const TArray2& s_mn, const TArray4& mnkl);

        // template <typename Tile, typename Policy>
        void init(const TArray2 &fock, const TArray2 &s_mn,
                  const TArray3 &Xmn) {

            auto fock_eig = tcc::array_ops::array_to_eigen(fock);
            auto s_mn_eig = tcc::array_ops::array_to_eigen(s_mn);

            Eigen::GeneralizedSelfAdjointEigenSolver<decltype(s_mn_eig)> es(
                    fock_eig, s_mn_eig);

            std::size_t n_frozen_core = trange1_engine_->get_nfrozen();
            std::size_t occupation = trange1_engine_->get_occ();
            std::size_t block_size = trange1_engine_->get_block_size();

            Eigen::VectorXd evals = es.eigenvalues().bottomRows(s_mn_eig.rows() - n_frozen_core);
            auto C_all = es.eigenvectors();

            // compute mo coefficient
            auto C_occ = C_all.block(0, n_frozen_core, s_mn_eig.rows(), occupation  - n_frozen_core);
            auto C_vir = C_all.rightCols(s_mn_eig.rows() - occupation);

            // compute mo blocking
            auto tr_0 = Xmn.trange().data().back();
            auto tr_occ = trange1_engine_->get_occ_tr1();
            auto tr_vir = trange1_engine_->get_vir_tr1();

            utility::print_par(fock.get_world(), "Block Size in MO     ", block_size, "\n");
            utility::print_par(fock.get_world(), "TiledRange1 Occupied ", tr_occ, "\n");
            utility::print_par(fock.get_world(), "TiledRange1 Virtual  ", tr_vir, "\n");

            // convert mo coefficient to TiledArray
            auto Ci = array_ops::eigen_to_array<TA::Tensor < double>> (fock.get_world(), C_occ, tr_0, tr_occ);

            auto Cv = array_ops::eigen_to_array<TA::Tensor < double>> (fock.get_world(), C_vir, tr_0, tr_vir);

            TArray3 Xmn_mo;
            Xmn_mo("X,i,a") = Xmn("X,mu,nu") * Ci("mu,i") * Cv("nu,a");
            // construct two electron mo (ia|jb)
            g_("i,a,j,b") = Xmn_mo("X,i,a") * Xmn_mo("X,j,b");

            // set energy
            orbital_energy_ = std::make_shared<Eigen::VectorXd>(std::move(evals));
        };


    private:
        // two electron mo (ia|jb)
        TArray4 g_;
        std::shared_ptr<Eigen::VectorXd> orbital_energy_;
        std::shared_ptr<tcc::TRange1Engine> trange1_engine_;
    };

}
#endif //TILECLUSTERCHEM_MP2_H
