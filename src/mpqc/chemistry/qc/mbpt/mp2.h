//
// Created by Chong Peng on 6/24/15.
//

#ifndef MPQC_MP2_H
#define MPQC_MP2_H

#include "../../../../../include/tiledarray.h"
#include "../../../../../common/namespaces.h"
#include <mpqc/chemistry/qc/integrals/molecular_integral.h>
#include <mpqc/chemistry/qc/scf/mo_build.h>
#include "../../../../../utility/trange1_engine.h"
#include "../../../../../utility/parallel_print.h"


using namespace mpqc;

namespace mpqc{
namespace mbpt{



    template<typename Tile, typename Policy>
    class MP2 {
    public:
        using TArray = TA::DistArray<Tile,Policy>;
        using MolecularIntegral = integrals::MolecularIntegral<Tile,Policy>;


      /// constructor using MO Integral with orbitals computed
        MP2(MolecularIntegral& mo_int, const Eigen::VectorXd& orbital_energy, const std::shared_ptr<TRange1Engine> tre)
                : mo_int_(mo_int), trange1_engine_(tre), orbital_energy_(std::make_shared<Eigen::VectorXd>(orbital_energy)) {}

      /// constructfor using MO Integral without orbitals computed
        MP2(MolecularIntegral& mo_int, const rapidjson::Document& in) : mo_int_(mo_int) {
            auto& ao_int = mo_int.atomic_integral();
            auto orbital_registry = mo_int.orbital_space();
            auto mol = mo_int.atomic_integral().molecule();
            int occ = mol.occupation(0)/2;
            Eigen::VectorXd orbital_energy;
            trange1_engine_ = closed_shell_obs_mo_build_eigen_solve(ao_int, *orbital_registry, orbital_energy, in, mol, occ);
            orbital_energy_ = std::make_shared<Eigen::VectorXd>(orbital_energy);
        }

        MP2() = default;

        double compute_df() {

            auto g_ijab = mo_int_.compute(L"<i j|G|a b>[df]");
            // compute mp2 energy
            double energy_mp2 = (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a"))).reduce(Mp2Energy(orbital_energy_, trange1_engine_->get_active_occ()));

            if (g_ijab.get_world().rank() == 0) {
                std::cout << "MP2 Energy With DF: " << energy_mp2 << std::endl;
            }

            return energy_mp2;
        }

        double compute_four_center() {

            auto g_ijab = mo_int_.compute(L"<i j|G|a b>");
            // compute mp2 energy
            double energy_mp2 = (g_ijab("i,j,a,b") * (2 * g_ijab("i,j,a,b") - g_ijab("i,j,b,a"))).reduce(Mp2Energy(orbital_energy_, trange1_engine_->get_active_occ()));

            if (g_ijab.get_world().rank() == 0) {
                std::cout << "MP2 Energy  " << energy_mp2 << std::endl;
            }

            return energy_mp2;
        }

        double compute(const rapidjson::Document& in){

            std::string method = in.HasMember("Method") ? in["Method"].GetString() : "df";

            double mp2_energy = 0.0;

            if(method == "four center"){
                mp2_energy = compute_four_center();
            }
            else if(method == "df"){
                mp2_energy = compute_df();
            }
            else{
                throw std::runtime_error("Wrong MP2 Method");
            }

            return mp2_energy;
        }

        const std::shared_ptr<Eigen::VectorXd> get_en() const {
            return orbital_energy_;
        }

    private:

        struct Mp2Energy {
            using result_type = double;
            using argument_type = Tile;

            std::shared_ptr<Eig::VectorXd> vec_;
            unsigned int n_occ_;

            Mp2Energy(std::shared_ptr<Eig::VectorXd> vec, int n_occ)
                    : vec_(std::move(vec)), n_occ_(n_occ) { }

            Mp2Energy(Mp2Energy const &) = default;

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
                auto stj = st[1];
                auto fnj = fn[1];
                auto sta = st[2];
                auto fna = fn[2];
                auto stb = st[3];
                auto fnb = fn[3];

                for (auto i = sti; i < fni; ++i) {
                    const auto e_i = vec[i];
                    for (auto j = stj; j < fnj; ++j) {
                        const auto e_ij = e_i + vec[j];
                        for (auto a = sta; a < fna; ++a) {
                            const auto e_ija = e_ij - vec[a + n_occ_];
                            for (auto b = stb; b < fnb; ++b, ++tile_idx) {
                                const auto e_iajb = e_ija - vec[b + n_occ_];
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

} // end of namespace mbpt
}//end of namespace mpqc
#endif //MPQC_MP2_H
