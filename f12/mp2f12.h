//
// Created by Chong Peng on 3/31/16.
//

#ifndef TILECLUSTERCHEM_MP2F12_H
#define TILECLUSTERCHEM_MP2F12_H

#include <string>

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../integrals/molecular_integral.h"
#include "../utility/trange1_engine.h"
//#include "../expression/formula.h"

namespace mpqc{
namespace f12{

class MP2F12 {

public:
    using Tile = TA::TensorD;
    using Policy = TA::SparsePolicy;
    using TArray = TA::DistArray<Tile, Policy>;
    using MolecularIntegral = integrals::MolecularIntegral<Tile,Policy>;

    MP2F12() = default;

    MP2F12(MolecularIntegral& mo_int, std::shared_ptr<TRange1Engine> tre, const Eigen::VectorXd& ens)
            : mo_int_(mo_int), tre_(tre), orbital_energy_(ens) {}

    void compute_mp2_f12_c();

private:

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
            for (auto i = st[0]; i < fn[0]; ++i) {
                for (auto j = st[1]; j < fn[1]; ++j) {
                    for (auto k = st[2]; k < fn[2]; ++k) {
                        for (auto l = st[3]; l < fn[3]; ++l, ++tile_idx) {
                            // iiii
                            if( (i==j) && (i==k) && (k==l)){
                                me += iiii*tile.data()[tile_idx];
                            }
                            // ijij
                            else if(i > j && k==i && l==j){
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

private:

    MolecularIntegral& mo_int_;
    std::shared_ptr<TRange1Engine> tre_;
    Eigen::VectorXd orbital_energy_;
};

}// end of namespce f12
} // end of namespace mpqc



#endif //TILECLUSTERCHEM_MP2F12_H
