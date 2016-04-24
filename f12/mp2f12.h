//
// Created by Chong Peng on 3/31/16.
//

#ifndef MPQC_MP2F12_H
#define MPQC_MP2F12_H

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
            : mo_int_(mo_int), tre_(tre), orbital_energy_(std::make_shared<Eigen::VectorXd>(ens)) {}

    double compute_mp2_f12_c_df();

    double compute_mp2_f12_c();

    double compute(const rapidjson::Document& in){

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

        return mp2_f12_energy;
    }

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

private:

    MolecularIntegral& mo_int_;
    std::shared_ptr<TRange1Engine> tre_;
    std::shared_ptr<Eigen::VectorXd> orbital_energy_;
};

}// end of namespce f12
} // end of namespace mpqc



#endif //MPQC_MP2F12_H
