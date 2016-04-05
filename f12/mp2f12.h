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

    TA::expressions::TsrExpr<TArray,true> mo_integral(const std::wstring& str){
        return mo_int_(str);
    };

private:

    MolecularIntegral& mo_int_;
    std::shared_ptr<TRange1Engine> tre_;
    Eigen::VectorXd orbital_energy_;
};

}// end of namespce f12
} // end of namespace mpqc



#endif //TILECLUSTERCHEM_MP2F12_H
