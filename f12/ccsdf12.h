//
// Created by Chong Peng on 4/12/16.
//

#ifndef TILECLUSTERCHEM_CCSDF12_H
#define TILECLUSTERCHEM_CCSDF12_H


#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../integrals/molecular_integral.h"
#include "../utility/trange1_engine.h"
#include "f12_intermediates.h"

namespace mpqc{
namespace f12{

template <typename Tile, typename Policy>
class CCSDF12 {
public:

    using TArray = TA::DistArray<Tile,Policy>;

    CCSDF12(integrals::MolecularIntegral<Tile, Policy> &mo_int_, const std::shared_ptr<TRange1Engine> &tre_,
            const std::shared_ptr<Eigen::VectorXd> &orbital_energy_, const TArray &t1_, const TArray &t2_)
            : mo_int_(mo_int_), tre_(tre_), orbital_energy_(orbital_energy_), t1_(t1_), t2_(t2_)
    { }

    void compute_c();

private:

    integrals::MolecularIntegral<Tile, Policy>& mo_int_;
    std::shared_ptr<TRange1Engine> tre_;
    std::shared_ptr<Eigen::VectorXd> orbital_energy_;

    TArray t1_; /// t1 amplitude
    TArray t2_; /// t2 amplitude

};

template <typename Tile, typename Policy>
void CCSDF12<Tile,Policy>::compute_c() {

    auto& mo_integral = mo_int_;

    // compute C_ijab
    TArray C_xyab = compute_C_ijab(mo_integral);

    // compute V_ijab
    TArray V_xyab = compute_V_xyab(mo_integral);

    // compute V_ijxy
    TArray  V_ijxy = compute_V_ijxy(mo_integral);

    // compute V_iaxy
    TArray V_iaxy = compute_V_iaxy(mo_integral);

    TArray V_bar_ijxy;
    V_bar_ijxy = V_ijxy;
    V_bar_ijxy("i,j,x,y") += 0.5*(V_xyab("x,y,a,b")+C_xyab("x,y,a,b"))*t2_("a,b,i,j");
    V_bar_ijxy("i,j,x,y") += V_iaxy("i,a,x,y")*t1_("j,a");
    V_bar_ijxy("i,j,x,y") += V_iaxy("j,a,x,y")*t1_("i,a");


}

}//end of namespace f12
}//end of namespace mpqc


#endif //TILECLUSTERCHEM_CCSDF12_H
