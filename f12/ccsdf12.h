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

template <typename Tile>
class CCSDF12 {
public:

    using Policy = TA::SparsePolicy;
    using TArray = TA::DistArray<Tile,Policy>;

    CCSDF12(integrals::MolecularIntegral<Tile, Policy> &mo_int_, const std::shared_ptr<TRange1Engine> &tre_,
            const std::shared_ptr<Eigen::VectorXd> &orbital_energy_, const TArray &t1_, const TArray &t2_)
            : mo_int_(mo_int_), tre_(tre_), orbital_energy_(orbital_energy_), t1_(t1_), t2_(t2_)
    { }

    /// standard approach, no approximation
    template<typename DirectArray>
    double compute_c(const DirectArray& darray);

private:

    integrals::MolecularIntegral<Tile, Policy>& mo_int_;
    std::shared_ptr<TRange1Engine> tre_;
    std::shared_ptr<Eigen::VectorXd> orbital_energy_;

    TArray t1_; /// t1 amplitude
    TArray t2_; /// t2 amplitude

};


template <typename Tile>
template <typename DirectArray>
double CCSDF12<Tile>::compute_c(const DirectArray& darray) {

    auto& mo_integral = mo_int_;
    auto& world = mo_integral.get_world();
    double E = 0.0;

    // create shape
    auto occ_tr1 = tre_->get_occ_tr1();
    TiledArray::TiledRange occ4_trange({occ_tr1,occ_tr1,occ_tr1,occ_tr1});
    auto ijij_ijji_shape = f12::make_ijij_ijji_shape(occ4_trange);

    // compute V_ijij_ijji
    TArray V_ijij_ijji = compute_V_ijij_ijji_df(mo_integral, ijij_ijji_shape);

    // VT2 contribution
    if(darray.is_initialized()){
        TArray tmp = compute_VT2_ijij_ijji_df(mo_integral, t2_, ijij_ijji_shape, darray);
        V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    }else{
        // compute C_ijab
        TArray C_ijab = compute_C_ijab_df(mo_integral);
        // compute V_ijab
        TArray V_ijab = compute_V_xyab_df(mo_integral);
        V_ijij_ijji("i1,j1,i2,j2") += ((V_ijab("i2,j2,a,b")+C_ijab("i2,j2,a,b"))*t2_("a,b,i1,j1")).set_shape(ijij_ijji_shape);
    }

    // compuate V_ijia
    TArray V_iaij = compute_V_iaxy_df(mo_integral);
    V_ijij_ijji("i1,j1,i2,j2") += V_iaij("i1,a,i2,j2")*t1_("a,j1");
    V_ijij_ijji("i1,j1,i2,j2") += V_iaij("j1,a,i2,j2")*t1_("a,i1");

    // V contribution to energy
    double E_v = V_ijij_ijji("i1,j1,i2,j2").reduce(f12::CLF12Energy<Tile>(1.0,2.5,-0.5));
    utility::print_par(world, "E_V: ", E_v, "\n");
    E += E_v;

    // compute X term
    TArray X_ijij_ijji = compute_X_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    // R_ipjq not needed
    mo_int_.registry().remove_formula(world, L"(i1 p|R|j1 q)[df]");

    auto Fij = mo_int_.compute(L"(i|F|j)[df]");
    auto Fij_eigen = array_ops::array_to_eigen(Fij);
    f12::convert_X_ijkl(X_ijij_ijji, Fij_eigen);

    double E_x = -X_ijij_ijji("i1,j1,i2,j2").reduce(f12::CLF12Energy<Tile>(0.25,0.4375,0.0625));

    utility::print_par(world, "E_X: ", E_x, "\n");
    E += E_x;

    // compute B term
    TArray B_ijij_ijji = compute_B_ijij_ijji_df(mo_integral, ijij_ijji_shape);
    double E_b = B_ijij_ijji("i1,j1,i2,j2").reduce(CLF12Energy<Tile>(0.25,0.4375,0.0625));
    utility::print_par(world, "E_B: ", E_b, "\n");
    E += E_b;

    utility::print_par(world, "E_F12: ", E, "\n");

    return E;
}

}//end of namespace f12
}//end of namespace mpqc


#endif //TILECLUSTERCHEM_CCSDF12_H
