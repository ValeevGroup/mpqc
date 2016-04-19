//
// Created by Chong Peng on 04/11/16.
//

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../integrals/molecular_integral.h"

namespace mpqc{
namespace f12{


/**
 * MP2-F12 C approach V term with Density Fitting, only ijij ijji part is computed
 * $V_{ij}^{ij}$  $V_{ij}^{ji}$
 * @param mo_integral reference to MolecularIntegral, has to use SparsePolicy
 * @param SparseShape that has ijij ijji shape
 * @return V(i1,j1,i2,j2)
 */
template<typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_V_ijij_ijji_df(
        integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral, TA::SparseShape<float> &shape)
{
    bool accurate_time = true;
    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    auto v_time0 = mpqc_time::now(world,accurate_time);

    TA::DistArray<Tile,TA::SparsePolicy> V_ijij_ijji;

    utility::print_par(world, "Compute V_ijij_ijji With DF \n" );
    {

        auto left = mo_integral(L"(Κ |GR|i2 i1)");
        auto middle = ao_integral(L"(Κ|GR|Λ)[inv]");
        auto right = mo_integral(L"(Λ |GR|j1 j2)");

        auto time0 = mpqc_time::now(world,accurate_time);

        V_ijij_ijji("i1,j1,i2,j2") = (left*middle*right).set_shape(shape);
        // all types of GR integral not needed
        mo_integral.remove_operation_all(world, L"GR");

        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term1 Time: ", time, " S\n");
    }


    {
        auto left = mo_integral(L"(i1 p|G|j1 q)[df]");
        auto right = mo_integral(L"(i2 p|R|j2 q)[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_ijij_ijji("i1,j1,i2,j2") -= (left*right).set_shape(shape);
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term2 Time: ", time, " S\n");
    }

    {
        auto left = mo_integral(L"(i1 m|G|j1 a')[df]");
        auto right = mo_integral(L"(i2 m|R|j2 a')[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        TA::DistArray<Tile,TA::SparsePolicy> tmp;
        tmp("i1,j1,i2,j2") = (left*right).set_shape(shape);
//    V_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|G|i1 a')[df]")*mo_integral(L"(j2 m|R|i2 a')[df]")).set_shape(shape);
        V_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
        V_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term3 Time: ", time, " S\n");
    }


    auto v_time1 = mpqc_time::now(world,accurate_time);
    auto v_time = mpqc_time::duration_in_s(v_time0,v_time1);
    utility::print_par(world,"V Term Total Time: ", v_time, " S\n");
    return V_ijij_ijji;
};


/**
 * MP2-F12 C approach X term, only ijij ijji part is computed
 * $X_{ij}^{ij}$  $X_{ij}^{ji}$
 * @param mo_integral reference to MolecularIntegral, has to use SparsePolicy
 * @param SparseShape that has ijij ijji shape
 * @return X(i1,j1,i2,j2)
 */
template <typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_X_ijij_ijji_df(
        integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral, TA::SparseShape<float> &ijij_ijji_shape)
{

    bool accurate_time = true;
    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    auto x_time0 = mpqc_time::now(world,accurate_time);

    TA::DistArray<Tile,TA::SparsePolicy> X_ijij_ijji;

    utility::print_par(world, "Compute X_ijij_ijji With DF \n" );
    {
        auto left = mo_integral(L"(Κ |R2|i1 i2)");
        auto middle = ao_integral(L"(Κ|R2|Λ)[inv]");
        auto right = mo_integral(L"(Λ |R2|j1 j2)");

        auto time0 = mpqc_time::now(world,accurate_time);
        X_ijij_ijji("i1,j1,i2,j2") = (left*middle*right).set_shape(ijij_ijji_shape);
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"X Term1 Time: ", time, " S\n");
    }


    {
        auto left = mo_integral(L"(i1 p|R|j1 q)[df]");
        auto right = mo_integral(L"(i2 p|R|j2 q)[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        X_ijij_ijji("i1,j1,i2,j2") -= (left*right).set_shape(ijij_ijji_shape);
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"X Term2 Time: ", time, " S\n");
    }


    {
        auto left = mo_integral(L"(i1 m|R|j1 a')[df]");
        auto right = mo_integral(L"(i2 m|R|j2 a')[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        TA::DistArray<Tile,TA::SparsePolicy> tmp;
        tmp("i1,j1,i2,j2") = (left*right).set_shape(ijij_ijji_shape);
//    X_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|R|i1 a')[df]")*mo_integral(L"(j2 m|R|i2 a')[df]")).set_shape(ijij_ijji_shape);
        X_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
        X_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"X Term3 Time: ", time, " S\n");
    }


    auto x_time1 = mpqc_time::now(world,accurate_time);
    auto x_time = mpqc_time::duration_in_s(x_time0,x_time1);
    utility::print_par(world,"X Term Total Time: ", x_time, " S\n");
    return X_ijij_ijji;
};


template <typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_B_ijij_ijji_df(
        integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral, TA::SparseShape<float> &ijij_ijji_shape)
{
    bool accurate_time = true;
    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    auto b_time0 = mpqc_time::now(world,accurate_time);

    TA::DistArray<Tile,TA::SparsePolicy> B_ijij_ijji;
    TA::DistArray<Tile,TA::SparsePolicy> tmp;

    utility::print_par(world, "Compute B_ijij_ijji With DF \n");

    {
        auto left = mo_integral(L"(Κ |dR2|i1 i2)");
        auto middle = ao_integral(L"(Κ|dR2|Λ)[inv]");
        auto right = mo_integral(L"(Λ |dR2|j1 j2)");

        auto time0 = mpqc_time::now(world,accurate_time);
        B_ijij_ijji("i1,j1,i2,j2") = (left*middle*right).set_shape(ijij_ijji_shape);
        mo_integral.remove_operation_all(world, L"dR2");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"B Term1 Time: ", time, " S\n");
    }

    {
        auto hJ = mo_integral(L"(P' | hJ | i2)[df]");
        auto left = mo_integral(L"(i1 P'|R2|j1 j2)[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        tmp("i1,j1,i2,j2") = (left*hJ).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(j1 P'|R2|i1 i2)[df]")*hJ("P',j2")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
        B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");

        mo_integral.remove_operation_all(world, L"R2");
        mo_integral.remove_operation_all(world, L"hJ");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"B Term2 Time: ", time, " S\n");
    }

    {
        auto left = mo_integral(L"(i1 P'|R|j1 Q')[df]");
        auto middle = mo_integral(L"(P'|K|R')[df]");
        auto right = mo_integral(L"(i2 R'|R|j2 Q')[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        tmp("i1,j1,i2,j2") = (left*middle*right).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(j2 R'|R|i2 Q')[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
        // AO R integral not needed
        mo_integral.atomic_integral().registry().remove_operation(world, L"R");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"B Term3 Time: ", time, " S\n");
    }

    {
        auto left = mo_integral(L"(i1 P'|R|j1 m)[df]");
        auto middle = mo_integral(L"(P'|F|R')[df]");
        auto right = mo_integral(L"(i2 R'|R|j2 m)[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        tmp("i1,j1,i2,j2") = (left*middle*right).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(j2 R'|R|i2 m)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"B Term4 Time: ", time, " S\n");
    }



    {
        auto left = mo_integral(L"(i1 m|R|j1 b')[df]");
        auto middle = mo_integral(L"(m|F|P')[df]");
        auto right = mo_integral(L"(i2 P'|R|j2 b')[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        tmp("i1,j1,i2,j2") = (2.0*left*middle*right).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(j2 P'|R|i2 b')[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

        // P' doesn't appear later
        mo_integral.registry().remove_orbital(world, L"P'");

        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"B Term5 Time: ", time, " S\n");
    }


    {
        auto left = mo_integral(L"(i1 p|R|j1 a)[df]");
        auto middle = mo_integral(L"(p|F|r)[df]");
        auto right = mo_integral(L"(i2 r|R|j2 a)[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        tmp("i1,j1,i2,j2") = (left*middle*right).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(j2 r|R|i2 a)[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"B Term6 Time: ", time, " S\n");
    }


    {
        auto left = mo_integral(L"(i1 m|R|j1 b')[df]");
        auto middle = mo_integral(L"(m|F|n)[df]");
        auto right = mo_integral(L"(i2 n|R|j2 b')[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        tmp("i1,j1,i2,j2") = (left*middle*right).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(j2 n|R|i2 b')[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
        B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"B Term7 Time: ", time, " S\n");
    }

    {
        auto left = mo_integral(L"(i1 p|R|j1 a)[df]");
        auto middle = mo_integral(L"(p|F|a')[df]");
        auto right = mo_integral(L"(j2 a|R|i2 a')[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        tmp("i1,j1,i2,j2") = (2.0*left*middle*right).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(i2 a|R|j2 a')[df]")).set_shape(ijij_ijji_shape);
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
        B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"B Term8 Time: ", time, " S\n");
    }


    auto b_time1 = mpqc_time::now(world,accurate_time);
    auto b_time = mpqc_time::duration_in_s(b_time0,b_time1);
    utility::print_par(world,"B Term Total Time: ", b_time, " S\n");
    return B_ijij_ijji;
};

/**
 * CC-F12 C approach V term
 * $V_{ia}^{xy}
 * @param mo_integral reference to MolecularIntegral
 * @return V("i,a,x,y")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_V_iaxy_df(integrals::MolecularIntegral <Tile, Policy> &mo_integral)
{

    auto& world = mo_integral.get_world();

    TA::DistArray<Tile,Policy> V_iaxy;
//    TA::DistArray<Tile,Policy> tmp;

    utility::print_par(world, "Compute V_iaxy With DF \n" );
    V_iaxy("i,a,k,l") = mo_integral(L"(Κ |GR|i k)")*mo_integral(L"(Κ|GR|Λ)[inv]")*mo_integral(L"(Λ |GR|a l)");

    V_iaxy("i,a,k,l") -= mo_integral(L"(i p|G|a q)[df]")*mo_integral(L"(k p|R|l q)[df]");
    V_iaxy("i,a,k,l") -= mo_integral(L"(i m|G|a a')[df]")*mo_integral(L"(k m|R|l a')[df]");
    V_iaxy("i,a,k,l") -= mo_integral(L"(a m|G|i a')[df]")*mo_integral(L"(l m|R|k a')[df]");
//    V_iaxy("i,a,k,l") -= tmp("i,a,k,l");
//    V_iaxy("i,a,k,l") -= tmp("a,i,l,k");

    return V_iaxy;

};


/**
 * CC-F12 C approach V term
 * $V_{xy}^{ab}$
 * @param mo_integral reference to MolecularIntegral
 * @return V("x,y,a,b")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_V_xyab_df(integrals::MolecularIntegral <Tile, Policy> &mo_integral)
{

    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();

    TA::DistArray<Tile,Policy> V_xyab;
    TA::DistArray<Tile,Policy> tmp;

    utility::print_par(world, "Compute V_xyab With DF \n" );
    V_xyab("i,j,a,b") = mo_integral(L"(Κ |GR|i a)")*ao_integral(L"(Κ|GR|Λ)[inv]")*mo_integral(L"(Λ |GR|j b)");

    V_xyab("i,j,a,b") -= mo_integral(L"(a p|G|b q)[df]")*mo_integral(L"(i p|R|j q)[df]");
//    tmp("i,j,a,b") = mo_integral(L"(a m|G|b a')[df]")*mo_integral(L"(i m|R|j a')[df]");
//    V_xyab("i,j,a,b") -= tmp("i,j,a,b");
//    V_xyab("i,j,a,b") -= tmp("j,i,b,a");

    return V_xyab;
};

/**
 * MP2-F12, CC-F12 C approach C term
 * $C_{ij}^{ab}
 * @param mo_integral reference to MolecularIntegral
 * @return C("i,j,a,b")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_C_ijab_df(integrals::MolecularIntegral <Tile, Policy> &mo_integral)
{
    auto& world = mo_integral.get_world();
    bool accurate_time = true;
    TA::DistArray<Tile,Policy> C_ijab;

    utility::print_par(world, "Compute C_ijab With DF \n" );
    auto time0 = mpqc_time::now(world,accurate_time);
    C_ijab("i,j,a,b") = mo_integral(L"(i a|R|j a')[df]")*mo_integral(L"(b|F|a')[df]");
    C_ijab("i,j,a,b") += C_ijab("j,i,b,a");
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"C Term Time: ", time, " S\n");
    return C_ijab;
};


/**
 * CC-F12 C approach VT term
 * $T_{ab}^{ij}$ * $V_{xy}^{ab}$
 * @param mo_integral reference to MolecularIntegral
 * @return V("i,j,j,i")
 */
template<typename Tile, typename Policy, typename DirectArray>
TA::DistArray<Tile,Policy> compute_VT2_ijij_ijji_df(integrals::MolecularIntegral <Tile, Policy> &mo_integral,
                                                    const TA::DistArray <Tile, Policy> &t2,
                                                    const TA::SparseShape<float> &ijij_ijji_shape,
                                                    DirectArray direct_array)
{

    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();

    TA::DistArray<Tile,Policy> V_ijji_ijji;
    TA::DistArray<Tile,Policy> tmp;

    utility::print_par(world, "Compute VT2_ijij_ijji With DF \n" );

    //C Term
    auto V_xyab = compute_C_ijab_df(mo_integral);

    V_xyab("i,j,a,b") += mo_integral(L"(Κ |GR|i a)")*ao_integral(L"(Κ|GR|Λ)[inv]")*mo_integral(L"(Λ |GR|j b)");
    mo_integral.registry().remove_operation(world,L"GR");

    V_ijji_ijji("i1,j1,i2,j2") = (t2("a,b,i1,j1")*V_xyab("i2,j2,a,b")).set_shape(ijij_ijji_shape);

    // clean V_xyab
    V_xyab = TA::DistArray<Tile,Policy>();

    auto Ca = mo_integral.orbital_space()->retrieve(OrbitalIndex(L"a")).array();
//    auto Cm = mo_integral.orbital_space()->retrieve(OrbitalIndex(L"m")).array();
    auto Cp = mo_integral.orbital_space()->retrieve(OrbitalIndex(L"p")).array();
//    auto Ca_prime = mo_integral.orbital_space()->retrieve(OrbitalIndex(L"a'")).array();
    // compuate intermediate U
    TA::DistArray<Tile,Policy> U;
    U("i,j,rho,mu") = (t2("a,b,i,j")*Ca("sigma,a")*Ca("nu,b"))*direct_array("rho, sigma, mu, nu");


    V_ijji_ijji("i1,j1,i2,j2") -= ((mo_integral(L"(i2 p|R|j2 q)[df]")*Cp("rho,p")*Cp("mu,q"))*U("i1,j1,rho,mu")).set_shape(ijij_ijji_shape);
//    tmp("i1,j1,i2,j2") = ((mo_integral(L"(i2 m|R|j2 a')[df]")*Cm("mu,m")*Ca_prime("nu,a'"))*U("i1,j1,mu,nu")).set_shape(ijij_ijji_shape);
//    V_ijji_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
//    V_ijji_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    return V_ijji_ijji;
};



} // end of namespace f12
} // end of namespace mpqc