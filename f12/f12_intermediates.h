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
 * \f$V_{ij}^{ij}\f$  \f$V_{ij}^{ji}\f$
 * @param mo_integral reference to MolecularIntegral, has to use SparsePolicy
 * @param shape SparseShape that has ijij ijji shape
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

    utility::print_par(world, "\nCompute V_ijij_ijji With DF \n" );
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
        auto left = mo_integral(L"<i1 j1|G|p q>[df]");
        auto right = mo_integral(L"<i2 j2|R|p q>[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_ijij_ijji("i1,j1,i2,j2") -= (left*right).set_shape(shape);
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term2 Time: ", time, " S\n");
    }

    {
        auto left = mo_integral(L"<i1 j1|G|m a'>[df]");
        auto right = mo_integral(L"<i2 j2|R|m a'>[df]");

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
 * \f$X_{ij}^{ij}\f$  \f$X_{ij}^{ji}\f$
 * @param mo_integral reference to MolecularIntegral, has to use SparsePolicy
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
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

    utility::print_par(world, "\nCompute X_ijij_ijji With DF \n" );
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
        auto left = mo_integral(L"<i1 j1|R|p q>[df]");
        auto right = mo_integral(L"<i2 j2|R|p q>[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        X_ijij_ijji("i1,j1,i2,j2") -= (left*right).set_shape(ijij_ijji_shape);
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"X Term2 Time: ", time, " S\n");
    }


    {
        auto left = mo_integral(L"<i1 j1|R|m a'>[df]");
        auto right = mo_integral(L"<i2 j2|R|m a'>[df]");

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


/**
 * MP2-F12 C approach B term, only ijij ijji part is computed
 * \f$B_{ij}^{ij}\f$  \f$B_{ij}^{ji}\f$
 * @param mo_integral reference to MolecularIntegral, has to use SparsePolicy
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return B(i1,j1,i2,j2)
 */
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

    utility::print_par(world, "\nCompute B_ijij_ijji With DF \n");

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
        auto hJ = mo_integral(L"<P' | hJ | i2>[df]");
        auto left = mo_integral(L"<i1 j1|R2|P' j2>[df]");

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
        auto left = mo_integral(L"<i1 j1|R|P' Q'>[df]");
        auto middle = mo_integral(L"<P'|K|R'>[df]");
        auto right = mo_integral(L"<i2 j2|R|R' Q'>[df]");

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
        auto left = mo_integral(L"<i1 j1|R|P' m>[df]");
        auto middle = mo_integral(L"<P'|F|R'>[df]");
        auto right = mo_integral(L"<i2 j2|R|R' m>[df]");

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
        auto left = mo_integral(L"<i1 j1|R|m b'>[df]");
        auto middle = mo_integral(L"<m|F|P'>[df]");
        auto right = mo_integral(L"<i2 j2|R|P' b'>[df]");

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
        auto left = mo_integral(L"<i1 j1|R|p a>[df]");
        auto middle = mo_integral(L"<p|F|r>[df]");
        auto right = mo_integral(L"<i2 j2|R|r a>[df]");

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
        auto left = mo_integral(L"<i1 j1|R|m b'>[df]");
        auto middle = mo_integral(L"<m|F|n>[df]");
        auto right = mo_integral(L"<i2 j2|R|n b'>[df]");

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
        auto left = mo_integral(L"<i1 j1|R|p a>[df]");
        auto middle = mo_integral(L"<p|F|a'>[df]");
        auto right = mo_integral(L"<i2 j2|R|a' a>[df]");

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
 * \f$V_{ia}^{xy}\f$
 * @param mo_integral reference to MolecularIntegral
 * @return V("i,a,x,y")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_V_iaxy_df(integrals::MolecularIntegral <Tile, Policy> &mo_integral)
{

    auto& world = mo_integral.get_world();
    bool accurate_time = true;
    TA::DistArray<Tile,Policy> V_iaxy;

    auto v_time0 = mpqc_time::now(world,accurate_time);

    utility::print_par(world, "\nCompute V_iaxy With DF \n" );
    {
        auto left = mo_integral(L"(Κ |GR|i k)");
        auto middle = mo_integral(L"(Κ|GR|Λ)[inv]");
        auto right = mo_integral(L"(Λ |GR|a l)");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_iaxy("i,a,k,l") = left*middle*right;
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term1 Time: ", time, " S\n");
    }

    {
        auto left = mo_integral(L"<i a|G|p q>[df]");
        auto right = mo_integral(L"<k l|R|p q>[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_iaxy("i,a,k,l") -= left*right;
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term2 Time: ", time, " S\n");
    }

    {
        auto left = mo_integral(L"<i a|G|m a'>[df]");
        auto right = mo_integral(L"<k l|R|m a'>[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_iaxy("i,a,k,l") -= left*right;
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term3 Time: ", time, " S\n");
    }

    {
        auto left = mo_integral(L"<a i|G|m a'>[df]");
        auto right = mo_integral(L"<l k|R|m a'>[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_iaxy("i,a,k,l") -= left*right;
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term4 Time: ", time, " S\n");
    }


    auto v_time1 = mpqc_time::now(world,accurate_time);
    auto v_time = mpqc_time::duration_in_s(v_time0,v_time1);
    utility::print_par(world,"V Term Total Time: ", v_time, " S\n");
    return V_iaxy;

};


/**
 * CC-F12 C approach V term
 * \f$V_{xy}^{ab}\f$
 * @param mo_integral reference to MolecularIntegral
 * @return V("x,y,a,b")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_V_xyab_df(integrals::MolecularIntegral <Tile, Policy> &mo_integral)
{

    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    bool accurate_time = true;

    auto v_time0 = mpqc_time::now(world,accurate_time);

    TA::DistArray<Tile,Policy> V_xyab;
    TA::DistArray<Tile,Policy> tmp;

    utility::print_par(world, "\nCompute V_xyab With DF \n" );

    {
        auto left = mo_integral(L"(Κ |GR|i a)");
        auto middle = ao_integral(L"(Κ|GR|Λ)[inv]");
        auto right = mo_integral(L"(Λ |GR|j b)");


        auto time0 = mpqc_time::now(world,accurate_time);
        V_xyab("i,j,a,b") = left*middle*right;
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term1 Time: ", time, " S\n");

    }

    {
        auto left = mo_integral(L"<a b|G|p q>[df]");
        auto right = mo_integral(L"<i j|R|p q>[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_xyab("i,j,a,b") -= left*right;
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"V Term2 Time: ", time, " S\n");

    }

//    tmp("i,j,a,b") = mo_integral(L"(a m|G|b a')[df]")*mo_integral(L"(i m|R|j a')[df]");
//    V_xyab("i,j,a,b") -= tmp("i,j,a,b");
//    V_xyab("i,j,a,b") -= tmp("j,i,b,a");

    auto v_time1 = mpqc_time::now(world,accurate_time);
    auto v_time = mpqc_time::duration_in_s(v_time0,v_time1);
    utility::print_par(world,"V Term Total Time: ", v_time, " S\n");

    return V_xyab;
};

/**
 * MP2-F12, CC-F12 C approach C term \f$C_{ij}^{ab} \f$
 * @param mo_integral reference to MolecularIntegral
 * @return C("i,j,a,b")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_C_ijab_df(integrals::MolecularIntegral <Tile, Policy> &mo_integral)
{
    auto& world = mo_integral.get_world();
    bool accurate_time = true;
    auto c_time0 = mpqc_time::now(world,accurate_time);
    TA::DistArray<Tile,Policy> C_ijab;

    utility::print_par(world, "\nCompute C_ijab With DF \n" );

    auto left = mo_integral(L"<i j|R|a a'>[df]");
    auto right = mo_integral(L"<a'|F|b>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    C_ijab("i,j,a,b") = left*right;
    C_ijab("i,j,a,b") += C_ijab("j,i,b,a");
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"C Term Time: ", time, " S\n");

    auto c_time = mpqc_time::duration_in_s(c_time0,time1);
    utility::print_par(world,"C Term Total Time: ", c_time, " S\n");
    return C_ijab;
};


/**
 * CC-F12 C approach VT2 term with direct integral
 * \f$T_{ab}^{ij} * (V_{xy}^{ab} + C_{xy}^{ab})\f$
 * @param mo_integral reference to MolecularIntegral
 * @param t2 t2 amplitude
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @param direct_array direct two electron integral \f$V_{\rho \sigma}^{\mu \nu}\f$
 * @return V("i1,j1,i2,j2")
 */
template<typename Tile, typename DirectArray>
TA::DistArray<Tile,TA::SparsePolicy> compute_VT2_ijij_ijji_df_direct(integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral,
                                                           const TA::DistArray <Tile, TA::SparsePolicy> &t2,
                                                           const TA::SparseShape<float> &ijij_ijji_shape,
                                                           DirectArray direct_array)
{

    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    bool accurate_time = true;

    TA::DistArray<Tile,TA::SparsePolicy> V_ijji_ijji;

    //C Term
    auto V_xyab = compute_C_ijab_df(mo_integral);

    utility::print_par(world, "\nCompute VT2_ijij_ijji With DF and Direct AO\n" );
    auto vt2_time0 = mpqc_time::now(world,accurate_time);
    {
        auto left = mo_integral(L"(Κ |GR|i a)");
        auto middle = ao_integral(L"(Κ|GR|Λ)[inv]");
        auto right = mo_integral(L"(Λ |GR|j b)");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_xyab("i,j,a,b") += left*middle*right;
        mo_integral.registry().remove_operation(world,L"GR");
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"VT2 Term1 Time: ", time, " S\n");
    }

    {
        auto right = t2("a,b,i1,j1");
        auto left = V_xyab("i2,j2,a,b");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_ijji_ijji("i1,j1,i2,j2") = (left*right).set_shape(ijij_ijji_shape);
        // clean V_xyab
        V_xyab = TA::DistArray<Tile,TA::SparsePolicy>();
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"VT2 Term2 Time: ", time, " S\n");
    }


    TA::DistArray<Tile,TA::SparsePolicy> U;
    auto Ca = mo_integral.orbital_space()->retrieve(OrbitalIndex(L"a")).array();
    auto Cp = mo_integral.orbital_space()->retrieve(OrbitalIndex(L"p")).array();
    {
        auto time0 = mpqc_time::now(world,accurate_time);

//    auto Cm = mo_integral.orbital_space()->retrieve(OrbitalIndex(L"m")).array();
//    auto Ca_prime = mo_integral.orbital_space()->retrieve(OrbitalIndex(L"a'")).array();
        // compuate intermediate U
        U("i,j,rho,mu") = (t2("a,b,i,j")*Ca("sigma,a")*Ca("nu,b"))*direct_array("rho, sigma, mu, nu");

        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"VT2 UTerm Time: ", time, " S\n");
    }


    {
        auto left = mo_integral(L"<i2 j2|R|p q>[df]");

        auto time0 = mpqc_time::now(world,accurate_time);
        V_ijji_ijji("i1,j1,i2,j2") -= ((left*Cp("rho,p")*Cp("mu,q"))*U("i1,j1,rho,mu")).set_shape(ijij_ijji_shape);
        auto time1 = mpqc_time::now(world,accurate_time);
        auto time = mpqc_time::duration_in_s(time0,time1);
        utility::print_par(world,"VT2 Term3 Time: ", time, " S\n");
    }
//    tmp("i1,j1,i2,j2") = ((mo_integral(L"(i2 m|R|j2 a')[df]")*Cm("mu,m")*Ca_prime("nu,a'"))*U("i1,j1,mu,nu")).set_shape(ijij_ijji_shape);
//    V_ijji_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
//    V_ijji_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    auto vt2_time1 = mpqc_time::now(world,accurate_time);
    auto vt2_time = mpqc_time::duration_in_s(vt2_time0,vt2_time1);
    utility::print_par(world,"VT2 Term Total Time: ", vt2_time, " S\n");

    return V_ijji_ijji;
};

/**
 * CC-F12 C approach VT2 term
 * \f$T_{ab}^{ij} * (V_{xy}^{ab} + C_{xy}^{ab})\f$
 * @param mo_integral reference to MolecularIntegral
 * @param t2 t2 amplitude
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return V("i1,j1,i2,j2")
 */

template<typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_VT2_ijij_ijji_df(
        integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral,
        const TA::DistArray <Tile, TA::SparsePolicy> &t2,
        const TA::SparseShape<float> &ijij_ijji_shape)
{
    auto& world = mo_integral.get_world();
    bool accurate_time = true;

    TA::DistArray<Tile,TA::SparsePolicy> V_ijij_ijji;

    // compute C_ijab
    TA::DistArray<Tile,TA::SparsePolicy> C_ijab = compute_C_ijab_df(mo_integral);

    // compute V_ijab
    TA::DistArray<Tile,TA::SparsePolicy> V_ijab = compute_V_xyab_df(mo_integral);

    auto vt2_time0 = mpqc_time::now(world,accurate_time);
    utility::print_par(world, "\nCompute VT2_ijij_ijji With DF\n" );
    V_ijij_ijji("i1,j1,i2,j2") = ((V_ijab("i2,j2,a,b")+C_ijab("i2,j2,a,b"))*t2("a,b,i1,j1")).set_shape(ijij_ijji_shape);

    auto vt2_time1 = mpqc_time::now(world,accurate_time);
    auto vt2_time = mpqc_time::duration_in_s(vt2_time0,vt2_time1);
    utility::print_par(world,"VT2 Term Total Time: ", vt2_time, " S\n");

    return V_ijij_ijji;
};



/**
 * CC-F12 C approach VT1 term
 * \f$T_{a}^{i} * V_{ia}^{xy}\f$
 * @param mo_integral reference to MolecularIntegral
 * @param t1 t1 amplitude
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return V("i1,j1,i2,j2")
 */

template<typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_VT1_ijij_ijji_df(
        integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral,
        const TA::DistArray <Tile, TA::SparsePolicy> &t1,
        const TA::SparseShape<float> &ijij_ijji_shape)
{
    auto& world = mo_integral.get_world();
    bool accurate_time = true;
    TA::DistArray<Tile,TA::SparsePolicy> V_ijij_ijji;
    TA::DistArray<Tile,TA::SparsePolicy> V_iaij = compute_V_iaxy_df(mo_integral);

    auto vt2_time0 = mpqc_time::now(world,accurate_time);
    utility::print_par(world, "\nCompute VT1_ijij_ijji With DF\n" );

    V_ijij_ijji("i1,j1,i2,j2") = V_iaij("i1,a,i2,j2")*t1("a,j1");
    V_ijij_ijji("i1,j1,i2,j2") += V_ijij_ijji("j1,i1,i2,j2");
//    V_ijij_ijji("i1,j1,i2,j2") += V_iaij("j1,a,i2,j2")*t1("a,i1");
    auto vt2_time1 = mpqc_time::now(world,accurate_time);
    auto vt2_time = mpqc_time::duration_in_s(vt2_time0,vt2_time1);
    utility::print_par(world,"VT1 Term Total Time: ", vt2_time, " S\n");

    return V_ijij_ijji;
};


} // end of namespace f12
} // end of namespace mpqc