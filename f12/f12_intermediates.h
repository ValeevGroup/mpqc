//
// Created by Chong Peng on 04/11/16.
//

#include "../include/tiledarray.h"
#include "../common/namespaces.h"
#include "../integrals/molecular_integral.h"

namespace mpqc{
namespace f12{

/**
 * MP2-F12 C approach V term, only ijij ijji part is computed
 * $V_{ij}^{ij}$  $V_{ij}^{ji}$
 * @param mo_integral reference to MolecularIntegral, has to use SparsePolicy
 * @param SparseShape that has ijij ijji shape
 * @return V(i1,j1,i2,j2)
 */
template<typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_V_ijij_ijji(integrals::MolecularIntegral<Tile, TA::SparsePolicy>& mo_integral, TA::SparseShape<float>& shape)
{
    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();

    TA::DistArray<Tile,TA::SparsePolicy> V_ijij_ijji;
    TA::DistArray<Tile,TA::SparsePolicy> tmp;

    utility::print_par(world, "Compute V_ijij_ijji With DF \n" );
    V_ijij_ijji("i1,j1,i2,j2") = (mo_integral(L"(Κ |GR|i2 i1)")*ao_integral(L"(Κ|GR|Λ)[inv]")*mo_integral(L"(Λ |GR|j1 j2)")).set_shape(shape);

    // all types of GR integral not needed
    mo_integral.remove_operation_all(world, L"GR");

    V_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|G|j1 q)[df]")*mo_integral(L"(i2 p|R|j2 q)[df]")).set_shape(shape);
    tmp("i1,j1,i2,j2") = (mo_integral(L"(i1 m|G|j1 a')[df]")*mo_integral(L"(i2 m|R|j2 a')[df]")).set_shape(shape);
//    V_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|G|i1 a')[df]")*mo_integral(L"(j2 m|R|i2 a')[df]")).set_shape(shape);
    V_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    V_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

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
TA::DistArray<Tile,TA::SparsePolicy> compute_X_ijij_ijji(integrals::MolecularIntegral<Tile, TA::SparsePolicy>& mo_integral, TA::SparseShape<float>& ijij_ijji_shape)
{

    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    TA::DistArray<Tile,TA::SparsePolicy> X_ijij_ijji;
    TA::DistArray<Tile,TA::SparsePolicy> tmp;

    utility::print_par(world, "Compute X_ijij_ijji With DF \n" );

    X_ijij_ijji("i1,j1,i2,j2") = (mo_integral(L"(Κ |R2|i1 i2)")*ao_integral(L"(Κ|R2|Λ)[inv]")*mo_integral(L"(Λ |R2|j1 j2)")).set_shape(ijij_ijji_shape);
    X_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 q)[df]")*mo_integral(L"(i2 p|R|j2 q)[df]")).set_shape(ijij_ijji_shape);
    tmp("i1,j1,i2,j2") = (mo_integral(L"(i1 m|R|j1 a')[df]")*mo_integral(L"(i2 m|R|j2 a')[df]")).set_shape(ijij_ijji_shape);
//    X_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|R|i1 a')[df]")*mo_integral(L"(j2 m|R|i2 a')[df]")).set_shape(ijij_ijji_shape);
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    return X_ijij_ijji;
};

/**
 * MP2-F12 C approach B term, only ijij ijji part is computed
 * $B_{ij}^{ij}$ $B_{ij}^{ji}$
 * @param mo_integral reference to MolecularIntegral, has to use SparsePolicy
 * @param SparseShape that has ijij ijji shape
 * @return B(i1,j1,i2,j2)
 */
template <typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_B_ijij_ijji_old(integrals::MolecularIntegral<Tile, TA::SparsePolicy>& mo_integral, TA::SparseShape<float>& ijij_ijji_shape)
{
    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    TA::DistArray<Tile,TA::SparsePolicy> B_ijij_ijji;
    TA::DistArray<Tile,TA::SparsePolicy> tmp;

    utility::print_par(world, "Compute B_ijij_ijji With DF \n");

    B_ijij_ijji("i1,j1,i2,j2") = (mo_integral(L"(Κ |dR2|i1 i2)")*ao_integral(L"(Κ|dR2|Λ)[inv]")*mo_integral(L"(Λ |dR2|j1 j2)")).set_shape(ijij_ijji_shape);

    mo_integral.remove_operation_all(world, L"dR2");
    auto hJ = mo_integral.compute(L"(P' | hJ | i)[df]");
    B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(i1 P'|R2|j1 j2)[df]")*hJ("P',i2")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(j1 P'|R2|i1 i2)[df]")*hJ("P',j2")).set_shape(ijij_ijji_shape);

    mo_integral.remove_operation_all(world, L"R2");
    mo_integral.remove_operation_all(world, L"hJ");

    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(i2 R'|R|j2 Q')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(j2 R'|R|i2 Q')[df]")).set_shape(ijij_ijji_shape);

    // AO R integral not needed
    mo_integral.atomic_integral().registry().remove_operation(world, L"R");

    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 P'|R|j1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(i2 R'|R|j2 m)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(j2 R'|R|i2 m)[df]")).set_shape(ijij_ijji_shape);


    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(i2 P'|R|j2 b')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(j2 P'|R|i2 b')[df]")).set_shape(ijij_ijji_shape);

    // P' doesn't appear later
    mo_integral.registry().remove_orbital(world, L"P'");

    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(i2 r|R|j2 a)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(j2 r|R|i2 a)[df]")).set_shape(ijij_ijji_shape);

    B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(i2 n|R|j2 b')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(j2 n|R|i2 b')[df]")).set_shape(ijij_ijji_shape);


    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(j2 a|R|i2 a')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(i2 a|R|j2 a')[df]")).set_shape(ijij_ijji_shape);

    return B_ijij_ijji;
};


template <typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_B_ijij_ijji(integrals::MolecularIntegral<Tile, TA::SparsePolicy>& mo_integral, TA::SparseShape<float>& ijij_ijji_shape)
{
    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    TA::DistArray<Tile,TA::SparsePolicy> B_ijij_ijji;
    TA::DistArray<Tile,TA::SparsePolicy> tmp;

    utility::print_par(world, "Compute B_ijij_ijji With DF \n");

    B_ijij_ijji("i1,j1,i2,j2") = (mo_integral(L"(Κ |dR2|i1 i2)")*ao_integral(L"(Κ|dR2|Λ)[inv]")*mo_integral(L"(Λ |dR2|j1 j2)")).set_shape(ijij_ijji_shape);

    mo_integral.remove_operation_all(world, L"dR2");
    auto hJ = mo_integral.compute(L"(P' | hJ | i)[df]");
    tmp("i1,j1,i2,j2") = (mo_integral(L"(i1 P'|R2|j1 j2)[df]")*hJ("P',i2")).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(j1 P'|R2|i1 i2)[df]")*hJ("P',j2")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");

    mo_integral.remove_operation_all(world, L"R2");
    mo_integral.remove_operation_all(world, L"hJ");

    tmp("i1,j1,i2,j2") = (mo_integral(L"(i1 P'|R|j1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(i2 R'|R|j2 Q')[df]")).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(j2 R'|R|i2 Q')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    // AO R integral not needed
    mo_integral.atomic_integral().registry().remove_operation(world, L"R");

    tmp("i1,j1,i2,j2") = (mo_integral(L"(i1 P'|R|j1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(i2 R'|R|j2 m)[df]")).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(j2 R'|R|i2 m)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");


    tmp("i1,j1,i2,j2") = (2.0*mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(i2 P'|R|j2 b')[df]")).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(j2 P'|R|i2 b')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    // P' doesn't appear later
    mo_integral.registry().remove_orbital(world, L"P'");

    tmp("i1,j1,i2,j2") = (mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(i2 r|R|j2 a)[df]")).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(j2 r|R|i2 a)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    tmp("i1,j1,i2,j2") = (mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(i2 n|R|j2 b')[df]")).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(j2 n|R|i2 b')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");


    tmp("i1,j1,i2,j2") = (2.0*mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(j2 a|R|i2 a')[df]")).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(i2 a|R|j2 a')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    return B_ijij_ijji;
};

/**
 * CC-F12 C approach V term
 * $V_{ij}^{xy}
 * @param mo_integral reference to MolecularIntegral
 * @return V("i,j,x,y")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_V_ijxy(integrals::MolecularIntegral<Tile, Policy>& mo_integral)
{

    auto& world = mo_integral.get_world();

    TA::DistArray<Tile,Policy> V_ijxy;

    utility::print_par(world, "Compute V_ijxy With DF \n" );
    V_ijxy("i,j,k,l") = mo_integral(L"(Κ |GR|i k)")*mo_integral(L"(Κ|GR|Λ)[inv]")*mo_integral(L"(Λ |GR|j l)");

    // all types of GR integral not needed
    mo_integral.remove_operation_all(world, L"GR");

    V_ijxy("i,j,k,l") -= mo_integral(L"(i p|G|j q)[df]")*mo_integral(L"(k p|R|l q)[df]");
    V_ijxy("i,j,k,l") -= mo_integral(L"(i m|G|j a')[df]")*mo_integral(L"(k m|R|l a')[df]");
    V_ijxy("i,j,k,l") -= mo_integral(L"(j m|G|i a')[df]")*mo_integral(L"(k m|R|l a')[df]");

    return V_ijxy;

};

/**
 * CC-F12 C approach V term
 * $V_{xy}^{ab}
 * @param mo_integral reference to MolecularIntegral
 * @return V("x,y,a,b")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_V_xyab(integrals::MolecularIntegral<Tile, Policy>& mo_integral)
{

    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();

    TA::DistArray<Tile,Policy> V_xyab;

    utility::print_par(world, "Compute V_xyab With DF \n" );
    V_xyab("i,j,a,b") = mo_integral(L"(Κ |GR|i a)")*ao_integral(L"(Κ|GR|Λ)[inv]")*mo_integral(L"(Λ |GR|j b)");

    // all types of GR integral not needed
    mo_integral.remove_operation_all(world, L"GR");

    V_xyab("i,j,a,b") -= mo_integral(L"(i p|G|j q)[df]")*mo_integral(L"(a p|R|b q)[df]");
    V_xyab("i,j,a,b") -= mo_integral(L"(i m|G|j a')[df]")*mo_integral(L"(a m|R|b a')[df]");
    V_xyab("i,j,a,b") -= mo_integral(L"(j m|G|i a')[df]")*mo_integral(L"(b m|R|a a')[df]");

    return V_xyab;
};

/**
 * CC-F12 C approach X term
 * $X_{xy}^{wz}
 * @param mo_integral reference to MolecularIntegral
 * @return V("x,y,w,z")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_X_xywz(integrals::MolecularIntegral<Tile, Policy>& mo_integral)
{
    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    TA::DistArray<Tile,Policy> X_xywz;

    utility::print_par(world, "Compute X_xywz With DF \n" );

    X_xywz("i1,j1,i2,j2") = mo_integral(L"(Κ |R2|i1 i2)")*ao_integral(L"(Κ|R2|Λ)[inv]")*mo_integral(L"(Λ |R2|j1 j2)");
    X_xywz("i1,j1,i2,j2") -= mo_integral(L"(i1 p|R|j1 q)[df]")*mo_integral(L"(i2 p|R|j2 q)[df]");
    X_xywz("i1,j1,i2,j2") -= mo_integral(L"(i1 m|R|j1 a')[df]")*mo_integral(L"(i2 m|R|j2 a')[df]");
    X_xywz("i1,j1,i2,j2") -= mo_integral(L"(j1 m|R|i1 a')[df]")*mo_integral(L"(j2 m|R|i2 a')[df]");

    return X_xywz;

};

/**
 * CC-F12 C approach B term
 * $B_{xy}^{wz}
 * @param mo_integral reference to MolecularIntegral
 * @return B("x,y,w,z")
 */
template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_B_xywz(integrals::MolecularIntegral<Tile, Policy>& mo_integral)
{

    auto& world = mo_integral.get_world();
    auto& ao_integral = mo_integral.atomic_integral();
    TA::DistArray<Tile,TA::SparsePolicy> B_xywz;

    utility::print_par(world, "Compute B_xywz With DF \n");

    B_xywz("i1,j1,i2,j2") = mo_integral(L"(Κ |dR2|i1 i2)")*ao_integral(L"(Κ|dR2|Λ)[inv]")*mo_integral(L"(Λ |dR2|j1 j2)");

    mo_integral.remove_operation_all(world, L"dR2");
    auto hJ = mo_integral.compute(L"(P' | hJ | i)[df]");
    B_xywz("i1,j1,i2,j2") += mo_integral(L"(i1 P'|R2|j1 j2)[df]")*hJ("P',i2");
    B_xywz("i1,j1,i2,j2") += mo_integral(L"(j1 P'|R2|i1 i2)[df]")*hJ("P',j2");

    mo_integral.remove_operation_all(world, L"R2");
    mo_integral.remove_operation_all(world, L"hJ");

    B_xywz("i1,j1,i2,j2") -= mo_integral(L"(i1 P'|R|j1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(i2 R'|R|j2 Q')[df]");
    B_xywz("i1,j1,i2,j2") -= mo_integral(L"(j1 P'|R|i1 Q')[df]")*mo_integral(L"(P'|K|R')[df]")*mo_integral(L"(j2 R'|R|i2 Q')[df]");

    // AO R integral not needed
    mo_integral.atomic_integral().registry().remove_operation(world, L"R");

    B_xywz("i1,j1,i2,j2") -= mo_integral(L"(i1 P'|R|j1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(i2 R'|R|j2 m)[df]");
    B_xywz("i1,j1,i2,j2") -= mo_integral(L"(j1 P'|R|i1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(j2 R'|R|i2 m)[df]");


    B_xywz("i1,j1,i2,j2") -= 2.0*mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(i2 P'|R|j2 b')[df]");
    B_xywz("i1,j1,i2,j2") -= 2.0*mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|P')[df]")*mo_integral(L"(j2 P'|R|i2 b')[df]");

    // P' doesn't appear later
    mo_integral.registry().remove_orbital(world, L"P'");

    B_xywz("i1,j1,i2,j2") -= mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(i2 r|R|j2 a)[df]");
    B_xywz("i1,j1,i2,j2") -= mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(j2 r|R|i2 a)[df]");

    B_xywz("i1,j1,i2,j2") += mo_integral(L"(i1 m|R|j1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(i2 n|R|j2 b')[df]");
    B_xywz("i1,j1,i2,j2") += mo_integral(L"(j1 m|R|i1 b')[df]")*mo_integral(L"(m|F|n)[df]")*mo_integral(L"(j2 n|R|i2 b')[df]");


    B_xywz("i1,j1,i2,j2") -= 2.0*mo_integral(L"(i1 p|R|j1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(j2 a|R|i2 a')[df]");
    B_xywz("i1,j1,i2,j2") -= 2.0*mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(i2 a|R|j2 a')[df]");

    return B_xywz;
};

} // end of namespace f12
} // end of namespace mpqc