//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC_DB_F12_INTERMEDIATES_H
#define MPQC_DB_F12_INTERMEDIATES_H


#include "../../../../../include/tiledarray.h"
#include "../../../../../common/namespaces.h"
#include <mpqc/chemistry/qc/f12/f12_intermediates.h>
#include <mpqc/chemistry/qc/integrals/molecular_integral.h>

/**
 * Dual Basis F12 Intermediates
 */

namespace mpqc{
namespace f12{

template<typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_V_ijij_ijji_db_df(
        integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral, TA::SparseShape<float> &shape)
{
  auto& world = mo_integral.get_world();
  bool accurate_time = mo_integral.accurate_time();
  auto& ao_integral = mo_integral.atomic_integral();
  auto v_time0 = mpqc_time::now(world,accurate_time);

  TA::DistArray<Tile,TA::SparsePolicy> V_ijij_ijji;

  utility::print_par(world, "\nCompute V_ijij_ijji With Dual Basis DF \n" );
  {

    auto left = mo_integral(L"(Κ |GR|i2 i1)");
    auto middle = ao_integral(L"(Κ|GR|Λ)[inv]");
    auto right = mo_integral(L"(Λ |GR|j1 j2)");

    auto time0 = mpqc_time::now(world,accurate_time);

    V_ijij_ijji("i1,j1,i2,j2") = (left*middle*right).set_shape(shape);

//    std::cout << V_ijij_ijji << std::endl;
    // all types of GR integral not needed
    mo_integral.remove_operation_all(world, L"GR");

    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term1 Time: ", time, " S\n");
  }


  {
    auto left = mo_integral(L"<i1 j1|G|m n>[df]");
    auto right = mo_integral(L"<i2 j2|R|m n>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_ijij_ijji("i1,j1,i2,j2") -= (left*right).set_shape(shape);
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term2 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i1 j1|G|m a>[df]");
    auto right = mo_integral(L"<i2 j2|R|m a>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    TA::DistArray<Tile,TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left*right).set_shape(shape);
    V_ijij_ijji("i1,j1,i2,j2") -= (tmp("i1,j1,i2,j2")).set_shape(shape);
    V_ijij_ijji("i1,j1,i2,j2") -= (tmp("j1,i1,j2,i2")).set_shape(shape);
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term3 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i1 j1|G|a b>[df]");
    auto right = mo_integral(L"<i2 j2|R|a b>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_ijij_ijji("i1,j1,i2,j2") -= (left*right).set_shape(shape);
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term4 Time: ", time, " S\n");
  }

//  std::cout << V_ijij_ijji << std::endl;

  {
    auto left = mo_integral(L"<i1 j1|G|m a'>[df]");
    auto right = mo_integral(L"<i2 j2|R|m a'>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    TA::DistArray<Tile,TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left*right).set_shape(shape);
//    std::cout << tmp << std::endl;
//    V_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 m|G|i1 a')[df]")*mo_integral(L"(j2 m|R|i2 a')[df]")).set_shape(shape);
    V_ijij_ijji("i1,j1,i2,j2") -= (tmp("i1,j1,i2,j2")).set_shape(shape);
    V_ijij_ijji("i1,j1,i2,j2") -= (tmp("j1,i1,j2,i2")).set_shape(shape);
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term5 Time: ", time, " S\n");
  }


  auto v_time1 = mpqc_time::now(world,accurate_time);
  auto v_time = mpqc_time::duration_in_s(v_time0,v_time1);
  utility::print_par(world,"V Term Total Time: ", v_time, " S\n");

//  std::cout << V_ijij_ijji << std::endl;
  return V_ijij_ijji;
};


template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_V_xyab_db_df(integrals::MolecularIntegral <Tile, Policy> &mo_integral)
{

  auto& world = mo_integral.get_world();
  auto& ao_integral = mo_integral.atomic_integral();
  bool accurate_time = mo_integral.accurate_time();

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
    auto left = mo_integral(L"<a b|G|c d>[df]");
    auto right = mo_integral(L"<i j|R|c d>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_xyab("i,j,a,b") -= left*right;
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term2 Time: ", time, " S\n");

  }

  {
    auto left = mo_integral(L"<m1 m2|G|a b>[df]");
    auto right = mo_integral(L"<i j|R|m1 m2>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_xyab("i,j,a,b") -= left*right;
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term3 Time: ", time, " S\n");

  }

  {
    auto left = mo_integral(L"<a b|G|m c>[df]");
    auto right = mo_integral(L"<i j|R|m c>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    tmp("i,j,a,b") = left*right;
    V_xyab("i,j,a,b") -= tmp("i,j,a,b");
    V_xyab("i,j,a,b") -= tmp("j,i,b,a");
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term4 Time: ", time, " S\n");

  }

  {
    auto right = mo_integral(L"<a b|G|m a'>[df]");
    auto left = mo_integral(L"<i j|R|m a'>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    tmp("i,j,a,b") = left*right;
    V_xyab("i,j,a,b") -= tmp("i,j,a,b");
    V_xyab("i,j,a,b") -= tmp("j,i,b,a");
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term5 Time: ", time, " S\n");
  }

  auto v_time1 = mpqc_time::now(world,accurate_time);
  auto v_time = mpqc_time::duration_in_s(v_time0,v_time1);
  utility::print_par(world,"V Term Total Time: ", v_time, " S\n");

  return V_xyab;
};


template<typename Tile, typename Policy>
TA::DistArray<Tile,Policy> compute_V_iaxy_db_df(integrals::MolecularIntegral <Tile, Policy> &mo_integral)
{

  auto& world = mo_integral.get_world();
  bool accurate_time = mo_integral.accurate_time();
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
    auto left = mo_integral(L"<i a|G|m1 m2>[df]");
    auto right = mo_integral(L"<k l|R|m1 m2>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_iaxy("i,a,k,l") -= left*right;
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term2 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i a|G|b c>[df]");
    auto right = mo_integral(L"<k l|R|b c>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_iaxy("i,a,k,l") -= left*right;
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term3 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i a|G|m b>[df]");
    auto right = mo_integral(L"<k l|R|m b>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_iaxy("i,a,k,l") -= left*right;
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term4 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i a|G|b m>[df]");
    auto right = mo_integral(L"<k l|R|b m>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_iaxy("i,a,k,l") -= left*right;
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term5 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i a|G|m a'>[df]");
    auto right = mo_integral(L"<k l|R|m a'>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_iaxy("i,a,k,l") -= left*right;
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term6 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i a|G|a' m>[df]");
    auto right = mo_integral(L"<k l|R|a' m>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_iaxy("i,a,k,l") -= left*right;
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term7 Time: ", time, " S\n");
  }


  auto v_time1 = mpqc_time::now(world,accurate_time);
  auto v_time = mpqc_time::duration_in_s(v_time0,v_time1);
  utility::print_par(world,"V Term Total Time: ", v_time, " S\n");
  return V_iaxy;

};




template <typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_X_ijij_ijji_db_df(
        integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral, TA::SparseShape<float> &ijij_ijji_shape)
{

  bool accurate_time = mo_integral.accurate_time();
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
    auto left = mo_integral(L"<i1 j1|R|m n>[df]");
    auto right = mo_integral(L"<i2 j2|R|m n>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") -= (left*right).set_shape(ijij_ijji_shape);
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"X Term2 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i1 j1|R|m a>[df]");
    auto right = mo_integral(L"<i2 j2|R|m a>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    TA::DistArray<Tile,TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left*right).set_shape(ijij_ijji_shape);
    X_ijij_ijji("i1,j1,i2,j2") -= (tmp("i1,j1,i2,j2")).set_shape(ijij_ijji_shape);
    X_ijij_ijji("i1,j1,i2,j2") -= (tmp("j1,i1,j2,i2")).set_shape(ijij_ijji_shape);
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term3 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i1 j1|R|a b>[df]");
    auto right = mo_integral(L"<i2 j2|R|a b>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") -= (left * right).set_shape(ijij_ijji_shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term4 Time: ", time, " S\n");
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
    utility::print_par(world,"X Term5 Time: ", time, " S\n");
  }


  auto x_time1 = mpqc_time::now(world,accurate_time);
  auto x_time = mpqc_time::duration_in_s(x_time0,x_time1);
  utility::print_par(world,"X Term Total Time: ", x_time, " S\n");

//    std::cout << X_ijij_ijji << std::endl;
  return X_ijij_ijji;
};


template <typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_B_ijij_ijji_db_df(
        integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral, TA::SparseShape<float> &ijij_ijji_shape)
{
  bool accurate_time = mo_integral.accurate_time();
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
    auto left1 = mo_integral(L"<i1 j1|R2|m j2>[df]");
    auto hJ1 = mo_integral(L"<m | hJ | i2>[df]");

    auto left2 = mo_integral(L"<i1 j1|R2|A' j2>[df]");
    auto hJ2 = mo_integral(L"<A' | hJ | i2>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    tmp("i1,j1,i2,j2") = (left1*hJ1 + left2*hJ2).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") += (mo_integral(L"(j1 P'|R2|i1 i2)[df]")*hJ("P',j2")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");

    mo_integral.remove_operation_all(world, L"R2");
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"B Term2 Time: ", time, " S\n");
  }

  {
    auto left1 = mo_integral(L"<i1 j1|R|m1 A'>[df]");
    auto middle1 = mo_integral(L"<m1 |K| m2>[df]");
    auto right1 = mo_integral(L"<i2 j2|R|m2 A'>[df]");

    auto left2 = mo_integral(L"<i1 j1|R|m1 A'>[df]");
    auto middle2 = mo_integral(L"<m1 |K| B'>[df]");
    auto right2 = mo_integral(L"<i2 j2|R|B' A'>[df]");

    auto left3 = mo_integral(L"<i1 j1|R|A' B'>[df]");
    auto middle3 = mo_integral(L"<A' |K| C'>[df]");
    auto right3 = mo_integral(L"<i2 j2|R|C' B'>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    tmp("i1,j1,i2,j2") = (left1*middle1*right1 + 2.0*left2*middle2*right2 + left3*middle3*right3).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    // AO R integral not needed
    mo_integral.atomic_integral().registry().remove_operation(world, L"R");
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"B Term3 Time: ", time, " S\n");
  }

  {
    auto left1 = mo_integral(L"<i1 j1|R|m1 m2>[df]");
    auto middle1 = mo_integral(L"<m1|hJ|m3>[df]");
    auto right1 = mo_integral(L"<i2 j2|R|m3 m2>[df]");

    auto left2 = mo_integral(L"<i1 j1|R|m1 A'>[df]");
    auto middle2 = mo_integral(L"<A'|hJ|m2>[df]");
    auto right2 = mo_integral(L"<i2 j2|R|m1 m2>[df]");

    auto left3 = mo_integral(L"<i1 j1|R|m A'>[df]");
    auto middle3 = mo_integral(L"<A'|hJ|B'>[df]");
    auto right3 = mo_integral(L"<i2 j2|R|m B'>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    tmp("i1,j1,i2,j2") = (left1*middle1*right1 + 2.0*left2*middle2*right2 + left3*middle3*right3).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 P'|R|i1 m)[df]")*mo_integral(L"(P'|F|R')[df]")*mo_integral(L"(j2 R'|R|i2 m)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    mo_integral.remove_operation_all(world, L"hJ");
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"B Term4 Time: ", time, " S\n");
  }



  {
    auto left1 = mo_integral(L"<i1 j1|R|m b'>[df]");
    auto middle1 = mo_integral(L"<m|F|n>[df]");
    auto right1 = mo_integral(L"<i2 j2|R|n b'>[df]");

    auto left2 = mo_integral(L"<i1 j1|R|m b'>[df]");
    auto middle2 = mo_integral(L"<m|F|A'>[df]");
    auto right2 = mo_integral(L"<i2 j2|R|A' b'>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    tmp("i1,j1,i2,j2") = (left1*middle1*right1 + 2.0*left2*middle2*right2).set_shape(ijij_ijji_shape);
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
    auto left1 = mo_integral(L"<i1 j1|R|m a>[df]");
    auto middle1 = mo_integral(L"<m|F|n>[df]");
    auto right1 = mo_integral(L"<i2 j2|R|n a>[df]");

    auto left2 = mo_integral(L"<i1 j1|R|m a>[df]");
    auto middle2 = mo_integral(L"<m|F|b>[df]");
    auto right2 = mo_integral(L"<i2 j2|R|b a>[df]");

    auto left3 = mo_integral(L"<i1 j1|R|b a>[df]");
    auto middle3 = mo_integral(L"<b|F|c>[df]");
    auto right3 = mo_integral(L"<i2 j2|R|c a>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    tmp("i1,j1,i2,j2") = (left1*middle1*right1 + 2*left2*middle2*right2 + left3*middle3*right3).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|r)[df]")*mo_integral(L"(j2 r|R|i2 a)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"B Term6 Time: ", time, " S\n");
  }

  {
    auto left1 = mo_integral(L"<i1 j1|R|m a>[df]");
    auto middle1 = mo_integral(L"<m|F|a'>[df]");
    auto right1 = mo_integral(L"<i2 j2|R|a' a>[df]");

    auto left2 = mo_integral(L"<i1 j1|R|b a>[df]");
    auto middle2 = mo_integral(L"<b|F|a'>[df]");
    auto right2 = mo_integral(L"<i2 j2|R|a' a>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    tmp("i1,j1,i2,j2") = (2.0*(left1*middle1*right1 + left2*middle2*right2)).set_shape(ijij_ijji_shape);
//    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*mo_integral(L"(j1 p|R|i1 a)[df]")*mo_integral(L"(p|F|a')[df]")*mo_integral(L"(i2 a|R|j2 a')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"B Term7 Time: ", time, " S\n");
  }


  auto b_time1 = mpqc_time::now(world,accurate_time);
  auto b_time = mpqc_time::duration_in_s(b_time0,b_time1);
  utility::print_par(world,"B Term Total Time: ", b_time, " S\n");

//    std::cout << B_ijij_ijji << std::endl;

  return B_ijij_ijji;
};



template<typename Tile>
TA::DistArray<Tile,TA::SparsePolicy> compute_VT2_ijij_ijji_db_df(
        integrals::MolecularIntegral <Tile, TA::SparsePolicy> &mo_integral,
        const TA::DistArray <Tile, TA::SparsePolicy> &t2,
        const TA::SparseShape<float> &ijij_ijji_shape)
{
  auto& world = mo_integral.get_world();
  bool accurate_time = mo_integral.accurate_time();

  TA::DistArray<Tile,TA::SparsePolicy> V_ijij_ijji;

  // compute C_ijab
  TA::DistArray<Tile,TA::SparsePolicy> C_ijab = compute_C_ijab_df(mo_integral);

  // compute V_ijab
  TA::DistArray<Tile,TA::SparsePolicy> V_ijab = compute_V_xyab_db_df(mo_integral);

  auto vt2_time0 = mpqc_time::now(world,accurate_time);
  utility::print_par(world, "\nCompute VT2_ijij_ijji With DF\n" );
  V_ijij_ijji("i1,j1,i2,j2") = ((V_ijab("i2,j2,a,b")+C_ijab("i2,j2,a,b"))*t2("a,b,i1,j1")).set_shape(ijij_ijji_shape);

  auto vt2_time1 = mpqc_time::now(world,accurate_time);
  auto vt2_time = mpqc_time::duration_in_s(vt2_time0,vt2_time1);
  utility::print_par(world,"VT2 Term Total Time: ", vt2_time, " S\n");

  return V_ijij_ijji;
};


}
}


#endif //MPQC_DB_F12_INTERMEDIATES_H
