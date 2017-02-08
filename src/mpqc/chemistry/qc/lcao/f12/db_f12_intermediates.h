//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DB_F12_INTERMEDIATES_H_
#define MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DB_F12_INTERMEDIATES_H_

#include <tiledarray.h>

#include "mpqc/chemistry/qc/lcao/f12/f12_intermediates.h"
#include "mpqc/chemistry/qc/lcao/factory/lcao_factory.h"

/**
 * Dual Basis F12 Intermediates
 */

namespace mpqc {
namespace lcao {
namespace f12 {

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_V_ijij_ijji_db_df(
    lcao::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &shape) {
  auto &world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto &ao_factory = lcao_factory.ao_factory();
  auto v_time0 = mpqc::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> V_ijij_ijji;

  utility::print_par(world, "\nCompute V_ijij_ijji With Dual Basis DF \n");
  {
    auto left = lcao_factory(L"(Κ |GR|i2 i1)");
    auto middle = ao_factory(L"(Κ|GR|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |GR|j1 j2)");

    auto time0 = mpqc::now(world, accurate_time);

    V_ijij_ijji("i1,j1,i2,j2") = (left * middle * right).set_shape(shape);
    V_ijij_ijji.truncate();

    //    std::cout << V_ijij_ijji << std::endl;
    // all types of GR integral not needed
    lcao_factory.purge_operator(world, L"GR");

    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|G|m n>[df]");
    auto right = lcao_factory(L"<i2 j2|R|m n>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_ijij_ijji("i1,j1,i2,j2") -= (left * right).set_shape(shape);
    V_ijij_ijji.truncate();
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<i j|R|m n>[df]");

  {
    auto left = lcao_factory(L"<i1 j1|G|m a>[df]");
    auto right = lcao_factory(L"<i2 j2|R|m a>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    TA::DistArray<Tile, TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left * right).set_shape(shape);
    tmp.truncate();
    V_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    V_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<i j|R|m a>[df]");

  {
    auto left = lcao_factory(L"<i1 j1|G|a b>[df]");
    auto right = lcao_factory(L"<i2 j2|R|a b>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_ijij_ijji("i1,j1,i2,j2") -= (left * right).set_shape(shape);
    V_ijij_ijji.truncate();
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term4 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<i j|R|a b>[df]");
  //  std::cout << V_ijij_ijji << std::endl;

  {
    auto left = lcao_factory(L"<i1 j1|G|m a'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|m a'>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    TA::DistArray<Tile, TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left * right).set_shape(shape);
    tmp.truncate();
    V_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    V_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term5 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<i j|R|m a'>[df]");
  auto v_time1 = mpqc::now(world, accurate_time);
  auto v_time = mpqc::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");

  //  std::cout << V_ijij_ijji << std::endl;
  return V_ijij_ijji;
};

template <typename Tile, typename Policy>
std::tuple<TA::DistArray<Tile, Policy>, TA::DistArray<Tile, Policy>>
compute_V_ixjy_ixyj_db_df(lcao::LCAOFactory<Tile, Policy> &lcao_factory,
                          bool cabs = true) {
  auto &world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  auto v_time0 = mpqc::now(world, accurate_time);

  TA::DistArray<Tile, Policy> V_ixjy, V_ixyj;

  {
    auto left = lcao_factory(L"<i x|GR|j y>[df]");
    auto right = lcao_factory(L"<i x|GR|y j>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_ixjy("i,x,j,y") = left;
    V_ixyj("i,x,y,j") = right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i x|G|m n>[df]");
    auto right1 = lcao_factory(L"<j y|R|m n>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_ixjy("i,x,j,y") -= left * right1;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);

    auto right2 = lcao_factory(L"<y j|R|m n>[df]");

    time0 = mpqc::now(world, accurate_time);
    V_ixyj("i,x,y,j") -= left * right2;
    time1 = mpqc::now(world, accurate_time);
    time += mpqc::duration_in_s(time0, time1);

    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<j y|R|m n>[df]");
  lcao_factory.purge_formula(world, L"<y j|R|m n>[df]");

  {
    auto left = lcao_factory(L"<i x|G|m a>[df]");
    auto right1 = lcao_factory(L"<j y|R|m a>[df]");
    auto right2 = lcao_factory(L"<y j|R|m a>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_ixjy("i,x,j,y") -= left * right1;
    V_ixyj("i,x,y,j") -= left * right2;

    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<j y|R|m a>[df]");
  lcao_factory.purge_formula(world, L"<y j|R|m a>[df]");

  {
    auto left = lcao_factory(L"<i x|G|a m>[df]");
    auto right1 = lcao_factory(L"<j y|R|a m>[df]");
    auto right2 = lcao_factory(L"<y j|R|a m>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_ixjy("i,x,j,y") -= left * right1;
    V_ixyj("i,x,y,j") -= left * right2;

    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term4 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<j y|R|a m>[df]");
  lcao_factory.purge_formula(world, L"<y j|R|a m>[df]");

  {
    auto left = lcao_factory(L"<i x|G|a b>[df]");
    auto right1 = lcao_factory(L"<j y|R|a b>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_ixjy("i,x,j,y") -= left * right1;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);

    auto right2 = lcao_factory(L"<y j|R|a b>[df]");
    time0 = mpqc::now(world, accurate_time);
    V_ixyj("i,x,y,j") -= left * right2;
    time1 = mpqc::now(world, accurate_time);
    time += mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term5 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<j y|R|a b>[df]");

  if (cabs) {
    {
      auto left = lcao_factory(L"<i x|G|m a'>[df]");
      auto right1 = lcao_factory(L"<j y|R|m a'>[df]");
      auto right2 = lcao_factory(L"<y j|R|m a'>[df]");

      auto time0 = mpqc::now(world, accurate_time);
      TA::DistArray<Tile, Policy> tmp;
      V_ixjy("i,x,j,y") -= left * right1;
      V_ixyj("i,x,y,j") -= left * right2;

      auto time1 = mpqc::now(world, accurate_time);
      auto time = mpqc::duration_in_s(time0, time1);
      utility::print_par(world, "V Term6 Time: ", time, " S\n");
    }
    lcao_factory.purge_formula(world, L"<j y|R|m a'>[df]");
    lcao_factory.purge_formula(world, L"<y j|R|m a'>[df]");

    {
      auto left = lcao_factory(L"<i x|G|a' m>[df]");
      auto right1 = lcao_factory(L"<j y|R|a' m>[df]");
      auto right2 = lcao_factory(L"<y j|R|a' m>[df]");

      auto time0 = mpqc::now(world, accurate_time);
      TA::DistArray<Tile, Policy> tmp;
      V_ixjy("i,x,j,y") -= left * right1;
      V_ixyj("i,x,y,j") -= left * right2;

      auto time1 = mpqc::now(world, accurate_time);
      auto time = mpqc::duration_in_s(time0, time1);
      utility::print_par(world, "V Term7 Time: ", time, " S\n");
    }
    lcao_factory.purge_formula(world, L"<j y|R|a' m>[df]");
    lcao_factory.purge_formula(world, L"<y j|R|a' m>[df]");
  }

  auto v_time1 = mpqc::now(world, accurate_time);
  auto v_time = mpqc::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");

  return std::make_tuple(V_ixjy, V_ixyj);
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_V_xyab_db_df(
    lcao::LCAOFactory<Tile, Policy> &lcao_factory) {
  auto &world = lcao_factory.world();
  auto &ao_factory = lcao_factory.ao_factory();
  bool accurate_time = lcao_factory.accurate_time();

  auto v_time0 = mpqc::now(world, accurate_time);

  TA::DistArray<Tile, Policy> V_xyab;
  TA::DistArray<Tile, Policy> tmp;

  utility::print_par(world, "\nCompute V_xyab With Dual Basis DF \n");

  {
    auto left = lcao_factory(L"(Κ |GR|i a)");
    auto middle = ao_factory(L"(Κ|GR|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |GR|j b)");

    auto time0 = mpqc::now(world, accurate_time);
    V_xyab("i,j,a,b") = left * middle * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<a b|G|c d>[df]");
    auto right = lcao_factory(L"<i j|R|c d>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_xyab("i,j,a,b") -= left * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<m1 m2|G|a b>[df]");
    auto right = lcao_factory(L"<i j|R|m1 m2>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_xyab("i,j,a,b") -= left * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<a b|G|m c>[df]");
    auto right = lcao_factory(L"<i j|R|m c>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    tmp("i,j,a,b") = left * right;
    V_xyab("i,j,a,b") -= tmp("i,j,a,b");
    V_xyab("i,j,a,b") -= tmp("j,i,b,a");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term4 Time: ", time, " S\n");
  }

  {
    auto right = lcao_factory(L"<a b|G|m a'>[df]");
    auto left = lcao_factory(L"<i j|R|m a'>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    tmp("i,j,a,b") = left * right;
    V_xyab("i,j,a,b") -= tmp("i,j,a,b");
    V_xyab("i,j,a,b") -= tmp("j,i,b,a");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term5 Time: ", time, " S\n");
  }

  auto v_time1 = mpqc::now(world, accurate_time);
  auto v_time = mpqc::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");

  return V_xyab;
}

template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_V_iaxy_db_df(
    lcao::LCAOFactory<Tile, Policy> &lcao_factory) {
  auto &world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  TA::DistArray<Tile, Policy> V_iaxy;

  auto v_time0 = mpqc::now(world, accurate_time);

  utility::print_par(world, "\nCompute V_iaxy With Dual Basis DF \n");
  {
    auto left = lcao_factory(L"(Κ |GR|i k)");
    auto middle = lcao_factory(L"(Κ|GR|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |GR|a l)");

    auto time0 = mpqc::now(world, accurate_time);
    V_iaxy("i,a,k,l") = left * middle * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|m1 m2>[df]");
    auto right = lcao_factory(L"<k l|R|m1 m2>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|b c>[df]");
    auto right = lcao_factory(L"<k l|R|b c>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|m b>[df]");
    auto right = lcao_factory(L"<k l|R|m b>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term4 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|b m>[df]");
    auto right = lcao_factory(L"<k l|R|b m>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term5 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|m a'>[df]");
    auto right = lcao_factory(L"<k l|R|m a'>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term6 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|a' m>[df]");
    auto right = lcao_factory(L"<k l|R|a' m>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term7 Time: ", time, " S\n");
  }

  auto v_time1 = mpqc::now(world, accurate_time);
  auto v_time = mpqc::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");
  return V_iaxy;
};

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_X_ijij_ijji_db_df(
    lcao::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &ijij_ijji_shape) {
  bool accurate_time = lcao_factory.accurate_time();
  auto &world = lcao_factory.world();
  auto &ao_factory = lcao_factory.ao_factory();
  auto x_time0 = mpqc::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> X_ijij_ijji;

  utility::print_par(world, "\nCompute X_ijij_ijji With Dual Basis DF \n");
  {
    auto left = lcao_factory(L"(Κ |R2|i1 i2)");
    auto middle = ao_factory(L"(Κ|R2|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |R2|j1 j2)");

    auto time0 = mpqc::now(world, accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") =
        (left * middle * right).set_shape(ijij_ijji_shape);
    X_ijij_ijji.truncate();
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "X Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|m n>[df]");
    auto right = lcao_factory(L"<i2 j2|R|m n>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") -= (left * right).set_shape(ijij_ijji_shape);
    X_ijij_ijji.truncate();
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "X Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|m a>[df]");
    auto right = lcao_factory(L"<i2 j2|R|m a>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    TA::DistArray<Tile, TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left * right).set_shape(ijij_ijji_shape);
    tmp.truncate();
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|a b>[df]");
    auto right = lcao_factory(L"<i2 j2|R|a b>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") -= (left * right).set_shape(ijij_ijji_shape);
    X_ijij_ijji.truncate();
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "V Term4 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|m a'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|m a'>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    TA::DistArray<Tile, TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left * right).set_shape(ijij_ijji_shape);
    tmp.truncate();
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "X Term5 Time: ", time, " S\n");
  }

  auto x_time1 = mpqc::now(world, accurate_time);
  auto x_time = mpqc::duration_in_s(x_time0, x_time1);
  utility::print_par(world, "X Term Total Time: ", x_time, " S\n");

  //    std::cout << X_ijij_ijji << std::endl;
  return X_ijij_ijji;
};

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_B_ijij_ijji_db_df(
    lcao::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &ijij_ijji_shape) {
  bool accurate_time = lcao_factory.accurate_time();
  auto &world = lcao_factory.world();
  auto &ao_factory = lcao_factory.ao_factory();
  auto b_time0 = mpqc::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> B_ijij_ijji;
  TA::DistArray<Tile, TA::SparsePolicy> tmp;

  utility::print_par(world, "\nCompute B_ijij_ijji C With Dual Basis DF \n");

  {
    auto left = lcao_factory(L"(Κ |dR2|i1 i2)");
    auto middle = ao_factory(L"(Κ|dR2|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |dR2|j1 j2)");

    auto time0 = mpqc::now(world, accurate_time);
    B_ijij_ijji("i1,j1,i2,j2") =
        (left * middle * right).set_shape(ijij_ijji_shape);
    B_ijij_ijji.truncate();
    lcao_factory.purge_operator(world, L"dR2");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "B Term1 Time: ", time, " S\n");
  }

  lcao_factory.purge_operator(world, L"dR2");

  {
    auto left1 = lcao_factory(L"<i1 j1|R2|m j2>[df]");
    auto hJ1 = lcao_factory(L"<m | hJ | i2>[df]");

    auto left2 = lcao_factory(L"<i1 j1|R2|A' j2>[df]");
    auto hJ2 = lcao_factory(L"<A' | hJ | i2>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left1 * hJ1).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") += (left2 * hJ2).set_shape(ijij_ijji_shape);
    tmp.truncate();
    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");

    lcao_factory.purge_operator(world, L"R2");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "B Term2 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<i1 j1|R2|m j2>[df]");
  lcao_factory.purge_formula(world, L"<i1 j1|R2|A' j2>[df]");

  {
    auto left1 = lcao_factory(L"<i1 j1|R|m1 A'>[df]");
    auto middle1 = lcao_factory(L"<m1 |K| m2>[df]");
    auto right1 = lcao_factory(L"<i2 j2|R|m2 A'>[df]");

    auto left2 = lcao_factory(L"<i1 j1|R|m1 A'>[df]");
    auto middle2 = lcao_factory(L"<m1 |K| B'>[df]");
    auto right2 = lcao_factory(L"<i2 j2|R|B' A'>[df]");

    auto left3 = lcao_factory(L"<i1 j1|R|A' B'>[df]");
    auto middle3 = lcao_factory(L"<A' |K| C'>[df]");
    auto right3 = lcao_factory(L"<i2 j2|R|C' B'>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left1 * middle1 * right1).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") +=
        (2.0 * left2 * middle2 * right2).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") += (left3 * middle3 * right3).set_shape(ijij_ijji_shape);
    tmp.truncate();
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    // AO R integral not needed
    lcao_factory.ao_factory().registry().purge_operator(world, L"R");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "B Term3 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<i1 j1|R|A' B'>[df]");

  {
    auto left1 = lcao_factory(L"<i1 j1|R|m1 m2>[df]");
    auto middle1 = lcao_factory(L"<m1|hJ|m3>[df]");
    auto right1 = lcao_factory(L"<i2 j2|R|m3 m2>[df]");

    auto left2 = lcao_factory(L"<i1 j1|R|m1 A'>[df]");
    auto middle2 = lcao_factory(L"<A'|hJ|m2>[df]");
    auto right2 = lcao_factory(L"<i2 j2|R|m1 m2>[df]");

    auto left3 = lcao_factory(L"<i1 j1|R|m A'>[df]");
    auto middle3 = lcao_factory(L"<A'|hJ|B'>[df]");
    auto right3 = lcao_factory(L"<i2 j2|R|m B'>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left1 * middle1 * right1).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") +=
        (2.0 * left2 * middle2 * right2).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") += (left3 * middle3 * right3).set_shape(ijij_ijji_shape);
    tmp.truncate();
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    lcao_factory.purge_operator(world, L"hJ");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "B Term4 Time: ", time, " S\n");
  }

  lcao_factory.purge_formula(world, L"<i1 j1|R|m1 A'>[df]");

  {
    auto left1 = lcao_factory(L"<i1 j1|R|m b'>[df]");
    auto middle1 = lcao_factory(L"<m|F|n>[df]");
    auto right1 = lcao_factory(L"<i2 j2|R|n b'>[df]");

    auto left2 = lcao_factory(L"<i1 j1|R|m b'>[df]");
    auto middle2 = lcao_factory(L"<m|F|A'>[df]");
    auto right2 = lcao_factory(L"<i2 j2|R|A' b'>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left1 * middle1 * right1).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") +=
        (2.0 * left2 * middle2 * right2).set_shape(ijij_ijji_shape);
    tmp.truncate();
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    // P' doesn't appear later
    lcao_factory.registry().purge_index(world, L"A'");

    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "B Term5 Time: ", time, " S\n");
  }

  {
    auto left1 = lcao_factory(L"<i1 j1|R|m a>[df]");
    auto middle1 = lcao_factory(L"<m|F|n>[df]");
    auto right1 = lcao_factory(L"<i2 j2|R|n a>[df]");

    auto left2 = lcao_factory(L"<i1 j1|R|m a>[df]");
    auto middle2 = lcao_factory(L"<m|F|b>[df]");
    auto right2 = lcao_factory(L"<i2 j2|R|b a>[df]");

    auto left3 = lcao_factory(L"<i1 j1|R|b a>[df]");
    auto middle3 = lcao_factory(L"<b|F|c>[df]");
    auto right3 = lcao_factory(L"<i2 j2|R|c a>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left1 * middle1 * right1).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") +=
        (2 * left2 * middle2 * right2).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") += (left3 * middle3 * right3).set_shape(ijij_ijji_shape);
    tmp.truncate();
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "B Term6 Time: ", time, " S\n");
  }

  {
    auto left1 = lcao_factory(L"<i1 j1|R|m a>[df]");
    auto middle1 = lcao_factory(L"<m|F|a'>[df]");
    auto right1 = lcao_factory(L"<i2 j2|R|a' a>[df]");

    auto left2 = lcao_factory(L"<i1 j1|R|b a>[df]");
    auto middle2 = lcao_factory(L"<b|F|a'>[df]");
    auto right2 = lcao_factory(L"<i2 j2|R|a' a>[df]");

    auto time0 = mpqc::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left1 * middle1 * right1).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") += (left2 * middle2 * right2).set_shape(ijij_ijji_shape);
    tmp.truncate();
    tmp("i1,j1,i2,j2") = 2.0 * tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc::now(world, accurate_time);
    auto time = mpqc::duration_in_s(time0, time1);
    utility::print_par(world, "B Term7 Time: ", time, " S\n");
  }

  auto b_time1 = mpqc::now(world, accurate_time);
  auto b_time = mpqc::duration_in_s(b_time0, b_time1);
  utility::print_par(world, "B Term Total Time: ", b_time, " S\n");

  //    std::cout << B_ijij_ijji << std::endl;

  return B_ijij_ijji;
}

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_VT2_ijij_ijji_db_df(
    lcao::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    const TA::DistArray<Tile, TA::SparsePolicy> &t2,
    const TA::SparseShape<float> &ijij_ijji_shape) {
  auto &world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();

  TA::DistArray<Tile, TA::SparsePolicy> V_ijij_ijji;

  // compute C_ijab
  TA::DistArray<Tile, TA::SparsePolicy> C_ijab =
      compute_C_ijab_df(lcao_factory);

  // compute V_ijab
  TA::DistArray<Tile, TA::SparsePolicy> V_ijab =
      compute_V_xyab_db_df(lcao_factory);

  auto vt2_time0 = mpqc::now(world, accurate_time);
  utility::print_par(world, "\nCompute VT2_ijij_ijji With Dual Basis DF\n");
  V_ijij_ijji("i1,j1,i2,j2") =
      ((V_ijab("i2,j2,a,b") + C_ijab("i2,j2,a,b")) * t2("a,b,i1,j1"))
          .set_shape(ijij_ijji_shape);
  V_ijij_ijji.truncate();

  auto vt2_time1 = mpqc::now(world, accurate_time);
  auto vt2_time = mpqc::duration_in_s(vt2_time0, vt2_time1);
  utility::print_par(world, "VT2 Term Total Time: ", vt2_time, " S\n");

  return V_ijij_ijji;
};

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_VT1_ijij_ijji_db_df(
    lcao::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    const TA::DistArray<Tile, TA::SparsePolicy> &t1,
    const TA::SparseShape<float> &ijij_ijji_shape) {
  auto &world = lcao_factory.world();
  bool accurate_time = lcao_factory.accurate_time();
  TA::DistArray<Tile, TA::SparsePolicy> V_ijij_ijji;
  TA::DistArray<Tile, TA::SparsePolicy> V_iaij =
      compute_V_iaxy_db_df(lcao_factory);

  auto vt1_time0 = mpqc::now(world, accurate_time);
  utility::print_par(world, "\nCompute VT1_ijij_ijji With Dual Basis DF\n");

  V_ijij_ijji("i1,j1,i2,j2") =
      (V_iaij("i1,a,i2,j2") * t1("a,j1")).set_shape(ijij_ijji_shape);
  V_ijij_ijji.truncate();
  V_ijij_ijji("i1,j1,i2,j2") += V_ijij_ijji("j1,i1,j2,i2");
  auto vt1_time1 = mpqc::now(world, accurate_time);
  auto vt1_time = mpqc::duration_in_s(vt1_time0, vt1_time1);
  utility::print_par(world, "VT1 Term Total Time: ", vt1_time, " S\n");

  return V_ijij_ijji;
};

}  // namespace f12
}  // namespace lcao
}  // namespace mpqc

#endif  // MPQC4_SRC_MPQC_CHEMISTRY_QC_F12_DB_F12_INTERMEDIATES_H_
