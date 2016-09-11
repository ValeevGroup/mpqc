//
// Created by Chong Peng on 04/11/16.
//

#ifndef MPQC_F12_INTERMEDIATES_H
#define MPQC_F12_INTERMEDIATES_H

#include "../../../../../common/namespaces.h"
#include "../../../../../include/tiledarray.h"
#include <mpqc/chemistry/qc/integrals/lcao_factory.h>
#include <mpqc/util/misc/string.h>

namespace mpqc {
namespace f12 {

/**
 * MP2-F12 C approach V term with Density Fitting, only ijij ijji part is
 * computed
 * \f$V_{ij}^{ij}\f$  \f$V_{ij}^{ji}\f$
 * @param lcao_factory reference to LCAOFactory, has to use SparsePolicy
 * @param shape SparseShape that has ijij ijji shape
 * @return V(i1,j1,i2,j2)
 */
template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_V_ijij_ijji_df(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &shape) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();
  auto &ao_integral = lcao_factory.atomic_integral();
  auto v_time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> V_ijij_ijji;

  utility::print_par(world, "\nCompute V_ijij_ijji With DF \n");
  {
    auto left = lcao_factory(L"(Κ |GR|i2 i1)");
    auto middle = ao_integral(L"(Κ|GR|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |GR|j1 j2)");

    auto time0 = mpqc_time::now(world, accurate_time);

    V_ijij_ijji("i1,j1,i2,j2") = (left * middle * right).set_shape(shape);
    // all types of GR integral not needed
    lcao_factory.purge_operator(world, L"GR");

    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|G|p q>[df]");
    auto right = lcao_factory(L"<i2 j2|R|p q>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_ijij_ijji("i1,j1,i2,j2") -= (left * right).set_shape(shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|G|m a'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|m a'>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    TA::DistArray<Tile, TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left * right).set_shape(shape);
    V_ijij_ijji("i1,j1,i2,j2") -= (tmp("i1,j1,i2,j2")).set_shape(shape);
    V_ijij_ijji("i1,j1,i2,j2") -= (tmp("j1,i1,j2,i2")).set_shape(shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  auto v_time1 = mpqc_time::now(world, accurate_time);
  auto v_time = mpqc_time::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");

  return V_ijij_ijji;
};

/**
 * MP2-F12 C approach V term, only ijij ijji part is computed
 * \f$V_{ij}^{ij}\f$  \f$V_{ij}^{ji}\f$
 * @param lcao_factory reference to LCAOFactory, has to use SparsePolicy
 * @param shape SparseShape that has ijij ijji shape
 * @return V(i1,j1,i2,j2)
 */
template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_V_ijij_ijji(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &shape) {
  bool accurate_time = lcao_factory.accurate_time();
  auto &world = lcao_factory.get_world();
  auto v_time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> V_ijij_ijji;

  utility::print_par(world, "\nCompute V_ijij_ijji Without DF \n");
  {
    auto left = lcao_factory(L"<i1 j1|GR|i2 j2>");

    auto time0 = mpqc_time::now(world, accurate_time);

    V_ijij_ijji("i1,j1,i2,j2") = (left).set_shape(shape);
    // all types of GR integral not needed
    lcao_factory.purge_operator(world, L"GR");

    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|G|p q>");
    auto right = lcao_factory(L"<i2 j2|R|p q>");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_ijij_ijji("i1,j1,i2,j2") -= (left * right).set_shape(shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|G|m a'>");
    auto right = lcao_factory(L"<i2 j2|R|m a'>");

    auto time0 = mpqc_time::now(world, accurate_time);
    TA::DistArray<Tile, TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left * right).set_shape(shape);
    //    V_ijij_ijji("i1,j1,i2,j2") -= (lcao_factory(L"(j1 m|G|i1
    //    a')[df]")*lcao_factory(L"(j2 m|R|i2 a')[df]")).set_shape(shape);
    V_ijij_ijji("i1,j1,i2,j2") -= (tmp("i1,j1,i2,j2")).set_shape(shape);
    V_ijij_ijji("i1,j1,i2,j2") -= (tmp("j1,i1,j2,i2")).set_shape(shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  auto v_time1 = mpqc_time::now(world, accurate_time);
  auto v_time = mpqc_time::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");
  return V_ijij_ijji;
};

/**
 * MP2-F12 C approach X term with DF, only ijij ijji part is computed
 * \f$X_{ij}^{ij}\f$  \f$X_{ij}^{ji}\f$
 * @param lcao_factory reference to LCAOFactory, has to use SparsePolicy
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return X(i1,j1,i2,j2)
 */
template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_X_ijij_ijji_df(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &ijij_ijji_shape) {
  bool accurate_time = lcao_factory.accurate_time();
  auto &world = lcao_factory.get_world();
  auto &ao_integral = lcao_factory.atomic_integral();
  auto x_time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> X_ijij_ijji;

  utility::print_par(world, "\nCompute X_ijij_ijji With DF \n");
  {
    auto left = lcao_factory(L"(Κ |R2|i1 i2)");
    auto middle = ao_integral(L"(Κ|R2|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |R2|j1 j2)");

    auto time0 = mpqc_time::now(world, accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") =
        (left * middle * right).set_shape(ijij_ijji_shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "X Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|p q>[df]");
    auto right = lcao_factory(L"<i2 j2|R|p q>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") -= (left * right).set_shape(ijij_ijji_shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "X Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|m a'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|m a'>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    TA::DistArray<Tile, TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left * right).set_shape(ijij_ijji_shape);
    //    X_ijij_ijji("i1,j1,i2,j2") -= (lcao_factory(L"(j1 m|R|i1
    //    a')[df]")*lcao_factory(L"(j2 m|R|i2
    //    a')[df]")).set_shape(ijij_ijji_shape);
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "X Term3 Time: ", time, " S\n");
  }

  auto x_time1 = mpqc_time::now(world, accurate_time);
  auto x_time = mpqc_time::duration_in_s(x_time0, x_time1);
  utility::print_par(world, "X Term Total Time: ", x_time, " S\n");

  //    std::cout << X_ijij_ijji << std::endl;
  return X_ijij_ijji;
};

/**
 * MP2-F12 C approach X term without DF, only ijij ijji part is computed
 * \f$X_{ij}^{ij}\f$  \f$X_{ij}^{ji}\f$
 * @param lcao_factory reference to LCAOFactory, has to use SparsePolicy
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return X(i1,j1,i2,j2)
 */
template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_X_ijij_ijji(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &ijij_ijji_shape) {
  bool accurate_time = lcao_factory.accurate_time();
  auto &world = lcao_factory.get_world();
  auto x_time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> X_ijij_ijji;

  utility::print_par(world, "\nCompute X_ijij_ijji Without DF \n");
  {
    auto left = lcao_factory(L"<i1 j1 |R2|i2 j2>");

    auto time0 = mpqc_time::now(world, accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") = (left).set_shape(ijij_ijji_shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "X Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|p q>");
    auto right = lcao_factory(L"<i2 j2|R|p q>");

    auto time0 = mpqc_time::now(world, accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") -= (left * right).set_shape(ijij_ijji_shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "X Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|m a'>");
    auto right = lcao_factory(L"<i2 j2|R|m a'>");

    auto time0 = mpqc_time::now(world, accurate_time);
    TA::DistArray<Tile, TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left * right).set_shape(ijij_ijji_shape);
    //    X_ijij_ijji("i1,j1,i2,j2") -= (lcao_factory(L"(j1 m|R|i1
    //    a')[df]")*lcao_factory(L"(j2 m|R|i2
    //    a')[df]")).set_shape(ijij_ijji_shape);
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    X_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "X Term3 Time: ", time, " S\n");
  }

  auto x_time1 = mpqc_time::now(world, accurate_time);
  auto x_time = mpqc_time::duration_in_s(x_time0, x_time1);
  utility::print_par(world, "X Term Total Time: ", x_time, " S\n");
  return X_ijij_ijji;
};

/**
 * MP2-F12 C approach B term, only ijij ijji part is computed
 * \f$B_{ij}^{ij}\f$  \f$B_{ij}^{ji}\f$
 * @param lcao_factory reference to LCAOFactory, has to use SparsePolicy
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return B(i1,j1,i2,j2)
 */
template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_B_ijij_ijji_df(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &ijij_ijji_shape) {
  bool accurate_time = lcao_factory.accurate_time();
  auto &world = lcao_factory.get_world();
  auto &ao_integral = lcao_factory.atomic_integral();
  auto b_time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> B_ijij_ijji;
  TA::DistArray<Tile, TA::SparsePolicy> tmp;

  utility::print_par(world, "\nCompute B_ijij_ijji With DF \n");

  {
    auto left = lcao_factory(L"(Κ |dR2|i1 i2)");
    auto middle = ao_integral(L"(Κ|dR2|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |dR2|j1 j2)");

    auto time0 = mpqc_time::now(world, accurate_time);
    B_ijij_ijji("i1,j1,i2,j2") =
        (left * middle * right).set_shape(ijij_ijji_shape);
    lcao_factory.purge_operator(world, L"dR2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term1 Time: ", time, " S\n");
  }

  {
    auto hJ = lcao_factory(L"<P' | hJ | i2>[df]");
    auto left = lcao_factory(L"<i1 j1|R2|P' j2>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * hJ).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") += (lcao_factory(L"(j1 P'|R2|i1
    //    i2)[df]")*hJ("P',j2")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");

    lcao_factory.purge_operator(world, L"R2");
    lcao_factory.purge_operator(world, L"hJ");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|P' Q'>[df]");
    auto middle = lcao_factory(L"<P'|K|R'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|R' Q'>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (lcao_factory(L"(j1 P'|R|i1
    //    Q')[df]")*lcao_factory(L"(P'|K|R')[df]")*lcao_factory(L"(j2 R'|R|i2
    //    Q')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    // AO R integral not needed
    lcao_factory.atomic_integral().registry().purge_operator(world, L"R");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term3 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|P' m>[df]");
    auto middle = lcao_factory(L"<P'|F|R'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|R' m>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (lcao_factory(L"(j1 P'|R|i1
    //    m)[df]")*lcao_factory(L"(P'|F|R')[df]")*lcao_factory(L"(j2 R'|R|i2
    //    m)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term4 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|m b'>[df]");
    auto middle = lcao_factory(L"<m|F|P'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|P' b'>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") =
        (2.0 * left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*lcao_factory(L"(j1 m|R|i1
    //    b')[df]")*lcao_factory(L"(m|F|P')[df]")*lcao_factory(L"(j2 P'|R|i2
    //    b')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    // P' doesn't appear later
    lcao_factory.registry().purge_index(world, L"P'");

    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term5 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|p a>[df]");
    auto middle = lcao_factory(L"<p|F|r>[df]");
    auto right = lcao_factory(L"<i2 j2|R|r a>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (lcao_factory(L"(j1 p|R|i1
    //    a)[df]")*lcao_factory(L"(p|F|r)[df]")*lcao_factory(L"(j2 r|R|i2
    //    a)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term6 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|m b'>[df]");
    auto middle = lcao_factory(L"<m|F|n>[df]");
    auto right = lcao_factory(L"<i2 j2|R|n b'>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") += (lcao_factory(L"(j1 m|R|i1
    //    b')[df]")*lcao_factory(L"(m|F|n)[df]")*lcao_factory(L"(j2 n|R|i2
    //    b')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term7 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|p a>[df]");
    auto middle = lcao_factory(L"<p|F|a'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|a' a>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") =
        (2.0 * left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*lcao_factory(L"(j1 p|R|i1
    //    a)[df]")*lcao_factory(L"(p|F|a')[df]")*lcao_factory(L"(i2 a|R|j2
    //    a')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term8 Time: ", time, " S\n");
  }

  auto b_time1 = mpqc_time::now(world, accurate_time);
  auto b_time = mpqc_time::duration_in_s(b_time0, b_time1);
  utility::print_par(world, "B Term Total Time: ", b_time, " S\n");

  //    std::cout << B_ijij_ijji << std::endl;

  return B_ijij_ijji;
};



/**
 * MP2-F12 D approach B term, only ijij ijji part is computed
 * \f$B_{ij}^{ij}\f$  \f$B_{ij}^{ji}\f$
 * @param lcao_factory reference to LCAOFactory, has to use SparsePolicy
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return B(i1,j1,i2,j2)
 */
template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_B_ijij_ijji_D_df(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &ijij_ijji_shape) {
  bool accurate_time = lcao_factory.accurate_time();
  auto &world = lcao_factory.get_world();
  auto &ao_integral = lcao_factory.atomic_integral();
  auto b_time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> B_ijij_ijji;
  TA::DistArray<Tile, TA::SparsePolicy> tmp;

  utility::print_par(world, "\nCompute B_ijij_ijji With DF \n");

  {
    auto left = lcao_factory(L"(Κ |dR2|i1 i2)");
    auto middle = ao_integral(L"(Κ|dR2|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |dR2|j1 j2)");

    auto time0 = mpqc_time::now(world, accurate_time);
    B_ijij_ijji("i1,j1,i2,j2") =
        (left * middle * right).set_shape(ijij_ijji_shape);

    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term1 Time: ", time, " S\n");
  }

  // operator dR2 no longer needed
  lcao_factory.purge_operator(world, L"dR2");

  {
    auto hJ = lcao_factory(L"<P' | hJ | i2>[df]");
    auto left = lcao_factory(L"<i1 j1|R2|P' j2>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * hJ).set_shape(ijij_ijji_shape);

    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");

    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term2 Time: ", time, " S\n");
  }
  lcao_factory.purge_operator(world, L"R2");
  lcao_factory.purge_operator(world, L"hJ");

  {
    auto left = lcao_factory(L"<i1 j1|R|P' q>[df]");
    auto middle = lcao_factory(L"<q|K|r>[df]");
    auto right = lcao_factory(L"<i2 j2|R|P' r>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);

    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term3 Time: ", time, " S\n");
  }
  // AO R integral not needed
  lcao_factory.atomic_integral().registry().purge_operator(world, L"R");

  {
    auto left = lcao_factory(L"<i1 j1|R|P' m>[df]");
    auto middle = lcao_factory(L"<P'|hJ|R'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|R' m>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);

    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term4 Time: ", time, " S\n");
  }

  // P' doesn't appear later
  lcao_factory.registry().purge_index(world, L"P'");

  {
    auto left = lcao_factory(L"<i1 j1|R|m p>[df]");
    auto middle = lcao_factory(L"<p|K|q>[df]");
    auto right = lcao_factory(L"<i2 j2|R|m q>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") =
        (left * middle * right).set_shape(ijij_ijji_shape);

    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");

    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term5 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|p a>[df]");
    auto middle = lcao_factory(L"<p|F|r>[df]");
    auto right = lcao_factory(L"<i2 j2|R|r a>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);

    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term6 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|b' m>[df]");
    auto middle = lcao_factory(L"<m|F|n>[df]");
    auto right = lcao_factory(L"<i2 j2|R|b' n>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);

    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term7 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|p a>[df]");
    auto middle = lcao_factory(L"<p|F|a'>[df]");
    auto right = lcao_factory(L"<i2 j2|R|a' a>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") =
        (2.0 * left * middle * right).set_shape(ijij_ijji_shape);

    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term8 Time: ", time, " S\n");
  }

  auto b_time1 = mpqc_time::now(world, accurate_time);
  auto b_time = mpqc_time::duration_in_s(b_time0, b_time1);
  utility::print_par(world, "B Term Total Time: ", b_time, " S\n");

  //    std::cout << B_ijij_ijji << std::endl;

  return B_ijij_ijji;
};

/**
 * MP2-F12 C approach B term without DF, only ijij ijji part is computed
 * \f$B_{ij}^{ij}\f$  \f$B_{ij}^{ji}\f$
 * @param lcao_factory reference to LCAOFactory, has to use SparsePolicy
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return B(i1,j1,i2,j2)
 */
template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_B_ijij_ijji(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    TA::SparseShape<float> &ijij_ijji_shape) {
  bool accurate_time = lcao_factory.accurate_time();
  auto &world = lcao_factory.get_world();
  auto b_time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, TA::SparsePolicy> B_ijij_ijji;
  TA::DistArray<Tile, TA::SparsePolicy> tmp;

  utility::print_par(world, "\nCompute B_ijij_ijji Without DF \n");

  {
    auto left = lcao_factory(L"<i1 j1 |dR2|i2 j2>");

    auto time0 = mpqc_time::now(world, accurate_time);
    B_ijij_ijji("i1,j1,i2,j2") = (left).set_shape(ijij_ijji_shape);
    lcao_factory.purge_operator(world, L"dR2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term1 Time: ", time, " S\n");
  }

  {
    auto hJ = lcao_factory(L"<P' | hJ | i2>");
    auto left = lcao_factory(L"<i1 j1|R2|P' j2>");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * hJ).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") += (lcao_factory(L"(j1 P'|R2|i1
    //    i2)[df]")*hJ("P',j2")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");

    lcao_factory.purge_operator(world, L"R2");
    lcao_factory.purge_operator(world, L"hJ");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|P' Q'>");
    auto middle = lcao_factory(L"<P'|K|R'>");
    auto right = lcao_factory(L"<i2 j2|R|R' Q'>");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (lcao_factory(L"(j1 P'|R|i1
    //    Q')[df]")*lcao_factory(L"(P'|K|R')[df]")*lcao_factory(L"(j2 R'|R|i2
    //    Q')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    // AO R integral not needed
    lcao_factory.atomic_integral().registry().purge_operator(world, L"R");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term3 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|P' m>");
    auto middle = lcao_factory(L"<P'|F|R'>");
    auto right = lcao_factory(L"<i2 j2|R|R' m>");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (lcao_factory(L"(j1 P'|R|i1
    //    m)[df]")*lcao_factory(L"(P'|F|R')[df]")*lcao_factory(L"(j2 R'|R|i2
    //    m)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term4 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|m b'>");
    auto middle = lcao_factory(L"<m|F|P'>");
    auto right = lcao_factory(L"<i2 j2|R|P' b'>");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") =
        (2.0 * left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*lcao_factory(L"(j1 m|R|i1
    //    b')[df]")*lcao_factory(L"(m|F|P')[df]")*lcao_factory(L"(j2 P'|R|i2
    //    b')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

    // P' doesn't appear later
    lcao_factory.registry().purge_index(world, L"P'");

    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term5 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|p a>");
    auto middle = lcao_factory(L"<p|F|r>");
    auto right = lcao_factory(L"<i2 j2|R|r a>");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (lcao_factory(L"(j1 p|R|i1
    //    a)[df]")*lcao_factory(L"(p|F|r)[df]")*lcao_factory(L"(j2 r|R|i2
    //    a)[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term6 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|m b'>");
    auto middle = lcao_factory(L"<m|F|n>");
    auto right = lcao_factory(L"<i2 j2|R|n b'>");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") = (left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") += (lcao_factory(L"(j1 m|R|i1
    //    b')[df]")*lcao_factory(L"(m|F|n)[df]")*lcao_factory(L"(j2 n|R|i2
    //    b')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") += tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term7 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i1 j1|R|p a>");
    auto middle = lcao_factory(L"<p|F|a'>");
    auto right = lcao_factory(L"<i2 j2|R|a' a>");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i1,j1,i2,j2") =
        (2.0 * left * middle * right).set_shape(ijij_ijji_shape);
    //    B_ijij_ijji("i1,j1,i2,j2") -= (2.0*lcao_factory(L"(j1 p|R|i1
    //    a)[df]")*lcao_factory(L"(p|F|a')[df]")*lcao_factory(L"(i2 a|R|j2
    //    a')[df]")).set_shape(ijij_ijji_shape);
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
    B_ijij_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "B Term8 Time: ", time, " S\n");
  }

  auto b_time1 = mpqc_time::now(world, accurate_time);
  auto b_time = mpqc_time::duration_in_s(b_time0, b_time1);
  utility::print_par(world, "B Term Total Time: ", b_time, " S\n");
  return B_ijij_ijji;
};

/**
 * CC-F12 C approach V term with DF
 * \f$V_{ia}^{xy}\f$
 * @param lcao_factory reference to LCAOFactory
 * @return V("i,a,x,y")
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_V_iaxy_df(
    integrals::LCAOFactory<Tile, Policy> &lcao_factory) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();
  TA::DistArray<Tile, Policy> V_iaxy;

  auto v_time0 = mpqc_time::now(world, accurate_time);

  utility::print_par(world, "\nCompute V_iaxy With DF \n");
  {
    auto left = lcao_factory(L"(Κ |GR|i k)");
    auto middle = lcao_factory(L"(Κ|GR|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |GR|a l)");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_iaxy("i,a,k,l") = left * middle * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|p q>[df]");
    auto right = lcao_factory(L"<k l|R|p q>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|m a'>[df]");
    auto right = lcao_factory(L"<k l|R|m a'>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<a i|G|m a'>[df]");
    auto right = lcao_factory(L"<l k|R|m a'>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term4 Time: ", time, " S\n");
  }

  auto v_time1 = mpqc_time::now(world, accurate_time);
  auto v_time = mpqc_time::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");
  return V_iaxy;
};

/**
 * CC-F12 C approach V term without DF
 * \f$V_{ia}^{xy}\f$
 * @param lcao_factory reference to LCAOFactory
 * @return V("i,a,x,y")
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_V_iaxy(
    integrals::LCAOFactory<Tile, Policy> &lcao_factory) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();
  TA::DistArray<Tile, Policy> V_iaxy;

  auto v_time0 = mpqc_time::now(world, accurate_time);

  utility::print_par(world, "\nCompute V_iaxy Without DF \n");
  {
    auto left = lcao_factory(L"<i a |GR|k l>");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_iaxy("i,a,k,l") = left;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|p q>");
    auto right = lcao_factory(L"<k l|R|p q>");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i a|G|m a'>");
    auto right = lcao_factory(L"<k l|R|m a'>");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<a i|G|m a'>");
    auto right = lcao_factory(L"<l k|R|m a'>");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_iaxy("i,a,k,l") -= left * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term4 Time: ", time, " S\n");
  }

  auto v_time1 = mpqc_time::now(world, accurate_time);
  auto v_time = mpqc_time::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");
  return V_iaxy;
};

/**
 * DF-based builder for the V intermediate with general indices,
 * \f$V_{xy}^{ab} \equiv R_{xy}^{\alpha \beta} g_{\alpha \beta}^{a b}\f$, for
 * the use in CC F12.
 * @param lcao_factory reference to LCAOFactory
 * @return V("x,y,a,b")
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_V_xyab_df(
    integrals::LCAOFactory<Tile, Policy> &lcao_factory) {
  auto &world = lcao_factory.get_world();
  auto &ao_integral = lcao_factory.atomic_integral();
  bool accurate_time = lcao_factory.accurate_time();

  auto v_time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, Policy> V_xyab;
  TA::DistArray<Tile, Policy> tmp;

  utility::print_par(world, "\nCompute V_xyab With DF \n");

  {
    auto left = lcao_factory(L"(Κ |GR|i a)");
    auto middle = ao_integral(L"(Κ|GR|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |GR|j b)");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_xyab("i,j,a,b") = left * middle * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto right = lcao_factory(L"<a b|G|p q>[df]");
    auto left = lcao_factory(L"<i j|R|p q>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_xyab("i,j,a,b") -= left * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  {
    auto right = lcao_factory(L"<a b|G|m a'>[df]");
    auto left = lcao_factory(L"<i j|R|m a'>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    tmp("i,j,a,b") = left * right;
    V_xyab("i,j,a,b") -= tmp("i,j,a,b");
    V_xyab("i,j,a,b") -= tmp("j,i,b,a");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term3 Time: ", time, " S\n");
  }

  auto v_time1 = mpqc_time::now(world, accurate_time);
  auto v_time = mpqc_time::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");

  return V_xyab;
};

/**
 * non-DF-based builder for the V intermediate with general indices,
 * \f$V_{xy}^{ab} \equiv R_{xy}^{\alpha \beta} g_{\alpha \beta}^{a b}\f$, for
 * the use in CC F12.
 * @param lcao_factory reference to LCAOFactory
 * @return V("x,y,a,b")
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_V_xyab(
    integrals::LCAOFactory<Tile, Policy> &lcao_factory) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();

  auto v_time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, Policy> V_xyab;
  TA::DistArray<Tile, Policy> tmp;

  utility::print_par(world, "\nCompute V_xyab Without DF \n");

  {
    auto left = lcao_factory(L"<i j|GR|a b>");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_xyab("i,j,a,b") = left;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term1 Time: ", time, " S\n");
  }

  {
    auto right = lcao_factory(L"<a b|G|p q>");
    auto left = lcao_factory(L"<i j|R|p q>");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_xyab("i,j,a,b") -= left * right;
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "V Term2 Time: ", time, " S\n");
  }

  {
    auto right = lcao_factory(L"<a b|G|m a'>");
    auto left = lcao_factory(L"<i j|R|m a'>");
    tmp("i,j,a,b") = left * right;
    V_xyab("i,j,a,b") -= tmp("i,j,a,b");
    V_xyab("i,j,a,b") -= tmp("j,i,b,a");
  }

  auto v_time1 = mpqc_time::now(world, accurate_time);
  auto v_time = mpqc_time::duration_in_s(v_time0, v_time1);
  utility::print_par(world, "V Term Total Time: ", v_time, " S\n");

  return V_xyab;
};

/**
 * MP2-F12, CC-F12 C approach C term \f$C_{ij}^{ab} \f$ with DF
 * @param lcao_factory reference to LCAOFactory
 * @return C("i,j,a,b")
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_C_ijab_df(
    integrals::LCAOFactory<Tile, Policy> &lcao_factory) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();
  auto c_time0 = mpqc_time::now(world, accurate_time);
  TA::DistArray<Tile, Policy> C_ijab;

  utility::print_par(world, "\nCompute C_ijab With DF \n");

  auto left = lcao_factory(L"<i j|R|a a'>[df]");
  auto right = lcao_factory(L"<a'|F|b>[df]");

  auto time0 = mpqc_time::now(world, accurate_time);
  C_ijab("i,j,a,b") = left * right;
  C_ijab("i,j,a,b") += C_ijab("j,i,b,a");
  auto time1 = mpqc_time::now(world, accurate_time);
  auto time = mpqc_time::duration_in_s(time0, time1);
  utility::print_par(world, "C Term Time: ", time, " S\n");

  auto c_time = mpqc_time::duration_in_s(c_time0, time1);
  utility::print_par(world, "C Term Total Time: ", c_time, " S\n");
  return C_ijab;
};

/**
 * MP2-F12, CC-F12 C approach C term \f$C_{ij}^{ab} \f$ without DF
 * @param lcao_factory reference to LCAOFactory
 * @return C("i,j,a,b")
 */
template <typename Tile, typename Policy>
TA::DistArray<Tile, Policy> compute_C_ijab(
    integrals::LCAOFactory<Tile, Policy> &lcao_factory) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();
  auto c_time0 = mpqc_time::now(world, accurate_time);
  TA::DistArray<Tile, Policy> C_ijab;

  utility::print_par(world, "\nCompute C_ijab Without DF \n");

  auto left = lcao_factory(L"<i j|R|a a'>");
  auto right = lcao_factory(L"<a'|F|b>");

  auto time0 = mpqc_time::now(world, accurate_time);
  C_ijab("i,j,a,b") = left * right;
  C_ijab("i,j,a,b") += C_ijab("j,i,b,a");
  auto time1 = mpqc_time::now(world, accurate_time);
  auto time = mpqc_time::duration_in_s(time0, time1);
  utility::print_par(world, "C Term Time: ", time, " S\n");

  auto c_time = mpqc_time::duration_in_s(c_time0, time1);
  utility::print_par(world, "C Term Total Time: ", c_time, " S\n");
  return C_ijab;
};

/**
 * CC-F12 C approach VT2 term with direct integral
 * \f$T_{ab}^{ij} * (V_{xy}^{ab} + C_{xy}^{ab})\f$
 * @param lcao_factory reference to LCAOFactory
 * @param t2 t2 amplitude
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @param direct_array direct two electron integral \f$V_{\rho \sigma}^{\mu
 * \nu}\f$
 * @return V("i1,j1,i2,j2")
 */
template <typename Tile, typename DirectArray>
TA::DistArray<Tile, TA::SparsePolicy> compute_VT2_ijij_ijji_df_direct(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    const TA::DistArray<Tile, TA::SparsePolicy> &t2,
    const TA::SparseShape<float> &ijij_ijji_shape, DirectArray direct_array) {
  auto &world = lcao_factory.get_world();
  auto &ao_integral = lcao_factory.atomic_integral();
  bool accurate_time = lcao_factory.accurate_time();

  TA::DistArray<Tile, TA::SparsePolicy> V_ijji_ijji;

  // C Term
  auto V_xyab = compute_C_ijab_df(lcao_factory);

  utility::print_par(world, "\nCompute VT2_ijij_ijji With DF and Direct AO\n");
  auto vt2_time0 = mpqc_time::now(world, accurate_time);
  {
    auto left = lcao_factory(L"(Κ |GR|i a)");
    auto middle = ao_integral(L"(Κ|GR|Λ)[inv]");
    auto right = lcao_factory(L"(Λ |GR|j b)");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_xyab("i,j,a,b") += left * middle * right;
    lcao_factory.registry().purge_operator(world, L"GR");
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "VT2 Term1 Time: ", time, " S\n");
  }

  {
    auto right = t2("a,b,i1,j1");
    auto left = V_xyab("i2,j2,a,b");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_ijji_ijji("i1,j1,i2,j2") = (left * right).set_shape(ijij_ijji_shape);
    // clean V_xyab
    V_xyab = TA::DistArray<Tile, TA::SparsePolicy>();
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "VT2 Term2 Time: ", time, " S\n");
  }

  TA::DistArray<Tile, TA::SparsePolicy> U;
  auto Ca = lcao_factory.orbital_space()->retrieve(OrbitalIndex(L"a")).array();
  auto Cp = lcao_factory.orbital_space()->retrieve(OrbitalIndex(L"p")).array();
  {
    auto time0 = mpqc_time::now(world, accurate_time);

    //    auto Cm =
    //    lcao_factory.orbital_space()->retrieve(OrbitalIndex(L"m")).array();
    //    auto Ca_prime =
    //    lcao_factory.orbital_space()->retrieve(OrbitalIndex(L"a'")).array();
    // compuate intermediate U
    U("i,j,rho,mu") = (t2("a,b,i,j") * Ca("sigma,a") * Ca("nu,b")) *
                      direct_array("rho, sigma, mu, nu");

    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "VT2 UTerm Time: ", time, " S\n");
  }

  {
    auto left = lcao_factory(L"<i2 j2|R|p q>[df]");

    auto time0 = mpqc_time::now(world, accurate_time);
    V_ijji_ijji("i1,j1,i2,j2") -=
        ((left * Cp("rho,p") * Cp("mu,q")) * U("i1,j1,rho,mu"))
            .set_shape(ijij_ijji_shape);
    auto time1 = mpqc_time::now(world, accurate_time);
    auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, "VT2 Term3 Time: ", time, " S\n");
  }
  //    tmp("i1,j1,i2,j2") = ((lcao_factory(L"(i2 m|R|j2
  //    a')[df]")*Cm("mu,m")*Ca_prime("nu,a'"))*U("i1,j1,mu,nu")).set_shape(ijij_ijji_shape);
  //    V_ijji_ijji("i1,j1,i2,j2") -= tmp("i1,j1,i2,j2");
  //    V_ijji_ijji("i1,j1,i2,j2") -= tmp("j1,i1,j2,i2");

  auto vt2_time1 = mpqc_time::now(world, accurate_time);
  auto vt2_time = mpqc_time::duration_in_s(vt2_time0, vt2_time1);
  utility::print_par(world, "VT2 Term Total Time: ", vt2_time, " S\n");

  return V_ijji_ijji;
};

/**
 * CC-F12 C approach VT2 term with DF
 * \f$T_{ab}^{ij} * (V_{xy}^{ab} + C_{xy}^{ab})\f$
 * @param lcao_factory reference to LCAOFactory
 * @param t2 t2 amplitude
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return V("i1,j1,i2,j2")
 */

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_VT2_ijij_ijji_df(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    const TA::DistArray<Tile, TA::SparsePolicy> &t2,
    const TA::SparseShape<float> &ijij_ijji_shape) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();

  TA::DistArray<Tile, TA::SparsePolicy> V_ijij_ijji;

  // compute C_ijab
  TA::DistArray<Tile, TA::SparsePolicy> C_ijab =
      compute_C_ijab_df(lcao_factory);

  // compute V_ijab
  TA::DistArray<Tile, TA::SparsePolicy> V_ijab =
      compute_V_xyab_df(lcao_factory);

  auto vt2_time0 = mpqc_time::now(world, accurate_time);
  utility::print_par(world, "\nCompute VT2_ijij_ijji With DF\n");
  V_ijij_ijji("i1,j1,i2,j2") =
      ((V_ijab("i2,j2,a,b") + C_ijab("i2,j2,a,b")) * t2("a,b,i1,j1"))
          .set_shape(ijij_ijji_shape);

  auto vt2_time1 = mpqc_time::now(world, accurate_time);
  auto vt2_time = mpqc_time::duration_in_s(vt2_time0, vt2_time1);
  utility::print_par(world, "VT2 Term Total Time: ", vt2_time, " S\n");

  return V_ijij_ijji;
};

/**
 * CC-F12 C approach VT2 term without DF
 * \f$T_{ab}^{ij} * (V_{xy}^{ab} + C_{xy}^{ab})\f$
 * @param lcao_factory reference to LCAOFactory
 * @param t2 t2 amplitude
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return V("i1,j1,i2,j2")
 */

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_VT2_ijij_ijji(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    const TA::DistArray<Tile, TA::SparsePolicy> &t2,
    const TA::SparseShape<float> &ijij_ijji_shape) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();

  TA::DistArray<Tile, TA::SparsePolicy> V_ijij_ijji;
  //    TA::DistArray<Tile,TA::SparsePolicy> tmp;

  // compute C_ijab
  TA::DistArray<Tile, TA::SparsePolicy> C_ijab = compute_C_ijab(lcao_factory);

  // compute V_ijab
  TA::DistArray<Tile, TA::SparsePolicy> V_ijab = compute_V_xyab(lcao_factory);

  auto vt2_time0 = mpqc_time::now(world, accurate_time);
  utility::print_par(world, "\nCompute VT2_ijij_ijji Without DF\n");
  V_ijij_ijji("i1,j1,i2,j2") =
      ((V_ijab("i2,j2,a,b") + C_ijab("i2,j2,a,b")) * t2("a,b,i1,j1"))
          .set_shape(ijij_ijji_shape);

  //    V_ijij_ijji("i1,j1,i2,j2") =
  //    ((V_ijab("i2,j2,a,b"))*t2("a,b,i1,j1")).set_shape(ijij_ijji_shape);
  //    std::cout << "VT2 Term " << std::endl;
  //    std::cout << V_ijij_ijji << std::endl;

  //    tmp("i1,j1,i2,j2") =
  //    ((C_ijab("i2,j2,a,b"))*t2("a,b,i1,j1")).set_shape(ijij_ijji_shape);
  //    std::cout << "CT2 Term " << std::endl;
  //    std::cout << tmp << std::endl;

  //    V_ijij_ijji("i1,j1,i2,j2") += tmp("i1,j1,i2,j2");

  auto vt2_time1 = mpqc_time::now(world, accurate_time);
  auto vt2_time = mpqc_time::duration_in_s(vt2_time0, vt2_time1);
  utility::print_par(world, "VT2 Term Total Time: ", vt2_time, " S\n");

  return V_ijij_ijji;
};

/**
 * CC-F12 C approach VT1 term with DF
 * \f$T_{a}^{i} * V_{ia}^{xy}\f$
 * @param lcao_factory reference to LCAOFactory
 * @param t1 t1 amplitude
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return V("i1,j1,i2,j2")
 */

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_VT1_ijij_ijji_df(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    const TA::DistArray<Tile, TA::SparsePolicy> &t1,
    const TA::SparseShape<float> &ijij_ijji_shape) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();
  TA::DistArray<Tile, TA::SparsePolicy> V_ijij_ijji;
  TA::DistArray<Tile, TA::SparsePolicy> V_iaij =
      compute_V_iaxy_df(lcao_factory);

  auto vt1_time0 = mpqc_time::now(world, accurate_time);
  utility::print_par(world, "\nCompute VT1_ijij_ijji With DF\n");

  V_ijij_ijji("i1,j1,i2,j2") =
      (V_iaij("i1,a,i2,j2") * t1("a,j1")).set_shape(ijij_ijji_shape);
  V_ijij_ijji("i1,j1,i2,j2") += V_ijij_ijji("j1,i1,j2,i2");
  auto vt1_time1 = mpqc_time::now(world, accurate_time);
  auto vt1_time = mpqc_time::duration_in_s(vt1_time0, vt1_time1);
  utility::print_par(world, "VT1 Term Total Time: ", vt1_time, " S\n");

  return V_ijij_ijji;
};

/**
 * CC-F12 C approach VT1 term without DF
 * \f$T_{a}^{i} * V_{ia}^{xy}\f$
 * @param lcao_factory reference to LCAOFactory
 * @param t1 t1 amplitude
 * @param ijij_ijji_shape SparseShape that has ijij ijji shape
 * @return V("i1,j1,i2,j2")
 */

template <typename Tile>
TA::DistArray<Tile, TA::SparsePolicy> compute_VT1_ijij_ijji(
    integrals::LCAOFactory<Tile, TA::SparsePolicy> &lcao_factory,
    const TA::DistArray<Tile, TA::SparsePolicy> &t1,
    const TA::SparseShape<float> &ijij_ijji_shape) {
  auto &world = lcao_factory.get_world();
  bool accurate_time = lcao_factory.accurate_time();
  TA::DistArray<Tile, TA::SparsePolicy> V_ijij_ijji;
  TA::DistArray<Tile, TA::SparsePolicy> V_iaij = compute_V_iaxy(lcao_factory);

  auto vt2_time0 = mpqc_time::now(world, accurate_time);
  utility::print_par(world, "\nCompute VT1_ijij_ijji Without DF\n");

  V_ijij_ijji("i1,j1,i2,j2") =
      (V_iaij("i1,a,i2,j2") * t1("a,j1")).set_shape(ijij_ijji_shape);

  //    std::cout << "VT1 First" << std::endl;
  //    std::cout << V_ijij_ijji << std::endl;

  V_ijij_ijji("i1,j1,i2,j2") += V_ijij_ijji("j1,i1,j2,i2");
  //    std::cout << "VT1 Second" << std::endl;
  //    std::cout << V_ijij_ijji << std::endl;

  auto vt2_time1 = mpqc_time::now(world, accurate_time);
  auto vt2_time = mpqc_time::duration_in_s(vt2_time0, vt2_time1);
  utility::print_par(world, "VT1 Term Total Time: ", vt2_time, " S\n");

  return V_ijij_ijji;
};

/**
 * the DF-based builder for the F12 intermediates V (or X). The V intermediate
 * is defined as
 * \f$V_{pq}^{rs} \equiv R_{pq}^{\alpha \beta} g_{\alpha \beta}^{r s}\f$ and
 * \f$V_{pq}^{sr} \equiv R_{pq}^{\alpha \beta} g_{\alpha \beta}^{r s}\f$,
 * and the X intermediate is obtained from these expressions by replacement \f$
 * g \to R \f$.
 * If \c p refers to the same space as \c q **or** \c r refers to the same space
 * as \c s , the two
 * tensors are equivalent and only the former is computed.
 *
 * @tparam String any string type (e.g., std::basic_string and char[])
 * @param target_str std::string, the only valid vales are "V" or "X"
 * @param lcao_factory reference to LCAOFactory
 * @param p an OrbitalIndex key (must be known to \c lcao_factory )
 * @param q an OrbitalIndex key (must be known to \c lcao_factory )
 * @param r an OrbitalIndex key (must be known to \c lcao_factory )
 * @param s an OrbitalIndex key (must be known to \c lcao_factory )
 * @param df if \c true , use density fitting (\c df=false is not yet supported)
 * @param cabs if \c false , skip the CABS contributions; the default value is
 * \c true
 * @return \c std::tuple with V("p,q,r,s") and V("p,q,s,r"); the latter
 *         is empty if \c p refers to the same space as \c q **or** \c r refers
 * to the same space as \c s .
 */
template <typename Tile, typename Policy, typename String>
std::tuple<TA::DistArray<Tile, Policy>, TA::DistArray<Tile, Policy>>
VX_pqrs_pqsr(const std::string &target_str,
             integrals::LCAOFactory<Tile, Policy> &lcao_factory,
             const String &p, const String &q, const String &r, const String &s,
             bool df = true, bool cabs = true) {
  using mpqc::utility::concat;
  using mpqc::utility::wconcat;
  using mpqc::utility::concatcm;
  using mpqc::utility::to_wstring;

  enum class Target { X, V };
  Target target;
  if (target_str == "V")
    target = Target::V;
  else if (target_str == "X")
    target = Target::X;
  else
    TA_USER_ASSERT(false, "invalid target type");

  TA_USER_ASSERT(df == true,
                 "non-DF-based generic VX build is not yet supported");
  const auto methodstr = df ? "[df]" : "";

  auto &world = lcao_factory.get_world();
  auto &ao_integral = lcao_factory.atomic_integral();
  const bool accurate_time = lcao_factory.accurate_time();

  const auto equiv_pq = OrbitalIndex(p) == OrbitalIndex(q);
  const auto equiv_rs = OrbitalIndex(r) == OrbitalIndex(s);
  // if !equiv(p,q) && !equiv(r,s), must compute both
  const auto need_pqsr = !equiv_pq && !equiv_rs;
  // if !equiv(p,q) || !equiv(r,s), particles are not equivalent =>
  // must include ma' and a'm contributions explicitly, rather than by
  // symmetrization
  const auto p1_equiv_p2 = equiv_pq && equiv_rs;

  const auto time0 = mpqc_time::now(world, accurate_time);

  TA::DistArray<Tile, Policy> A_pqrs, A_pqsr;

  if (need_pqsr) {
    utility::print_par(world, "Compute ", target_str, "(", p, ",", q, ",", r,
                       ",", s, ") and ", target_str, "(", p, ",", q, ",", s,
                       ",", r, ") DF=", std::to_string(df), "\n");
  } else {
    utility::print_par(world, "Compute ", target_str, "(", p, ",", q, ",", r,
                       ",", s, ") DF=", std::to_string(df), "\n");
  }

  {
    const char *opstr = (target == Target::V) ? "GR" : "R2";
    const auto left_pr =
        lcao_factory(wconcat(L"(Κ |", opstr, "|", p, " ", r, ")"));
    const auto middle = ao_integral(wconcat(L"(Κ|", opstr, L"|Λ)[inv]"));
    const auto right_qs =
        lcao_factory(wconcat(L"(Λ |", opstr, "|", q, " ", s, ")"));

    const auto time0 = mpqc_time::now(world, accurate_time);
    A_pqrs(concatcm(p, q, r, s)) = left_pr * middle * right_qs;
    if (need_pqsr) {
      const auto left_ps =
          lcao_factory(wconcat(L"(Κ |", opstr, "|", p, " ", s, ")"));
      const auto right_qr =
          lcao_factory(wconcat(L"(Λ |", opstr, "|", q, " ", r, ")"));
      A_pqsr(concatcm(p, q, s, r)) = left_ps * middle * right_qr;
    }
    // remove all generated integrals as they are likely not needed any longer
    lcao_factory.purge_operator(world, to_wstring(opstr));
    const auto time1 = mpqc_time::now(world, accurate_time);
    const auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, target_str, " Term1 Time: ", time, " S\n");
  }

  {
    const auto rightopstr = (target == Target::V) ? "G" : "R";

    const auto left =
        lcao_factory(wconcat("<", p, " ", q, "|R|p0 q0>", methodstr));
    const auto right_rs = lcao_factory(
        wconcat("<", r, " ", s, "|", rightopstr, "|p0 q0>", methodstr));

    const auto time0 = mpqc_time::now(world, accurate_time);
    A_pqrs(concatcm(p, q, r, s)) -= left * right_rs;
    if (need_pqsr) {
      const auto right_sr = lcao_factory(
          wconcat("<", s, " ", r, "|", rightopstr, "|p0 q0>", methodstr));
      A_pqsr(concatcm(p, q, s, r)) -= left * right_sr;
    }
    const auto time1 = mpqc_time::now(world, accurate_time);
    const auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, target_str, " Term2 Time: ", time, " S\n");
  }

  if (cabs) {
    const auto rightopstr = (target == Target::V) ? "G" : "R";

    const auto pqRmA =
        lcao_factory(wconcat("<", p, " ", q, "|R|m0 a'0>", methodstr));
    // Y == R (if need X) or G (if need V)
    const auto rsYmA = lcao_factory(
        wconcat("<", r, " ", s, "|", rightopstr, "|m0 a'0>", methodstr));
    const auto qpRmA =
        lcao_factory(wconcat("<", q, " ", p, "|R|m0 a'0>", methodstr));
    const auto srYmA = lcao_factory(
        wconcat("<", s, " ", r, "|", rightopstr, "|m0 a'0>", methodstr));

    const auto time0 = mpqc_time::now(world, accurate_time);
    {
      const auto pqrs_str = concatcm(p, q, r, s);
      if (p1_equiv_p2) {  // compute by symmetrization
        TA::DistArray<Tile, Policy> tmp;
        tmp(pqrs_str) = pqRmA * rsYmA;
        A_pqrs(pqrs_str) -= tmp(pqrs_str);
        A_pqrs(pqrs_str) -= tmp(concatcm(q, p, s, r));
      } else {
        A_pqrs(pqrs_str) -= pqRmA * rsYmA;
        A_pqrs(pqrs_str) -= qpRmA * srYmA;
      }
    }

    if (need_pqsr) {
      const auto pqsr_str = concatcm(p, q, s, r);
      if (p1_equiv_p2) {  // compute by symmetrization
        TA::DistArray<Tile, Policy> tmp;
        tmp(pqsr_str) = pqRmA * srYmA;
        A_pqsr(pqsr_str) -= tmp(pqsr_str);
        A_pqsr(pqsr_str) -= tmp(concatcm(q, p, r, s));
      } else {
        A_pqsr(pqsr_str) -= pqRmA * srYmA;
        A_pqsr(pqsr_str) -= qpRmA * rsYmA;
      }
    }

    const auto time1 = mpqc_time::now(world, accurate_time);
    const auto time = mpqc_time::duration_in_s(time0, time1);
    utility::print_par(world, target_str, " Term3 Time: ", time, " S\n");
  }

  const auto time1 = mpqc_time::now(world, accurate_time);
  const auto time = mpqc_time::duration_in_s(time0, time1);
  utility::print_par(world, target_str, " Total Time: ", time, " S\n");

  //  std::cout << "A_pqrs:\n" << A_pqrs << std::endl;
  //  if (need_pqsr) std::cout << "A_pqsr:\n" << A_pqsr << std::endl;

  return std::make_tuple(A_pqrs, A_pqsr);
}

}  // end of namespace f12
}  // end of namespace mpqc

#endif  // MPQC_F12_INTERMEDIATES_H
