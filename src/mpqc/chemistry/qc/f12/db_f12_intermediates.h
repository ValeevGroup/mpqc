//
// Created by Chong Peng on 7/11/16.
//

#ifndef MPQC_DB_F12_INTERMEDIATES_H
#define MPQC_DB_F12_INTERMEDIATES_H


#include "../../../../../include/tiledarray.h"
#include "../../../../../common/namespaces.h"
#include <mpqc/chemistry/qc/integrals/molecular_integral.h>

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

    std::cout << V_ijij_ijji << std::endl;
    // all types of GR integral not needed
    mo_integral.remove_operation_all(world, L"GR");

    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term1 Time: ", time, " S\n");
  }


  {
    auto left = mo_integral(L"<i1 j1|G|i3 j3>[df]");
    auto right = mo_integral(L"<i2 j2|R|i3 j3>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    V_ijij_ijji("i1,j1,i2,j2") -= (left*right).set_shape(shape);
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"V Term2 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i1 j1|G|i3 a>[df]");
    auto right = mo_integral(L"<i2 j2|R|i3 a>[df]");

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

  std::cout << V_ijij_ijji << std::endl;

  {
    auto left = mo_integral(L"<i1 j1|G|m a'>[df]");
    auto right = mo_integral(L"<i2 j2|R|m a'>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    TA::DistArray<Tile,TA::SparsePolicy> tmp;
    tmp("i1,j1,i2,j2") = (left*right).set_shape(shape);
    std::cout << tmp << std::endl;
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

  std::cout << V_ijij_ijji << std::endl;
  return V_ijij_ijji;
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
    auto left = mo_integral(L"<i1 j1|R|i3 j3>[df]");
    auto right = mo_integral(L"<i2 j2|R|i3 j3>[df]");

    auto time0 = mpqc_time::now(world,accurate_time);
    X_ijij_ijji("i1,j1,i2,j2") -= (left*right).set_shape(ijij_ijji_shape);
    auto time1 = mpqc_time::now(world,accurate_time);
    auto time = mpqc_time::duration_in_s(time0,time1);
    utility::print_par(world,"X Term2 Time: ", time, " S\n");
  }

  {
    auto left = mo_integral(L"<i1 j1|R|i3 a>[df]");
    auto right = mo_integral(L"<i2 j2|R|i3 a>[df]");

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

}
}


#endif //MPQC_DB_F12_INTERMEDIATES_H
