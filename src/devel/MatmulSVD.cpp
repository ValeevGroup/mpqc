/*
 *  This file is a part of TiledArray.
 *  Copyright (C) 2014  Virginia Tech
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "TileSVD.hpp"
#include <Eigen/Dense>
#include <chemistry/qc/libint2/libint2.h>
#include <util/madness/init.h>
#include <util/madness/world.h>
#include <chemistry/qc/basis/tiledbasisset.hpp>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/integralenginepool.hpp>
#include <chemistry/qc/basis/taskintegrals.hpp>
#include <chemistry/qc/scf/cldfgengine.hpp>
#include <math/elemental/eigensolver.hpp>
#include <chemistry/qc/lcao/soad.h>
#include <utility>
#include <iostream>

using namespace sc;
using namespace mpqc;
using namespace mpqc::TA;

using Matrix = TiledArray::Array<double, 2>;
using Eri3 = TiledArray::Array<double, 3>;
using Eri4 = TiledArray::Array<double, 4>;
using onebpool = IntegralEnginePool<Ref<OneBodyInt> >;
using e3pool = IntegralEnginePool<Ref<TwoBodyThreeCenterInt> >;
using e2pool = IntegralEnginePool<Ref<TwoBodyTwoCenterInt> >;
using e4pool = IntegralEnginePool<Ref<TwoBodyInt> >;
using mat_pmap_iter = decltype(std::declval<Matrix>().get_pmap()->begin());
using eri3_pmap_iter = decltype(std::declval<Eri3>().get_pmap()->begin());
using eri4_pmap_iter = decltype(std::declval<Eri4>().get_pmap()->begin());

void mat_task(mat_pmap_iter it, Matrix array, double cut,
              long* t_mat_elems,
              long* t_svd_elems, double* t_mat_diff,
              std::shared_ptr<onebpool> pool){

  using PoolPtrType = typename std::pointer_traits<std::shared_ptr<onebpool>>::element_type;
  typename PoolPtrType::engine_type engine =
      pool->instance();

  typename Matrix::value_type tile(array.trange().make_tile_range(*it));
  get_integrals(tile, engine);

  std::size_t size0 = tile.range().size()[0];
  std::size_t size1 = tile.range().size()[1];

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
  Emap = TiledArray::eigen_map(tile, size0, size1);

  Eigen::MatrixXd Mat = Emap;
  SVDTile SVDMat(Mat,cut);

  *t_mat_elems += Mat.rows() * Mat.cols();
  *t_svd_elems += SVDMat.nelements();

  double mat_diff = Eigen::MatrixXd(Mat -
     SVDMat.U_contract() * SVDMat.Vt()).lpNorm<Eigen::Infinity>();
  *t_mat_diff += mat_diff;
}

void eri2_task(mat_pmap_iter it, Matrix array, double cut,
              long* t_mat_elems,
              long* t_svd_elems, double* t_mat_diff,
              std::shared_ptr<e2pool> pool){

  using PoolPtrType = typename std::pointer_traits<std::shared_ptr<e2pool>>::element_type;
  typename PoolPtrType::engine_type engine =
      pool->instance();

  typename Matrix::value_type tile(array.trange().make_tile_range(*it));
  get_integrals(tile, engine);

  std::size_t size0 = tile.range().size()[0];
  std::size_t size1 = tile.range().size()[1];

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
  Emap = TiledArray::eigen_map(tile, size0, size1);

  Eigen::MatrixXd Mat = Emap;
  SVDTile SVDMat(Mat,cut);

  *t_mat_elems += Mat.rows() * Mat.cols();
  *t_svd_elems += SVDMat.nelements();

  double mat_diff = Eigen::MatrixXd(Mat -
     SVDMat.U_contract() * SVDMat.Vt()).lpNorm<Eigen::Infinity>();
  *t_mat_diff += mat_diff;
}

void eri3_task(eri3_pmap_iter it, Eri3 array, double cut,
              long* t_mat_elems,
              long* t_svd_elems, double* t_mat_diff,
              std::shared_ptr<e3pool> pool){

  using PoolPtrType = typename std::pointer_traits<std::shared_ptr<e3pool>>::element_type;
  typename PoolPtrType::engine_type engine =
      pool->instance();

  typename Matrix::value_type tile(array.trange().make_tile_range(*it));
  get_integrals(tile, engine);

  std::size_t size0 = tile.range().size()[0];
  std::size_t size1 = tile.range().size()[1];
  std::size_t size2 = tile.range().size()[2];

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
  Emap = TiledArray::eigen_map(tile, size0*size1, size2);

  Eigen::MatrixXd Mat = Emap;
  SVDTile SVDMat(Mat,cut);

  *t_mat_elems += Mat.rows() * Mat.cols();
  *t_svd_elems += SVDMat.nelements();

  double mat_diff = Eigen::MatrixXd(Mat -
     SVDMat.U_contract() * SVDMat.Vt()).lpNorm<Eigen::Infinity>();
  *t_mat_diff += mat_diff;
}

#if 1
void eri4_task(eri4_pmap_iter it, Eri4 array, double cut,
              long* t_mat_elems,
              long* t_svd_elems, double* t_mat_diff,
              std::shared_ptr<e4pool> pool){

  using PoolPtrType = typename std::pointer_traits<std::shared_ptr<e4pool>>::element_type;
  typename PoolPtrType::engine_type engine =
      pool->instance();

  typename Matrix::value_type tile(array.trange().make_tile_range(*it));
  get_integrals(tile, engine);

  std::size_t size0 = tile.range().size()[0];
  std::size_t size1 = tile.range().size()[1];
  std::size_t size2 = tile.range().size()[2];
  std::size_t size3 = tile.range().size()[3];

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
  Emap = TiledArray::eigen_map(tile, size0*size1, size2*size3);

  Eigen::MatrixXd Mat = Emap;
  SVDTile SVDMat(Mat,cut);

  *t_mat_elems += Mat.rows() * Mat.cols();
  *t_svd_elems += SVDMat.nelements();

  double mat_diff = Eigen::MatrixXd(Mat -
     SVDMat.U_contract() * SVDMat.Vt()).lpNorm<Eigen::Infinity>();
  *t_mat_diff += mat_diff;
}
#endif

void percent_compression(Matrix &mat, std::shared_ptr<onebpool> pool, double cut){

  long total_mat_elems = 0;
  long total_svd_elems = 0;
  double total_mat_diff = 0;
  long ntiles = 0;
  double svd_time0 = madness::wall_time();
  for(auto it = mat.get_pmap()->begin(); it != mat.get_pmap()->end(); ++it){
    mat.get_world().taskq.add(&mat_task, it, mat, cut, &total_mat_elems,
                              &total_svd_elems, &total_mat_diff, pool);
    ++ntiles;
  }
  mat.get_world().gop.fence();
  double svd_time1 = madness::wall_time();

  std::cout << "\t\tTotal Time = " <<
            svd_time1 - svd_time0  << "s" << std::endl;
  std::cout << "\t\tTotal Matrix Elements = " << total_mat_elems << std::endl;
  std::cout << "\t\tTotal SVD Elements = " << total_svd_elems << std::endl;
  std::cout << "\t\tAverage tile diff = " << total_mat_diff/ntiles << std::endl;
  std::cout << "\t\tPercentage saved = " << (1 -
               double(total_svd_elems)/double(total_mat_elems)) * 100
            << std::endl;
}

void percent_compression(Matrix &mat, std::shared_ptr<e2pool> pool, double cut){

  long total_mat_elems = 0;
  long total_svd_elems = 0;
  double total_mat_diff = 0;
  long ntiles = 0;
  double svd_time0 = madness::wall_time();
  for(auto it = mat.get_pmap()->begin(); it != mat.get_pmap()->end(); ++it){
    mat.get_world().taskq.add(&eri2_task, it, mat, cut, &total_mat_elems,
                              &total_svd_elems, &total_mat_diff, pool);
    ++ntiles;
  }
  mat.get_world().gop.fence();
  double svd_time1 = madness::wall_time();


  std::cout << "\t\tTotal Time = " <<
            svd_time1 - svd_time0  << "s" << std::endl;
  std::cout << "\t\tTotal Matrix Elements = " << total_mat_elems << std::endl;
  std::cout << "\t\tTotal SVD Elements = " << total_svd_elems << std::endl;
  std::cout << "\t\tAverage tile diff = " << total_mat_diff/ntiles << std::endl;
  std::cout << "\t\tPercentage saved = " << (1 -
               double(total_svd_elems)/double(total_mat_elems)) * 100
            << std::endl;
}

void percent_compression(Eri3 &mat, std::shared_ptr<e3pool> pool, double cut){

  long total_mat_elems = 0;
  long total_svd_elems = 0;
  double total_mat_diff = 0;
  long ntiles = 0;
  double svd_time0 = madness::wall_time();
  for(auto it = mat.get_pmap()->begin(); it != mat.get_pmap()->end(); ++it){
    mat.get_world().taskq.add(&eri3_task, it, mat, cut, &total_mat_elems,
                              &total_svd_elems, &total_mat_diff, pool);
    ++ntiles;
  }
  mat.get_world().gop.fence();
  double svd_time1 = madness::wall_time();

  std::cout << "\t\tNumber of tiles = " << ntiles << " in "
            << svd_time1 - svd_time0 << "s" << std::endl;
  std::cout << "\t\tTotal Matrix Elements = " << total_mat_elems << std::endl;
  std::cout << "\t\tTotal SVD Elements = " << total_svd_elems << std::endl;
  std::cout << "\t\tAverage tile diff = " << total_mat_diff/ntiles << std::endl;
  std::cout << "\t\tPercentage saved = " << (1 -
               double(total_svd_elems)/double(total_mat_elems)) * 100
            << std::endl;
}

#if 1
void percent_compression(Eri4 &mat, std::shared_ptr<e4pool> pool, double cut){

  long total_mat_elems = 0;
  long total_svd_elems = 0;
  double total_mat_diff = 0;
  long ntiles = 0;
  double svd_time0 = madness::wall_time();
  for(auto it = mat.get_pmap()->begin(); it != mat.get_pmap()->end(); ++it){
    mat.get_world().taskq.add(&eri4_task, it, mat, cut, &total_mat_elems,
                              &total_svd_elems, &total_mat_diff, pool);
    ++ntiles;
  }
  mat.get_world().gop.fence();
  double svd_time1 = madness::wall_time();

  std::cout << "\t\tNumber of tiles = " << ntiles << " in "
            << svd_time1 - svd_time0 << "s" << std::endl;
  std::cout << "\t\tTotal Matrix Elements = " << total_mat_elems << std::endl;
  std::cout << "\t\tTotal SVD Elements = " << total_svd_elems << std::endl;
  std::cout << "\t\tAverage tile diff = " << total_mat_diff/ntiles << std::endl;
  std::cout << "\t\tPercentage saved = " << (1 -
               double(total_svd_elems)/double(total_mat_elems)) * 100
            << std::endl;
}
#endif

int main(int argc, char** argv){
  sc::ExEnv::init(argc, argv);
  mpqc::MADNESSRuntime::initialize();

  Ref<World> world = new World();

  const char *input = "./benzene_trimer.kv";
  Ref<KeyVal> kv = new ParsedKeyVal(input);
  Ref<TiledBasisSet> tbs; tbs << kv->describedclassvalue("basis");
  Ref<TiledBasisSet> dftbs; dftbs << kv->describedclassvalue("dfbasis");

  tbs->print();
  dftbs->print();

  Ref<IntegralLibint2> ints_fac;
  ints_fac << kv->describedclassvalue("integrals");

  ints_fac->set_basis(tbs);

  // Local scope for one body ints
  {
    // Integral engine typedef

    std::array<TiledArray::TiledRange1, 2> ob_array =
                                        {{tbs->trange1(), tbs->trange1()}};
    TiledArray::TiledRange ob_trange(ob_array.begin(), ob_array.end());

    // S Range
    Matrix S(*world->madworld(), ob_trange);
    auto S_pool = std::make_shared<onebpool>(ints_fac->overlap()->clone());

    // Compute
    std::cout << "\nS\n";
    std::cout << "\t1e-4" << std::endl;
    percent_compression(S, S_pool, 1e-4);
    std::cout << "\t1e-5" << std::endl;
    percent_compression(S, S_pool, 1e-5);
    std::cout << "\t1e-6" << std::endl;
    percent_compression(S, S_pool, 1e-6);
    std::cout << "\t1e-7" << std::endl;
    percent_compression(S, S_pool, 1e-7);
    std::cout << "\n" << std::endl;

    world->madworld()->gop.fence();

    // H Range
    Matrix H(*world->madworld(), ob_trange);
    auto H_pool = std::make_shared<onebpool>(ints_fac->hcore()->clone());
    world->madworld()->gop.fence();


    std::cout << "\nH\n";
    std::cout << "\t1e-4" << std::endl;
    percent_compression(H, H_pool, 1e-4);
    std::cout << "\t1e-5" << std::endl;
    percent_compression(H, H_pool, 1e-5);
    std::cout << "\t1e-6" << std::endl;
    percent_compression(H, H_pool, 1e-6);
    std::cout << "\t1e-7" << std::endl;
    percent_compression(H, H_pool, 1e-7);
    std::cout << "\n" << std::endl;
  }

  // Density fitting integrals
  {
    // Set basis and grab a clone of the engine we need
    ints_fac->set_basis(tbs, tbs, dftbs);
    auto eri3_clone = ints_fac->electron_repulsion3()->clone();

    // Make an integral engine pool out of our clone
    auto eri3_ptr = std::make_shared<e3pool>(eri3_clone);

    std::array<TiledArray::TiledRange1,3> ec_array = {tbs->trange1(), tbs->trange1(),
                                                      dftbs->trange1()};

    TiledArray::TiledRange ec_range(ec_array.begin(), ec_array.end());

    // Using the df_ints as temporary storage for the twobody three center ints
    TiledArray::Array<double,3> df_ints_(*world->madworld(), ec_range);
    world->madworld()->gop.fence();

    std::cout << "\nEri3\n";
    std::cout << "\t1e-4" << std::endl;
    percent_compression(df_ints_, eri3_ptr, 1e-4);
    std::cout << "\t1e-5" << std::endl;
    percent_compression(df_ints_, eri3_ptr, 1e-5);
    std::cout << "\t1e-6" << std::endl;
    percent_compression(df_ints_, eri3_ptr, 1e-6);
    std::cout << "\t1e-7" << std::endl;
    percent_compression(df_ints_, eri3_ptr, 1e-7);
    std::cout << "\n" << std::endl;


    // Set basis and grab a clone of the engine we need
    ints_fac->set_basis(dftbs, dftbs);
    auto eri2_clone = ints_fac->electron_repulsion2()->clone();

    // Make an integral engine pool out of our clone
    auto eri2_ptr = std::make_shared<e2pool>(eri2_clone);

    std::array<TiledArray::TiledRange1,2> et_array = {dftbs->trange1(),
                                                      dftbs->trange1()};

    TiledArray::TiledRange et_range(et_array.begin(), et_array.end());

    // Using the df_ints as temporary storage for the twobody three center ints
    TiledArray::Array<double,2> two_ints(*world->madworld(), et_range);
    world->madworld()->gop.fence();

    std::cout << "\nEri2\n";
    std::cout << "\t1e-4" << std::endl;
    percent_compression(two_ints, eri2_ptr, 1e-4);
    std::cout << "\t1e-5" << std::endl;
    percent_compression(two_ints, eri2_ptr, 1e-5);
    std::cout << "\t1e-6" << std::endl;
    percent_compression(two_ints, eri2_ptr, 1e-6);
    std::cout << "\t1e-7" << std::endl;
    percent_compression(two_ints, eri2_ptr, 1e-7);
    std::cout << "\n" << std::endl;
  }

  // Four center ints
  {
    // Set basis and grab a clone of the engine we need
    ints_fac->set_basis(tbs);
    auto eri_clone = ints_fac->electron_repulsion()->clone();

    // Make an integral engine pool out of our clone
    auto eri_ptr = std::make_shared<e4pool>(eri_clone);

    std::array<TiledArray::TiledRange1, 4> e4_array = {tbs->trange1(), tbs->trange1(),
                                                       tbs->trange1(), tbs->trange1()};

    TiledArray::TiledRange e4_range(e4_array.begin(), e4_array.end());

    // Using the df_ints as temporary storage for the twobody three center ints
    TiledArray::Array<double,4> four_ints(*world->madworld(), e4_range);
    world->madworld()->gop.fence();

    std::cout << "\nFour Center\n";
    std::cout << "\t1e-4" << std::endl;
    percent_compression(four_ints, eri_ptr, 1e-4);
    std::cout << "\t1e-5" << std::endl;
    percent_compression(four_ints, eri_ptr, 1e-5);
    std::cout << "\t1e-6" << std::endl;
    percent_compression(four_ints, eri_ptr, 1e-6);
    std::cout << "\t1e-7" << std::endl;
    percent_compression(four_ints, eri_ptr, 1e-7);
    std::cout << "\n" << std::endl;

  }
  madness::finalize();
  return 0;
}

