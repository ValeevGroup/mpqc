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
#include <algorithm>

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
using bpmap_iter = decltype(std::declval<TiledArray::detail::BlockedPmap>().begin());

template<typename T>
void mat_task(const bpmap_iter first, const bpmap_iter last,
              TiledArray::TiledRange trange,
              double cut,
              long* t_mat_elems,
              long* t_svd_elems,
              double* t_mat_diff,
              T pool);

template<>
void mat_task(bpmap_iter first, bpmap_iter last, TiledArray::TiledRange trange, double cut,
              long* t_mat_elems,
              long* t_svd_elems, double* t_mat_diff,
              std::shared_ptr<onebpool> pool){

  using PoolPtrType = typename std::pointer_traits<std::shared_ptr<onebpool>>::element_type;
  typename PoolPtrType::engine_type engine =
      pool->instance();

  for(auto j = first; j != last; ++j){
    typename Matrix::value_type  tile = trange.make_tile_range(*j);
    get_integrals(tile, engine);

    std::size_t size0 = tile.range().size()[0];
    std::size_t size1 = tile.range().size()[1];

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Emap = TiledArray::eigen_map(tile, size0, size1);

    Eigen::MatrixXd Mat = Emap;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Mat);
    qr.setThreshold(cut);

    auto rank = qr.rank();
    Eigen::MatrixXd Q = Eigen::MatrixXd(qr.householderQ()).leftCols(rank);
    Eigen::MatrixXd R = qr.matrixR().topLeftCorner(rank, Mat.cols()).triangularView<Eigen::Upper>();

    Eigen::MatrixXd permuted_cols = qr.colsPermutation();

    Eigen::MatrixXd FinalVersion = Q * R * permuted_cols.transpose();

    *t_mat_elems += Mat.rows() * Mat.cols();
    *t_svd_elems += Q.size() + R.size();

    double mat_diff = Eigen::MatrixXd(Mat - FinalVersion).lpNorm<Eigen::Infinity>();
    *t_mat_diff += mat_diff;
  }
}

template<>
void mat_task(bpmap_iter first, bpmap_iter last, TiledArray::TiledRange trange, double cut,
              long* t_mat_elems,
              long* t_svd_elems, double* t_mat_diff,
              std::shared_ptr<e2pool> pool){

  using PoolPtrType = typename std::pointer_traits<std::shared_ptr<e2pool>>::element_type;
  typename PoolPtrType::engine_type engine =
      pool->instance();

  for(auto j = first; j != last; ++j){
    typename Matrix::value_type  tile = trange.make_tile_range(*j);
    get_integrals(tile, engine);

    std::size_t size0 = tile.range().size()[0];
    std::size_t size1 = tile.range().size()[1];

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Emap = TiledArray::eigen_map(tile, size0, size1);

    Eigen::MatrixXd Mat = Emap;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Mat);
    qr.setThreshold(cut);

    auto rank = qr.rank();
    Eigen::MatrixXd Q = Eigen::MatrixXd(qr.householderQ()).leftCols(rank);
    Eigen::MatrixXd R = qr.matrixR().topLeftCorner(rank, Mat.cols()).triangularView<Eigen::Upper>();

    Eigen::MatrixXd permuted_cols = qr.colsPermutation();

    Eigen::MatrixXd FinalVersion = Q * R * permuted_cols.transpose();

    *t_mat_elems += Mat.rows() * Mat.cols();
    *t_svd_elems += Q.size() + R.size();

    double mat_diff = Eigen::MatrixXd(Mat - FinalVersion).lpNorm<Eigen::Infinity>();
    *t_mat_diff += mat_diff;
  }
}

template <>
void mat_task(bpmap_iter first, bpmap_iter last, TiledArray::TiledRange trange, double cut,
              long* t_mat_elems,
              long* t_svd_elems, double* t_mat_diff,
              std::shared_ptr<e3pool> pool){

  using PoolPtrType = typename std::pointer_traits<std::shared_ptr<e3pool>>::element_type;
  typename PoolPtrType::engine_type engine =
      pool->instance();

  for(auto j = first; j != last; ++j){
    typename Matrix::value_type  tile = trange.make_tile_range(*j);
    get_integrals(tile, engine);

    std::size_t size0 = tile.range().size()[0];
    std::size_t size1 = tile.range().size()[1];
    std::size_t size2 = tile.range().size()[2];

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Emap = TiledArray::eigen_map(tile, size0*size1, size2);

    Eigen::MatrixXd Mat = Emap;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Mat);
    qr.setThreshold(cut);

    auto rank = qr.rank();
    Eigen::MatrixXd Q = Eigen::MatrixXd(qr.householderQ()).leftCols(rank);
    Eigen::MatrixXd R = qr.matrixR().topLeftCorner(rank, Mat.cols()).triangularView<Eigen::Upper>();

    Eigen::MatrixXd permuted_cols = qr.colsPermutation();

    Eigen::MatrixXd FinalVersion = Q * R * permuted_cols.transpose();

    *t_mat_elems += Mat.rows() * Mat.cols();
    *t_svd_elems += Q.size() + R.size();

    double mat_diff = Eigen::MatrixXd(Mat - FinalVersion).lpNorm<Eigen::Infinity>();
    *t_mat_diff += mat_diff;
  }
}

template<>
void mat_task(bpmap_iter first, bpmap_iter last, TiledArray::TiledRange trange, double cut,
              long* t_mat_elems,
              long* t_svd_elems,
              double* t_mat_diff,
              std::shared_ptr<e4pool> pool){

  using PoolPtrType = typename std::pointer_traits<std::shared_ptr<e4pool>>::element_type;
  typename PoolPtrType::engine_type engine =
      pool->instance();

  for(auto j = first; j != last; ++j){
    typename Matrix::value_type  tile = trange.make_tile_range(*j);
    get_integrals(tile, engine);

    std::size_t size0 = tile.range().size()[0];
    std::size_t size1 = tile.range().size()[1];
    std::size_t size2 = tile.range().size()[2];
    std::size_t size3 = tile.range().size()[3];

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Emap = TiledArray::eigen_map(tile, size0*size1, size2*size3);

    Eigen::MatrixXd Mat = Emap;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(Mat);
    qr.setThreshold(cut);

    auto rank = qr.rank();
    Eigen::MatrixXd Q = Eigen::MatrixXd(qr.householderQ()).leftCols(rank);
    Eigen::MatrixXd R = qr.matrixR().topLeftCorner(rank, Mat.cols()).triangularView<Eigen::Upper>();

    Eigen::MatrixXd permuted_cols = qr.colsPermutation();

    Eigen::MatrixXd FinalVersion = Q * R * permuted_cols.transpose();

    *t_mat_elems += Mat.rows() * Mat.cols();
    *t_svd_elems += Q.size() + R.size();

    double mat_diff = Eigen::MatrixXd(Mat - FinalVersion).lpNorm<Eigen::Infinity>();
    *t_mat_diff += mat_diff;
  }
}

template<typename T>
void percent_compression(TiledArray::TiledRange &trange, T pool, double cut, madness::World &world){

  long total_mat_elems = 0;
  long total_svd_elems = 0;
  double total_mat_diff = 0;
  long ntiles = trange.tiles().volume();
  double svd_time0 = madness::wall_time();

  TiledArray::detail::BlockedPmap bpmap(world, ntiles);

  auto first = bpmap.begin();
  auto end = bpmap.end();

  auto nthreads = madness::ThreadPool::size();
  // Attempt to divide the number of tasks into nthread^2 chuncks
  auto ideal_guess = bpmap.local_size()/(std::pow(nthreads,2));
  // Also attmpt to divide the number into nthread chuncks
  auto ok_guess = bpmap.local_size()/nthreads;

  // If the ideal block_size guess fails try the ok block_size guess
  auto block_size_guess = (ideal_guess == 0) ? ok_guess : ideal_guess;

  // But if block_size_guess went to 0 then set to at least 1;
  long block_size = std::max(std::size_t(block_size_guess), std::size_t(1));

  auto last = first;
  while(first != end){
    last = first;
    std::advance(last, std::min(block_size, std::distance(first, end)));

    world.taskq.add(&mat_task<T>, first, last, trange, cut, &total_mat_elems,
                                &total_svd_elems, &total_mat_diff, pool);
    first = last;
  }

  double svd_time1 = madness::wall_time();

  if(world.rank() == 0){
    std::cout << "\t\tTotal Time = " <<
              svd_time1 - svd_time0  << "s\n";
    std::cout << "\t\tTotal Matrix Elements = " << total_mat_elems << "\n";
    std::cout << "\t\tTotal SVD Elements = " << total_svd_elems <<  "\n";
    std::cout << "\t\tAverage tile diff = " << total_mat_diff/ntiles << "\n";
    std::cout << "\t\tPercentage saved = " << (1 -
                 double(total_svd_elems)/double(total_mat_elems)) * 100
              << "\n";
  }
}


int main(int argc, char** argv){
  sc::ExEnv::init(argc, argv);
  mpqc::MADNESSRuntime::initialize();

  int debug = (argc >= 2) ? atol(argv[1]) : 0;

  Ref<World> world = new World();

  const char *input = "./benzene_trimer.kv";
  Ref<KeyVal> kv = new ParsedKeyVal(input);
  Ref<TiledBasisSet> tbs; tbs << kv->describedclassvalue("basis");
  Ref<TiledBasisSet> dftbs; dftbs << kv->describedclassvalue("dfbasis");

  if(world->madworld()->rank() == 0){
    tbs->print();
    dftbs->print();
  }

  if(world->madworld()->rank() == 0 && debug != 0){
    volatile int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
      sleep(5);
  }

  Ref<IntegralLibint2> ints_fac;
  ints_fac << kv->describedclassvalue("integrals");

  ints_fac->set_basis(tbs);

  std::array<double, 6> cuts = {1e-4,1e-5, 1e-6, 1e-7, 1e-8, 1e-9};

  // Local scope for one body ints
  {
    // Integral engine typedef

    std::array<TiledArray::TiledRange1, 2> ob_array =
                                        {{tbs->trange1(), tbs->trange1()}};
    TiledArray::TiledRange ob_trange(ob_array.begin(), ob_array.end());

    // S Range
    //Matrix S(*world->madworld(), ob_trange);
    auto S_pool = std::make_shared<onebpool>(ints_fac->overlap()->clone());

    if(world->madworld()->rank() == 0){
      std::cout << "\nS\n";
    }
    for(auto cut : cuts){
      if(world->madworld()->rank() == 0)
        std::cout << "\n\t" << cut << std::endl;
      percent_compression(ob_trange, S_pool, cut, *world->madworld());
    }

#if 1
    // H Range
    //Matrix H(*world->madworld(), ob_trange);
    auto H_pool = std::make_shared<onebpool>(ints_fac->hcore()->clone());

    if(world->madworld()->rank() == 0){
      std::cout << "\nH\n";
    }
    for(auto cut : cuts){
      if(world->madworld()->rank() == 0)
        std::cout << "\n\t" << cut << std::endl;
      percent_compression(ob_trange, H_pool, cut, *world->madworld());
    }
#endif
  }

  // Density fitting integrals
#if 1
  if(true){
    // Set basis and grab a clone of the engine we need
    ints_fac->set_basis(tbs, tbs, dftbs);
    auto eri3_clone = ints_fac->electron_repulsion3()->clone();

    // Make an integral engine pool out of our clone
    auto eri3_ptr = std::make_shared<e3pool>(eri3_clone);

    std::array<TiledArray::TiledRange1,3> ec_array = {tbs->trange1(), tbs->trange1(),
                                                      dftbs->trange1()};

    TiledArray::TiledRange ec_range(ec_array.begin(), ec_array.end());

    if(world->madworld()->rank() == 0){
      std::cout << "\nEri3\n";
    }
    for(auto cut : cuts){
      if(world->madworld()->rank() == 0)
        std::cout << "\n\t" << cut << std::endl;
      percent_compression(ec_range, eri3_ptr, cut, *world->madworld());
    }

    // Set basis and grab a clone of the engine we need
    ints_fac->set_basis(dftbs, dftbs);
    auto eri2_clone = ints_fac->electron_repulsion2()->clone();

    // Make an integral engine pool out of our clone
    auto eri2_ptr = std::make_shared<e2pool>(eri2_clone);
    std::array<TiledArray::TiledRange1,2> et_array = {dftbs->trange1(),
                                                      dftbs->trange1()};
    TiledArray::TiledRange et_range(et_array.begin(), et_array.end());

    if(world->madworld()->rank() == 0){
      std::cout << "\nEri2\n";
    }
    for(auto cut : cuts){
      if(world->madworld()->rank() == 0)
        std::cout << "\n\t" << cut << std::endl;
      percent_compression(et_range, eri2_ptr, cut, *world->madworld());
    }
  }
#endif

#if 1
  // Four center ints
  if(true){
    // Set basis and grab a clone of the engine we need
    ints_fac->set_basis(tbs);
    auto eri_clone = ints_fac->electron_repulsion()->clone();

    // Make an integral engine pool out of our clone
    auto eri_ptr = std::make_shared<e4pool>(eri_clone);

    std::array<TiledArray::TiledRange1, 4> e4_array = {tbs->trange1(), tbs->trange1(),
                                                       tbs->trange1(), tbs->trange1()};

    TiledArray::TiledRange e4_range(e4_array.begin(), e4_array.end());


    if(world->madworld()->rank() == 0){
      std::cout << "\nEri4\n";
    }
    for(auto cut : cuts){
      if(world->madworld()->rank() == 0)
        std::cout << "\n\t" << cut << std::endl;
      percent_compression(e4_range, eri_ptr, cut, *world->madworld());
    }
  }
#endif
  mpqc::MADNESSRuntime::finalize();

  return 0;
}

