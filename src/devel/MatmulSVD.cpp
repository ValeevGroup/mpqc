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
#include <iostream>

using namespace sc;
using namespace mpqc;
using namespace mpqc::TA;

using Matrix = TiledArray::Array<double, 2, TiledArray::Tensor<double> >;
using DistMatrix = elem::DistMatrix<double>;

void mat_task(typename Matrix::iterator it, double cut,
               long* t_mat_elems,
               long* t_svd_elems, double* t_mat_diff){
  typename Matrix::value_type tile = *it;

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

void eri3_task(TiledArray::Array<double,3>::iterator it, double cut,
               long* t_mat_elems,
               long* t_svd_elems, double* t_mat_diff){
  typename TiledArray::Array<double,3>::value_type tile = *it;

  std::size_t size0 = tile.range().size()[0];
  std::size_t size1 = tile.range().size()[1];
  std::size_t size2 = tile.range().size()[2];

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
  Emap = TiledArray::eigen_map(tile, size0 * size1, size2);

  Eigen::MatrixXd Mat = Emap;
  SVDTile SVDMat(Mat,cut);

  *t_mat_elems += Mat.rows() * Mat.cols();
  *t_svd_elems += SVDMat.nelements();

  double mat_diff = Eigen::MatrixXd(Mat -
     SVDMat.U_contract() * SVDMat.Vt()).lpNorm<Eigen::Infinity>();
  *t_mat_diff += mat_diff;
}

void percent_compression(Matrix &mat, double cut){

  long total_mat_elems = 0;
  long total_svd_elems = 0;
  double total_mat_diff = 0;
  long ntiles = 0;
  double svd_time0 = madness::wall_time();
  for(typename Matrix::iterator it = mat.begin(); it != mat.end(); ++it){
    mat.get_world().taskq.add(&mat_task, it, cut, &total_mat_elems,
                              &total_svd_elems, &total_mat_diff);
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

void percent_compression(TiledArray::Array<double ,3> &mat, double cut){

  long total_mat_elems = 0;
  long total_svd_elems = 0;
  double total_mat_diff = 0;
  long ntiles = 0;
  double svd_time0 = madness::wall_time();
  for(typename TiledArray::Array<double,3>::iterator it = mat.begin(); it != mat.end(); ++it){
    mat.get_world().taskq.add(&eri3_task, it, cut, &total_mat_elems,
                              &total_svd_elems, &total_mat_diff);
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

  {
  // Get pools
  using onebpool = IntegralEnginePool<Ref<OneBodyInt> >;
  auto S_pool = std::make_shared<onebpool>(ints_fac->overlap()->clone());
  auto H_pool = std::make_shared<onebpool>(ints_fac->hcore()->clone());

  // Compute ints
  Matrix S = Integrals(*world->madworld(), S_pool, tbs);
  world->madworld()->gop.fence();
  Matrix H = Integrals(*world->madworld(), H_pool, tbs);
  world->madworld()->gop.fence();

  std::cout << "\nS\n";
  std::cout << "\t1e-4" << std::endl;
  percent_compression(S, 1e-4);
  std::cout << "\t1e-5" << std::endl;
  percent_compression(S, 1e-5);
  std::cout << "\t1e-6" << std::endl;
  percent_compression(S, 1e-6);
  std::cout << "\t1e-7" << std::endl;
  percent_compression(S, 1e-7);
  std::cout << "\n" << std::endl;

  std::cout << "\nH\n";
  std::cout << "\t1e-4" << std::endl;
  percent_compression(H, 1e-4);
  std::cout << "\t1e-5" << std::endl;
  percent_compression(H, 1e-5);
  std::cout << "\t1e-6" << std::endl;
  percent_compression(H, 1e-6);
  std::cout << "\t1e-7" << std::endl;
  percent_compression(H, 1e-7);
  std::cout << "\n" << std::endl;
  }

  // Set basis and grab a clone of the engine we need
  ints_fac->set_basis(tbs, tbs, dftbs);
  auto eri3_clone = ints_fac->electron_repulsion3()->clone();

  // Make an integral engine pool out of our clone
  using eri3pool = IntegralEnginePool<sc::Ref<sc::TwoBodyThreeCenterInt> >;
  auto eri3_ptr = std::make_shared<eri3pool>(eri3_clone);

  // Using the df_ints as temporary storage for the twobody three center ints
  TiledArray::Array<double,3>
      df_ints_ =  Integrals(*world->madworld(), eri3_ptr, tbs, dftbs);
  world->madworld()->gop.fence();

  std::cout << "\nEri3\n";
  std::cout << "\t1e-4" << std::endl;
  percent_compression(df_ints_, 1e-4);
  std::cout << "\t1e-5" << std::endl;
  percent_compression(df_ints_, 1e-5);
  std::cout << "\t1e-6" << std::endl;
  percent_compression(df_ints_, 1e-6);
  std::cout << "\t1e-7" << std::endl;
  percent_compression(df_ints_, 1e-7);
  std::cout << "\n" << std::endl;

#if 0
  std::array<TiledArray::TiledRange1, 2>
        blocking{{tbs->trange1(), tbs->trange1()}};

  TiledArray::TiledRange trange(blocking.begin(), blocking.end());

  // Intialize Density matrix
  Matrix dens(*world->madworld(), trange);
  dens.set_all_local(0.0);
  world->madworld()->gop.fence();

  Ref<AssignedKeyVal> akv = new AssignedKeyVal();
  akv->assign("molecule", tbs->molecule().pointer());
  akv->assign("basis", tbs.pointer());
  akv->assign("dfbasis", dftbs.pointer());
  akv->assign("world", world.pointer());
  akv->assign("integrals", ints_fac.pointer());

  // Construct a density guess based on a SOAD object
  using Soad = sc::SuperpositionOfAtomicDensities;
  sc::Ref<Soad> guess = new Soad(sc::Ref<sc::KeyVal>(akv));


  ClDFGEngine GC(static_cast<Ref<KeyVal> >(akv));
  GC.set_densities({&dens});
  world->madworld()->gop.fence();

  // Make Fock matrix
  Matrix F = H("i,j") + GC("i,j");
  world->madworld()->gop.fence();

  auto occ = tbs->molecule()->total_Z();

  Matrix C = eigensolver_occ_Coeff(F,S,occ);
  dens("mu,nu") = C("mu,i") * C("nu,i");
  GC.set_coefficients({&C});

  std::cout << "\nD\n";
  std::cout << "\t1e-4" << std::endl;
  percent_compression(dens, 1e-4);
  std::cout << "\t1e-5" << std::endl;
  percent_compression(dens, 1e-5);
  std::cout << "\t1e-6" << std::endl;
  percent_compression(dens, 1e-6);
  std::cout << "\t1e-7" << std::endl;
  percent_compression(dens, 1e-7);
  std::cout << "\n" << std::endl;

  std::cout << "\nC\n";
  std::cout << "\t1e-4" << std::endl;
  percent_compression(C, 1e-4);
  std::cout << "\t1e-5" << std::endl;
  percent_compression(C, 1e-5);
  std::cout << "\t1e-6" << std::endl;
  percent_compression(C, 1e-6);
  std::cout << "\t1e-7" << std::endl;
  percent_compression(C, 1e-7);
  std::cout << "\n" << std::endl;


  std::cout << "\nF\n";
  std::cout << "\t1e-4" << std::endl;
  percent_compression(F, 1e-4);
  std::cout << "\t1e-5" << std::endl;
  percent_compression(F, 1e-5);
  std::cout << "\t1e-6" << std::endl;
  percent_compression(F, 1e-6);
  std::cout << "\t1e-7" << std::endl;
  percent_compression(F, 1e-7);
  std::cout << "\n" << std::endl;
#endif

  madness::finalize();
  return 0;
}

