#include "lr_tile.h"
#include <iostream>
#include <vector>
#include "../include/tiledarray.h"

int main(int argc, char **argv) {
  madness::World &world = madness::initialize(argc, argv);

  std::vector<unsigned int> blocking;
  for (auto i : {0, 5, 10}) {
    blocking.push_back(i);
  }

  std::vector<TiledArray::TiledRange1> blocking2(
      2, TiledArray::TiledRange1(blocking.begin(), blocking.end()));

  TiledArray::TiledRange trange(blocking2.begin(), blocking2.end());

  TiledArray::Array<double, 2, LRTile<double>> A(world, trange);

  auto it = A.begin();
  auto end = A.end();

  for (; it != end; ++it) {
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(5, 5);
    Q = Q.transpose() * Q;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Q);
    Eigen::MatrixXd C = es.eigenvectors();
    Eigen::VectorXd V = es.eigenvalues();

    std::for_each(V.data(), V.data() + 4, [](double &x) { x = 0.0; });
    Q = C * V.asDiagonal() * C.transpose();

    auto tile_trange = A.trange().make_tile_range(it.ordinal());
    using tile_value_type = TiledArray::Array
        <double, 2, LRTile<double>>::value_type;
    tile_value_type tile(std::move(tile_trange), Q);

    *it = tile;
  }

  TiledArray::Array<double, 2, LRTile<double>> TB;
  TB("i,j") = A("i,j") + A("i,j");
  world.gop.fence();
  std::cout << "TB = \n" << TB << std::endl;

  TB("i,j") = A("i,k") * A("k,j");
  world.gop.fence();
  std::cout << "TB = \n" << TB << std::endl;

  madness::finalize();
  return 0;
}
