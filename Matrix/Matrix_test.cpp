#include "lr_tile.h"
#include <iostream>
#include <vector>
#include "../include/tiledarray.h"

int main(int argc, char** argv) {
  madness::World& world = madness::initialize(argc, argv);

  // Test that multiplication works
  Eigen::MatrixXd Z = Eigen::MatrixXd::Random(1000, 1000);
  Eigen::MatrixXd Q = Eigen::MatrixXd::Random(1000, 1000);

  LRTile<double> B(Z);
  LRTile<double> C(Q);
  LRTile<double> D = B.mult(C);
  std::cout << "Does multiply work (1:yes,0:no)? "
            << D.matrixLR().isApprox(Z * Q) << "\n";

  LRTile<double> E = B.add(C);
  std::cout << "Does add work (1:yes,0:no)? " << E.matrixLR().isApprox(Z + Q)
            << "\n";

  std::cout << "Size of added LRTile = " << E.size()
            << " size of full matrix = " << E.matrixLR().size() << "\n";

  E.compress();
  std::cout << "Size of added LRTile after compression = " << E.size()
            << " size of full matrix = " << E.matrixLR().size() << "\n";

  std::vector<unsigned int> blocking;
  for (auto i : {0, 2, 4}) {
    blocking.push_back(i);
  }

  std::vector<TiledArray::TiledRange1> blocking2(
      2, TiledArray::TiledRange1(blocking.begin(), blocking.end()));

  TiledArray::TiledRange trange(blocking2.begin(), blocking2.end());

  TiledArray::Array<double, 2, LRTile<double>> A(world, trange);

  madness::finalize();
  return 0;
}
