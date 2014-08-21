#include "lr_tile.h"
#include <iostream>
#include <vector>
#include "../include/tiledarray.h"

int main(int argc, char **argv) {
  madness::World &world = madness::initialize(argc, argv);

  // Test that multiplication works
  Eigen::MatrixXd Z = Eigen::MatrixXd::Random(1000, 1000);
  Eigen::MatrixXd Q = Eigen::MatrixXd::Random(1000, 1000);

  Z = Z.transpose() * Z;
  Q = Q.transpose() * Q;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Z);
  Eigen::MatrixXd C = es.eigenvectors();
  Eigen::VectorXd V = es.eigenvalues();

  std::for_each(V.data(), V.data() + 990, [](double &x) { x = 0.0; });
  Z = C * V.asDiagonal() * C.transpose();

  es.compute(Q);
  C = es.eigenvectors();
  V = es.eigenvalues();

  std::for_each(V.data(), V.data() + 990, [](double &x) { x = 0.0; });
  Q = C * V.asDiagonal() * C.transpose();

  LRTile<double> B(Z);
  std::cout << "Rank of B = " << B.rank() << std::endl;
  LRTile<double> F(Q);
  std::cout << "Rank of F = " << F.rank() << std::endl;
  LRTile<double> D = B.mult(F);
  std::cout << "Does multiply work (1:yes,0:no)? "
            << D.matrixLR().isApprox(Z * Q) << "\n";

  LRTile<double> E = B.add(F);
  std::cout << "Does add work (1:yes,0:no)? " << E.matrixLR().isApprox(Z + Q)
            << "\n";

  std::cout << "Size of added LRTile = " << E.size()
            << " size of full matrix = " << E.matrixLR().size() << "\n";

  E.compress();
  std::cout << "Size of added LRTile after compression = " << E.size()
            << " size of full matrix = " << E.matrixLR().size() << "\n";

  std::cout << "Does compressed add work (1:yes,0:no)? "
            << E.matrixLR().isApprox(Z + Q) << "\n";

  Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(E.matrixLR());
  std::cout << "The QRP rank of E is " << qr.rank() << std::endl;

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

    TiledArray::Array<double, 2, LRTile<double>>::value_type tile(Q);

    *it = tile;
  }
  std::cout << "A = \n" << A << std::endl;

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
