#include "lr_tile.h"
#include "../include/tiledarray.h"

int main(int argc, char** argv){
  madness::World &world = madness::initialize(argc, argv);

  // Test that multiplication works
  Eigen::MatrixXd Z = Eigen::MatrixXd::Random(1000,1000);
  Eigen::MatrixXd Q = Eigen::MatrixXd::Random(1000,1000);
  if(Z.isApprox(Q)){
    std::cout << "Z and Q were the same matrix\n";
  }
  LRTile<double> B(Z);
  LRTile<double> C(Q);
  LRTile<double> D = B.mult(C);
  std::cout << "Does multiply work (1:yes,0:no)? " << D.matrixLR().isApprox(Z*Q) << "\n";

  std::vector<unsigned int> blocking;
  for(auto i : {0,2,4}){
    blocking.push_back(i);
  }

  std::vector<TiledArray::TiledRange1> blocking2(2,
      TiledArray::TiledRange1(blocking.begin(), blocking.end()));

  TiledArray::TiledRange
    trange(blocking2.begin(), blocking2.end());

  TiledArray::Array<double, 2, LRTile<double> > A(world, trange);

  madness::finalize();
  return 0;
}
