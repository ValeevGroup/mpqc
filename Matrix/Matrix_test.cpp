#include "lr_tile.h"
#include "../include/tiledarray.h"

int main(int argc, char** argv){
  madness::World &world = madness::initialize(argc, argv);

  Eigen::MatrixXd Z = Eigen::MatrixXd::Random(10,10);
  LRTile<double> B(Z);

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
