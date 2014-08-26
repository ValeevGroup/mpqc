#include <iostream>
#include <vector>
#include <chrono>
#include <string>

#include "include/eigen.h"
#include "include/tbb.h"
#include "include/tiledarray.h"

#include "matrix/lr_tile.h"

Eigen::MatrixXd read_matrix(const std::string &filename);

template <typename T>
void purify(TiledArray::Array<double, 2, T> &D,
            const TiledArray::Array<double, 2, T> &S) {
  D("i,j") = D("i,k") * S("k,l") * D("l,j");
}

TiledArray::Array<double, 2, LRTile<double>>
make_lr_array(madness::World &,
              TiledArray::TiledRange &,
              const Eigen::MatrixXd &);

template<typename T, typename LR>
bool check_equal(const TiledArray::Array<double, 2, T> &Full,
                 const TiledArray::Array<double, 2, LR> &Low){
  auto fit = Full.begin();
  auto fend = Full.end();
  auto LRit = Low.begin();

  for(; fit != fend; ++fit, ++LRit){
    Eigen::MatrixXd LRmat = LRit->get().matrixLR();
    auto range = fit->get.range();
  }

  return false;
}

int main(int argc, char **argv) {
  std::string Sfile;
  std::string Dfile;
  std::string Ffile;
  unsigned long blocksize;
  if (argc < 5) {
    std::cout << "Please provide input files for S, D and F, and blocksize"
              << std::endl;
    return 0;
  } else {
    Sfile = argv[1];
    Dfile = argv[2];
    Ffile = argv[3];
    blocksize = std::stoul(argv[4]);
  }
  madness::World &world = madness::initialize(argc, argv);
  world.rank();

  Eigen::MatrixXd ES = read_matrix(Sfile);
  Eigen::MatrixXd ED = read_matrix(Dfile);
  Eigen::MatrixXd EF = read_matrix(Ffile);

  const auto nbasis = ES.rows();

  std::vector<unsigned int> blocking;
  auto i = 0;
  for (; i < nbasis; i += blocksize) {
    blocking.emplace_back(i);
  }
  blocking.emplace_back(nbasis);
  for (auto elem : blocking) {
    std::cout << elem << " ";
  }
  std::cout << "\n";


  std::vector<TiledArray::TiledRange1> blocking2(
      2, TiledArray::TiledRange1(blocking.begin(), blocking.end()));

  TiledArray::TiledRange trange(blocking2.begin(), blocking2.end());

  TiledArray::Array<double, 2> S = TiledArray::eigen_to_array
      <TiledArray::Array<double, 2>>(world, trange, ES);
  TiledArray::Array<double, 2> D = TiledArray::eigen_to_array
      <TiledArray::Array<double, 2>>(world, trange, ED);
  TiledArray::Array<double, 2> F = TiledArray::eigen_to_array
      <TiledArray::Array<double, 2>>(world, trange, EF);


  std::cout << "S = \n" << S << std::endl;
  std::cout << "F = \n" << F << std::endl;
  std::cout << "D = \n" << D << std::endl;

  purify(D, S);


  /************************************************
   * Perform LR version
   ************************************************/

  TiledArray::Array<double, 2, LRTile<double>> LR_S
      = make_lr_array(world, trange, ES);
  TiledArray::Array<double, 2, LRTile<double>> LR_D
      = make_lr_array(world, trange, ED);
  TiledArray::Array<double, 2, LRTile<double>> LR_F
      = make_lr_array(world, trange, EF);

  std::cout << "LR_S = \n" << LR_S << std::endl;
  std::cout << "LR_F = \n" << LR_F << std::endl;
  std::cout << "LR_D = \n" << LR_D << std::endl;

  purify(LR_D, LR_S);

  std::cout << "Are the matrices approximately equal? " << check_equal(D, LR_D) << std::endl;


  world.gop.fence();
  madness::finalize();
  return 0;
}

Eigen::MatrixXd read_matrix(const std::string &filename) {
  int cols = 0, rows = 0;

  // Read numbers from file into buffer.
  std::ifstream infile;
  infile.open(filename);
  if (!infile.eof()) {
    std::string line;
    std::getline(infile, line);
    std::stringstream stream(line);
    stream >> cols;
    stream >> rows;
  }

  Eigen::MatrixXd out_mat(rows, cols);

  int current_row = 0;

  while (!infile.eof()) {
    std::string line;
    std::getline(infile, line);

    std::stringstream stream(line);
    auto i = 0;
    for (; i < cols && !stream.eof(); ++i) {
      stream >> out_mat(current_row, i);
    }
    ++current_row;

    if (i != cols) {
      std::cout << "EIGEN MATRIX READ FAILED for cols" << std::endl;
    }
  }
  if (current_row != rows) {
    std::cout << "EIGEN MATRIX READ FAILED for rows" << std::endl;
  }

  infile.close();

  return out_mat;
}

TiledArray::Array<double, 2, LRTile<double>>
make_lr_array(madness::World &world,
              TiledArray::TiledRange &trange,
              const Eigen::MatrixXd &mat) {
  TiledArray::Array<double, 2, LRTile<double>> A(world, trange);
  auto i = A.begin();
  auto end = A.end();
  for (; i != end; ++i) {
    auto range = trange.make_tile_range(i.ordinal());
    Eigen::MatrixXd mat_block = mat.block(range.start()[0], range.start()[1],
                                          range.size()[0], range.size()[1]);

    TiledArray::Array
        <double, 2, LRTile<double>>::value_type tile(range, mat_block);

    *i = tile;
  }
  return A;
}
