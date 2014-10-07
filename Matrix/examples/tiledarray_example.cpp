#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <iomanip>

#include "../../include/eigen.h"
#include "../../include/tbb.h"
#include "../../include/tiledarray.h"

#include "../tile_pimpl.h"

Eigen::MatrixXd read_matrix(const std::string &filename);

TiledArray::Array<double, 2, TilePimpl<double>, TiledArray::SparsePolicy>
make_lr_array(madness::World &, TiledArray::TiledRange &,
              const Eigen::MatrixXd &);

TiledArray::Array<double, 2, TiledArray::Tensor<double>,
                  TiledArray::SparsePolicy>
make_f_array(madness::World &, TiledArray::TiledRange &,
             const Eigen::MatrixXd &);

void compress(TiledArray::Array<double, 2, TilePimpl<double>,
                                TiledArray::SparsePolicy> &A) {
  for(auto i = A.begin(); i != A.end(); ++i){
     i->get().compress();
  }
}

template <typename T, typename LR>
bool check_equal(
    const TiledArray::Array<double, 2, T, TiledArray::SparsePolicy> &Full,
    const TiledArray::Array<double, 2, LR, TiledArray::SparsePolicy> &Low);

TiledArray::TiledRange create_trange(int nbasis, int blocksize) {
    std::vector<unsigned int> blocking;
    for (auto i = 0; i < nbasis; i += blocksize) {
        blocking.emplace_back(i);
    }

    blocking.emplace_back(nbasis);

    std::vector<TiledArray::TiledRange1> blocking2(
        2, TiledArray::TiledRange1(blocking.begin(), blocking.end()));

    return TiledArray::TiledRange{blocking2.begin(), blocking2.end()};
}

int main(int argc, char **argv) {
    std::string Sfile;
    unsigned long blocksize;
    if (argc < 3) {
        std::cout << "Please provide an input files for the array and a "
                     "blocksize\n";
        return 0;
    } else {
        Sfile = argv[1];
        blocksize = std::stoul(argv[2]);
    }

    TiledArray::SparseShape<float>::threshold(1e-10);
    madness::World &world = madness::initialize(argc, argv);
    Eigen::MatrixXd ES = read_matrix(Sfile);

    auto trange = create_trange(ES.rows(), blocksize);
    auto LR_S = make_lr_array(world, trange, ES);
    auto S = make_f_array(world, trange, ES);
   // if (check_equal(S, LR_S)) {
   // i     std::cout << "Tiles started off equal" << std::endl;
   // i }
    //compress(LR_S);
    for (auto i = 0; i < 10; ++i) {
        std::cout << "Iteration " << i << std::endl;
        S("i,j") = S("i,k") * S("k,j");
        //LR_S("i,j") = LR_S("i,k") * LR_S("k,j");
        //check_equal(S, LR_S);
        //compress(LR_S);
    }

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

TiledArray::Array<double, 2, TilePimpl<double>, TiledArray::SparsePolicy>
make_lr_array(madness::World &world, TiledArray::TiledRange &trange,
              const Eigen::MatrixXd &mat) {
    // Shape tensor
    TiledArray::Tensor<float> shape_tensor(trange.tiles(), 0);
    for (auto i = 0ul; i != trange.tiles().volume(); ++i) {
        auto range = trange.make_tile_range(i);
        shape_tensor[i]
            = mat.block(range.start()[0], range.start()[1], range.size()[0],
                        range.size()[1]).lpNorm<2>();
    }

    TiledArray::SparseShape<float> shape(world, shape_tensor, trange);

    TiledArray::Array<double, 2, TilePimpl<double>, TiledArray::SparsePolicy> A(
        world, trange, shape);

    for (auto i = A.begin(); i != A.end(); ++i) {
        auto range = trange.make_tile_range(i.ordinal());
        decltype(A)::value_type tile{};

        Eigen::MatrixXd mat_block
            = mat.block(range.start()[0], range.start()[1], range.size()[0],
                        range.size()[1]);

        Eigen::MatrixXd L, R;

        auto tile_is_full_rank
            = algebra::Decompose_Matrix(mat_block, L, R, 1e-07);

        if (!tile_is_full_rank) {
            tile = TilePimpl<double>{std::move(range),
                                     TileVariant<double>{LowRankTile<double>{
                                         std::move(L), std::move(R)}},
                                     1e-07};
        } else {
            tile = TilePimpl<double>{
                std::move(range),
                TileVariant<double>{FullRankTile<double>{mat_block}}, 1e-07};
        }

        *i = std::move(tile);
    }

    return A;
}

TiledArray::Array<double, 2, TiledArray::Tensor<double>,
                  TiledArray::SparsePolicy>
make_f_array(madness::World &world, TiledArray::TiledRange &trange,
             Eigen::MatrixXd const &mat) {
    // Shape tensor
    TiledArray::Tensor<float> shape_tensor(trange.tiles(), 0);
    for (auto i = 0ul; i != trange.tiles().volume(); ++i) {
        auto range = trange.make_tile_range(i);
        shape_tensor[i]
            = mat.block(range.start()[0], range.start()[1], range.size()[0],
                        range.size()[1]).lpNorm<2>();
    }

    TiledArray::SparseShape<float> shape(world, shape_tensor, trange);

    TiledArray::Array<double, 2, TiledArray::Tensor<double>,
                      TiledArray::SparsePolicy> A(world, trange, shape);

    for (auto i = A.begin(); i != A.end(); ++i) {
        auto range = trange.make_tile_range(i.ordinal());
        decltype(A)::value_type tile{range};

        TiledArray::eigen_map(tile, range.size()[0], range.size()[1])
            = mat.block(range.start()[0], range.start()[1], range.size()[0],
                        range.size()[1]);
        *i = std::move(tile);
    }
    return A;
}

template <typename T, typename LR>
bool check_equal(
    const TiledArray::Array<double, 2, T, TiledArray::SparsePolicy> &Full,
    const TiledArray::Array<double, 2, LR, TiledArray::SparsePolicy> &Low) {
    auto fit = Full.begin();
    auto fend = Full.end();
    auto LRit = Low.begin();

    bool same = true;
    for (; fit != fend; ++fit, ++LRit) {
        if (fit->get().empty()) {
        } else {
            auto range = fit->get().range();
            Eigen::MatrixXd Fmat = TiledArray::eigen_map(
                fit->get(), range.size().at(0), range.size().at(1));
            Eigen::MatrixXd LRmat = LRit->get().tile().matrix();
            auto inner_same = ((Fmat - LRmat).lpNorm<2>() < 1e-06);
            if (inner_same == false) {
                std::cout << "\n\tTile = (" << fit.index()[0] << ","
                          << fit.index()[1] << ")"
                          << " norm of diff = " << (Fmat - LRmat).lpNorm<2>()
                          << std::endl;
                same = inner_same;
            }
        }
    }

    return same;
}
