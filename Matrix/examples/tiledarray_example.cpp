#include <iostream>
#include <vector>
#include <chrono>
#include <string>

#include "../../include/eigen.h"
#include "../../include/tbb.h"
#include "../../include/tiledarray.h"

#include "../tile_pimpl.h"

Eigen::MatrixXd read_matrix(const std::string &filename);

TiledArray::Array<double, 2, TilePimpl<double>>
make_lr_array(madness::World &, TiledArray::TiledRange &,
              const Eigen::MatrixXd &);

template <typename T, typename LR>
bool check_equal(const TiledArray::Array<double, 2, T> &Full,
                 const TiledArray::Array<double, 2, LR> &Low);


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
    madness::World &world = madness::initialize(argc, argv);

    Eigen::MatrixXd ES = read_matrix(Sfile);

    const auto nbasis = ES.rows();

    std::vector<unsigned int> blocking;
    auto i = 0;
    for (; i < nbasis; i += blocksize) {
        blocking.emplace_back(i);
    }

    blocking.emplace_back(nbasis);

    std::vector<TiledArray::TiledRange1> blocking2(
        2, TiledArray::TiledRange1(blocking.begin(), blocking.end()));

    TiledArray::TiledRange trange(blocking2.begin(), blocking2.end());

    TiledArray::Array<double, 2> S
        = TiledArray::eigen_to_array<TiledArray::Array<double, 2>>(world,
                                                                   trange, ES);

    TiledArray::Array<double, 2, TilePimpl<double>> LR_S
        = make_lr_array(world, trange, ES);


    std::cout << "\nChecking arrays for approximate equality. . . . ";
    bool passed_check = check_equal(S, LR_S);
    if (!passed_check) {
        std::cout << "Arrays were not equal!";
    } else {
        std::cout << "Ok!";
    }
    std::cout << "\n\n";

    world.gop.fence();
    auto lr_time = madness::wall_time();
    LR_S("i,j") = LR_S("i,k") * LR_S("k,j");

    world.gop.fence();
    lr_time = madness::wall_time() - lr_time;
    std::cout << "LR time trace decrease was " << lr_time << " s\n";

    auto full_time = madness::wall_time();
    S("i,j") = S("i,k") * S("k,j");
    world.gop.fence();
    full_time = madness::wall_time() - full_time;
    std::cout << "full time trace decrease was " << full_time << " s\n";

    std::cout << "Checking arrays for approximate equality after trace "
                 "decrease like operations . . . . ";
    passed_check = check_equal(S, LR_S);
    if (!passed_check) {
        std::cout << "Arrays were not equal!";
    } else {
        std::cout << "Ok!";
    }
    std::cout << "\n\n";

    world.gop.fence();
    lr_time = madness::wall_time();
    LR_S("i,j") = 2 * LR_S("i,j") - LR_S("i,k") * LR_S("k,j");

    world.gop.fence();
    lr_time = madness::wall_time() - lr_time;
    std::cout << "LR time trace increase was " << lr_time << " s\n";

    full_time = madness::wall_time();
    S("i,j") = 2 * S("i,j") - S("i,k") * S("k,j");
    world.gop.fence();
    full_time = madness::wall_time() - full_time;
    std::cout << "full time trace increase was " << full_time << " s\n";

    std::cout << "Checking arrays for approximate equality after trace "
                 "increase like operation . . . . ";
    passed_check = check_equal(S, LR_S);
    if (!passed_check) {
        std::cout << "Arrays were not equal!";
    } else {
        std::cout << "Ok!";
    }
    std::cout << "\n\n";

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

TiledArray::Array<double, 2, TilePimpl<double>>
make_lr_array(madness::World &world, TiledArray::TiledRange &trange,
              const Eigen::MatrixXd &mat) {
    TiledArray::Array<double, 2, TilePimpl<double>> A(world, trange);
    auto i = A.begin();
    auto end = A.end();
    for (; i != end; ++i) {
        auto range = trange.make_tile_range(i.ordinal());
        Eigen::MatrixXd mat_block
            = mat.block(range.start()[0], range.start()[1], range.size()[0],
                        range.size()[1]);

        TiledArray::Array<double, 2, TilePimpl<double>>::value_type tile{};

        Eigen::MatrixXd L, R;

        auto tile_is_full_rank
            = algebra::Decompose_Matrix(mat_block, L, R, 1e-07);

        auto norm = mat_block.lpNorm<2>();
        if (!tile_is_full_rank && norm > 1e-15) {
            tile = TilePimpl<double>{std::move(range),
                                     TileVariant<double>{LowRankTile<double>{
                                         std::move(L), std::move(R)}},
                                     1e-07};
        } else {
            tile = TilePimpl<double>{
                std::move(range),
                TileVariant<double>{FullRankTile<double>{mat_block}}, 1e-07};
        }

        *i = tile;
    }
    return A;
}

template <typename T, typename LR>
bool check_equal(const TiledArray::Array<double, 2, T> &Full,
                 const TiledArray::Array<double, 2, LR> &Low) {
    auto fit = Full.begin();
    auto fend = Full.end();
    auto LRit = Low.begin();

    bool same = true;

    for (; fit != fend; ++fit, ++LRit) {
        Eigen::MatrixXd LRmat = LRit->get().tile().matrix();
        auto range = fit->get().range();
        Eigen::MatrixXd Fmat = TiledArray::eigen_map(
            fit->get(), range.size()[0], range.size()[1]);
        if (LRmat.size() == 0) {
            double norm = Fmat.lpNorm<2>();
            std::cout << "\n\tLRmat for tile(" << fit.index()[0] << ","
                      << fit.index()[1]
                      << ") was empty Fmat has norm = " << norm << " ";
            if (norm > 1e-16) {
                same = false;
            }
        } else {
            auto inner_same = ((Fmat - LRmat).lpNorm<2>() < 1e-06);
            if (inner_same == false) {
                std::cout << "\n\tTile = (" << fit.index()[0] << ","
                          << fit.index()[1] << ")"
                          << "\n\t\t2 norm of diff = "
                          << (Fmat - LRmat).lpNorm<2>() << std::endl;
                same = inner_same;
            }
        }
    }

    return same;
}
