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
             const Eigen::MatrixXd &, TiledArray::Array<double, 2> &dense);

template <typename T, typename LR>
bool check_equal(
    const TiledArray::Array<double, 2, T, TiledArray::SparsePolicy> &Full,
    const TiledArray::Array<double, 2, T, TiledArray::DensePolicy> &Full_dense,
    const TiledArray::Array<double, 2, LR, TiledArray::SparsePolicy> &Low);


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
    TiledArray::SparseShape<float>::threshold(1e-15);
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

    // auto SS = TiledArray::eigen_to_array<TiledArray::Array<double, 2>>(
    //     world, trange, ES);
    // auto S = make_f_array(world, trange, ES, SS);
    auto LR_S = make_lr_array(world, trange, ES);

    // S("i,j") = S("i,k") * S("k,j");
    world.gop.fence();

#if 0
    SS("i,j") = SS("i,k") * SS("k,j");
    std::cout << "\ttile\t\test\t\treal_norm\n";
    for (auto i = SS.begin(); i != SS.end(); ++i) {
        auto tile = i->get();
        auto range = tile.range();
        auto norm_dense
            = Eigen::MatrixXd(TiledArray::eigen_map(tile, range.size().at(0),
                                                    range.size().at(1))).norm();

        auto shape = S.get_shape();
        auto gemmh = TiledArray::math::GemmHelper(
            madness::cblas::CBLAS_TRANSPOSE::NoTrans,
            madness::cblas::CBLAS_TRANSPOSE::NoTrans, 2, 2, 2);
        auto squared_shape = shape.gemm(shape, 1, gemmh);

        std::cout << "\t" << i.ordinal() << "\t\t"
                  << 96 * 96 * squared_shape[i.ordinal()] << "\t\t"
                  << norm_dense << std::endl;
    }
    world.gop.fence();
    S("i,j") = S("i,k") * S("k,j");
    for (auto i = S.begin(); i != S.end(); ++i) {
        auto tile = i->get();
        auto range = tile.range();
        auto norm_dense
            = Eigen::MatrixXd(TiledArray::eigen_map(tile, range.size().at(0),
                                                    range.size().at(1))).norm();

        std::cout << i.ordinal() << "\t" << norm_dense << std::endl;
    }
#endif
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

    Eigen::MatrixXd Lbasis;
    Eigen::MatrixXd Rbasis;
    Eigen::MatrixXd Tile1;
    Eigen::MatrixXd L1, R1;

    for (auto i = A.begin(); i != A.end(); ++i) {
        auto range = trange.make_tile_range(i.ordinal());
        decltype(A)::value_type tile{};

        Eigen::MatrixXd mat_block
            = mat.block(range.start()[0], range.start()[1], range.size()[0],
                        range.size()[1]);

        Eigen::MatrixXd L, R;

        auto tile_is_full_rank
            = algebra::Decompose_Matrix(mat_block, L, R, 1e-7);

        if (!tile_is_full_rank) {
            if (i.index()[0] == 0 && i.index()[1] == 3) {
                Tile1 = mat_block;
                L1 = L;
                R1 = R;
            }
            if (i.index()[0] == 0) {
                if (Lbasis.size() == 0) {
                    Lbasis = mat_block * mat_block.transpose();
                } else {
                    Lbasis += mat_block * mat_block.transpose();
                }
            }
            if (i.index()[1] == 3) {
                if (Rbasis.size() == 0) {
                    Rbasis = mat_block.transpose() * mat_block;
                } else {
                    Rbasis += mat_block.transpose() * mat_block;
                }
            }

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

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Lbasis);
    auto evals = es.eigenvalues();
    auto rankL = 0;
    for (auto i = 0; i < evals.size(); ++i) {
        evals[i] = evals[i] * evals[i];
        if (evals[i] >= 1e-14) {
            ++rankL;
        }
    }
    std::reverse(evals.data(), evals.data() + evals.size());
    std::cout << "Rank 0 = " << rankL << std::endl;

    Eigen::MatrixXd basisr0 = es.eigenvectors().rightCols(rankL);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> esR(Rbasis);
    auto evalsR = esR.eigenvalues();
    auto evalsRc = evalsR;
    auto rankR = 0;
    for (auto i = 0; i < evalsR.size(); ++i) {
        evalsR[i] = evalsR[i] * evalsR[i];
        if (evalsR[i] >= 1e-14) {
            ++rankR;
        } else {
            evalsRc[i] = 0;
        }
    }
    std::reverse(evalsR.data(), evalsR.data() + evalsR.size());
    std::cout << "Rank 1 = " << rankR << std::endl;

    Eigen::MatrixXd basisc1 = esR.eigenvectors().rightCols(rankR);
    std::cout << "Correct basis norm should be low = "
              << (Rbasis - basisc1
                  * Eigen::MatrixXd(evalsRc.asDiagonal())
                        .bottomRightCorner(rankR, rankR)
                  * basisc1.transpose()).norm() << std::endl;

    std::cout << "Dim L1 = " << L1.rows() << "x" << L1.cols() << std::endl;
    std::cout << "Dim R1 = " << R1.rows() << "x" << R1.cols() << std::endl;
    std::cout << "Diff between tile and low rank version = "
              << (Tile1 - L1 * R1).norm() << std::endl;

    Eigen::MatrixXd hversion = basisr0.transpose() * L1 * R1 * basisc1;
    std::cout << "Hversion dims = " << hversion.rows() << "x" << hversion.cols()
              << std::endl;
    std::cout << "Diff between tile and hversion = "
              << (Tile1 - basisr0 * hversion * basisc1.transpose()).norm()
              << std::endl;

    return A;
}

TiledArray::Array<double, 2, TiledArray::Tensor<double>,
                  TiledArray::SparsePolicy>
make_f_array(madness::World &world, TiledArray::TiledRange &trange,
             const Eigen::MatrixXd &mat, TiledArray::Array<double, 2> &dense) {
    // Shape tensor
    TiledArray::Tensor<float> shape_tensor(trange.tiles(), 0);
    for (auto i = 0ul; i != trange.tiles().volume(); ++i) {
        auto range = trange.make_tile_range(i);
        shape_tensor[i]
            = mat.block(range.start()[0], range.start()[1], range.size()[0],
                        range.size()[1]).lpNorm<2>();
    }
    auto gemmh = TiledArray::math::GemmHelper(
        madness::cblas::CBLAS_TRANSPOSE::NoTrans,
        madness::cblas::CBLAS_TRANSPOSE::NoTrans, 2, 2, 2);

    auto gemmed_shape_tensor = shape_tensor.add(shape_tensor, 1.0);

    TiledArray::SparseShape<float> shape(world, shape_tensor, trange);

    auto gemmed_norm_shape_tensor = shape.data().gemm(shape.data(), 1.0, gemmh);
    for (auto i = 0; i < gemmed_norm_shape_tensor.size(); ++i) {
        // std::cout << " tile " << i << " estimated squared norm = " <<
        // gemmed_norm_shape_tensor[i] << std::endl;
    }

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

    for (auto i = dense.begin(); i != dense.end(); ++i) {
        if (A.is_zero(i.ordinal())) {
            for (auto j = 0; j < i->get().range().volume(); ++j) {
                i->get()[j] = 0;
            }
        }
        auto range = i->get().range();
    }
    return A;
}

template <typename T, typename LR>
bool check_equal(
    const TiledArray::Array<double, 2, T, TiledArray::SparsePolicy> &Full,
    const TiledArray::Array<double, 2, T, TiledArray::DensePolicy> &Full_dense,
    const TiledArray::Array<double, 2, LR, TiledArray::SparsePolicy> &Low) {
    auto fit = Full.begin();
    auto fend = Full.end();
    auto LRit = Low.begin();

    bool same = true;
    bool empty_tile = false;

    std::cout << "\nTile\tnorm(zeroed dense)\tnorm(sparse)\tempty\n";
    for (auto dit = Full_dense.begin(); dit != Full_dense.end(); ++dit) {
        auto tile = Full_dense.find(dit.ordinal()).get();
        auto range = tile.range();
        auto norm_dense
            = Eigen::MatrixXd(TiledArray::eigen_map(tile, range.size().at(0),
                                                    range.size().at(1)))
                  .lpNorm<2>();
        bool empty = true;
        if (!Full.is_zero(dit.ordinal())) {
            empty = Full.find(dit.ordinal()).get().empty();
        }
        auto norm_sparse = Full.get_shape().data()[dit.ordinal()];
        std::cout << std::setprecision(16);
        std::cout << dit.ordinal() << "\t" << double(norm_dense) << "\t"
                  << norm_sparse * range.volume() << "\t" << empty << "\n";
    }
    fit = Full.begin();
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
                          << "\n\t\t2 norm of diff = "
                          << (Fmat - LRmat).lpNorm<2>() << std::endl;
                same = inner_same;
            }
        }
    }
    if (!empty_tile) {
        std::cout << "We tried to access 0 empty tiles . . .  ";
    }


    return same;
}
