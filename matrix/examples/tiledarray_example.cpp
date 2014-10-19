#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <string>
#include <iomanip>

#include "../../include/eigen.h"
#include "../../include/tiledarray.h"

#include "../tile_pimpl.h"

Eigen::MatrixXd read_matrix(const std::string &filename);

TiledArray::Array<double, 2, TilePimpl<double>>
make_lr_array(madness::World &, TiledArray::TiledRange const &,
              Eigen::MatrixXd const &);

TiledArray::Array<double, 2, TilePimpl<double>> make_f_array(madness::World &,
                                          TiledArray::TiledRange const &);

template <typename Array>
Array copy_to_lr(Array const &in, int rank) {
    Array out(in.get_world(), in.trange());
    using T = double;

    auto in_it = in.begin();
    for (auto it = out.begin(); it != out.end(); ++it, ++in_it) {
        auto ctile = in_it->get();
        typename Array::value_type tile{ctile.range()};
        if (it.index()[0] != it.index()[1]) {
            auto mat = ctile.tile().matrix();
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU
                                                       | Eigen::ComputeThinV);
            Eigen::MatrixXd diag = svd.singularValues().asDiagonal();
            Eigen::MatrixXd L = svd.matrixU().leftCols(rank)
                                * diag.block(0, 0, rank, rank);
            Eigen::MatrixXd R = svd.matrixV().transpose().topRows(rank);
            tile = TilePimpl<T>{ctile.range(), TileVariant<T>{LowRankTile<T>{
                                                   std::move(L), std::move(R)}},
                                1e-07};
        } else {
            tile = ctile.clone();
        }
        *it = tile;
    }

    return out;
}

struct norm_diff_info {
    double sum_squares_ = 0;
    double sum_ = 0;
    double min_ = 1;
    double max_ = 0;
    long ntiles_ = 0;

    void update(double norm) {
        max_ = (norm > max_) ? norm : max_;
        min_ = (norm < min_) ? norm : min_;
        ++ntiles_;
        sum_ += norm;
        sum_squares_ += norm * norm;
    }

    double max() { return max_; }
    double min() { return min_; }
    double avg() { return sum_ / double(ntiles_); }
    double F_norm() { return std::sqrt(sum_squares_); }
};

template <typename T>
norm_diff_info compute_norms_full(TiledArray::Array<double, 2, T> const &Full,
                                  TiledArray::Array<double, 2, T> const &Low);

template <typename T, typename LR>
norm_diff_info compute_norms(TiledArray::Array<double, 2, T> const &Full,
                             TiledArray::Array<double, 2, LR> const &Low);

int main(int argc, char **argv) {
    unsigned long matsize;
    unsigned long blocksize;
    if (argc < 3) {
        std::cout << "Please provide an input matrix size for the array and a "
                     "blocksize\n";
        return 0;
    } else {
        matsize = std::stoul(argv[1]);
        blocksize = std::stoul(argv[2]);
    }

    madness::World &world = madness::initialize(argc, argv);


    std::vector<unsigned int> blocking;
    auto i = 0;
    for (; i < matsize; i += blocksize) {
        blocking.emplace_back(i);
    }
    blocking.emplace_back(matsize);

    std::vector<TiledArray::TiledRange1> blocking2(
        2, TiledArray::TiledRange1(blocking.begin(), blocking.end()));

    TiledArray::TiledRange trange(blocking2.begin(), blocking2.end());

    auto S = make_f_array(world, trange);
    auto Sc = copy_to_lr(S, blocksize/5);
    TiledArray::Array<double, 2, TilePimpl<double>> C(world,trange);
    TiledArray::Array<double, 2, TilePimpl<double>> Cc(world,trange);

    world.gop.fence();
    double full_time = madness::wall_time();
    C("i,j") = S("i,k") * S("k,j");
    world.gop.fence();
    full_time = madness::wall_time() - full_time;
    
    double lr_time = madness::wall_time();
    Cc("i,j") = Sc("i,k") * Sc("k,j");
    lr_time = madness::wall_time() - lr_time;

    std::cout << "Full time took " << full_time << " lr time took " << lr_time << std::endl;
    auto diff = compute_norms_full(S, Sc);
    std::cout << "\tF_norm of diff is " << diff.F_norm() << "\n"
              << "\tavg for tiles is " << diff.avg() << "\n"
              << "\tmax tile diff is " << diff.max() << "\n"
              << "\tmin tile diff is " << diff.min() << std::endl;


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
make_lr_array(madness::World &world, TiledArray::TiledRange const &trange,
              Eigen::MatrixXd const &mat) {
    TiledArray::Array<double, 2, TilePimpl<double>> A(world, trange);
    for (auto i = A.begin(); i != A.end(); ++i) {
        auto range = trange.make_tile_range(i.ordinal());
        decltype(A)::value_type tile{};

        Eigen::MatrixXd mat_block
            = mat.block(range.start()[0], range.start()[1], range.size()[0],
                        range.size()[1]);

        auto norm = mat_block.norm();


        auto tile_is_full_rank = true; // Assume Full Rank
        Eigen::MatrixXd L, R;
        if (norm >= 1e-16) { // If large norm check rank
            tile_is_full_rank
                = algebra::Decompose_Matrix(mat_block, L, R, 1e-16);
        } else { // If small norm not full rank
            tile_is_full_rank = false;
        }

        if (!tile_is_full_rank) {
            if (norm >= 1e-16) {
                tile = TilePimpl<double>{
                    std::move(range), TileVariant<double>{LowRankTile<double>{
                                          std::move(L), std::move(R)}},
                    1e-16};
            } else {
                tile = TilePimpl<double>{
                    std::move(range),
                    TileVariant<double>{LowRankTile<double>{true}}, 1e-16};
            }
        } else {
            tile = TilePimpl<double>{
                std::move(range),
                TileVariant<double>{FullRankTile<double>{mat_block}}, 1e-16};
        }
        *i = std::move(tile);
    }

    return A;
}

TiledArray::Array<double, 2, TilePimpl<double>> make_f_array(madness::World &world,
                                          TiledArray::TiledRange const &trange){
    using TiledArray::eigen_map;

    TiledArray::Array<double, 2, TilePimpl<double>> A(world, trange);

    for (auto i = A.begin(); i != A.end(); ++i) {

        const auto range = trange.make_tile_range(i.ordinal());
        const auto &size = range.size();
        const auto &start = range.start();

        decltype(A)::value_type tile{range};
        tile = TilePimpl<double>{
            range, TileVariant<double>{FullRankTile<double>{
                       Eigen::MatrixXd::Random(size[0], size[1])}}};

        *i = std::move(tile);
    }

    return A;
}

template <typename T, typename LR>
norm_diff_info compute_norms(TiledArray::Array<double, 2, T> const &Full,
                             TiledArray::Array<double, 2, LR> const &Low) {

    using TiledArray::eigen_map;

    norm_diff_info diff;

    const auto fend = Full.end();
    auto LRit = Low.begin();
    for (auto fit = Full.begin(); fit != fend; ++fit, ++LRit) {
        auto size = fit->get().range().size();

        auto Fmat = eigen_map(fit->get(), size[0], size[1]);
        auto LRmat = LRit->get().tile().matrix();

        if (LRit->get().tile().iszero()) {
            diff.update(Fmat.norm());
        } else {
            diff.update((Fmat - LRmat).norm());
        }
    }

    return diff;
}

template <typename T>
norm_diff_info compute_norms_full(TiledArray::Array<double, 2, T> const &Full,
                                  TiledArray::Array<double, 2, T> const &Low) {

    using TiledArray::eigen_map;
    norm_diff_info diff;

    const auto fend = Full.end();
    auto LRit = Low.begin();
    for (auto fit = Full.begin(); fit != fend; ++fit, ++LRit) {
        const auto &size = fit->get().range().size();

        auto Fmat = fit->get().tile().matrix();
        auto LRmat = LRit->get().tile().matrix();

        diff.update((Fmat - LRmat).norm());
    }

    return diff;
}
