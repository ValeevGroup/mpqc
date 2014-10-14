#include <iostream>
#include <vector>
#include <chrono>
#include <random>
#include <string>
#include <iomanip>

#include "../../include/eigen.h"
#include "../../include/tbb.h"
#include "../../include/tiledarray.h"

#include "../tile_pimpl.h"

Eigen::MatrixXd read_matrix(const std::string &filename);

TiledArray::Array<double, 2, TilePimpl<double>>
make_lr_array(madness::World &, TiledArray::TiledRange const &,
              Eigen::MatrixXd const &);

TiledArray::Array<double, 2> make_f_array(madness::World &,
                                          TiledArray::TiledRange const &,
                                          Eigen::MatrixXd const &);

template <typename Array>
Array create_noise(Array const &in, double epsilon) {
    std::default_random_engine generator;
    std::uniform_real_distribution<double> distribution(-epsilon, epsilon);

    Array out(in.get_world(), in.trange());

    auto in_it = in.begin();
    for (auto it = out.begin(); it != out.end(); ++it, ++in_it) {
        auto ctile = in_it->get();
        typename Array::value_type tile{ctile.range()};
        for (auto i = 0ul; i < tile.size(); ++i) {
            *(tile.data() + i) = *(ctile.data() + i) + distribution(generator);
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

template <typename Array>
void create_DiagDom_mat(Array &m) {
    auto counter = 10;
    for (auto it = m.begin(); it != m.end(); ++it) {
        typename Array::value_type tile{
            m.trange().make_tile_range(it.ordinal())};
        if (it.index()[0] == it.index()[1]) {
            for (auto i = 0ul; i < tile.size(); ++i) {
                *(tile.data() + i) = 100;
            }
        } else {
            for (auto i = 0ul; i < tile.size(); ++i) {
                *(tile.data() + i) = 1.0 / double(it.index()[0] * counter
                                                  + it.index()[1] * counter);
            }
            counter += 1;
        }
        *it = tile;
    }
}

template <typename T>
norm_diff_info compute_norms_full(TiledArray::Array<double, 2, T> const &Full,
                                  TiledArray::Array<double, 2, T> const &Low);

template <typename T, typename LR>
norm_diff_info compute_norms(TiledArray::Array<double, 2, T> const &Full,
                             TiledArray::Array<double, 2, LR> const &Low);

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

    // auto LR_S = make_f_array(world, trange, ES);

    // Hmat from chain of H's not important though.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ev(ES);
    auto evals = ev.eigenvalues();
    evals(0) = 20;
    for (auto i = 1; i < evals.size(); ++i) {
        evals(i) = 10.0 / double(i * i);
    }
    auto S = TiledArray::eigen_to_array<TiledArray::Array<double, 2>>(
        world, trange, ES);

    for (auto noise = 1e-16; noise < 1e-9; noise *= 10) {
        world.gop.fence();
        auto S_copy = create_noise(S, noise);
        std::cout << "Noise range per element is " << -noise << " - " << noise
                  << std::endl;

        for (auto i = 1; i < 10; ++i) {
            world.gop.fence();
            std::cout << "\tStarting Iteration " << i << std::endl;

            S("i,j") = S("i,k") * S("k,j");
            S_copy("i,j") = S_copy("i,k") * S_copy("k,j");
            auto diff = compute_norms_full(S, S_copy);
            std::cout << "\tF_norm of diff is " << diff.F_norm() << "\n"
                      << "\tavg for tiles is " << diff.avg() << "\n"
                      << "\tmax tile diff is " << diff.max() << "\n"
                      << "\tmin tile diff is " << diff.min() << std::endl;

            ES = TiledArray::array_to_eigen(S);
            Eigen::MatrixXd ESC = TiledArray::array_to_eigen(S_copy);

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ev(ES);
            auto min_ev_S = ev.eigenvalues().minCoeff();
            auto max_ev_S = ev.eigenvalues().maxCoeff();

            ev.compute(ESC);
            auto min_ev_SC = ev.eigenvalues().minCoeff();
            auto max_ev_SC = ev.eigenvalues().maxCoeff();

            //        LR_S("i,j") = LR_S("i,k") * LR_S("k,j");

            std::cout << "\tmax eigen value of S is " << max_ev_S << "\n"
                      << "\tmin eigen value of S is " << min_ev_S << "\n"
                      << "\tmax eigen value of S_copy is " << max_ev_SC << "\n"
                      << "\tmin eigen value of S_copy is " << min_ev_SC << "\n"
                      << std::endl;
        }
    }

    // S("i,j") = S("i,j") - S_copy("i,j");
    world.gop.fence();
    // std::cout << S << std::endl;

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

TiledArray::Array<double, 2> make_f_array(madness::World &world,
                                          TiledArray::TiledRange const &trange,
                                          Eigen::MatrixXd const &mat) {
    using TiledArray::eigen_map;

    TiledArray::Array<double, 2> A(world, trange);

    for (auto i = A.begin(); i != A.end(); ++i) {

        const auto range = trange.make_tile_range(i.ordinal());
        const auto &size = range.size();
        const auto &start = range.start();

        decltype(A)::value_type tile{range};

        auto mat_block = mat.block(start[0], start[1], size[0], size[1]);
        auto norm = mat_block.norm();

        Eigen::MatrixXd L, R;
        if (norm <= 1e-14) {
            auto tstart = tile.data();
            const auto end = tile.data() + tile.size();
            std::for_each(tstart, end, [](double &x) { x = 0; });
        } else if (!algebra::Decompose_Matrix(mat_block, L, R, 1e-14)) {
            eigen_map(tile, size[0], size[1]) = L * R;
        } else {
            eigen_map(tile, size[0], size[1]) = mat_block;
        }

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

        auto Fmat = eigen_map(fit->get(), size[0], size[1]);
        auto LRmat = eigen_map(LRit->get(), size[0], size[1]);

        diff.update((Fmat - LRmat).norm());
    }

    return diff;
}
