#include <iostream>
#include <vector>
#include <chrono>
#include <string>

#include "include/eigen.h"
#include "include/tbb.h"
#include "include/tiledarray.h"

#include "matrix/lr_tile.h"

Eigen::MatrixXd read_matrix(const std::string &filename);

template<typename T>
void check_data_sparsity(const TiledArray::Array<double, 2, T> &A){
    TiledArray::Array<double, 2, LRTile<double>> B(A.get_world(), A.trange());

    auto it_A = A.begin();
    auto it_B = B.begin();
    auto end = B.end();
    for(; it_B != end; ++it_B, ++it_A) {
        TiledArray::Array<double, 2, LRTile<double>>::value_type tile(
            it_A->get().range(), TiledArray::eigen_map(it_A->get()), true, 1e-03);
        *it_B = tile;
    }

    int total_tiles = 0;
    int full_tiles = 0;
    for (const auto it : B) {
        full_tiles += int(it.get().is_full());
        ++total_tiles;
    }
    std::cout << "Percentage of full tiles in Mat = " << double(full_tiles)
                                                  / double(total_tiles) << "\n";
}

template <typename T>
double compute_trace(const TiledArray::Array<double, 2, T> &A) {
    double trace = 0.0;
    for (auto it = A.begin(); it != A.end(); ++it) {
        auto pos = it.index();
        if (pos[0] == pos[1]) {
            const auto range = it->get().range();
            trace += TiledArray::eigen_map(it->get(), range.size()[0],
                                           range.size()[1]).trace();
        }
    }
    return trace;
}

template <>
double compute_trace(const TiledArray::Array<double, 2, LRTile<double>> &A) {
    double trace = 0.0;
    for (auto it = A.begin(); it != A.end(); ++it) {
        auto pos = it.index();
        if (pos[0] == pos[1]) {
            trace += it->get().matrixLR().trace();
        }
    }
    return trace;
}

template <typename T>
void purify(TiledArray::Array<double, 2, T> &D,
            const TiledArray::Array<double, 2, T> &S) {

    TiledArray::Array<double, 2, T> DS(D.get_world(), D.trange());
    DS("i,j") = D("i,k") * S("k,j");
    double trace = compute_trace(DS);
    TiledArray::Array<double, 2, T> D2(D.get_world(), D.trange());

    // TODO add occ
    while (std::abs(trace - 2.0) > 1e-05) {
        // Compute D^2 just capture expression
        auto D2 = DS("i,k") * D("k,j");

        if (trace < 2) { // Raise trace
            D("i,j") = 2 * D("i,j") - D2;
        } else { // Lower trace
            D("i,j") = D2;
        }

        DS("i,j") = D("i,k") * S("k,j");

        trace = compute_trace(DS);
    }
}

TiledArray::Array<double, 2, LRTile<double>>
make_lr_array(madness::World &, TiledArray::TiledRange &,
              const Eigen::MatrixXd &);

template <typename T, typename LR>
bool check_equal(const TiledArray::Array<double, 2, T> &Full,
                 const TiledArray::Array<double, 2, LR> &Low) {
    auto fit = Full.begin();
    auto fend = Full.end();
    auto LRit = Low.begin();

    bool same = true;

    for (; fit != fend && same; ++fit, ++LRit) {
        Eigen::MatrixXd LRmat = LRit->get().matrixLR();
        auto range = fit->get().range();
        Eigen::MatrixXd Fmat = TiledArray::eigen_map(
            fit->get(), range.size()[0], range.size()[1]);
        same = ((Fmat - LRmat).lpNorm<2>() < 1e-06);
        if (same == false) {
            std::cout << "Tile = " << fit.ordinal()
                      << " 2 norm of diff = " << (Fmat - LRmat).lpNorm<2>()
                      << std::endl;
        }
    }

    return same;
}

int main(int argc, char **argv) {
    std::string Sfile;
    std::string Dfile;
    std::string Ffile;
    unsigned long blocksize;
    if (argc < 5) {
        std::cout << "Please provide input files for S, D and F, and "
                     "blocksize\n";
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

    /************************************************
     * Perform LR version
     ************************************************/

    TiledArray::Array<double, 2, LRTile<double>> LR_S
        = make_lr_array(world, trange, ES);
    TiledArray::Array<double, 2, LRTile<double>> LR_D
        = make_lr_array(world, trange, ED);
    TiledArray::Array<double, 2, LRTile<double>> LR_F
        = make_lr_array(world, trange, EF);

    bool passed_check = (check_equal(S, LR_S) == true)
                            ? ((check_equal(D, LR_D) == true)
                                   ? ((check_equal(F, LR_F) == true) ?: false)
                                   : false)
                            : false;

    if (!passed_check) {
        std::cout << "Arrays were not equal!\n";
    }

    int total_tiles = 0;
    int full_tiles = 0;
    for (const auto it : LR_S) {
        full_tiles += int(it.get().is_full());
        ++total_tiles;
    }
    std::cout << "Percentage of full tiles in S = " << double(full_tiles)
                                                  / double(total_tiles) << "\n";


    std::cout << "\n";
    world.gop.fence();
    auto lr_time = madness::wall_time();
    // purify(LR_D, LR_S);
    LR_D("i,j") = 2 * LR_D("i,j") - LR_D("i,k") * LR_S("k,l") * LR_D("l,j");
    //LR_D("i,j") = LR_S("i,k") * LR_S("k,j");
    world.gop.fence();
    lr_time = madness::wall_time() - lr_time;
    std::cout << "LR time was " << lr_time << " s\n";

    world.gop.fence();
    auto full_time = madness::wall_time();
    // purify(D, S);
    D("i,j") = 2 * D("i,j") - D("i,k") * S("k,l") * D("l,j");
    //D("i,j") = S("i,k") * S("k,j");
    world.gop.fence();
    full_time = madness::wall_time() - full_time;
    std::cout << "Full time was " << full_time << " s\n";

    passed_check = (check_equal(D, LR_D) == true) ? true : false;

    full_tiles = 0;
    total_tiles = 0;
    for (const auto it : LR_D) {
        full_tiles += int(it.get().is_full());
        ++total_tiles;
    }
    std::cout << "Percentage of full tiles after = " << double(full_tiles)
                                                  / double(total_tiles) << "\n";

    if (!passed_check) {
        std::cout << "Arrays were not equal!\n";
    } else {
        std::cout << "All Arrays were equal!\n";
    }

    check_data_sparsity(D);

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
make_lr_array(madness::World &world, TiledArray::TiledRange &trange,
              const Eigen::MatrixXd &mat) {
    TiledArray::Array<double, 2, LRTile<double>> A(world, trange);
    auto i = A.begin();
    auto end = A.end();
    for (; i != end; ++i) {
        auto range = trange.make_tile_range(i.ordinal());
        Eigen::MatrixXd mat_block
            = mat.block(range.start()[0], range.start()[1], range.size()[0],
                        range.size()[1]);

        TiledArray::Array<double, 2, LRTile<double>>::value_type tile(
            range, mat_block, true, 1e-03);

        *i = tile;
    }
    return A;
}
