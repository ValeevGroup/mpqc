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
void check_data_sparsity(const TiledArray::Array<double, 2, T> &A) {
    double recompress = madness::wall_time();
    TiledArray::Array<double, 2, LRTile<double>> B(A.get_world(), A.trange());

    auto it_A = A.begin();
    auto it_B = B.begin();
    auto end = B.end();
    for (; it_B != end; ++it_B, ++it_A) {
        TiledArray::Array<double, 2, LRTile<double>>::value_type tile(
            it_A->get().range(), TiledArray::eigen_map(it_A->get()), true,
            1e-08);
        *it_B = tile;
    }
    recompress = madness::wall_time() - recompress;

    int total_tiles = 0;
    int full_tiles = 0;
    for (const auto it : B) {
        full_tiles += int(it.get().is_full());
        ++total_tiles;
    }
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

    // TODO_TCC add occ
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

Eigen::MatrixXd
create_plottable_array(const TiledArray::Array<double, 2, LRTile<double>> &A);

template <typename T, typename LR>
bool check_equal(const TiledArray::Array<double, 2, T> &Full,
                 const TiledArray::Array<double, 2, LR> &Low) {
    auto fit = Full.begin();
    auto fend = Full.end();
    auto LRit = Low.begin();

    bool same = true;

    for (; fit != fend; ++fit, ++LRit) {
        Eigen::MatrixXd LRmat = LRit->get().matrixLR();
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
                          << "\n\t\t2 norm of diff = " << (Fmat - LRmat).lpNorm
                                                          <2>() << std::endl;
                same = inner_same;
            }
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

    Eigen::MatrixXd ES = read_matrix(Sfile);
    //Eigen::MatrixXd ED = read_matrix(Dfile);
    //EiEigen::MatrixXd EF = read_matrix(Ffile);

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
    //TiledArray::Array<double, 2> D = TiledArray::eigen_to_array
    //Ti    <TiledArray::Array<double, 2>>(world, trange, ED);
    //TiTiledArray::Array<double, 2> F = TiledArray::eigen_to_array
    //Ti    <TiledArray::Array<double, 2>>(world, trange, EF);

    /************************************************
     * Perform LR version
     ************************************************/

    TiledArray::Array<double, 2, LRTile<double>> LR_S
        = make_lr_array(world, trange, ES);
    //TiledArray::Array<double, 2, LRTile<double>> LR_D
    //Ti    = make_lr_array(world, trange, ED);
    //TiTiledArray::Array<double, 2, LRTile<double>> LR_F
    //Ti    = make_lr_array(world, trange, EF);


    std::cout << "\nChecking arrays for approximate equality.\n";
    std::cout << "Checking Overlap . . . ";
    bool passed_check = check_equal(S, LR_S);
    if (!passed_check) {
        std::cout << "Overlaps were not equal!";
    } else {
        std::cout << "Ok!";
    }
  //  std::cout << "\nChecking Density . . . ";
  //    passed_check = check_equal(D, LR_D);
  //    if (!passed_check) {
  //        std::cout << "Densities were not equal!";
  //    } else {
  //        std::cout << "Ok!";
  //    }
  //    std::cout << "\n";
  //    std::cout << "Checking Fock . . . ";
  //    passed_check = check_equal(F, LR_F);
  //    if (!passed_check) {
  //        std::cout << "Focks were not equal!";
  //    } else {
  //        std::cout << "Ok!";
  //    }
    std::cout << "\n";


    world.gop.fence();

    auto lr_time = madness::wall_time();
    // purify(LR_D, LR_S);
    // LR_D("i,j") = 2 * LR_D("i,j") - LR_D("i,k") * LR_S("k,l") * LR_D("l,j");
    // LR_D("i,j") = LR_S("i,k") * LR_S("k,j");
    world.gop.fence();
    lr_time = madness::wall_time() - lr_time;
    std::cout << "LR time was " << lr_time << " s\n";

    world.gop.fence();
    auto full_time = madness::wall_time();
    // purify(D, S);
    // D("i,j") = 2 * D("i,j") - D("i,k") * S("k,l") * D("l,j");
    // D("i,j") = S("i,k") * S("k,j");
    world.gop.fence();
    full_time = madness::wall_time() - full_time;
    std::cout << "Full time was " << full_time << " s\n";

  //  std::cout << "\nChecking arrays for approximate equality.\n";
  //    std::cout << "Checking Overlap . . . ";
  //    passed_check = check_equal(S, LR_S);
  //    if (!passed_check) {
  //        std::cout << "Overlaps were not equal!";
  //    } else {
  //        std::cout << "Ok!";
  //    }
  //    std::cout << "\nChecking Density . . . ";
  //    passed_check = check_equal(D, LR_D);
  //    if (!passed_check) {
  //        std::cout << "Densities were not equal!";
  //    } else {
  //        std::cout << "Ok!";
  //    }
  //    std::cout << "\n";
  //    std::cout << "Checking Fock . . . ";
  //    passed_check = check_equal(F, LR_F);
  //    if (!passed_check) {
  //        std::cout << "Focks were not equal!";
  //    } else {
  //        std::cout << "Ok!";
  //    }
  //    std::cout << "\n";



    auto Eig_LR_F = create_plottable_array(LR_S);
    std::ofstream matrix_out;
    matrix_out.open("LR_out_test.dat");
    matrix_out << Eig_LR_F << std::endl;
    matrix_out.close();

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
            range, mat_block, true, 1e-09);

        *i = tile;
    }
    return A;
}

Eigen::MatrixXd
create_plottable_array(const TiledArray::Array<double, 2, LRTile<double>> &A) {
    // Construct the Eigen matrix
    Eigen::MatrixXd matrix(A.trange().elements().size()[0],
                           A.trange().elements().size()[1]);
    for (auto it = A.begin(); it != A.end(); ++it) {
        auto const &tile = it->get();
        auto L = tile.matrixL();
        std::for_each(L.data(), L.data() + L.size(), [](double &x){ x = 1;});
        auto range = tile.range();
        if (tile.is_full()) {
            matrix.block(range.start()[0], range.start()[1], range.size()[0],
                         range.size()[1]) = L;
        } else {
            auto R = tile.matrixR();
            std::for_each(R.data(), R.data() + R.size(), [](double &x){ x = 1;});
            Eigen::MatrixXd temp = Eigen::MatrixXd::Zero(L.rows(), R.cols());
            temp.topRows(tile.rank()) = R;
            temp.leftCols(tile.rank()) = L;
            matrix.block(range.start()[0], range.start()[1], range.size()[0],
                         range.size()[1]) = temp;
        }
    }
    return matrix;
}
