#include <memory>
#include <fstream>
#include <algorithm>

#include "include/tbb.h"
#include "include/libint.h"
#include "include/tiledarray.h"
#include "include/btas.h"

#include "utility/make_array.h"

#include "molecule/atom.h"
#include "molecule/cluster.h"
#include "molecule/molecule.h"
#include "molecule/clustering_functions.h"

#include "basis/atom_basisset.h"
#include "basis/basis_set.h"
#include "basis/cluster_shells.h"
#include "basis/basis.h"

#include "integrals/ta_compute_functors.h"
#include "integrals/integral_engine_pool.h"
#include "integrals/task_integrals.h"
#include "integrals/sparse_task_integrals.h"

#include "purification/purification_devel.h"
#include "purification/sqrt_inv.h"

using namespace tcc;

// Bring in Eigen::Matrix
using Matrix
    = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

Matrix compute_1body_ints(const std::vector<libint2::Shell> &shells);
Matrix compute_3center(const std::vector<libint2::Shell> &shells);

molecule::Molecule read_xyz(std::ifstream &f) {
    // Get number of atoms.
    unsigned long natoms = 0;
    f >> natoms;

    constexpr auto ang_to_bohr = 1.0 / 0.52917721092;

    std::string line;
    std::vector<molecule::Clusterable> clusterables;
    while (std::getline(f, line)) {
        if (!line.empty()) {
            std::stringstream ss(line);
            std::string atom = "";
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            ss >> atom;
            ss >> x;
            ss >> y;
            ss >> z;
            x *= ang_to_bohr;
            y *= ang_to_bohr;
            z *= ang_to_bohr;
            if (atom == "H") {
                clusterables.emplace_back(molecule::Atom({x, y, z}, 1, 1));
            } else if (atom == "O") {
                clusterables.emplace_back(molecule::Atom({x, y, z}, 16, 8));
            }
        }
    }
    return molecule::Molecule{std::move(clusterables)};
}

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_file = "";
    int nclusters = 0;
    if (argc == 4) {
        mol_file = argv[1];
        basis_file = argv[2];
        nclusters = std::stoi(argv[3]);
    } else {
        std::cout << "input is $./program mol_file basis_file "
                     "nclusters \n";
        return 0;
    }

    std::ifstream molecule_file(mol_file);
    auto mol = read_xyz(molecule_file);

    basis::BasisSet bs{basis_file};

    if (world.rank() == 0) {
        std::cout << mol.nelements() << " elements with " << nclusters
                  << " clusters" << std::endl;
    }

    std::vector<std::shared_ptr<molecule::Cluster>> clusters;
    clusters.reserve(nclusters);

    for (auto &&cluster : mol.attach_H_and_kmeans(nclusters)) {
        clusters.push_back(
            std::make_shared<molecule::Cluster>(std::move(cluster)));
    }

    basis::Basis basis{bs.create_basis(clusters)};

    auto max_am = basis.max_am();
    auto max_nprim = basis.max_nprim();

    std::cout << "Checking 1 body ints.\n";
    libint2::init();
    libint2::OneBodyEngine overlap(libint2::OneBodyEngine::overlap, max_nprim,
                                   max_am, 0);

    auto overlap_pool = integrals::make_pool(overlap);
    auto TA_overlap = integrals::SparseIntegrals(
        world, overlap_pool, utility::make_array(basis, basis),
        integrals::compute_functors::TaTileFunctor<double>{});

    world.gop.fence();
    TiledArray::Array<double, 2> Dense(world, TA_overlap.trange());

    auto dim = TA_overlap.trange().elements().size()[0];
    Matrix TCC_overlap = Matrix::Zero(dim, dim);
    for (std::size_t i = 0; i < TA_overlap.size(); ++i) {
        if (!TA_overlap.is_zero(i)) {
            TiledArray::Tensor<double> tensor = TA_overlap.find(i);
            auto range = tensor.range();
            TCC_overlap.block(range.start()[0], range.start()[1],
                              range.size()[0], range.size()[1])
                = TiledArray::eigen_map(tensor, range.size()[0],
                                        range.size()[1]);
        }
    }

    std::vector<libint2::Shell> flattened_shells;
    for (auto const &cluster : basis.cluster_shells()) {
        auto shells = cluster.flattened_shells();
        auto end = flattened_shells.end();
        flattened_shells.insert(end, shells.begin(), shells.end());
    }
    auto simple_overlap = compute_1body_ints(flattened_shells);

    auto norm = (simple_overlap - TCC_overlap).lpNorm<2>();
    std::cout << "Norm of TA - Simple = " << norm << std::endl;
    std::cout << "\tNorm/size = " << norm/simple_overlap.size() << std::endl;
    if (norm >= 1e-6 && simple_overlap.cols() <= 20) {
        std::cout << "TCC overlap = \n" << TCC_overlap << std::endl;
        std::cout << "Simple overlap = \n" << simple_overlap << std::endl;
    }

    std::cout << "Checking 3 center 2 body ints.\n";
    libint2::TwoBodyEngine<libint2::Coulomb> eri(max_nprim, max_am);
    auto eri_pool = integrals::make_pool(eri);
    auto TA_eri3 = integrals::SparseIntegrals(
        world, eri_pool, utility::make_array(basis, basis, basis),
        integrals::compute_functors::TaTileFunctor<double>{});

    world.gop.fence();

    auto eri3_dim = TA_eri3.trange().elements().size()[0];
    Matrix TCC_eri3 = Matrix::Zero(eri3_dim, eri3_dim * eri3_dim);
    for (std::size_t i = 0; i < TA_eri3.size(); ++i) {
        if (!TA_eri3.is_zero(i)) {
            TiledArray::Tensor<double> tensor = TA_eri3.find(i);
            auto range = tensor.range();
            auto const &size = range.size();
            auto const &start = range.start();

            const auto data = tensor.data();
            auto counter = 0;
            for (auto i = start[0]; i < start[0] + size[0]; ++i) {
                for (auto j = start[1]; j < start[1] + size[1]; ++j) {
                    for (auto k = start[2]; k < start[2] + size[2];
                         ++k, ++counter) {
                        TCC_eri3(i, eri3_dim * j + k) = *(data + counter);
                    }
                }
            }
        }
    }

    auto simple_eri3 = compute_3center(flattened_shells);
    auto eri3_norm = (simple_eri3 - TCC_eri3).lpNorm<2>();
    std::cout << "Norm diff for eri3 was = " << eri3_norm << std::endl;
    std::cout << "\tNorm/size = " << eri3_norm/simple_eri3.size() << std::endl;
    if (eri3_norm >= 1e-7 && simple_eri3.rows() < 11
        && simple_eri3.cols() < 101) {
        std::cout << "TCC eri3 = \n" << TCC_eri3 << std::endl;
        std::cout << "Simple eri3 = \n" << simple_eri3 << std::endl;
        std::cout << "Diff = \n" << Matrix(simple_eri3 - TCC_eri3) << std::endl;
    }


    libint2::cleanup();
    madness::finalize();
    return 0;
}

Matrix compute_1body_ints(const std::vector<libint2::Shell> &shells) {
    auto n = 0;
    auto max_am = 0;
    auto max_nprim = 0ul;
    for (auto const &shell : shells) {
        n += shell.size();
        max_nprim = std::max(max_nprim, shell.nprim());
        for (auto c : shell.contr) {
            max_am = std::max(c.l, max_am);
        }
    }

    Matrix result = Matrix::Zero(n, n);

    // construct the overlap integrals engine
    libint2::OneBodyEngine engine(libint2::OneBodyEngine::overlap, max_nprim,
                                  max_am, 0);

    std::vector<size_t> shell2bf;
    shell2bf.reserve(shells.size());
    auto nbf = 0;
    for (auto const &shell : shells) {
        shell2bf.push_back(nbf);
        nbf += shell.size();
    }

    // loop over unique shell pairs, {s1,s2} such that s1 >= s2
    // this is due to the permutational symmetry of the real integrals over
    // Hermitian operators: (1|2) = (2|1)
    for (auto s1 = 0; s1 != shells.size(); ++s1) {

        auto bf1 = shell2bf[s1]; // first basis function in this shell
        auto n1 = shells[s1].size();

        for (auto s2 = 0; s2 < shells.size(); ++s2) {

            auto bf2 = shell2bf[s2];
            auto n2 = shells[s2].size();

            // compute shell pair; return is the pointer to the buffer
            const auto *buf = engine.compute(shells[s1], shells[s2]);

            // "map" buffer to a const Eigen Matrix, and copy it to the
            // corresponding blocks of the result
            Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
            result.block(bf1, bf2, n1, n2) = buf_mat;
        }
    }

    return result;
}

Matrix compute_3center(const std::vector<libint2::Shell> &shells) {
    auto n = 0;
    auto max_am = 0;
    auto max_nprim = 0ul;
    for (auto const &shell : shells) {
        n += shell.size();
        max_nprim = std::max(max_nprim, shell.nprim());
        for (auto c : shell.contr) {
            max_am = std::max(c.l, max_am);
        }
    }

    Matrix result = Matrix::Zero(n, n * n);

    // construct the overlap integrals engine
    libint2::TwoBodyEngine<libint2::Coulomb> engine(max_nprim, max_am, 0);


    libint2::Shell unit{
        {0.0}, // exponent
        {{0, false, {1.0}}},
        {{0.0, 0.0, 0.0}} // placed at origin
    };
    unit.renorm();

    auto s1_start = 0;
    for (auto s1 = 0; s1 != shells.size(); ++s1) {
        auto n1 = shells[s1].size();

        auto s2_start = 0;
        for (auto s2 = 0; s2 != shells.size(); ++s2) {
            auto n2 = shells[s2].size();

            auto s3_start = 0;
            for (auto s3 = 0; s3 != shells.size(); ++s3) {
                auto n3 = shells[s3].size();

                const auto *buf_1234
                    = engine.compute(shells[s1], unit, shells[s2], shells[s3]);

                auto counter = 0;
                for (auto i = s1_start; i < s1_start + n1; ++i) {
                    for (auto j = s2_start; j < s2_start + n2; ++j) {
                        for (auto k = s3_start; k < s3_start + n3;
                             ++k, ++counter) {
                            result(i, n * j + k) = *(buf_1234 + counter);
                        }
                    }
                }

                s3_start += n3;
            }
            s2_start += n2;
        }
        s1_start += n1;
    }
    return result;
}
