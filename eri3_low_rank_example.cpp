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

#include "integrals/btas_to_ta_tensor.h"
#include "integrals/btas_to_low_rank_tensor.h"
#include "integrals/ta_tensor_to_low_rank_tensor.h"
#include "integrals/integral_engine_pool.h"
#include "integrals/dense_task_integrals.h"
#include "integrals/sparse_task_integrals.h"

#include "purification/purification_devel.h"
#include "purification/sqrt_inv.h"

using namespace tcc;


template <typename Policy>
double array_diff(TiledArray::Array<double, 3, TiledArray::Tensor<double>,
                                    Policy> const &other,
                  TiledArray::Array<double, 3, tensor::TilePimpl<double>,
                                    TiledArray::SparsePolicy> const &lr) {
    double out = 0;
    auto const &pmap_ptr = other.get_pmap();
    const auto end = pmap_ptr->end();
    auto it_lr = lr.get_pmap()->begin();

    for (auto it = pmap_ptr->begin(); it != end; ++it, ++it_lr) {
        const auto range = other.trange().make_tile_range(*it);
        const auto rows = range.size()[0];
        const auto cols = range.size()[1] * range.size()[2];
        Eigen::MatrixXd s_tile = Eigen::MatrixXd::Zero(rows, cols);
        Eigen::MatrixXd lr_tile = Eigen::MatrixXd::Zero(rows, cols);

        if (!other.is_zero(*it)) {
            const auto tensor = other.find(*it).get();
            s_tile = TiledArray::eigen_map(tensor, rows, cols);
        }
        if (!lr.is_zero(*it_lr)) {
            tensor::TilePimpl<double> tile_lr = lr.find(*it_lr);
            lr_tile = tile_lr.tile().matrix();
        }
        const auto norm = (s_tile - lr_tile).lpNorm<2>();
        out += norm * norm;
    }

    // sum outs from all nodes
    other.get_world().gop.sum(&out, 1);
    other.get_world().gop.fence();
    out = std::sqrt(out) / other.trange().elements().volume();
    return out;
}


int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    std::string df_basis_name = "";
    int nclusters = 0;
    int debug = 0;
    if (argc >= 5) {
        mol_file = argv[1];
        basis_name = argv[2];
        df_basis_name = argv[3];
        nclusters = std::stoi(argv[4]);
    } else {
        std::cout << "input is $./program mol_file basis_file df_basis_file "
                     "nclusters \n";
        return 0;
    }
    if (argc == 6) {
        debug = std::stoi(argv[5]);
    }

    auto mol = molecule::read_xyz(mol_file);

    basis::BasisSet bs{basis_name};
    basis::BasisSet df_bs{df_basis_name};

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
    basis::Basis df_basis{df_bs.create_basis(clusters)};

    auto max_am = std::max(basis.max_am(), df_basis.max_am());
    auto max_nprim = std::max(basis.max_nprim(), df_basis.max_nprim());

    libint2::init();
    world.gop.fence();

    libint2::TwoBodyEngine<libint2::Coulomb> eri{max_nprim,
                                                 static_cast<int>(max_am)};


    auto eri_pool = integrals::make_pool(std::move(eri));

    world.gop.fence();
    auto eri2_dense = integrals::DenseIntegrals(
        world, eri_pool, utility::make_array(df_basis, df_basis),
        integrals::compute_functors::BtasToTaTensor{});

    auto eri2_inv_sqrt = pure::inverse_sqrt(eri2_dense);
    {
        decltype(eri2_dense) ident;
        ident("i,j") = eri2_inv_sqrt("i,k") * eri2_dense("k,l")
                       * eri2_inv_sqrt("l,j");
        decltype(eri2_dense) zero;
        zero("i,j") = pure::create_diagonal_matrix(eri2_dense, 1.0)("i,j")
                      - ident("i,j");
        auto norm = zero("i,j").norm().get();
        if (world.rank() == 0) {
            std::cout << "Norm of S^{-1/2} * S * S^{-1/2} - ident = " << norm
                      << std::endl;
        }
    }

    auto eri3_dense = integrals::DenseIntegrals(
        world, eri_pool, utility::make_array(df_basis, basis, basis),
        integrals::compute_functors::BtasToTaTensor{});

    eri3_dense("X,a,b") = eri2_inv_sqrt("X,P") * eri3_dense("P,a,b");

    if (world.rank() == 0) {
        std::cout << "\nConverting to low rank block sparse array "
                     "from dense" << std::endl;
    }
    {
        auto eri3 = TiledArray::to_sparse(eri3_dense);
        auto eri3_low_rank = TiledArray::conversion::to_new_tile_type(
            eri3, integrals::compute_functors::TaToLowRankTensor<3>{1e-8});

        if (world.rank() == 0) {
            std::cout << "Eri3 sparse low rank from dense storage ";
        }

        array_storage(eri3_low_rank);
        double dvlr_eri3_diff = array_diff(eri3_dense, eri3_low_rank);
        double svlr_eri3_diff = array_diff(eri3, eri3_low_rank);
        if (world.rank() == 0) {
            std::cout << "Block sparse low rank error with dense array in the "
                         "contracted product = " << dvlr_eri3_diff << std::endl;
            std::cout << "Block sparse low rank error with block sparse array "
                         "in the contracted product = " << dvlr_eri3_diff
                      << std::endl;
        }
    }

    const auto thresh = TiledArray::SparseShape<float>::threshold();
    if (world.rank() == 0) {
        std::cout << "Finally the sparse shape threshold was = " << thresh
                  << std::endl;
    }

    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}
