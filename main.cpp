#include <memory>
#include <algorithm>

#include "include/tbb.h"
#include "include/libint.h"
#include "include/tiledarray.h"
#include "include/btas.h"

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

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    basis::BasisSet bs{"3-21G_basis_G94.txt"};

    std::vector<molecule::Clusterable> clusterables;
    for (auto i = 0; i < 500; ++i) {
        clusterables.emplace_back(molecule::Atom{{0, 0, i}, 1, 1});
    }

    molecule::Molecule mol{std::move(clusterables)};

    auto cluster_func = molecule::clustering::kmeans{42};
    std::vector<std::shared_ptr<molecule::Cluster>> clusters;
    clusters.reserve(20);

    auto cluster_num = 0;
    for (auto &&cluster : mol.cluster_molecule(cluster_func, 20)) {
        std::cout << "cluster " << cluster_num << " has " << cluster.nelements()
                  << " atoms\n";
        ++cluster_num;
        for (auto i = cluster.begin(); i != cluster.end(); ++i) {
            std::cout << "\t" << i->center().transpose() << std::endl;
        }
        clusters.push_back(
            std::make_shared<molecule::Cluster>(std::move(cluster)));
    }

    basis::Basis basis{bs.create_basis(clusters)};

    auto max_nprim = 0ul;
    auto max_am = 0ul;

    for (auto const &cluster : basis.cluster_shells()) {
        auto const &shell_vec = cluster.flattened_shells();

        auto temp_nprim = std::max_element(shell_vec.begin(), shell_vec.end(),
                                           [](libint2::Shell const &s,
                                              libint2::Shell const &t) {
                                               return t.nprim() > s.nprim();
                                           })->nprim();

        auto temp_am = cluster.max_am();

        max_nprim = (temp_nprim > max_nprim) ? temp_nprim : max_nprim;
        max_am = (temp_am > max_am) ? temp_am : max_am;
    }

    libint2::OneBodyEngine overlap{libint2::OneBodyEngine::overlap, max_nprim,
                                   static_cast<int>(max_am)};

    auto overlap_pool = integrals::make_pool(std::move(overlap));
    auto S = integrals::Integrals(
        world, overlap_pool, basis,
        integrals::compute_functors::TaTileFunctor<double>{});

    libint2::OneBodyEngine kinetic{libint2::OneBodyEngine::kinetic, max_nprim,
                                   static_cast<int>(max_am)};

    auto kinetic_pool = integrals::make_pool(std::move(kinetic));
    auto T = integrals::Integrals(
        world, std::move(kinetic_pool), basis,
        integrals::compute_functors::TaTileFunctor<double>{});

    libint2::OneBodyEngine nuclear{libint2::OneBodyEngine::nuclear, max_nprim,
                                   static_cast<int>(max_am)};

    auto nuclear_pool = integrals::make_pool(std::move(nuclear));
    auto V = integrals::Integrals(
        world, std::move(nuclear_pool), basis,
        integrals::compute_functors::TaTileFunctor<double>{});


    decltype(S) H;
    H("i,j") = V("i,j") + T("i,j");

    // std::cout << "S = \n" << S << std::endl;
    // std::cout << "H = \n" << H << std::endl;

    //    auto sparse_S = integrals::SparseIntegrals(
    //        world, overlap_pool, basis,
    //        integrals::compute_functors::TaTileFunctor<double>{});
    //
    //    auto num = 0;
    //    for (auto it = sparse_S.begin(); it != sparse_S.end(); ++it) {
    //        auto elements = it->get().range().volume();
    //        std::cout << "Tile is not sparse " << it.ordinal() << " has "
    //                  << elements << " elements " << std::endl;
    //        ++num;
    //    }
    //    std::cout << "Given should have " << 100 << " tiles, but only had " <<
    //    num
    //              << std::endl;
    //
    //    auto elements = sparse_S.elements().volume();
    //    auto dim = std::sqrt(elements);
    //    std::cout << "Matrix has " << elements << " elements give a dimension
    //    of "
    //              << dim << std::endl;
    //
    //    TiledArray::Array<double, 2, TiledArray::Tensor<double>,
    //                      TiledArray::SparsePolicy> S2;
    //    world.gop.fence();
    //    S2("i,j") = sparse_S("i,k") * sparse_S("k,j");
    //    world.gop.fence();
    //
    //    std::cout << "Going to print S2" << std::endl;
    //    num = 0;
    //    for (auto it = S2.begin(); it != S2.end(); ++it) {
    //        auto elements = it->get().range().volume();
    //        std::cout << "Tile is not sparse " << it.ordinal() << " has "
    //                  << elements << " elements " << std::endl;
    //
    //        ++num;
    //    }
    //    std::cout << "Given should have " << 100 << " tiles, but only had " <<
    //    num
    //              << std::endl;

    auto SsqrtInv = pure::inverse_sqrt(S);
    auto eig_S = TiledArray::array_to_eigen(S);
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(eig_S);
    auto real_inv = es.operatorInverseSqrt();
    // std::cout << "Mine = \n" << TiledArray::array_to_eigen(SsqrtInv) <<
    // std::endl;
    // std::cout << "Correct = \n" << real_inv << std::endl;
    auto diff = real_inv - TiledArray::array_to_eigen(SsqrtInv);
    std::cout << "Norm Diff between correct and mine is " << diff.lpNorm<2>()
              << std::endl;

    //  auto P = pure::purifier()(H, S, 10);
    //  std::cout << "First density = \n" << P << std::endl;
    //  decltype(S) PS; PS("i,j") = P("i,k") * S("k,j");
    //  std::cout << "PS trace = " << PS("i,j").trace().get() << std::endl;
    //    auto X = create_eval_scaled_guess(H,S);
    //    std::cout << "X = \n" << X << std::endl;

    world.gop.fence();
    return 0;
}
