#include <memory>
#include <algorithm>

#include "include/tbb.h"
#include "include/libint.h"
#include "include/tiledarray.h"
#include "include/btas.h"

#include "molecule/atom.h"
#include "molecule/cluster.h"

#include "basis/atom_basisset.h"
#include "basis/basis_set.h"
#include "basis/cluster_shells.h"
#include "basis/basis.h"

#include "integrals/integral_engine_pool.h"
#include "integrals/task_integrals.h"
#include "integrals/low_tile_functors.h"

#include "purification/purification_devel.h"

using namespace tcc;

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    basis::BasisSet bs{"3-21G_basis_G94.txt"};
    std::cout << "Basis set is " << bs << std::endl;
    molecule::Atom h1{{0, 0, 0}, 1, 1};
    molecule::Atom h2{{0, 0, 1}, 1, 1};

    auto cluster = std::make_shared<molecule::Cluster>();
    cluster->add_clusterable(std::move(h1));
    cluster->add_clusterable(std::move(h2));

    basis::Basis basis{bs.create_basis({cluster})};

    std::cout << basis << std::endl;

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
        world, std::move(overlap_pool), basis,
        integrals::compute_functors::LowRankTileFunctor<double>{1e-8});

    libint2::OneBodyEngine kinetic{libint2::OneBodyEngine::kinetic, max_nprim,
                                   static_cast<int>(max_am)};

    auto kinetic_pool = integrals::make_pool(std::move(kinetic));
    auto T = integrals::Integrals(
        world, std::move(kinetic_pool), basis,
        integrals::compute_functors::LowRankTileFunctor<double>{1e-8});

    libint2::OneBodyEngine nuclear{libint2::OneBodyEngine::nuclear, max_nprim,
                                   static_cast<int>(max_am)};

    auto nuclear_pool = integrals::make_pool(std::move(nuclear));
    auto V = integrals::Integrals(
        world, std::move(nuclear_pool), basis,
        integrals::compute_functors::LowRankTileFunctor<double>{1e-8});

    decltype(S) H(world, S.trange());
    H("i,j") = V("i,j") + T("i,j");

    for(auto it = S.begin(); it != S.end(); ++it){
        std::cout << "\nS Matrix = \n" << it->get().tile().matrix() << std::endl;
    }

    for(auto it = T.begin(); it != T.end(); ++it){
        std::cout << "\nT Matrix = \n" << it->get().tile().matrix() << std::endl;
    }

    for(auto it = V.begin(); it != V.end(); ++it){
        std::cout << "\nV Matrix = \n" << it->get().tile().matrix() << std::endl;
    }

    for(auto it = H.begin(); it != H.end(); ++it){
        std::cout << "\nH Matrix = \n" << it->get().tile().matrix() << std::endl;
    }

    pure::purifier()(H,2);

    world.gop.fence();
    return 0;
}
