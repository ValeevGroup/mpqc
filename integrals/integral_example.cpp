#include <memory>
#include <algorithm>

#include "../include/tbb.h"
#include "../include/libint.h"
#include "../include/tiledarray.h"
#include "../include/btas.h"

#include "../molecule/atom.h"
#include "../molecule/cluster.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "integral_engine_pool.h"
#include "task_integrals.h"

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

        auto temp_nprim = std::max_element(
            shell_vec.begin(), shell_vec.end(),
            [](libint2::Shell const &s, libint2::Shell const
               &t) { return t.nprim() > s.nprim(); })->nprim();

        auto temp_am = cluster.max_am();

        max_nprim = (temp_nprim > max_nprim) ? temp_nprim : max_nprim;
        max_am = (temp_am > max_am) ? temp_am : max_am;
    }

    libint2::TwoBodyEngine
        <libint2::Coulomb> coulomb_ints{max_nprim, static_cast<int>(max_am)};
    auto eri_pool = integrals::make_pool(std::move(coulomb_ints));

    auto a = integrals::Integrals(world, std::move(eri_pool), basis);
    std::cout << a << std::endl;
    return 0;
}
