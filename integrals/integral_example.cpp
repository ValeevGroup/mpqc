#include <memory>

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
    molecule::Atom h1{{0, 0, 0}, 1, 1};
    molecule::Atom h2{{0, 0, 1}, 1, 1};

    auto cluster = std::make_shared<molecule::Cluster>();
    cluster->add_clusterable(std::move(h1));
    cluster->add_clusterable(std::move(h2));

    basis::Basis basis{bs.create_basis({cluster})};

    libint2::TwoBodyEngine<libint2::Coulomb> coulomb_ints;
    auto eri_pool = integrals::make_pool(std::move(coulomb_ints));

    auto a = integrals::Integrals(world, std::move(eri_pool), basis);
    std::cout << a << std::endl;
    return 0;
}
