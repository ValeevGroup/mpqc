#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../include/libint.h"
#include "../include/tiledarray.h"

#include "../utility/make_array.h"

#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "../integrals/integral_engine_pool.h"
#include "../integrals/task_integrals.h"
#include "../integrals/make_engine.h"

#include <memory>

using namespace tcc;

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    int nclusters = 0;
    double threshold = 1e-13;
    if (argc == 4) {
        mol_file = argv[1];
        basis_name = argv[2];
        nclusters = std::stoi(argv[3]);
    } else {
        std::cout << "input is $./program mol_file basis_file nclusters ";
        return 0;
    }
    TiledArray::SparseShape<float>::threshold(threshold);

    auto mol = molecule::read_xyz(mol_file);
    auto clusters = molecule::attach_hydrogens_kmeans(mol, nclusters);

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::BasisSet bs{basis_name};
    basis::Basis basis{bs.create_basis(clusters)};
    std::cout.rdbuf(cout_sbuf);

    libint2::init();

    auto eri_pool = integrals::make_pool(integrals::make_2body(basis));

    auto ta_pass_through = [](TA::TensorD &&ten){
        return TA::TensorD(std::move(ten));
    };

    auto eri2 = mpqc_ints::TaskInts<DnPolicy>(world, eri_pool, 
            utility::make_array(basis, basis),
            ta_pass_through);

    auto eri2_norm = eri2("i,j").norm(world).get();

    if(world.rank() == 0){
        std::cout << "All done norm of ints was " << eri2_norm << std::endl;
    }

    libint2::cleanup();
    madness::finalize();
    return 0;
}

