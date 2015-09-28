#include <memory>
// #include <fstream>
// #include <algorithm>

#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../include/libint.h"
#include "../include/tiledarray.h"
#include "../include/btas.h"

#include "../utility/make_array.h"
// #include "../utility/parallel_print.h"
// #include "../utility/array_storage.h"
// #include "../utility/ta_helpers.h"
// #include "../utility/time.h"

#include "../molecule/atom.h"
#include "../molecule/cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/cluster_shells.h"
#include "../basis/basis.h"

#include "../integrals/btas_to_ta_tensor.h"
#include "../integrals/integral_engine_pool.h"
#include "../integrals/sparse_task_integrals.h"
#include "../integrals/make_engine.h"

#include "../tensor/tcc_tile.h"
#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_algebra.h"
#include "../tensor/decomposed_tensor_unary.h"

// #include "../ta_routines/sqrt_inv.h"
// #include "../ta_routines/inverse.h"
// #include <madness/world/array_addons.h>

using namespace tcc;

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    int nclusters = 0;
    double threshold = 1e-13;
    if (argc >= 4) {
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

    struct convert_2d {
        double cut_;
        convert_2d(double thresh) : cut_{thresh} {}
        using TileType = tensor::Tile<tensor::DecomposedTensor<double>>;
        TileType operator()(tensor::ShallowTensor<2> const &bt) {
            auto range = bt.range();

            auto const extent = range.extent();
            const auto X = extent[0];
            const auto i = extent[1];
            auto local_range = TA::Range{X, i};

            auto tensor = TA::Tensor<double>(local_range);
            const auto b_data = bt.tensor().data();
            const auto size = bt.tensor().size();
            std::copy(b_data, b_data + size, tensor.data());

            auto dense
                  = tensor::DecomposedTensor<double>(cut_, std::move(tensor));

            auto test = tensor::algebra::two_way_decomposition(dense);
            if (!test.empty()) {
                dense = std::move(test);
            }

            return tensor::Tile<tensor::DecomposedTensor<double>>(
                  range, std::move(dense));
        }
    };

    auto eri2 = integrals::BlockSparseIntegrals(
          world, eri_pool, utility::make_array(basis, basis),
          convert_2d(1e-8));

    auto eri2_norm = eri2("i,j").norm(world).get();

    if(world.rank() == 0){
        std::cout << "All done norm of ints was " << eri2_norm << std::endl;
    }

    libint2::cleanup();
    madness::finalize();
    return 0;
}
