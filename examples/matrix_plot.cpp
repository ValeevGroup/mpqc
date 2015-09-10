#include <memory>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include "../include/libint.h"
#include "../include/tiledarray.h"
#include "../include/btas.h"

#include "../utility/make_array.h"
#include "../utility/parallel_print.h"
#include "../utility/array_storage.h"
#include "../utility/time.h"
#include "../utility/ta_helpers.h"

#include "../tensor/conversions/tile_pimpl_to_ta_tensor.h"

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
#include "../integrals/btas_to_low_rank_tensor.h"
#include "../integrals/make_engine.h"
#include "../integrals/ta_tensor_to_low_rank_tensor.h"
#include "../integrals/integral_engine_pool.h"
#include "../integrals/sparse_task_integrals.h"
#include "../integrals/dense_task_integrals.h"

#include "../density/sqrt_inv.h"
#include "../density/purification_devel.h"

using namespace tcc;
namespace ints = integrals;

int main(int argc, char **argv) {
    auto &world = madness::initialize(argc, argv);
    if(argc != 6){
        std::cout << "Bad input" << std::endl;
    }
    std::string mol_file = argv[1];
    std::string out_file = argv[2];
    std::string basis_name = argv[3];
    int bs_nclusters = std::stoi(argv[4]);
    auto low_rank_threshold = std::stod(argv[5]);

    auto mol = molecule::read_xyz(mol_file);

    auto bs_clusters = molecule::attach_hydrogens_kmeans(mol, bs_nclusters);

    basis::BasisSet bs{basis_name};

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::Basis basis{bs.create_basis(bs_clusters)};
    std::cout.rdbuf(cout_sbuf);

    libint2::init();
    auto eri_pool = ints::make_pool(ints::make_2body(basis));

    auto dfbasis_array = utility::make_array(basis, basis);
    auto eri2
        = BlockSparseIntegrals(world, eri_pool, dfbasis_array,
                               integrals::compute_functors::BtasToTaTensor{});

    auto eri2_lr = TA::to_new_tile_type(
        eri2,
        integrals::compute_functors::TaToLowRankTensor<2>(low_rank_threshold));

    auto const &dims = eri2_lr.trange().elements().size();
    RowMatrixXd eig_lr = RowMatrixXd::Zero(dims[0], dims[1]);

    for(auto it = eri2_lr.begin(); it != eri2_lr.end(); ++it){
        auto const &tile = it->get();
        auto const &range = tile.range();
        auto const &start = range.start();
        auto const &size = range.extent();

        RowMatrixXd eig_block = RowMatrixXd::Zero(size[0], size[1]);
        if(tile.isFull()){
            eig_block = RowMatrixXd(tile.tile().ftile().matrix());
        }else {
            auto const &L = tile.tile().lrtile().matrixL();
            auto const &R = tile.tile().lrtile().matrixR();

            eig_block.topRows(R.rows()) = R;
            eig_block.leftCols(L.cols()) = L;
        }

        eig_lr.block(start[0], start[1], size[0], size[1]) = eig_block;
    }

    std::ofstream out_stream(out_file.c_str());
    out_stream << std::setprecision(15);
    out_stream << std::fixed;
    out_stream << eig_lr << std::endl;
    out_stream.close();

    madness::finalize();
    return 0;
}
