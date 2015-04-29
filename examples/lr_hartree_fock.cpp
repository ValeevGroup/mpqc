#include <memory>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <chrono>

#include "../include/tbb.h"
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

#include "../purification/sqrt_inv.h"
#include "../purification/purification_devel.h"

using namespace tcc;
namespace ints = integrals;

void main_print_clusters(
      std::vector<std::shared_ptr<molecule::Cluster>> const &bs,
      std::ostream &os) {
    std::vector<std::vector<molecule::Atom>> clusters;
    auto total_atoms = 0ul;
    for (auto const &cluster : bs) {
        std::vector<molecule::Atom> temp;
        for (auto atom : molecule::collapse_to_atoms(*cluster)) {
            temp.push_back(std::move(atom));
            ++total_atoms;
        }
        clusters.push_back(std::move(temp));
    }

    os << total_atoms << std::endl;
    os << "Whole molecule" << std::endl;

    for (auto const &cluster : clusters) {
        for (auto const &atom : cluster) {
            auto center = 0.52917721092 * atom.center();
            os << atom.charge() << " " << center[0] << " " << center[1] << " "
               << center[2] << std::endl;
        }
    }
    os << std::endl;
    auto counter = 0ul;
    for (auto const &cluster : clusters) {
        os << cluster.size() << std::endl;
        os << "Cluster " << counter++ << std::endl;
        for (auto const &atom : cluster) {
            auto center = 0.52917721092 * atom.center();
            os << atom.charge() << " " << center[0] << " " << center[1] << " "
               << center[2] << std::endl;
        }
        os << std::endl;
    }
}

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);

    std::string mol_file = (argc >= 2) ? argv[1] : "";
    int bs_nclusters = (argc >= 3) ? std::stoi(argv[2]) : 0;
    int dfbs_nclusters = (argc >= 4) ? std::stoi(argv[3]) : 0;

    if (mol_file.empty() || 0 == bs_nclusters || 0 == dfbs_nclusters) {
        std::cout << "input is $./program mol_file "
                     "bs_cluster dfbs_clusters "
                     "basis_set(cc-pvdz) df_basis_set(cc-pvdz-ri) "
                     "sparse_threshold(1e-11) "
                     "low_rank_threshhold(1e-7) print_cluster_xyz(true)\n";
        return 0;
    }

    std::string basis_name = (argc >= 5) ? argv[4] : "cc-pvdz";
    std::string df_basis_name = (argc >= 6) ? argv[5] : "cc-pvdz-ri";

    auto threshold = (argc >= 7) ? std::stod(argv[6]) : 1e-11;
    auto low_rank_threshold = (argc >= 8) ? std::stod(argv[7]) : 1e-7;
    bool print_clusters = (argc >= 9) ? std::stod(argv[8]) : true;

    if (world.rank() == 0) {
        std::cout << "Mol file is " << mol_file << std::endl;
        std::cout << "basis is " << basis_name << std::endl;
        std::cout << "df basis is " << df_basis_name << std::endl;
        std::cout << "Using " << bs_nclusters << " bs clusters" << std::endl;
        std::cout << "Using " << dfbs_nclusters << " dfbs clusters"
                  << std::endl;
        std::cout << "low rank threshhold is " << low_rank_threshold
                  << std::endl;
        if (print_clusters) {
            std::cout << "Printing clusters to clusters_bs.xyz and "
                         "cluster_dfbs.xyz." << std::endl;
        }
    }

    TiledArray::SparseShape<float>::threshold(threshold);
    utility::print_par(world, "Sparse threshold is ",
                       TiledArray::SparseShape<float>::threshold(), "\n");

    auto mol = molecule::read_xyz(mol_file);
    auto charge = 0;
    auto occupation = mol.occupation(charge);
    auto repulsion_energy = mol.nuclear_repulsion();

    utility::print_par(world, "Nuclear repulsion_energy = ", repulsion_energy,
                       "\n");

    auto bs_clusters = molecule::attach_hydrogens_kmeans(mol, bs_nclusters);
    auto dfbs_clusters = molecule::attach_hydrogens_kmeans(mol, dfbs_nclusters);

    if (world.rank() == 0) {
        if (print_clusters) {
            std::cout << "Printing clusters\n";
            std::ofstream bs_cluster_file("clusters_bs.xyz");
            main_print_clusters(bs_clusters, bs_cluster_file);
            bs_cluster_file.close();
            std::ofstream dfbs_cluster_file("clusters_dfbs.xyz");
            main_print_clusters(dfbs_clusters, dfbs_cluster_file);
            dfbs_cluster_file.close();
        }
    }
    world.gop.fence();

    basis::BasisSet bs{basis_name};
    basis::BasisSet df_bs{df_basis_name};

    std::streambuf *cout_sbuf = std::cout.rdbuf(); // Silence libint printing.
    std::ofstream fout("/dev/null");
    std::cout.rdbuf(fout.rdbuf());
    basis::Basis basis{bs.create_basis(bs_clusters)};
    basis::Basis df_basis{df_bs.create_basis(dfbs_clusters)};
    std::cout.rdbuf(cout_sbuf);


    libint2::init();

    // Make a btas to decomposed tensor function
    struct convert_2d {
        double cut_;
        convert_2d(double thresh) : cut_{thresh} {}
        using TileType = tensor::Tile<tensor::DecomposedTensor<double>>;
        TileType operator()(tensor::ShallowTensor<2> const &bt) {
            auto range = bt.range();

            auto const &extent = range.size();
            const auto i = extent[0];
            const auto j = extent[1];
            auto local_range = TA::Range{i, j};

            auto tensor = TA::Tensor<double>(local_range);
            const auto b_data = bt.tensor().data();
            const auto size = bt.tensor().size();
            std::copy(b_data, b_data + size, tensor.data());

            auto dense
                  = tensor::DecomposedTensor<double>(cut_, std::move(tensor));

            return tensor::Tile<tensor::DecomposedTensor<double>>(
                  range, std::move(dense));
        }
    };

    auto bs_basis_array = utility::make_array(basis, basis);

    // Compute overlap.
    auto overlap_pool = ints::make_pool(ints::make_1body("overlap", basis));
    auto S = BlockSparseIntegrals(world, overlap_pool, bs_basis_array,
                                  convert_2d(low_rank_threshold));

    auto to_ta = [](tensor::Tile<tensor::DecomposedTensor<double>> const &t) {
        auto tensor = tensor::algebra::combine(t.tile());
        auto range = t.range();
        return TA::Tensor<double>(range, tensor.data());
    };
    auto S_TA = TA::to_new_tile_type(S, to_ta);

    // Invert overlap
    utility::print_par(world, "\nComputing overlap inverse\n");
    auto S_inv_sqrt = pure::inverse_sqrt(S_TA);

    // Compute T
    auto kinetic_pool = ints::make_pool(ints::make_1body("kinetic", basis));
    auto T = BlockSparseIntegrals(world, kinetic_pool, bs_basis_array,
                                  convert_2d(low_rank_threshold));

    /* // Compute V */
    auto nuclear_pool = ints::make_pool(ints::make_1body("nuclear", basis));
    auto V = BlockSparseIntegrals(world, nuclear_pool, bs_basis_array,
                                  convert_2d(low_rank_threshold));

    /* // Compute Hcore */
    utility::print_par(world, "Computing Hcore\n");

    auto T_TA = TA::to_new_tile_type(T, to_ta);
    auto V_TA = TA::to_new_tile_type(V, to_ta);
    decltype(V_TA) H_TA;
    H_TA("i,j") = T_TA("i,j") + V_TA("i,j");
    world.gop.fence();
    /* utility::print_size_info(H, "Hcore"); */

#if 1
    // Super werid race condition
    decltype(V) H;
    H("i,j") = T("i,j") + V("i,j");
    world.gop.fence();
    std::cout << "\nH ranges" << std::endl;
    for(auto fut : H){
        std::cout << fut.get().range() << std::endl;
    }
    std::cout << "\nH_TA ranges" << std::endl;
    for(auto fut : H_TA){
        std::cout << fut.get().range() << std::endl;
    }
    auto H_race = TA::to_new_tile_type(H,to_ta);
    std::cout << "\nH_race ranges" << std::endl;
    for(auto fut : H_race){
        std::cout << fut.get().range() << std::endl;
    }
    std::cout << "H_TA " << H_TA << "\n" << std::endl;
    std::cout << "H_race " << H_race << "\n" << std::endl;

    // If ranges bad will fail here.
    decltype(H_TA) diff;
    diff("i,j") = H_race("i,j") - H_TA("i,j");
    double norm_diff = diff("i,j").norm();
    std::cout << "Diff norm = " << norm_diff << std::endl;
    world.gop.fence();
#endif

    /* // Compute intial density */
    utility::print_par(world, "\nComputing Density\n");
    auto purifier = pure::make_orthogonal_tr_reset_pure(S_inv_sqrt);
    auto D_TA = purifier(H_TA, occupation);
    /* utility::print_size_info(D_TA, "D initial"); */

    /* auto to_decomp = [=](TA::Tensor<double> const &t) { */
    /*     auto range = t.range(); */

    /*     auto const &extent = range.size(); */
    /*     const auto i = extent[0]; */
    /*     const auto j = extent[1]; */
    /*     auto local_range = TA::Range{i, j}; */

    /*     auto tensor = TA::Tensor<double>(local_range, t.data()); */
    /*     auto dense = tensor::DecomposedTensor<double>(low_rank_threshold, */
    /*                                                   std::move(tensor)); */

    /*     return tensor::Tile<tensor::DecomposedTensor<double>>(range, */
    /*                                                           std::move(dense)); */
    /* }; */
    /* auto D = TA::to_new_tile_type(D_TA, to_decomp); */

    /* // Begin Two electron integrals section. */
    /* auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));
     */

    /* // Computing Eri2 */
    /* utility::print_par(world, "\n"); */
    /* auto eri2 = time_and_print_block_sparse( */
    /*       world, eri_pool, utility::make_array(df_basis, df_basis), */
    /*       integrals::compute_functors::BtasToTaTensor{}, "Eri2
     * Integrals");
     */

    /* // Computing the sqrt inverse of Eri2 */
    /* utility::print_par(world, "\nComputing eri2 sqrt Inverse\n"); */
    /* auto inv_timer */
    /*       = tcc_time::make_timer([&]() { return pure::inverse_sqrt(eri2);
     * });
     */
    /* auto eri2_sqrt_inv = inv_timer.apply(); */
    /* utility::print_par(world, "Eri2 inverse computation time = ", */
    /*                    inv_timer.time(), "\n"); */
    /* utility::print_size_info(eri2_sqrt_inv, "Eri2 sqrt inverse"); */

    /*  * Start using Low rank arrays */
    /* // Compute center integrals */
    /* utility::print_par(world, "\n"); */
    /* auto Xab = time_and_print_block_sparse( */
    /*       world, eri_pool, utility::make_array(df_basis, basis, basis),
     */
    /*       integrals::compute_functors::BtasToTaTensor{}, "Eri3
     * integrals");
     */


    /* Xab("X,i,j") = eri2_sqrt_inv("X,P") * Xab("P,i,j"); */
    /* Xab.truncate(); */

    /* { */
    /*     auto Xab_lr = TA::to_new_tile_type( */
    /*           Xab, integrals::compute_functors::TaToLowRankTensor<3>{ */
    /*                      low_rank_threshold}); */

    /*     utility::print_size_info(Xab_lr, "Xab lr"); */
    /* } */

    /* decltype(D) J, K, F; */
    /* J("i,j") = (Xab("X,a,b") * D("a,b")) * Xab("X,i,j"); */
    /* K("i,j") = (Xab("X,i,b") * D("b,a")) * Xab("X,a,j"); */
    /* F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j"); */

    /* D = purifier(F, occupation); */
    /* auto energy = D("i,j").dot(F("i,j") + H("i,j"), world).get(); */

    /* utility::print_par(world, "\nStarting SCF iterations\n"); */
    /* auto diis = TiledArray::DIIS<decltype(D)>{3, 7}; */
    /* auto iter = 1; */
    /* decltype(F) Ferror; */
    /* auto error = 1.0; */
    /* const auto volume = double(F.trange().elements().volume()); */
    /* while (error >= 1e-12 && iter <= 35) { */
    /*     auto t0 = tcc_time::now(); */

    /*     J("i,j") = (Xab("X,a,b") * D("a,b")) * Xab("X,i,j"); */
    /*     J.truncate(); */

    /*     decltype(Xab) X_temp; */
    /*     X_temp("X,a,i") = Xab("X,i,b") * D("b,a"); */

    /*     auto X_temp_lr = TA::to_new_tile_type( */
    /*           X_temp, integrals::compute_functors::TaToLowRankTensor<3>{
     */
    /*                         low_rank_threshold}); */

    /*     utility::print_size_info(X_temp_lr, "X_temp lr"); */
    /*     utility::print_par(world, "\n"); */

    /*     K("i,j") = X_temp("X,a,i") * Xab("X,a,j"); */

    /*     F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j"); */
    /*     Ferror("i,j") = F("i,k") * D("k,l") * S("l,j") */
    /*                     - S("i,k") * D("k,l") * F("l,j"); */
    /*     Ferror.truncate(); */
    /*     error = Ferror("i,j").norm().get() / volume; */
    /*     diis.extrapolate(F, Ferror); */

    /*     D = purifier(F, occupation); */
    /*     energy = D("i,j").dot(F("i,j") + H("i,j"), world).get(); */
    /*     world.gop.fence(); */

    /*     auto t1 = tcc_time::now(); */

    /*     auto time = tcc_time::duration_in_s(t0, t1); */
    /*     /1* auto ktime = tcc_time::duration_in_s(k0, k1); *1/ */

    /*     /1* utility::print_par(world, "Iteration: ", iter++, " has energy
     * ",
     * *1/ */
    /*     /1*                    std::setprecision(11), energy, " with
     * error ",
     */
    /*      * error, *1/ */
    /*     /1*                    " in ", time, " s with K time ", ktime,
     * "\n");
     * *1/ */
    /*     utility::print_par(world, "Iteration: ", iter++, " has energy ",
     */
    /*                        std::setprecision(11), energy, " with error ",
     * error, */
    /*                        " in ", time, " s \n"); */
    /* } */

    /* utility::print_par(world, "Final energy = ", std::setprecision(11),
     */
    /*                    energy + repulsion_energy, "\n"); */

    /* world.gop.fence(); */
    libint2::cleanup();
    madness::finalize();
    return 0;
}
