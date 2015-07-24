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
#include "../utility/parallel_break_point.h"
#include "../utility/array_storage.h"
#include "../utility/time.h"
#include "../utility/json_input.h"

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
#include "../integrals/make_engine.h"
#include "../integrals/integral_engine_pool.h"
#include "../integrals/sparse_task_integrals.h"

#include "../scf/soad.h"
#include "../scf/diagonalize_for_coffs.hpp"

#include "../ta_routines/array_to_eigen.h"

using namespace tcc;
namespace ints = integrals;


void
main_print_clusters(std::vector<std::shared_ptr<molecule::Cluster>> const &bs,
                    std::ostream &os);

int try_main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);

    rapidjson::Document in;
    parse_input(argc, argv, in);

    if (!in.HasMember("xyz file") || !in.HasMember("number of bs clusters")
        || !in.HasMember("number of dfbs clusters")) {
        if (world.rank() == 0) {
            std::cout << "At a minimum your input file must provide\n";
            std::cout << "\"xyz file\", which is path to an xyz input\n";
            std::cout << "\"number of bs clusters\", which is the number of "
                         "clusters in the obs\n";
            std::cout << "\"number of dfbs clusters\", which is the number of "
                         "clusters in the dfbs\n";
        }
    }
    // Get necessary info
    std::string mol_file = in["xyz file"].GetString();
    int bs_nclusters = in["number of bs clusters"].GetInt();
    int dfbs_nclusters = in["number of dfbs clusters"].GetInt();
    int occ_nclusters = in.HasMember("number of occupied clusters")
                              ? in["number of occupied clusters"].GetInt()
                              : 1;

    // Get basis info
    std::string basis_name = in.HasMember("basis") ? in["basis"].GetString()
                                                   : "cc-pvdz";
    std::string df_basis_name = in.HasMember("df basis")
                                      ? in["df basis"].GetString()
                                      : "cc-pvdz-ri";

    // Get thresh info
    auto threshold = in.HasMember("block sparse threshold")
                           ? in["block sparse threshold"].GetDouble()
                           : 1e-13;
    auto low_rank_threshold = in.HasMember("low rank threshold")
                                    ? in["low rank threshold"].GetDouble()
                                    : 1e-8;

    // Get printing info
    bool print_clusters = in.HasMember("print clusters")
                                ? in["print clusters"].GetBool()
                                : false;

    // Use Chol Vectors? 
    bool use_chol_vectors = in.HasMember("use cholesky vectors")
                                ? in["use cholesky vectors"].GetBool()
                                : false;

    volatile int debug
          = in.HasMember("debug break") ? in["debug break"].GetInt() : 0;
    utility::parallal_break_point(world, debug);

    if (world.rank() == 0) {
        std::cout << "Mol file is " << mol_file << std::endl;
        std::cout << "basis is " << basis_name << std::endl;
        std::cout << "df basis is " << df_basis_name << std::endl;
        std::cout << "Using " << bs_nclusters << " bs clusters" << std::endl;
        std::cout << "Using " << dfbs_nclusters << " dfbs clusters"
                  << std::endl;
        std::cout << "Using " << occ_nclusters << " occupied clusters"
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

    world.gop.fence();
    if (world.rank() == 0) {
        if (print_clusters) {
            std::string obs_file = in.HasMember("basis clusters file")
                                         ? in["basis clusters file"].GetString()
                                         : "clusters_bs.xyz";
            std::string dfbs_file = in.HasMember("df basis clusters file")
                          ? in["df basis clusters file"].GetString()
                          : "clusters_dfbs.xyz";

            std::cout << "Printing clusters\n";
            std::cout << "\tobs clusters to " << obs_file
                      << ": over ride with \"basis clusters file\" keyword\n";
            std::cout
                  << "\tdfbs clusters to " << dfbs_file
                  << ": over ride with \"df basis clusters file\" keyword\n";

            std::ofstream bs_cluster_file(obs_file);
            main_print_clusters(bs_clusters, bs_cluster_file);
            bs_cluster_file.close();
            std::ofstream dfbs_cluster_file(dfbs_file);
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

    if (world.rank() == 0) {
        std::cout << "Basis trange " << std::endl;
        TA::TiledRange1 bs_range = basis.create_flattend_trange1();
        std::cout << bs_range << std::endl;
        TA::TiledRange1 dfbs_range = df_basis.create_flattend_trange1();
        std::cout << "DF Basis trange " << std::endl;
        std::cout << dfbs_range << std::endl;
    }

    libint2::init();

    // Make a btas to decomposed tensor function
    struct convert_2d {
        double cut_;
        convert_2d(double thresh) : cut_{thresh} {}
        using TileType = tensor::Tile<tensor::DecomposedTensor<double>>;
        TileType operator()(tensor::ShallowTensor<2> const &bt) {
            auto range = bt.range();

            auto const extent = range.extent();
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
    decltype(V) H;
    H("i,j") = T("i,j") + V("i,j");
    world.gop.fence();

    auto H_TA = TA::to_new_tile_type(H, to_ta);

    auto to_decomp = [=](TA::Tensor<double> const &t) {
        auto range = t.range();

        auto const extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];
        auto local_range = TA::Range{i, j};

        auto tensor = TA::Tensor<double>(local_range, t.data());
        auto dense = tensor::DecomposedTensor<double>(low_rank_threshold,
                                                      std::move(tensor));

        return tensor::Tile<tensor::DecomposedTensor<double>>(range,
                                                              std::move(dense));
    };

    /* // Begin Two electron integrals section. */
    auto eri_pool = ints::make_pool(ints::make_2body(basis, df_basis));

    auto to_decomp_with_decompose = [=](TA::Tensor<double> const &t) {
        auto range = t.range();

        auto const extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];
        auto local_range = TA::Range{i, j};

        auto tensor = TA::Tensor<double>(local_range, t.data());
        auto dense = tensor::DecomposedTensor<double>(low_rank_threshold,
                                                      std::move(tensor));

        auto test = tensor::algebra::two_way_decomposition(dense);
        if (!test.empty()) {
            dense = std::move(test);
        }

        return tensor::Tile<tensor::DecomposedTensor<double>>(range,
                                                              std::move(dense));
    };
    /* // Computing Eri2 */
    utility::print_par(world, "Starting 2 Center Integrals\n");
    TA::Array<double, 2, tensor::Tile<tensor::DecomposedTensor<double>>,
              TA::SparsePolicy> V_inv_oh;
    {
        utility::print_par(world, "\tStarting V\n");
        auto eri2 = ints::BlockSparseIntegrals(
              world, eri_pool, utility::make_array(df_basis, df_basis),
              integrals::compute_functors::BtasToTaTensor{});

        decltype(eri2) L_inv_TA;
        {
            utility::print_par(world, "\tReplicating to Eigen\n");
            auto eig_E2 = array_ops::array_to_eigen(eri2);
            Eig::LLT<decltype(eig_E2)> llt(eig_E2);
            eig_E2 = llt.matrixL();
            decltype(eig_E2) eig_L_inv = eig_E2.inverse();
            utility::print_par(world, "\tConverting back to TA\n");
            L_inv_TA = array_ops::eigen_to_array<TA::Tensor<double>>(
                  world, eig_L_inv, eri2.trange().data()[0],
                  eri2.trange().data()[1]);
        }

        utility::print_par(world, "\tDecomposing Tiles\n");
        V_inv_oh = TA::to_new_tile_type(L_inv_TA, to_decomp_with_decompose);
        utility::print_size_info(V_inv_oh, "V^{-1/2}");
    }
    decltype(V_inv_oh)::wait_for_lazy_cleanup(world, 2);


    struct convert_3d {
        double cut_;
        convert_3d(double thresh) : cut_{thresh} {}
        using TileType = tensor::Tile<tensor::DecomposedTensor<double>>;
        TileType operator()(tensor::ShallowTensor<3> const &bt) {
            auto range = bt.range();

            auto const extent = range.extent();
            const auto X = extent[0];
            const auto i = extent[1];
            const auto j = extent[2];
            auto local_range = TA::Range{X, i, j};

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

    // Compute center integrals
    utility::print_par(world, "\nStarting 3 Center Integrals\n");
    auto E0 = tcc_time::now();
    auto Xab = ints::BlockSparseIntegrals(
          world, eri_pool, utility::make_array(df_basis, basis, basis),
          convert_3d(low_rank_threshold));
    auto E1 = tcc_time::now();
    auto etime = tcc_time::duration_in_s(E0, E1);
    utility::print_par(world, "Time to compute 3 center integrals ", etime,
                       " s\n");
    utility::print_size_info(Xab, "E");
    decltype(Xab)::wait_for_lazy_cleanup(world, 60);

    // Make B tensor
    auto B0 = tcc_time::now();
    Xab("X,i,j") = V_inv_oh("X,P") * Xab("P,i,j");
    Xab.truncate();
    auto B1 = tcc_time::now();
    auto btime = tcc_time::duration_in_s(B0, B1);
    utility::print_par(world, "\nTime to compute B ", btime, " s\n");
    utility::print_size_info(Xab, "B Tensor");
    decltype(Xab)::wait_for_lazy_cleanup(world, 60);

    decltype(H) F;
    utility::print_par(world, "\nStarting SOAD guess");
    auto soad0 = tcc_time::now();
    F = ints::scf::fock_from_minimal_v_oh(world, basis, df_basis, eri_pool, H,
                                          V_inv_oh, Xab, bs_clusters,
                                          low_rank_threshold * 100,
                                          convert_3d(low_rank_threshold));
    auto soad1 = tcc_time::now();
    auto soad_time = tcc_time::duration_in_s(soad0, soad1);
    utility::print_par(world, "\nSoad time ", soad_time, " s\n");
    decltype(F)::wait_for_lazy_cleanup(world);

    utility::print_par(world, "\nConverting Fock to TA::Tensor...\n");
    auto F_TA = TA::to_new_tile_type(F, to_ta);
    world.gop.fence();

    auto n_occ = occupation / 2;
    auto tr_i = scf::tr_occupied(occ_nclusters, n_occ);
    utility::print_par(world, "Computing MO coeffs...\n");
    auto Coeffs_TA = scf::Coeffs_from_fock(F_TA, S_TA, tr_i, n_occ, use_chol_vectors);
    utility::print_par(world, "Converting Coeffs to Decomp Form...\n");
    auto Coeffs = TA::to_new_tile_type(Coeffs_TA, to_decomp);

    decltype(Coeffs_TA) D_TA;
    utility::print_par(world, "Forming Density...\n");
    D_TA("i,j") = Coeffs_TA("i,a") * Coeffs_TA("j,a");

    utility::print_par(world, "Computing Initial energy...\n");
    auto energy = D_TA("i,j").dot(F_TA("i,j") + H_TA("i,j"), world).get();
    utility::print_par(world, "Initial energy = ", energy + repulsion_energy,
                       "\n");

    auto D = to_new_tile_type(D_TA, to_decomp);

    decltype(D) J, K, SC, Kocc;
    utility::print_par(world, "\nStarting SCF iterations\n");
    auto diis = TiledArray::DIIS<decltype(D_TA)>(1);
    auto iter = 1;
    decltype(F_TA) Ferror;
    auto error = 1.0;
    auto old_e = energy;
    auto delta_e = energy;
    const auto volume = double(F.trange().elements().volume());
    decltype(Xab) W;
    double time;
    double ktime, jtime;
    while ((error >= 1e-13 || std::abs(delta_e) >= 1e-12) && iter <= 35) {
        utility::print_par(world, "Iteration: ", iter, "\n");
        auto t0 = tcc_time::now();
        D = to_new_tile_type(D_TA, to_decomp);

        utility::print_par(world, "\tStarting W...  ");
        auto w0 = tcc_time::now();
        W("X,a,i") = Xab("X,a,b") * Coeffs("b,i");
        auto w1 = tcc_time::now();
        auto wtime = tcc_time::duration_in_s(w0, w1);
        utility::print_par(world, wtime, " s\n");

        utility::print_par(world, "\tStarting Coulomb...  ");
        auto j0 = tcc_time::now();
        J("i,j") = Xab("X,i,j") * (W("X,a,i") * Coeffs("a,i"));
        auto j1 = tcc_time::now();
        jtime = tcc_time::duration_in_s(j0, j1);
        utility::print_par(world, jtime, " s\n");

        utility::print_par(world, "\tStarting Exchange... ");
        auto k0 = tcc_time::now();
        W("X,i,a") = W("X,a,i");
        K("a,j") = W("X,i,a") * (W("X,i,b") * Coeffs("b,j"));
        Kocc("i,j") = Coeffs("a,i") * K("a,j");
        SC("a,i") = S("a,k") * Coeffs("k,i");
        K("a,b") = SC("a,i") * K("b,i") + K("a,i") * SC("b,i") - 
            SC("a,i") * Kocc("i,j") * SC("b,j");
        auto k1 = tcc_time::now();
        ktime = tcc_time::duration_in_s(k0, k1);
        utility::print_par(world, ktime, " s\n");


        F("i,j") = H("i,j") + 2 * J("i,j") - K("i,j");
        F_TA = TA::to_new_tile_type(F, to_ta);
        // to make energy variational need to use same density for F and energy!
        energy = D_TA("i,j").dot(F_TA("i,j") + H_TA("i,j"), world).get();

        Ferror("i,j") = F_TA("i,k") * D_TA("k,l") * S_TA("l,j")
                        - S_TA("i,k") * D_TA("k,l") * F_TA("l,j");

        error = Ferror("i,j").norm().get() / volume;
        diis.extrapolate(F_TA, Ferror);

        Coeffs_TA = scf::Coeffs_from_fock(F_TA, S_TA, tr_i, n_occ, use_chol_vectors);
        Coeffs = TA::to_new_tile_type(Coeffs_TA, to_decomp);
        D_TA("i,j") = Coeffs_TA("i,a") * Coeffs_TA("j,a");

        delta_e = energy - old_e;
        old_e = energy;

        auto t1 = tcc_time::now();
        time = tcc_time::duration_in_s(t0, t1);

        utility::print_par(world, "\tHas energy ", std::setprecision(17),
                           energy + repulsion_energy, " with error ", error,
                           " in ", time, " s \n");
        utility::print_par(world, "\t\tDelta e = ", delta_e, "\n");
        if (iter == 5) {
            utility::print_par(world, "\n\nTemporary Info:\n");
            utility::print_size_info(W, "Exch Temp");
            utility::print_par(world, "\n\n");
        }
        ++iter;
    }

    utility::print_par(world, "\nFinal energy = ", std::setprecision(17),
                       energy + repulsion_energy, "\n");


    world.gop.fence();
    libint2::cleanup();
    madness::finalize();
    return 0;
}

int main(int argc, char** argv){
    try{
        try_main(argc, argv);
    } catch(const madness::MadnessException &e){
        std::cout << "Madness Exception Says " << e.what() << std::endl;
    } catch(const TiledArray::Exception &e){
        std::cout << "TA Exception Says " << e.what() << std::endl;
    } catch(const std::exception &e){
        std::cout << "std Exception Says " << e.what() << std::endl;
    } catch(...){
        std::cout << "Caught unknown exception" << std::endl;
    }
    return 0;
}


void
main_print_clusters(std::vector<std::shared_ptr<molecule::Cluster>> const &bs,
                    std::ostream &os) {

    // Collect atoms
    std::vector<molecule::Atom> atoms;
    for (auto const &shared_cluster : bs) {
        auto current_atoms = molecule::collapse_to_atoms(*shared_cluster);
        atoms.insert(atoms.end(), current_atoms.begin(), current_atoms.end());
    }

    // Print all atoms
    os << atoms.size() << std::endl;
    os << "Entire Molecule" << std::endl;
    for (auto const &atom : atoms) {
        os << atom.xyz_string() << std::endl;
    }
    os << std::endl;

    // Print clusters
    for (auto const &shared_cluster : bs) {
        os << *shared_cluster << std::endl;
    }
}
