#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../include/libint.h"
#include "../include/tiledarray.h"

#include "../utility/make_array.h"
#include "../clustering/kmeans.h"

#include "../molecule/atom.h"
#include "../molecule/atom_based_cluster.h"
#include "../molecule/molecule.h"
#include "../molecule/clustering_functions.h"
#include "../molecule/make_clusters.h"

#include <rapidjson/document.h>
#include <rapidjson/stringbuffer.h>
#include <rapidjson/writer.h>

#include "../utility/json_input.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/basis.h"

#include "../integrals/integrals.h"

#include "../utility/time.h"
#include "../utility/array_info.h"
#include "../ta_routines/array_to_eigen.h"
#include "../ta_routines/minimize_storage.h"

#include "../scf/diagonalize_for_coffs.hpp"
#include "../scf/soad.h"
#include "../scf/orbital_localization.h"
#include "../scf/clusterd_coeffs.h"
#include "../scf/cadf_helper_functions.h"
#include "../scf/scf.h"

#include "../scf/cadf_builder.h"
#include "../scf/cadf_builder_forced_shape.h"


#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_nonintrusive_interface.h"
#include "../tensor/mpqc_tile.h"
#include "../tensor/tensor_transforms.h"

#include <memory>

using namespace mpqc;
namespace ints = mpqc::integrals;

bool tensor::detail::recompress = true;

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::cout << std::setprecision(15);

    rapidjson::Document in;
    parse_input(argc, argv, in);

    double ta_threshold = in.HasMember("sparse threshold")
                                ? in["sparse threshold"].GetDouble()
                                : 1e-11;

    TiledArray::SparseShape<float>::threshold(ta_threshold);
    std::cout << "TA::Sparse Threshold: " << ta_threshold << std::endl;

    if (in.HasMember("engine precision")) {
        integrals::detail::integral_engine_precision
              = in["engine precision"].GetDouble();
    }

    std::string mol_file;
    if (!in.HasMember("xyz file")) {
        if (world.rank() == 0) {
            std::cout << "Did not detect xyz file in input.\n";
        }
        madness::finalize();
        return 0;
    } else {
        mol_file = in["xyz file"].GetString();
    }

    auto atom_mol = molecule::read_xyz(mol_file);

    molecule::Molecule clustered_mol;
    if (in.HasMember("cluster by atom") && in["cluster by atom"].GetBool()) {
        clustered_mol = std::move(atom_mol);
    } else if (in.HasMember("number of clusters")) {
        auto nclusters = in["number of clusters"].GetInt();
        clustered_mol
              = molecule::attach_hydrogens_and_kmeans(atom_mol.clusterables(),
                                                      nclusters);
    } else {
        if (world.rank() == 0) {
            std::cout << "Did not specify the number of clusters\n";
        }
        madness::finalize();
        return 0;
    }

    auto repulsion_energy = clustered_mol.nuclear_repulsion();
    std::cout << "Nuclear Repulsion Energy: " << repulsion_energy << std::endl;
    auto occ = clustered_mol.occupation(0);


    std::string basis_name = in["obs name"].GetString();
    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));
    std::cout << "Basis has " << basis.nfunctions() << " functions\n";

    auto df_clustered_mol = clustered_mol;

    std::string df_basis_name = in["df name"].GetString();
    basis::BasisSet dfbs(df_basis_name);
    basis::Basis df_basis(dfbs.get_cluster_shells(df_clustered_mol));

    libint2::init();
    const auto bs_array = utility::make_array(basis, basis);

    // Overlap ints
    auto overlap_e = ints::make_1body_shr_pool("overlap", basis, clustered_mol);
    auto S = ints::sparse_integrals(world, overlap_e, bs_array);

    // Overlap ints
    auto kinetic_e = ints::make_1body_shr_pool("kinetic", basis, clustered_mol);
    auto T = ints::sparse_integrals(world, kinetic_e, bs_array);

    auto nuclear_e = ints::make_1body_shr_pool("nuclear", basis, clustered_mol);
    auto V = ints::sparse_integrals(world, nuclear_e, bs_array);

    decltype(T) H;
    H("i,j") = T("i,j") + V("i,j");

    const auto dfbs_array = utility::make_array(df_basis, df_basis);
    auto eri_e = ints::make_2body_shr_pool(df_basis, basis);

    auto three_c_array = utility::make_array(df_basis, basis, basis);
    decltype(H) M = ints::sparse_integrals(world, eri_e, dfbs_array);

    const auto schwarz_thresh = in.HasMember("schwarz threshold")
                                      ? in["schwarz threshold"].GetDouble()
                                      : 1e-12;
    if (world.rank() == 0) {
        std::cout << "Schwarz Threshold: " << schwarz_thresh << std::endl;
    }
    auto sbuilder = ints::init_schwarz_screen(schwarz_thresh);
    auto shr_screen = std::make_shared<ints::SchwarzScreen>(
          sbuilder(world, eri_e, df_basis, basis));

    auto eri3 = ints::direct_sparse_integrals(world, eri_e, three_c_array,
                                              shr_screen);

    // Trying my hand at making C for CADF.
    std::unordered_map<int, TA::TensorD> tiles;
    DArray<3, TA::TensorD, SpPolicy> C_df;
    {
        // OBS by atom
        auto sort_atoms_in_mol = false;
        std::vector<molecule::AtomBasedClusterable> atoms;
        std::unordered_map<std::size_t, std::size_t> obs_atom_to_cluster_map;
        auto cluster_ind = 0;
        auto atom_ind = 0;
        for (auto const &cluster : clustered_mol.clusterables()) {
            for (auto atom : cluster.atoms()) {
                atoms.push_back(std::move(atom));
                obs_atom_to_cluster_map[atom_ind] = cluster_ind;
                ++atom_ind;
            }
            ++cluster_ind;
        }
        basis::Basis obs_by_atom(bs.get_cluster_shells(
              molecule::Molecule(std::move(atoms), sort_atoms_in_mol)));

        // DFBS by atom
        atoms.clear(); // just in case the move didn't clear it
        for (auto const &cluster : df_clustered_mol.clusterables()) {
            for (auto atom : cluster.atoms()) {
                atoms.push_back(std::move(atom));
            }
        }
        basis::Basis df_by_atom(dfbs.get_cluster_shells(
              molecule::Molecule(std::move(atoms), sort_atoms_in_mol)));

        const auto dfbs_array = utility::make_array(df_by_atom, df_by_atom);
        auto eri_e = ints::make_2body_shr_pool(df_by_atom, obs_by_atom);

        auto three_c_array
              = utility::make_array(df_by_atom, obs_by_atom, obs_by_atom);
        decltype(H) M = ints::sparse_integrals(world, eri_e, dfbs_array);

        const auto schwarz_thresh = in.HasMember("schwarz threshold")
                                          ? in["schwarz threshold"].GetDouble()
                                          : 1e-12;
        auto sbuilder = ints::init_schwarz_screen(schwarz_thresh);
        auto shr_screen = std::make_shared<ints::SchwarzScreen>(
              sbuilder(world, eri_e, df_basis, basis));

        auto eri3 = ints::direct_sparse_integrals(world, eri_e, three_c_array,
                                                  shr_screen);

        auto const &eri3_tiles = eri3.array().trange().tiles();
        auto const &M_tiles = M.trange().tiles();

        auto same_center_tile = [&](unsigned long ind) {

            auto eri3_idx = std::array<unsigned long, 3>{{ind, ind, ind}};
            auto eri3_ord = eri3_tiles.ordinal(eri3_idx);
            TA::TensorD eri3_tile = eri3.array().find(eri3_ord).get();
            auto eri3_extent = eri3_tile.range().extent();
            MatrixD eri3_eig = TA::eigen_map(eri3_tile, eri3_extent[0],
                                             eri3_extent[1] * eri3_extent[2]);

            auto M_idx = std::array<unsigned long, 2>{{ind, ind}};
            auto M_ord = M_tiles.ordinal(M_idx);
            auto M_tile = M.find(M_ord).get();
            auto M_extent = M_tile.range().extent();
            MatrixD M_eig = TA::eigen_map(M_tile, M_extent[0], M_extent[1]);

            TA::TensorD out_tile(eri3_tile.range());

            MatrixD out_eig = M_eig.inverse() * eri3_eig;
            TA::eigen_map(out_tile, eri3_extent[0],
                          eri3_extent[1] * eri3_extent[2]) = out_eig;

            return std::make_pair<int, TA::TensorD>(eri3_ord,
                                                    std::move(out_tile));
        };

        // Loop over number of tiles
        for (auto i = 0; i < obs_by_atom.nclusters(); ++i) {

            // Deal with same tile mu nu.
            auto ord_tile_pair = same_center_tile(i);
            tiles[ord_tile_pair.first] = ord_tile_pair.second;


            for (auto j = 0; j < obs_by_atom.nclusters(); ++j) {
                if (j == i) continue;

                unsigned long il = i;
                unsigned long jl = j;

                auto eri3_idx_0 = std::array<unsigned long, 3>{{il, il, jl}};
                auto eri3_ord_0 = eri3_tiles.ordinal(eri3_idx_0);
                auto eri3_idx_1 = std::array<unsigned long, 3>{{jl, il, jl}};
                auto eri3_ord_1 = eri3_tiles.ordinal(eri3_idx_1);

                auto eri3_range_0
                      = eri3.array().trange().make_tile_range(eri3_idx_0);
                auto eri3_range_1
                      = eri3.array().trange().make_tile_range(eri3_idx_1);

                auto eri3_extent_0 = eri3_range_0.extent();
                auto eri3_extent_1 = eri3_range_1.extent();

                MatrixD eri3_eig_0;
                if (!eri3.array().is_zero(eri3_idx_0)) {
                    TA::TensorD eri3_tile_0
                          = eri3.array().find(eri3_ord_0).get();

                    eri3_eig_0
                          = TA::eigen_map(eri3_tile_0, eri3_extent_0[0],
                                          eri3_extent_0[1] * eri3_extent_0[2]);
                } else {
                    eri3_eig_0.resize(eri3_extent_0[0],
                                      eri3_extent_0[1] * eri3_extent_0[2]);
                    eri3_eig_0.setZero();
                }

                MatrixD eri3_eig_1;
                if (!eri3.array().is_zero(eri3_idx_1)) {
                    TA::TensorD eri3_tile_1
                          = eri3.array().find(eri3_ord_1).get();

                    eri3_eig_1
                          = TA::eigen_map(eri3_tile_1, eri3_extent_1[0],
                                          eri3_extent_1[1] * eri3_extent_1[2]);
                } else {
                    eri3_eig_1.resize(eri3_extent_1[0],
                                      eri3_extent_1[1] * eri3_extent_1[2]);
                    eri3_eig_1.setZero();
                }


                auto M_idx_11 = std::array<unsigned long, 2>{{il, il}};
                auto M_idx_12 = std::array<unsigned long, 2>{{il, jl}};
                auto M_idx_21 = std::array<unsigned long, 2>{{jl, il}};
                auto M_idx_22 = std::array<unsigned long, 2>{{jl, jl}};

                auto M_ord_11 = M_tiles.ordinal(M_idx_11);
                auto M_ord_12 = M_tiles.ordinal(M_idx_12);
                auto M_ord_21 = M_tiles.ordinal(M_idx_21);
                auto M_ord_22 = M_tiles.ordinal(M_idx_22);

                auto M_tile_11 = M.find(M_ord_11).get();
                auto M_tile_12 = M.find(M_ord_12).get();
                auto M_tile_21 = M.find(M_ord_21).get();
                auto M_tile_22 = M.find(M_ord_22).get();

                auto M_extent_11 = M_tile_11.range().extent();
                auto M_extent_12 = M_tile_12.range().extent();
                auto M_extent_21 = M_tile_21.range().extent();
                auto M_extent_22 = M_tile_22.range().extent();

                auto rows = M_extent_11[0] + M_extent_21[0];
                auto cols = M_extent_11[1] + M_extent_12[1];

                MatrixD M_combo(rows, cols);
                // Write M11
                M_combo.block(0, 0, M_extent_11[0], M_extent_11[1])
                      = TA::eigen_map(M_tile_11, M_extent_11[0],
                                      M_extent_11[1]);

                // Write M12
                M_combo.block(0, M_extent_11[1], M_extent_12[0], M_extent_12[1])
                      = TA::eigen_map(M_tile_12, M_extent_12[0],
                                      M_extent_12[1]);

                // Write M21
                M_combo.block(M_extent_11[0], 0, M_extent_21[0], M_extent_21[1])
                      = TA::eigen_map(M_tile_21, M_extent_21[0],
                                      M_extent_21[1]);

                // Write M22
                M_combo.block(M_extent_11[0], M_extent_11[1], M_extent_22[0],
                              M_extent_22[1])
                      = TA::eigen_map(M_tile_22, M_extent_22[0],
                                      M_extent_22[1]);

                MatrixD M_combo_inv = M_combo.inverse();

                // Doing the block wise GEMM by hand for now.
                auto block11
                      = M_combo_inv.block(0, 0, M_extent_11[0], M_extent_11[1]);
                auto block12
                      = M_combo_inv.block(0, M_extent_11[1], M_extent_12[0],
                                          M_extent_12[1]);
                auto block21
                      = M_combo_inv.block(M_extent_11[0], 0, M_extent_21[0],
                                          M_extent_21[1]);
                auto block22
                      = M_combo_inv.block(M_extent_11[0], M_extent_11[1],
                                          M_extent_22[0], M_extent_22[1]);

                MatrixD C0 = block11 * eri3_eig_0 + block12 * eri3_eig_1;
                MatrixD C1 = block21 * eri3_eig_0 + block22 * eri3_eig_1;

                TA::TensorD out_tile_0(eri3_range_0);
                TA::TensorD out_tile_1(eri3_range_1);


                TA::eigen_map(out_tile_0, eri3_extent_0[0],
                              eri3_extent_0[1] * eri3_extent_0[2]) = C0;

                TA::eigen_map(out_tile_1, eri3_extent_1[0],
                              eri3_extent_1[1] * eri3_extent_1[2]) = C1;

                if (out_tile_0.norm() / out_tile_0.range().volume()
                    > TA::SparseShape<float>::threshold()) {
                    tiles[eri3_ord_0] = out_tile_0;
                }
                if (out_tile_1.norm() / out_tile_1.range().volume()
                    > TA::SparseShape<float>::threshold()) {
                    tiles[eri3_ord_1] = out_tile_1;
                }
            }
        }

        TA::Tensor<float> C_shape(eri3.array().trange().tiles(), 0.0);
        for (auto &pair : tiles) {
            *(C_shape.data() + pair.first) = pair.second.norm();
        }

        TA::SparseShape<float> shape(C_shape, eri3.array().trange());
        DArray<3, TA::TensorD, SpPolicy> C_df_(world, eri3.array().trange(),
                                               shape);

        for (auto &pair : tiles) {
            C_df_.set(pair.first, pair.second);
        }
        world.gop.fence();
        C_df_.truncate();

        std::cout << "C df sparsity (by atom) = "
                  << C_df_.get_shape().sparsity() << std::endl;

        auto by_atom_trange = integrals::detail::create_trange(three_c_array);
        auto by_cluster_trange = integrals::detail::create_trange(
              utility::make_array(df_basis, basis, basis));
        std::cout << "Trange by atom: " << by_atom_trange << std::endl;
        std::cout << "Trange by cluster: " << by_cluster_trange << std::endl;

        C_df = scf::reblock_from_atoms(C_df_, obs_atom_to_cluster_map,
                                       obs_atom_to_cluster_map,
                                       by_cluster_trange);

        std::cout << "C df sparsity (by cluster) = "
                  << C_df.get_shape().sparsity() << std::endl;
    }

    {
        decltype(M) M_inv_oh;
        auto M_eig = TA::array_to_eigen(M);
        Eig::SelfAdjointEigenSolver<MatrixD> es(M_eig);
        auto M_eig_inv_oh = es.operatorInverseSqrt();

        M_inv_oh = array_ops::eigen_to_array<TA::TensorD>(
              world, M_eig_inv_oh, M.trange().data()[0], M.trange().data()[1]);

        auto F_soad = scf::fock_from_soad_low_mem(world, clustered_mol, basis,
                                                  eri_e, H);

        auto multi_pool
              = ints::make_1body_shr_pool("emultipole2", basis, clustered_mol);

        auto r_xyz = ints::sparse_xyz_integrals(world, multi_pool, bs_array);

        const auto clr_threshold = in.HasMember("clr threshold")
                                         ? in["clr threshold"].GetDouble()
                                         : 1e-6;
        if (world.rank() == 0) {
            std::cout << "CLR threshold = " << clr_threshold << std::endl;
        }

        auto deri3 = ints::direct_sparse_integrals(
              world, eri_e, three_c_array, shr_screen,
              tensor::TaToDecompTensor(clr_threshold));

        auto dC_df = TA::to_new_tile_type(C_df, tensor::TaToDecompTensor(
                                                      clr_threshold));

        // Don't CLR compress M
        auto dM = TA::to_new_tile_type(
              M, tensor::TaToDecompTensor(clr_threshold, false));

        decltype(dC_df) dG_df;
        world.gop.fence();
        auto old_compress = tensor::detail::recompress;
        tensor::detail::recompress = true;
        dG_df("X, mu, nu")
              = (deri3("X, mu, nu") - 0.5 * dM("X,Y") * dC_df("Y,mu,nu"));
        ta_routines::minimize_storage(dG_df, clr_threshold);
        world.gop.fence();
        tensor::detail::recompress = old_compress;

        auto g_store = utility::array_storage(dG_df);
        if (world.rank() == 0) {
            std::cout << "G_df storage = \n"
                      << "\tFull    " << g_store[0] << "\n"
                      << "\tSparse  " << g_store[1] << "\n"
                      << "\tCLR     " << g_store[2] << "\n" << std::flush;
        }
        world.gop.fence();


        scf::ClosedShellSCF hf;

        double TcutC = 0.0;
        if (in.HasMember("TcutC")) {
            TcutC = in["TcutC"].GetDouble();
        }

        if (in.HasMember("use forced shape")
            && in["use forced shape"].GetBool()) {

            if (world.rank() == 0) {
                std::cout << "Using forced shape build\n";
            }

            scf::CADFForcedShapeFockBuilder<decltype(eri3)> cadf_builder(
                  M, eri3, dC_df, dG_df, clr_threshold);

            if (in.HasMember("coulomb sparse threshold")) {
                cadf_builder.set_J_sparse_thresh(
                      in["coulomb sparse threshold"].GetDouble());

                if (world.rank() == 0) {
                    std::cout << "Changing J threshold to "
                              << in["coulomb sparse threshold"].GetDouble()
                              << std::endl;
                }
            }

            hf = scf::ClosedShellSCF(H, S, occ / 2, repulsion_energy,
                                     r_xyz, std::move(cadf_builder), F_soad, TcutC);
        } else {
            scf::CADFFockBuilder<decltype(eri3)> cadf_builder(
                  M, eri3, dC_df, dG_df, clr_threshold);

            if (in.HasMember("coulomb sparse threshold")) {
                cadf_builder.set_J_sparse_thresh(
                      in["coulomb sparse threshold"].GetDouble());

                if (world.rank() == 0) {
                    std::cout << "Changing J threshold to "
                              << in["coulomb sparse threshold"].GetDouble()
                              << std::endl;
                }
            }

            hf = scf::ClosedShellSCF(H, S, occ / 2, repulsion_energy,
                                     r_xyz, std::move(cadf_builder), F_soad, TcutC);
        }

        world.gop.fence();

        hf.solve(in["scf max iter"].GetInt(),
                 in["scf convergence threshold"].GetDouble());
    }

    madness::finalize();
    return 0;
}
