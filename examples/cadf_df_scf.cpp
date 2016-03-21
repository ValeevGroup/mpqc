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
#include <rapidjson/prettywriter.h>

#include "../utility/json_handling.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/basis.h"

#include "../integrals/integrals.h"

#include "../utility/time.h"
#include "../utility/array_info.h"

#include "../ta_routines/array_to_eigen.h"
#include "../ta_routines/minimize_storage.h"
#include "../ta_routines/diagonal_array.h"
#include "../ta_routines/cholesky_inverse.h"

#include "../scf/diagonalize_for_coffs.hpp"
#include "../scf/soad.h"
#include "../scf/orbital_localization.h"
#include "../scf/clusterd_coeffs.h"
#include "../scf/cadf_helper_functions.h"
#include "../scf/cadf_fitting_coeffs.h"
#include "../scf/scf.h"

#include "../scf/eigen_solve_density_builder.h"
#include "../scf/purification_density_build.h"

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
    json::parse_input(argc, argv, in);

    auto owriter = json::init_json_writer(in);
    auto &out_doc = owriter->doc();

    double ta_threshold = in.HasMember("sparse threshold")
                                ? in["sparse threshold"].GetDouble()
                                : 1e-11;

    out_doc.AddMember("ta sparse threshold", ta_threshold,
                      out_doc.GetAllocator());

    TiledArray::SparseShape<float>::threshold(ta_threshold);
    if (world.rank() == 0) {
        std::cout << "TA::Sparse Threshold: " << ta_threshold << std::endl;
    }

    if (in.HasMember("engine precision")) {
        integrals::detail::integral_engine_precision
              = in["engine precision"].GetDouble();
    }

    out_doc.AddMember("integral engine precision",
                      integrals::detail::integral_engine_precision,
                      out_doc.GetAllocator());

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

    out_doc.AddMember("number of atoms", int(atom_mol.atoms().size()),
                      out_doc.GetAllocator());

    auto nclusters = 0;
    molecule::Molecule clustered_mol;
    if (in.HasMember("cluster by atom") && in["cluster by atom"].GetBool()) {
        clustered_mol = std::move(atom_mol);
        nclusters = clustered_mol.nclusters();
    } else if (in.HasMember("number of clusters")) {
        nclusters = in["number of clusters"].GetInt();
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

    out_doc.AddMember("number of obs clusters", clustered_mol.nclusters(),
                      out_doc.GetAllocator());

    auto repulsion_energy = clustered_mol.nuclear_repulsion();
    if (world.rank() == 0) {
        std::cout << "Nuclear Repulsion Energy: " << repulsion_energy
                  << std::endl;
    }
    auto occ = clustered_mol.occupation(0);


    std::string basis_name = in["obs name"].GetString();
    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));
    if (world.rank() == 0) {
        std::cout << "Basis has " << basis.nfunctions() << " functions\n";
    }

    out_doc.AddMember("obs basis",
                      rapidjson::Value(basis_name.c_str(), basis_name.size(),
                                       out_doc.GetAllocator()),
                      out_doc.GetAllocator());

    out_doc.AddMember("number of obs basis functions", basis.nfunctions(),
                      out_doc.GetAllocator());

    auto df_clustered_mol = clustered_mol;

    out_doc.AddMember("number of dfbs clusters", df_clustered_mol.nclusters(),
                      out_doc.GetAllocator());

    std::string df_basis_name = in["df name"].GetString();
    basis::BasisSet dfbs(df_basis_name);
    basis::Basis df_basis(dfbs.get_cluster_shells(df_clustered_mol));
    if (world.rank() == 0) {
        std::cout << "DF Basis has " << df_basis.nfunctions() << " functions\n";
    }

    out_doc.AddMember("df basis", rapidjson::Value(df_basis_name.c_str(),
                                                   df_basis_name.size(),
                                                   out_doc.GetAllocator()),
                      out_doc.GetAllocator());

    out_doc.AddMember("number of df basis functions", df_basis.nfunctions(),
                      out_doc.GetAllocator());

    libint2::init();
    const auto bs_array = utility::make_array(basis, basis);

    // Overlap ints
    auto s0 = mpqc_time::fenced_now(world);
    auto overlap_e = ints::make_1body_shr_pool("overlap", basis, clustered_mol);
    auto S = ints::sparse_integrals(world, overlap_e, bs_array);
    auto s1 = mpqc_time::fenced_now(world);
    auto stime = mpqc_time::duration_in_s(s0, s1);
    if (world.rank() == 0) {
        std::cout << "Overlap time: " << stime << std::endl;
    }
    out_doc.AddMember("overlap time", stime, out_doc.GetAllocator());

    auto h0 = mpqc_time::fenced_now(world);
    auto kinetic_e = ints::make_1body_shr_pool("kinetic", basis, clustered_mol);
    auto T = ints::sparse_integrals(world, kinetic_e, bs_array);

    auto nuclear_e = ints::make_1body_shr_pool("nuclear", basis, clustered_mol);
    auto V = ints::sparse_integrals(world, nuclear_e, bs_array);

    decltype(T) H;
    H("i,j") = T("i,j") + V("i,j");
    auto h1 = mpqc_time::fenced_now(world);
    auto htime = mpqc_time::duration_in_s(h0, h1);
    if (world.rank() == 0) {
        std::cout << "Hcore time: " << stime << std::endl;
    }
    out_doc.AddMember("hcore time", htime, out_doc.GetAllocator());

    const auto dfbs_array = utility::make_array(df_basis, df_basis);
    auto eri_e = ints::make_2body_shr_pool(df_basis, basis);

    auto three_c_array = utility::make_array(df_basis, basis, basis);
    auto m0 = mpqc_time::fenced_now(world);
    decltype(H) M = ints::sparse_integrals(world, eri_e, dfbs_array);
    auto m1 = mpqc_time::fenced_now(world);
    auto mtime = mpqc_time::duration_in_s(m0, m1);
    if (world.rank() == 0) {
        std::cout << "Metric time: " << mtime << std::endl;
    }
    out_doc.AddMember("metric time", mtime, out_doc.GetAllocator());

    const auto schwarz_thresh = in.HasMember("schwarz threshold")
                                      ? in["schwarz threshold"].GetDouble()
                                      : 1e-12;
    if (world.rank() == 0) {
        std::cout << "Schwarz Threshold: " << schwarz_thresh << std::endl;
    }
    out_doc.AddMember("schwarz thresh", schwarz_thresh, out_doc.GetAllocator());

    auto ss0 = mpqc_time::fenced_now(world);
    auto sbuilder = ints::init_schwarz_screen(schwarz_thresh);
    auto shr_screen = std::make_shared<ints::SchwarzScreen>(
          sbuilder(world, eri_e, df_basis, basis));
    auto ss1 = mpqc_time::fenced_now(world);
    auto sstime = mpqc_time::duration_in_s(ss0, ss1);
    if (world.rank() == 0) {
        std::cout << "Screener time: " << sstime << std::endl;
    }

    auto cdf0 = mpqc_time::fenced_now(world);
    std::unordered_map<std::size_t, std::size_t> obs_atom_to_cluster_map;
    std::unordered_map<std::size_t, std::size_t> dfbs_atom_to_cluster_map;
    auto C_df_ = scf::compute_atomic_fitting_coeffs(world, clustered_mol,
                                                    df_clustered_mol, bs, dfbs,
                                                    obs_atom_to_cluster_map,
                                                    dfbs_atom_to_cluster_map);

    auto cdf1 = mpqc_time::fenced_now(world);
    auto cdftime = mpqc_time::duration_in_s(cdf0, cdf1);
    if (world.rank() == 0) {
        std::cout << "Cdf time: " << cdftime << std::endl;
    }

    out_doc.AddMember("Cdf build time", cdftime, out_doc.GetAllocator());

    auto array_storage_Cdf = utility::array_storage(C_df_);
    if (world.rank() == 0) {
        std::cout << "C_df by atom storage = \n"
                  << "\tFull   " << array_storage_Cdf[0] << "\n"
                  << "\tSparse " << array_storage_Cdf[1] << "\n"
                  << "\tCLR    " << array_storage_Cdf[2] << "\n";
    }
    out_doc.AddMember("Cdf By Atom Full Storage", array_storage_Cdf[0],
                      out_doc.GetAllocator());

    out_doc.AddMember("Cdf By Atom Sparse Storage", array_storage_Cdf[1],
                      out_doc.GetAllocator());

    out_doc.AddMember("Cdf By Atom CLR Storage", array_storage_Cdf[2],
                      out_doc.GetAllocator());

    auto by_cluster_trange = integrals::detail::create_trange(
          utility::make_array(df_basis, basis, basis));


    decltype(C_df_) C_df;
    if (in.HasMember("cluster by atom") && in["cluster by atom"].GetBool()) {
        C_df = C_df_;
    } else {
        auto reblock0 = mpqc_time::fenced_now(world);
        C_df = scf::reblock_from_atoms(C_df_, obs_atom_to_cluster_map,
                                       dfbs_atom_to_cluster_map,
                                       by_cluster_trange);

        auto reblock1 = mpqc_time::fenced_now(world);
        auto reblock_time = mpqc_time::duration_in_s(reblock0, reblock1);
        if (world.rank() == 0) {
            std::cout << "Reblock C_df time " << reblock_time << std::endl;
        }

        out_doc.AddMember("Cdf reblock time", reblock_time,
                          out_doc.GetAllocator());
    }

    array_storage_Cdf = utility::array_storage(C_df);
    if (world.rank() == 0) {
        std::cout << "C_df storage = \n"
                  << "\tFull   " << array_storage_Cdf[0] << "\n"
                  << "\tSparse " << array_storage_Cdf[1] << "\n"
                  << "\tCLR    " << array_storage_Cdf[2] << "\n";
    }
    out_doc.AddMember("Cdf Full Storage", array_storage_Cdf[0],
                      out_doc.GetAllocator());

    out_doc.AddMember("Cdf Sparse Storage", array_storage_Cdf[1],
                      out_doc.GetAllocator());

    out_doc.AddMember("Cdf CLR Storage", array_storage_Cdf[2],
                      out_doc.GetAllocator());

    // Begin scf
    auto soad0 = mpqc_time::fenced_now(world);
    decltype(S) F_soad;
    if (in.HasMember("Soad Method")
        && std::string(in["Soad Method"].GetString()) == "low memory") {
        F_soad = scf::fock_from_soad_low_mem(world, clustered_mol, basis, eri_e,
                                             H);
        out_doc.AddMember("Soad Method", "Low memory", out_doc.GetAllocator());
    } else {
        F_soad = scf::fock_from_soad(world, clustered_mol, basis, eri_e, H);
        out_doc.AddMember("Soad Method", "High memory", out_doc.GetAllocator());
    }
    auto soad1 = mpqc_time::fenced_now(world);
    auto soadtime = mpqc_time::duration_in_s(soad0, soad1);
    if (world.rank() == 0) {
        std::cout << "SOAD time: " << soadtime << std::endl;
    }
    out_doc.AddMember("Soad Time", soadtime, out_doc.GetAllocator());

    auto rxyz0 = mpqc_time::fenced_now(world);
    auto multi_pool
          = ints::make_1body_shr_pool("emultipole2", basis, clustered_mol);

    auto r_xyz = ints::sparse_xyz_integrals(world, multi_pool, bs_array);
    auto rxyz1 = mpqc_time::fenced_now(world);
    auto rxyztime = mpqc_time::duration_in_s(rxyz0, rxyz1);
    if (world.rank() == 0) {
        std::cout << "dipole integrals time: " << rxyztime << std::endl;
    }

    const auto clr_threshold = in.HasMember("clr threshold")
                                     ? in["clr threshold"].GetDouble()
                                     : 1e-6;
    if (world.rank() == 0) {
        std::cout << "CLR threshold = " << clr_threshold << std::endl;
    }

    out_doc.AddMember("CLR threshold", clr_threshold, out_doc.GetAllocator());

    auto deri3 = ints::direct_sparse_integrals(
          world, eri_e, three_c_array, shr_screen,
          tensor::TaToDecompTensor(clr_threshold));

    auto dC_df
          = TA::to_new_tile_type(C_df, tensor::TaToDecompTensor(clr_threshold));

    // Don't CLR compress M
    auto dM = TA::to_new_tile_type(M, tensor::TaToDecompTensor(clr_threshold,
                                                               false));

    decltype(dC_df) dG_df;
    auto G_df0 = mpqc_time::fenced_now(world);
    auto old_compress = tensor::detail::recompress;
    tensor::detail::recompress = true;
    dG_df("X, mu, nu")
          = (deri3("X, mu, nu") - 0.5 * dM("X,Y") * dC_df("Y,mu,nu"));
    ta_routines::minimize_storage(dG_df, clr_threshold);
    auto G_df1 = mpqc_time::fenced_now(world);
    tensor::detail::recompress = old_compress;
    auto G_df_time = mpqc_time::duration_in_s(rxyz0, rxyz1);
    if (world.rank() == 0) {
        std::cout << "G_df time: " << G_df_time << std::endl;
    }
    out_doc.AddMember("G_df time", G_df_time, out_doc.GetAllocator());

    auto g_store = utility::array_storage(dG_df);
    if (world.rank() == 0) {
        std::cout << "G_df storage = \n"
                  << "\tFull    " << g_store[0] << "\n"
                  << "\tSparse  " << g_store[1] << "\n"
                  << "\tCLR     " << g_store[2] << "\n" << std::flush;
    }
    world.gop.fence();
    out_doc.AddMember("G_df Full Storage", g_store[0], out_doc.GetAllocator());

    out_doc.AddMember("G_df Sparse Storage", g_store[1],
                      out_doc.GetAllocator());

    out_doc.AddMember("G_df CLR Storage", g_store[2], out_doc.GetAllocator());

    double TcutC = 0.0;
    if (in.HasMember("TcutC")) {
        TcutC = in["TcutC"].GetDouble();
    }

    out_doc.AddMember("TcutC", TcutC, out_doc.GetAllocator());

    auto localize = true;
    auto ebuilder = scf::ESolveDensityBuilder(S, r_xyz, occ / 2, nclusters,
                                              TcutC, "inverse sqrt", localize);

    std::unique_ptr<scf::DensityBuilder> d_builder
          = make_unique<decltype(ebuilder)>(std::move(ebuilder));

    auto eri3 = ints::direct_sparse_integrals(world, eri_e, three_c_array,
                                              shr_screen);

    std::unique_ptr<scf::FockBuilder> f_builder;
    if (in.HasMember("use forced shape") && in["use forced shape"].GetBool()) {

        if (world.rank() == 0) {
            std::cout << "Using forced shape build\n";
        }

        scf::CADFForcedShapeFockBuilder<decltype(eri3)> forced_shape(
              M, eri3, dC_df, dG_df, clr_threshold);

        if (in.HasMember("coulomb sparse threshold")) {
            forced_shape.set_J_sparse_thresh(
                  in["coulomb sparse threshold"].GetDouble());
        }
        f_builder = std::unique_ptr<scf::FockBuilder>(
              make_unique<decltype(forced_shape)>(std::move(forced_shape)));

        out_doc.AddMember("Using forced shape", true, out_doc.GetAllocator());

    } else {
        scf::CADFFockBuilder<decltype(eri3)> cadf_builder(M, eri3, dC_df, dG_df,
                                                          clr_threshold);

        if (in.HasMember("coulomb sparse threshold")) {
            cadf_builder.set_J_sparse_thresh(
                  in["coulomb sparse threshold"].GetDouble());
        }
        f_builder = std::unique_ptr<scf::FockBuilder>(
              make_unique<decltype(cadf_builder)>(std::move(cadf_builder)));

        out_doc.AddMember("Using forced shape", false, out_doc.GetAllocator());
    }

    auto hf = scf::ClosedShellSCF(H, S, repulsion_energy, std::move(f_builder),
                                  std::move(d_builder), F_soad);

    world.gop.fence();

    auto converged = hf.solve(in["scf max iter"].GetInt(),
                              in["scf convergence threshold"].GetDouble());

    out_doc.AddMember("SCF Converged", converged, out_doc.GetAllocator());
    out_doc.AddMember("SCF", hf.results(out_doc), out_doc.GetAllocator());

    if (in.HasMember("compute mp2") && in["compute mp2"].GetBool()) {
        decltype(S) D = hf.density();
        decltype(S) Ci = hf.coefficents();
        decltype(S) F = hf.fock();
        decltype(S) Cv;
        {
            auto F_eig = array_ops::array_to_eigen(F);
            auto S_eig = array_ops::array_to_eigen(S);

            auto mp2_vir = basis.nfunctions() - occ / 2;
            Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                       S_eig);
            decltype(S_eig) Ceigv = es.eigenvectors().rightCols(mp2_vir);
            auto tr_vir = scf::tr_occupied(clustered_mol.nclusters(), mp2_vir);

            Cv = array_ops::eigen_to_array<TA::TensorD>(
                  world, Ceigv, S.trange().data()[0], tr_vir);
        }

        auto I = array_ops::create_diagonal_matrix(S, 1.0);
        decltype(S) Q;
        Q("i,j") = I("i,j") - D("i,k") * S("k,j");
        Q.truncate();
        Q = Cv;

        decltype(S) Focc, Fbar, Sbar;
        Fbar("i,j") = Q("k,i") * F("k,l") * Q("l,j");
        Sbar("i,j") = Q("k,i") * S("k,l") * Q("l,j");
        Focc("i,j") = Ci("mu, i") * F("mu, nu") * Ci("nu, j");
        Fbar.truncate();
        Sbar.truncate();
        Focc.truncate();

        decltype(C_df) W;
        W("Y, i, sig") = ((C_df("Y, mu, nu") * Ci("nu, i")) * Q("mu, sig"));
        W.truncate();

        auto W_store = utility::array_storage(W);
        if (world.rank() == 0) {
            std::cout << "W sizes:\n"
                      << "\tFull   = " << W_store[0] << "\n"
                      << "\tSparse = " << W_store[1] << "\n"
                      << "\tCLR    = " << W_store[2] << "\n";
        }

        decltype(C_df) G;
        G("i, s, j, r") = W("X,i,s") * (M("X,Y") * W("Y, j, r"));
        G.truncate();

        auto G_store = utility::array_storage(G);
        if (world.rank() == 0) {
            std::cout << "G sizes:\n"
                      << "\tFull   = " << G_store[0] << "\n"
                      << "\tSparse = " << G_store[1] << "\n"
                      << "\tCLR    = " << G_store[2] << "\n";
        }

        // Begin MP2
        Eig::VectorXd F_occ_diag, F_pao_diag;
        {
            auto F_occ_eig = TA::array_to_eigen(Focc);
            auto F_pao_eig = TA::array_to_eigen(Fbar);
            F_occ_diag = F_occ_eig.diagonal();
            F_pao_diag = F_pao_eig.diagonal();
        }

        auto r_to_T = [&](TA::TensorD &t) {
            const auto start = t.range().lobound();
            const auto end = t.range().upbound();
            const auto extent = t.range().extent();

            for (auto i = 0; i < extent[0]; ++i) {
                const auto f_ii = F_occ_diag[i + start[0]];
                const auto i_ind = i * extent[1] * extent[2] * extent[3];

                for (auto r = 0; r < extent[1]; ++r) {
                    const auto f_rr = -F_pao_diag[r + start[1]];
                    const auto ir_ind = i_ind + r * extent[2] * extent[3];

                    for (auto j = 0; j < extent[2]; ++j) {
                        const auto f_jj = F_occ_diag[j + start[2]];
                        const auto irj_ind = ir_ind + j * extent[3];

                        for (auto s = 0; s < extent[3]; ++s) {
                            const auto f_ss = -F_pao_diag[s + start[3]];
                            const auto ind = irj_ind + s;

                            auto ptr = t.data() + ind;
                            auto val = *ptr;
                            *ptr = val / (f_ii + f_jj + f_ss + f_rr);
                        }
                    }
                }
            }

            return t.norm();
        };

        decltype(G) T, R;
        R("i, p, j, q") = G("i,p,j,q");
        R.truncate();

        double norm = R("i,p,j,q").norm();

        TA::foreach_inplace(R, r_to_T);
        T("i,p,j,q") = R("i,p,j,q");
        T.truncate();

        double energy
              = TA::dot(2 * G("i,p,j,q") - G("i, q, j, p"), T("i,p,j,q"));

        if (world.rank() == 0) {
            std::cout << "Energy after initialization: " << energy << std::endl;
            std::cout << "R norm after initialization: " << norm << std::endl;
        }

        for (auto i = 1; i <= 30 && norm > 1e-7; ++i) {
            R("i, p, j, q") = G("i, p, j, q")
                              + Fbar("p, r") * T("i, r, j, s") * Sbar("s, q")
                              + Sbar("p, r") * T("i, r, j, s") * Fbar("s, q")
                              - Sbar("p, r")
                                * (Focc("i, k") * T("k, r, j, s")
                                   + T("i, r, k, s") * Focc("k, j"))
                                * Sbar("s, q");
            R.truncate();

            norm = R("i,p,j,q").norm();

            TA::foreach_inplace(R, r_to_T);
            T("i, p, j, q") = T("i, p, j, q") + R("i, p, j, q");
            T.truncate();

            energy = TA::dot(2 * G("i,p,j,q") - G("i, q, j, p"), T("i,p,j,q"));

            if (world.rank() == 0) {
                std::cout << i << " Energy: " << energy << std::endl;
                std::cout << i << " R norm: " << norm << std::endl;
            }
        }
    }

    if (world.rank() == 0) {
        if (in.HasMember("output file")) {
            std::string out_file_name = in["output file"].GetString();
            std::ofstream out_file(out_file_name, std::ofstream::out);
            owriter->finalize_and_print(out_file);
        } else {
            owriter->finalize_and_print(std::cout);
        }
    }

    madness::finalize();
    return 0;
}
