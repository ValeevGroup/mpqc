#include "../common/namespaces.h"
#include "../common/typedefs.h"

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
#include "../integrals/atomic_integral.h"

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
#include "../scf/linear_cadf_builder.h"

#include "../tensor/tensor_transforms.h"

#include "../f12/f12_utility.h"

using namespace mpqc;
namespace ints = mpqc::integrals;

// Causes general recompression in clr based tensor expressions
bool tensor::detail::recompress = false;

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

  tensor::detail::recompress =
      in.HasMember("recompress") ? in["recompress"].GetBool() : false;

  out_doc.AddMember("ta sparse threshold", ta_threshold,
                    out_doc.GetAllocator());

  TiledArray::SparseShape<float>::threshold(ta_threshold);
  if (world.rank() == 0) {
    std::cout << "TA::Sparse Threshold: " << ta_threshold << std::endl;
  }

  if (in.HasMember("engine precision")) {
    integrals::detail::integral_engine_precision =
        in["engine precision"].GetDouble();
  }

  out_doc.AddMember("integral engine precision",
                    integrals::detail::integral_engine_precision,
                    out_doc.GetAllocator());

  std::string mol_file_name;
  if (!in.HasMember("xyz file")) {
    if (world.rank() == 0) {
      std::cout << "Did not detect xyz file in input.\n";
    }
    madness::finalize();
    return 0;
  } else {
    mol_file_name = in["xyz file"].GetString();
  }

  std::fstream mol_file(mol_file_name);
  auto atom_mol = molecule::Molecule(mol_file, false);
  mol_file.close();

  out_doc.AddMember("number of atoms", int(atom_mol.atoms().size()),
                    out_doc.GetAllocator());

  auto nclusters = 0;
  molecule::Molecule clustered_mol;
  if (in.HasMember("cluster by atom") && in["cluster by atom"].GetBool()) {
    clustered_mol = std::move(atom_mol);
    nclusters = clustered_mol.nclusters();
  } else if (in.HasMember("number of clusters")) {
    nclusters = in["number of clusters"].GetInt();
    clustered_mol = molecule::attach_hydrogens_and_kmeans(
        atom_mol.clusterables(), nclusters);
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
    std::cout << "Nuclear Repulsion Energy: " << repulsion_energy << std::endl;
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

  out_doc.AddMember(
      "df basis", rapidjson::Value(df_basis_name.c_str(), df_basis_name.size(),
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

  decltype(H) Mj = ints::sparse_integrals(world, eri_e, dfbs_array);

  auto m1 = mpqc_time::fenced_now(world);
  auto mtime = mpqc_time::duration_in_s(m0, m1);
  if (world.rank() == 0) {
    std::cout << "Metric time: " << mtime << std::endl;
  }
  out_doc.AddMember("metric time", mtime, out_doc.GetAllocator());
  auto Mj_store = utility::array_storage(Mj);

  if (world.rank() == 0) {
    std::cout << "(X|Y) storage:\n"
              << "\tFull   " << Mj_store[0] << "\n"
              << "\tSparse " << Mj_store[1] << "\n"
              << "\tCLR    " << Mj_store[2] << "\n";
  }

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

  decltype(Mj) Mk;
  TA::DistArray<TA::TensorD, SpPolicy> C_df_;
  auto cdf0 = mpqc_time::fenced_now(world);
  std::unordered_map<std::size_t, std::size_t> obs_atom_to_cluster_map;
  std::unordered_map<std::size_t, std::size_t> dfbs_atom_to_cluster_map;
  if (in.HasMember("Metric")) {
    if (in["Metric"].GetInt() == 0) {
      Mk = Mj;
      C_df_ = scf::compute_atomic_fitting_coeffs(
          world, clustered_mol, df_clustered_mol, bs, dfbs, eri_e,
          obs_atom_to_cluster_map, dfbs_atom_to_cluster_map);
    } else if (in["Metric"].GetInt() == 1) {
      auto overlap_eM =
          ints::make_1body_shr_pool("overlap", df_basis, clustered_mol);

      Mk = ints::sparse_integrals(world, overlap_eM, dfbs_array);

      C_df_ = scf::compute_atomic_fitting_coeffs(
          world, clustered_mol, df_clustered_mol, bs, dfbs, overlap_eM,
          obs_atom_to_cluster_map, dfbs_atom_to_cluster_map);
// Needs to be fixed to work with new engine. 
   //  } else if (in["Metric"].GetInt() == 2) {
   //    auto cgtg_e = ints::make_2body_cGTG_shr_pool(df_basis);

   //    Mk = ints::sparse_integrals(world, cgtg_e, dfbs_array);

   //    C_df_ = scf::compute_atomic_fitting_coeffs(
   //        world, clustered_mol, df_clustered_mol, bs, dfbs, cgtg_e,
   //        obs_atom_to_cluster_map, dfbs_atom_to_cluster_map);
   //  } else if (in["Metric"].GetInt() == 3) {
   //    auto cgtgC_e = ints::make_2body_cGTG_C_shr_pool(df_basis);

   //    Mk = ints::sparse_integrals(world, cgtgC_e, dfbs_array);

   //    C_df_ = scf::compute_atomic_fitting_coeffs(
   //        world, clustered_mol, df_clustered_mol, bs, dfbs, cgtgC_e,
   //        obs_atom_to_cluster_map, dfbs_atom_to_cluster_map);
    } else {
      Mk = Mj;
      C_df_ = scf::compute_atomic_fitting_coeffs(
          world, clustered_mol, df_clustered_mol, bs, dfbs, eri_e,
          obs_atom_to_cluster_map, dfbs_atom_to_cluster_map);
    }
    if (world.rank() == 0) {
      if (in["Metric"].GetInt() == 0) {
        std::cout << "Using Coulomb Metric for Cdf." << std::endl;
      } else if (in["Metric"].GetInt() == 1) {
        std::cout << "Using Overlap Metric for Cdf." << std::endl;
      } else if (in["Metric"].GetInt() == 2) {
        std::cout << "Using cGTG Metric for Cdf." << std::endl;
      } else if (in["Metric"].GetInt() == 3) {
        std::cout << "Using cGTG_times Coulomb Metric for Cdf." << std::endl;
      } else {
        std::cout << "Didn't detect an integer between 0 and 3, using Coulomb "
                     "Metric for Cdf." << std::endl;
      }
    }
  } else {
    C_df_ = scf::compute_atomic_fitting_coeffs(
        world, clustered_mol, df_clustered_mol, bs, dfbs, eri_e,
        obs_atom_to_cluster_map, dfbs_atom_to_cluster_map);
  }
  if (!in.HasMember("test m") || in["test m"].GetBool() == false) {
    Mk = Mj;
  }

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
                                   dfbs_atom_to_cluster_map, by_cluster_trange);

    auto reblock1 = mpqc_time::fenced_now(world);
    auto reblock_time = mpqc_time::duration_in_s(reblock0, reblock1);
    if (world.rank() == 0) {
      std::cout << "Reblock C_df time " << reblock_time << std::endl;
    }

    out_doc.AddMember("Cdf reblock time", reblock_time, out_doc.GetAllocator());
  }

  // Begin scf
  auto soad0 = mpqc_time::fenced_now(world);
  decltype(S) F_soad =
      scf::fock_from_soad(world, clustered_mol, basis, eri_e, H);
  auto soad1 = mpqc_time::fenced_now(world);
  auto soadtime = mpqc_time::duration_in_s(soad0, soad1);
  if (world.rank() == 0) {
    std::cout << "SOAD time: " << soadtime << std::endl;
  }
  out_doc.AddMember("Soad Time", soadtime, out_doc.GetAllocator());

  auto rxyz0 = mpqc_time::fenced_now(world);
  auto multi_pool =
      ints::make_1body_shr_pool("emultipole1", basis, clustered_mol);

  auto r_xyz = ints::sparse_xyz_integrals(world, multi_pool, bs_array);
  auto rxyz1 = mpqc_time::fenced_now(world);
  auto rxyztime = mpqc_time::duration_in_s(rxyz0, rxyz1);
  if (world.rank() == 0) {
    std::cout << "dipole integrals time: " << rxyztime << std::endl;
  }

  const auto clr_threshold =
      in.HasMember("clr threshold") ? in["clr threshold"].GetDouble() : 1e-6;
  if (world.rank() == 0) {
    std::cout << "CLR threshold = " << clr_threshold << std::endl;
  }

  out_doc.AddMember("CLR threshold", clr_threshold, out_doc.GetAllocator());

  auto deri3 =
      ints::direct_sparse_integrals(world, eri_e, three_c_array, shr_screen,
                                    tensor::TaToDecompTensor(clr_threshold));

  auto dC_df =
      TA::to_new_tile_type(C_df, tensor::TaToDecompTensor(clr_threshold));

  array_storage_Cdf = utility::array_storage(dC_df);
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

  auto dMj =
      TA::to_new_tile_type(Mj, tensor::TaToDecompTensor(clr_threshold, false));

  auto dMk =
      TA::to_new_tile_type(Mk, tensor::TaToDecompTensor(clr_threshold, false));

  const auto occ_nclusters =
      in.HasMember("occ nclusters") ? in["occ nclusters"].GetInt() : nclusters;

  out_doc.AddMember("number occupied clusters", occ_nclusters,
                    out_doc.GetAllocator());

  auto localize = true;
  auto ebuilder = scf::ESolveDensityBuilder(S, r_xyz, occ / 2, occ_nclusters,
                                            0.0, "cholesky inverse", localize);

  std::unique_ptr<scf::DensityBuilder> d_builder =
      make_unique<decltype(ebuilder)>(std::move(ebuilder));

  std::unique_ptr<scf::FockBuilder> f_builder;

  bool force_shape =
      in.HasMember("force shape") ? in["force shape"].GetBool() : false;

  double force_thresh = in.HasMember("force threshold")
                            ? in["force threshold"].GetDouble()
                            : 1e-7;

  double mo_thresh =
      in.HasMember("mo threshold") ? in["mo threshold"].GetDouble() : 0.0;
  if (world.rank() == 0) {
    std::cout << "MO Threshold = " << mo_thresh << std::endl;
    if (force_shape) {
      std::cout << "Force Threshold = " << force_thresh << std::endl;
    }
  }

  if (in.HasMember("stored integrals") &&
      in["stored integrals"].GetBool() == true) {
    auto e0 = mpqc_time::fenced_now(world);
    auto deri3s =
        ints::sparse_integrals(world, eri_e, three_c_array, shr_screen,
                               tensor::TaToDecompTensor(clr_threshold));
    auto e1 = mpqc_time::fenced_now(world);
    auto etime = mpqc_time::duration_in_s(e0, e1);
    if (world.rank() == 0) {
      std::cout << "3 center time: " << etime << std::endl;
    }
    auto e_store = utility::array_storage(deri3s);
    if (world.rank() == 0) {
      std::cout << "E storage:"
                << "\n\tFull     = " << e_store[0]
                << "\n\tSparse   = " << e_store[1]
                << "\n\tCLR      = " << e_store[2] << std::endl;
    }
    out_doc.AddMember("E time", etime, out_doc.GetAllocator());
    out_doc.AddMember("E Full Storage", e_store[0], out_doc.GetAllocator());

    out_doc.AddMember("E Sparse Storage", e_store[1], out_doc.GetAllocator());

    out_doc.AddMember("E CLR Storage", e_store[2], out_doc.GetAllocator());

    scf::ONCADFFockBuilder<decltype(deri3s)> test_scf(
        dMj, dMk, deri3s, dC_df, clr_threshold, force_thresh, mo_thresh,
        force_shape);
    f_builder = make_unique<decltype(test_scf)>(std::move(test_scf));

    out_doc.AddMember("Stored integrals", true, out_doc.GetAllocator());
  } else {
    scf::ONCADFFockBuilder<decltype(deri3)> test_scf(
        dMj, dMk, deri3, dC_df, clr_threshold, force_thresh, mo_thresh,
        force_shape);
    f_builder = make_unique<decltype(test_scf)>(std::move(test_scf));
    out_doc.AddMember("Stored integrals", false, out_doc.GetAllocator());
  }

  auto hf = scf::ClosedShellSCF(H, S, repulsion_energy, std::move(f_builder),
                                std::move(d_builder), F_soad);

  world.gop.fence();

  auto converged = hf.solve(in["scf max iter"].GetInt(),
                            in["scf convergence threshold"].GetDouble());

  out_doc.AddMember("SCF Converged", converged, out_doc.GetAllocator());
  out_doc.AddMember("SCF", hf.results(out_doc), out_doc.GetAllocator());

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
