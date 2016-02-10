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

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/basis.h"

#include "../integrals/integrals.h"

#include "../utility/time.h"
// #include "../utility/array_storage.h"
#include "../ta_routines/array_to_eigen.h"

#include "../scf/diagonalize_for_coffs.hpp"
#include "../scf/soad.h"
#include "../scf/orbital_localization.h"
#include "../scf/clusterd_coeffs.h"
#include "../scf/cadf_helper_functions.h"

#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_nonintrusive_interface.h"
#include "../tensor/mpqc_tile.h"
#include "../tensor/tensor_transforms.h"

#include <memory>

using namespace mpqc;
namespace ints = mpqc::integrals;

bool tensor::detail::recompress = true;

std::array<double, 3>
storage_for_array(TA::DistArray<TA::TensorD, SpPolicy> const &a) {
    std::atomic<long> full_a(0.0);
    std::atomic<long> sparse_a(0.0);

    auto task_f = [&](int ord) {
        auto const &trange = a.trange();
        auto size = trange.make_tile_range(ord).volume();
        full_a += size;

        if (!a.is_zero(ord)) {
            sparse_a += size;
        }
    };

    auto &pmap = *a.get_pmap();
    for (auto it : pmap) {
        a.get_world().taskq.add(task_f, it);
    }
    a.get_world().gop.fence();

    double full = full_a * 1e-9 * 8;
    double sparse = sparse_a * 1e-9 * 8;

    return {full, sparse, 0.0};
}

template <typename Tile>
std::array<double, 3>
storage_for_array(TA::DistArray<Tile, SpPolicy> const &a) {
    std::atomic<long> full_a(0.0);
    std::atomic<long> sparse_a(0.0);
    std::atomic<long> clr_a(0.0);

    auto task_f = [&](int ord) {
        auto const &trange = a.trange();
        auto size = trange.make_tile_range(ord).volume();
        full_a += size;

        if (!a.is_zero(ord)) {
            sparse_a += size;

            tensor::Tile<tensor::DecomposedTensor<double>> wrapper_tile
                  = a.find(ord).get();

            auto const &tile = wrapper_tile.tile();

            if (tile.ndecomp() == 1) {
                clr_a += size;
            } else {
                auto nelems = 0;
                for (auto const &t : tile.tensors()) {
                    nelems += t.range().volume();
                }

                clr_a += nelems;
            }
        }
    };

    auto &pmap = *a.get_pmap();
    for (auto it : pmap) {
        a.get_world().taskq.add(task_f, it);
    }
    a.get_world().gop.fence();

    double full = full_a * 1e-9 * 8;
    double sparse = sparse_a * 1e-9 * 8;
    double clr = clr_a * 1e-9 * 8;

    return {full, sparse, clr};
}

struct JobSummary {
    // Basis and cluster info
    std::string obs;
    std::string dfbs;

    int num_obs = 0;
    int num_dfbs = 0;

    int obs_nfuncs = 0;
    int dfbs_nfuncs = 0;

    // Times
    double avg_w_time = 0;
    double avg_occ_ri_k_time = 0;
    double avg_k_time = 0;
    double avg_loc_time = 0;

    // Storages
    double w_fully_dense_storage = 0;
    double avg_w_sparse_only_storage = 0;
    double avg_w_sparse_clr_storage = 0;

    double eri3_fully_dense_storage = 0;
    double eri3_sparse_only_storage = 0;
    double eri3_sparse_clr_storage = 0;

    double B_fully_dense_storage = 0;
    double B_sparse_only_storage = 0;
    double B_sparse_clr_storage = 0;

    // Thresholds
    double ta_sparse_threshold = 0;
    double clr_threshold = 0;
    double schwarz_threshold = 0;

    // Computables
    double hf_energy = 0;
    double hf_energy_diff_error = 0;
    double hf_grad_div_volume = 0;
    double exchange_energy = 0;

    double dipole_vector = 0;
    double dipole_x = 0;
    double dipole_y = 0;
    double dipole_z = 0;

    double mp2_energy = 0;

    // Bools
    bool did_mp2 = false;
    bool computed_dipole = false;

    void print_report() const {
        std::cout << "Job Summary:\n";

        std::cout << "number obs clusters, ";
        std::cout << "number dfbs clusters, ";

        std::cout << "obs basis, ";
        std::cout << "dfbs basis, ";

        std::cout << "obs nfuctions, ";
        std::cout << "dfbs nfuctions, ";

        std::cout << "ta sparse threshold, ";
        std::cout << "schwarz threshold, ";
        std::cout << "clr threshold, ";

        std::cout << "avg w time, ";
        std::cout << "avg k time, ";
        std::cout << "avg occ-ri k time, ";

        std::cout << "eri3 fully dense storage, ";
        std::cout << "eri3 sparse only, ";
        std::cout << "eri3 sparse + clr storage, ";

        std::cout << "B fully dense storage, ";
        std::cout << "B sparse only, ";
        std::cout << "B sparse + clr storage, ";

        std::cout << "w fully dense storage, ";
        std::cout << "w sparse only, ";
        std::cout << "w sparse + clr storage, ";

        std::cout << "hf energy, ";
        std::cout << "error change in energy, ";
        std::cout << "error rms div volume, ";
        std::cout << "exchange energy(iter 0), ";

        std::cout << "computed dipole, ";

        std::cout << "dipole x, ";
        std::cout << "dipole y, ";
        std::cout << "dipole z, ";
        std::cout << "dipole vec, ";

        std::cout << "computed mp2, ";
        std::cout << "mp2 energy\n";

        std::cout << num_obs << ", ";
        std::cout << num_dfbs << ", ";

        std::cout << obs << ", ";
        std::cout << dfbs << ", ";

        std::cout << obs_nfuncs << ", ";
        std::cout << dfbs_nfuncs << ", ";

        std::cout << ta_sparse_threshold << ", ";
        std::cout << schwarz_threshold << ", ";
        std::cout << clr_threshold << ", ";


        std::cout << avg_w_time << ", ";
        std::cout << avg_k_time << ", ";
        std::cout << avg_occ_ri_k_time << ", ";

        std::cout << eri3_fully_dense_storage << ", ";
        std::cout << eri3_sparse_only_storage << ", ";
        std::cout << eri3_sparse_clr_storage << ", ";

        std::cout << B_fully_dense_storage << ", ";
        std::cout << B_sparse_only_storage << ", ";
        std::cout << B_sparse_clr_storage << ", ";

        std::cout << w_fully_dense_storage << ", ";
        std::cout << avg_w_sparse_only_storage << ", ";
        std::cout << avg_w_sparse_clr_storage << ", ";

        std::cout << hf_energy << ", ";
        std::cout << hf_energy_diff_error << ", ";
        std::cout << hf_grad_div_volume << ", ";
        std::cout << exchange_energy << ", ";

        std::cout << computed_dipole << ", ";
        std::cout << dipole_x << ", ";
        std::cout << dipole_y << ", ";
        std::cout << dipole_z << ", ";
        std::cout << dipole_vector << ", ";

        std::cout << did_mp2 << ", ";
        std::cout << mp2_energy << "\n";
    }
};

class ThreeCenterScf {
  private:
    using array_type = DArray<2, TA::TensorD, SpPolicy>;
    using dtile = tensor::Tile<tensor::DecomposedTensor<double>>;
    using darray_type = DArray<2, dtile, SpPolicy>;

    array_type H_;
    array_type S_;

    array_type F_;
    array_type K_;
    array_type D_;

    array_type V_inv_oh_;

    darray_type dH_;
    darray_type dC_;

    std::vector<array_type> r_xyz_ints_;
    TiledArray::DIIS<array_type> diis_;

    std::vector<double> k_times_;
    std::vector<double> occ_k_times_;
    std::vector<double> j_times_;
    std::vector<double> w_times_;
    std::vector<double> scf_times_;
    std::vector<double> w_sparse_store_;
    std::vector<double> w_sparse_clr_store_;
    std::vector<double> localization_times_;

    int64_t occ_;
    double repulsion_;
    double final_rms_error_;
    double final_ediff_error_;
    double final_energy_;
    double w_full_storage_;

    double recompress_w_time_;
    double clr_w_no_recompress_;

    double clr_thresh_;
    double exchange_energy_;

    void compute_density(int64_t occ) {
        auto F_eig = array_ops::array_to_eigen(F_);
        auto S_eig = array_ops::array_to_eigen(S_);

        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);
        decltype(S_eig) C = es.eigenvectors().leftCols(occ);
        auto tr_ao = S_.trange().data()[0];

        auto tr_occ = scf::tr_occupied(1, occ);
        auto C_TA = array_ops::eigen_to_array<TA::TensorD>(H_.get_world(), C,
                                                           tr_ao, tr_occ);

        dC_ = TA::to_new_tile_type(C_TA, tensor::TaToDecompTensor(clr_thresh_,
                                                                  false));

        decltype(S_eig) D = C * C.transpose();

        D_ = array_ops::eigen_to_array<TA::TensorD>(H_.get_world(), D, tr_ao,
                                                    tr_ao);
    }

    template <typename Integral, typename Array>
    void form_fock(Integral const &eri3, Array const &C_df, Array const &G_df) {
        auto &world = F_.get_world();

        auto get_time = [&]() {
            world.gop.fence();
            return utility::time::now();
        };

        using tp = decltype(get_time());

        auto calc_time = [](tp const &a, tp const &b) {
            return utility::time::duration_in_s(a, b);
        };

        array_type J;
        auto j0 = get_time();
        J("mu, nu") = eri3("Y,mu,nu")
                      * (V_inv_oh_("Y,Z")
                         * (V_inv_oh_("Z,X") * (eri3("X,r,s") * D_("r,s"))));
        auto j1 = get_time();
        j_times_.push_back(calc_time(j0, j1));

        darray_type dD_ = TA::to_new_tile_type(
              D_, tensor::TaToDecompTensor(clr_thresh_, false));

        using GTile = tensor::Tile<tensor::DecomposedTensor<double>>;
        DArray<3, GTile, SpPolicy> D_df;
        DArray<3, GTile, SpPolicy> F_df;

        bool occ_ri = false;
        if (!occ_ri) {
            auto d0 = get_time();
            D_df("X, mu, i") = C_df("X, mu, nu") * dC_("nu, i");

            auto storage = storage_for_array(D_df);
            std::cout << "D_df full = " << storage[0] << std::endl;
            std::cout << "D_df Sparse = " << storage[1] << std::endl;
            std::cout << "D_df CLR = " << storage[2] << std::endl;

            if (clr_thresh_ != 0 && tensor::detail::recompress) {
                TA::foreach_inplace(
                      D_df, [](tensor::Tile<tensor::DecomposedTensor<double>> &
                                     t_tile) {
                          auto &t = t_tile.tile();
                          auto input_norm = norm(t);

                          auto compressed_norm = input_norm;
                          if (t.cut() != 0.0) {
                              if (t.ndecomp() == 1) {
                                  auto test = tensor::algebra::
                                        two_way_decomposition(t);
                                  if (!test.empty()) {
                                      t = test;
                                      compressed_norm = norm(t);
                                  }
                              } else {
                                  tensor::algebra::recompress(t);
                                  compressed_norm = norm(t);
                              }
                          }

                          // Both are always larger than or equal to the real
                          // norm.
                          return std::min(input_norm, compressed_norm);
                      });
            } else {
                D_df.truncate();
            }

            D_df.get_world().gop.fence();
            storage = storage_for_array(D_df);
            std::cout << "D_df full = " << storage[0] << std::endl;
            std::cout << "D_df Sparse = " << storage[1] << std::endl;
            std::cout << "D_df CLR = " << storage[2] << std::endl;

            D_df("X, i, mu") = D_df("X, mu, i");
            auto d1 = get_time();
            w_times_.push_back(calc_time(d0, d1));

            F_df("X, i, mu") = G_df("X, mu, nu") * dC_("nu, i");

            auto k0 = get_time();
            darray_type dL;
            dL("mu, nu") = D_df("X, i, mu") * F_df("X, i, nu");


            auto L = TA::to_new_tile_type(dL, tensor::DecompToTaTensor{});

            K_("mu, nu") = L("mu, nu") + L("nu, mu");
            auto k1 = get_time();
            k_times_.push_back(calc_time(k0, k1));

            F_("i,j") = H_("i,j") + 2 * J("i,j") - K_("i,j");
        } else {
            auto d0 = get_time();
            D_df("X, nu, i") = C_df("X, nu, mu") * dC_("mu, i");
            D_df.truncate();
            D_df("X, i, nu") = D_df("X, nu, i");
            D_df("X, i, j") = D_df("X, i, nu") * dC_("nu,j");

            auto storage = storage_for_array(D_df);
            std::cout << "D_df full = " << storage[0] << std::endl;
            std::cout << "D_df Sparse = " << storage[1] << std::endl;
            std::cout << "D_df CLR = " << storage[2] << std::endl;

            if (clr_thresh_ != 0 && tensor::detail::recompress) {
                TA::foreach_inplace(
                      D_df, [](tensor::Tile<tensor::DecomposedTensor<double>> &
                                     t_tile) {
                          auto &t = t_tile.tile();
                          auto input_norm = norm(t);

                          auto compressed_norm = input_norm;
                          if (t.cut() != 0.0) {
                              if (t.ndecomp() == 1) {
                                  auto test = tensor::algebra::
                                        two_way_decomposition(t);
                                  if (!test.empty()) {
                                      t = test;
                                      compressed_norm = norm(t);
                                  }
                              } else {
                                  tensor::algebra::recompress(t);
                                  compressed_norm = norm(t);
                              }
                          }

                          // Both are always larger than or equal to the real
                          // norm.
                          return std::min(input_norm, compressed_norm);
                      });
            } else {
                D_df.truncate();
            }

            D_df.get_world().gop.fence();
            storage = storage_for_array(D_df);
            std::cout << "D_df full = " << storage[0] << std::endl;
            std::cout << "D_df Sparse = " << storage[1] << std::endl;
            std::cout << "D_df CLR = " << storage[2] << std::endl;

            F_df("X, mu, i") = G_df("X, mu, nu") * dC_("nu, i");
            F_df("X, i, mu") = F_df("X, mu, i");

            auto d1 = get_time();
            w_times_.push_back(calc_time(d0, d1));

            auto k0 = get_time();
            darray_type dL;
            dL("j, nu") = D_df("X, i, j") * F_df("X, i, nu");
            dL("nu,j") = dL("j,nu");

            auto dS = TA::to_new_tile_type(
                  S_, tensor::TaToDecompTensor(clr_thresh_, false));
            decltype(dL) dLij, Sc;
            dLij("i,j") = dC_("mu, i") * dL("mu, j");
            Sc("mu, j") = dS("mu, lam") * dC_("lam, j");
            dL("mu, nu") = Sc("mu, j") * dL("nu, j") + dL("mu, j") * Sc("nu,j")
                           - (Sc("mu,i") * dLij("i,j")) * Sc("nu,j");

            auto L = TA::to_new_tile_type(dL, tensor::DecompToTaTensor{});

            K_("mu, nu") = L("mu, nu") + L("nu, mu");
            auto k1 = get_time();
            k_times_.push_back(calc_time(k0, k1));

            F_("i,j") = H_("i,j") + 2 * J("i,j") - K_("i,j");
        }
    }

  public:
    ThreeCenterScf(array_type const &H, array_type const &F_guess,
                   array_type const &S, array_type const &V_inv_oh,
                   std::vector<array_type> const &rxyz, int64_t occ, double rep,
                   double clr_thresh)
            : H_(H),
              F_(F_guess),
              S_(S),
              V_inv_oh_(V_inv_oh),
              r_xyz_ints_(rxyz),
              occ_(occ),
              repulsion_(rep),
              clr_thresh_(clr_thresh) {

        compute_density(occ_);
    }

    template <typename Integral, typename Array>
    void solve(int64_t max_iters, double thresh, Integral const &eri3,
               Array const &C_df, Array const &G_df) {
        auto iter = 0;
        auto error = std::numeric_limits<double>::max();
        auto rms_error = std::numeric_limits<double>::max();
        auto old_energy = 0.0;
        const double volume = F_.trange().elements().volume();

        while (iter < max_iters && (thresh < error || thresh < rms_error)
               && (thresh / 100.0 < error && thresh / 100.0 < rms_error)) {
            auto s0 = mpqc_time::now();
            F_.get_world().gop.fence();
            form_fock(eri3, C_df, G_df);

            auto current_energy = energy();
            error = std::abs(old_energy - current_energy);
            old_energy = current_energy;

            if (iter == 0) {
                exchange_energy_
                      = D_("i,j").dot(K_("i,j"), F_.get_world()).get();
            }

            array_type Grad;
            Grad("i,j") = F_("i,k") * D_("k,l") * S_("l,j")
                          - S_("i,k") * D_("k,l") * F_("l,j");

            rms_error = Grad("i,j").norm().get() / volume;

            diis_.extrapolate(F_, Grad);

            // Lastly update density
            compute_density(occ_);

            F_.get_world().gop.fence();
            auto s1 = mpqc_time::now();
            scf_times_.push_back(mpqc_time::duration_in_s(s0, s1));


            std::cout << "Iteration: " << (iter + 1)
                      << " energy: " << old_energy << " error: " << error
                      << " RMS error: " << rms_error;
            std::cout
                  << "\n\tD df time: " << w_times_.back()
                  // << "\n\tW recompress time: " << recompress_w_time_
                  << "\n\tJ time: " << j_times_.back()
                  << "\n\tK time: " << k_times_.back() << "\n\titer time: "
                  << scf_times_.back()
                  // << "\n\tW sparse only storage: " << w_sparse_store_.back()
                  // << "\n\tW sparse clr(no recompress): " <<
                  // clr_w_no_recompress_
                  // << "\n\tW sparse clr storage: " <<
                  // w_sparse_clr_store_.back()
                  // << "\n\tLocalization time: " << localization_times_.back()
                  << std::endl;
            // if (iter == 0) {
            //     std::cout << "\tExchange energy = " << exchange_energy_
            //               << std::endl;
            // }


            ++iter;
        }

        final_rms_error_ = rms_error;
        final_ediff_error_ = error;
        final_energy_ = old_energy;
    }

    double energy() {
        return repulsion_
               + D_("i,j").dot(F_("i,j") + H_("i,j"), D_.get_world()).get();
    }

    array_type fock_matrix() const { return F_; }

    JobSummary init_summary() {
        auto avg = [](std::vector<double> const &v) {
            auto sum = 0.0;
            for (auto const &e : v) {
                sum += e;
            }
            return sum / v.size();
        };

        JobSummary js;

        // Times
        js.avg_w_time = avg(w_times_);
        js.avg_k_time = avg(k_times_);
        js.avg_occ_ri_k_time = avg(occ_k_times_);
        js.avg_loc_time = avg(localization_times_);

        // Storages
        js.avg_w_sparse_only_storage = avg(w_sparse_store_);
        js.avg_w_sparse_clr_storage = avg(w_sparse_clr_store_);
        js.w_fully_dense_storage = w_full_storage_;

        // Thresholds
        js.ta_sparse_threshold = TA::SparseShape<float>::threshold();

        // Energies
        js.hf_energy = final_energy_;
        js.hf_energy_diff_error = final_ediff_error_;
        js.hf_grad_div_volume = final_rms_error_;
        js.exchange_energy = exchange_energy_;

        return js;
    }
};


int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    std::string df_basis_name = "";
    int nclusters = 0;
    std::cout << std::setprecision(15);
    double threshold = 1e-12;
    double clr_threshold = 1e-8;
    std::string eri3_storage_method;
    std::string recompression_method;
    double mp2_mem_thresh = 60;
    if (argc >= 9) {
        mol_file = argv[1];
        basis_name = argv[2];
        df_basis_name = argv[3];
        nclusters = std::stoi(argv[4]);
        threshold = std::stod(argv[5]);
        clr_threshold = std::stod(argv[6]);
        eri3_storage_method = argv[7];
        recompression_method = argv[8];
    }
    if (argc == 10) {
        mp2_mem_thresh = std::stod(argv[9]);
    }
    if (argc >= 11 || argc < 9) {
        std::cout << "input is $./program mol_file basis_file df_basis_file "
                     "nclusters sparse_thresh clr_thresh eri3_method(direct, "
                     "stored) recompression_method(compress, nocompress)";
        return 0;
    }

    if (recompression_method == "nocompress") {
        std::cout << "Not using rounded addition." << std::endl;
        tensor::detail::recompress = false;
    } else {
        std::cout << "Using rounded addition." << std::endl;
    }
    world.gop.fence();

    bool use_direct_ints;
    if (eri3_storage_method == "direct") {
        use_direct_ints = true;
        std::cout << "Using direct 3 center integrals" << std::endl;
    } else if (eri3_storage_method == "stored") {
        use_direct_ints = false;
        std::cout << "Using stored 3 center integrals" << std::endl;
    } else {
        std::cout << "Integral storage method must be either direct or stored"
                  << std::endl;
        return 0;
    }

    TiledArray::SparseShape<float>::threshold(threshold);
    std::cout << "TA::Sparse Threshold: " << threshold << std::endl;

    auto clustered_mol = molecule::attach_hydrogens_and_kmeans(
          molecule::read_xyz(mol_file).clusterables(), nclusters);

    auto repulsion_energy = clustered_mol.nuclear_repulsion();
    std::cout << "Nuclear Repulsion Energy: " << repulsion_energy << std::endl;
    auto occ = clustered_mol.occupation(0);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));
    std::cout << "Basis has " << basis.nfunctions() << " functions\n";

    // auto df_nclusters = std::max(nclusters / 2, 1);
    // auto df_clustered_mol = molecule::attach_hydrogens_and_kmeans(
    //       molecule::read_xyz(mol_file).clusterables(), df_nclusters);
    auto df_clustered_mol = clustered_mol;


    std::cout << "\nDF clustering:\n";
    for (auto const &c : df_clustered_mol.clusterables()) {
        auto atoms = c.atoms();
        std::cout << atoms.size() << "\n\n";
        for (auto const &atom : c.atoms()) {
            std::cout << atom.xyz_string(true) << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

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

        const auto schwarz_thresh = 1e-12;
        auto sbuilder = ints::init_schwarz_screen(schwarz_thresh);
        auto shr_screen = std::make_shared<ints::SchwarzScreen>(
              sbuilder(world, eri_e, df_by_atom, obs_by_atom));

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

        std::cout << "C df sparsity (by atom) = " << C_df_.get_shape().sparsity()
                  << std::endl;

        auto by_atom_trange = integrals::detail::create_trange(three_c_array);
        auto by_cluster_trange = integrals::detail::create_trange(
              utility::make_array(df_basis, basis, basis));
        std::cout << "Trange by atom: " << by_atom_trange << std::endl;
        std::cout << "Trange by cluster: " << by_cluster_trange << std::endl;

        C_df = scf::reblock_from_atoms(C_df_, obs_atom_to_cluster_map,
                                obs_atom_to_cluster_map, by_cluster_trange);

        std::cout << "C df sparsity (by cluster) = " << C_df.get_shape().sparsity()
                  << std::endl;
    }

#if 1
    // Continue with scf
    {
        const auto dfbs_array = utility::make_array(df_basis, df_basis);
        auto eri_e = ints::make_2body_shr_pool(df_basis, basis);

        auto three_c_array = utility::make_array(df_basis, basis, basis);
        decltype(H) M = ints::sparse_integrals(world, eri_e, dfbs_array);

        const auto schwarz_thresh = 1e-14;
        std::cout << "Schwarz Threshold: " << schwarz_thresh << std::endl;
        auto sbuilder = ints::init_schwarz_screen(schwarz_thresh);
        auto shr_screen = std::make_shared<ints::SchwarzScreen>(
              sbuilder(world, eri_e, df_basis, basis));

        auto eri3 = ints::direct_sparse_integrals(world, eri_e, three_c_array,
                                                  shr_screen);

        // Make M_inv
        decltype(M) M_inv_oh;
        auto M_eig = TA::array_to_eigen(M);
        Eig::SelfAdjointEigenSolver<MatrixD> es(M_eig);
        auto M_eig_inv_oh = es.operatorInverseSqrt();

        M_inv_oh = array_ops::eigen_to_array<TA::TensorD>(
              world, M_eig_inv_oh, M.trange().data()[0], M.trange().data()[1]);

        auto deri3 = ints::direct_sparse_integrals(
              world, eri_e, three_c_array, shr_screen,
              tensor::TaToDecompTensor(clr_threshold));

        auto dC_df = TA::to_new_tile_type(C_df, tensor::TaToDecompTensor(
                                                      clr_threshold));

        // Don't CLR compress M
        auto dM = TA::to_new_tile_type(
              M, tensor::TaToDecompTensor(clr_threshold, false));

        world.gop.fence();
        auto old_compress = tensor::detail::recompress;
        tensor::detail::recompress = true;

        decltype(dC_df) dG_df;
        dG_df("X,i,j") = deri3("X,i,j") - 0.5 * dM("X,Y") * dC_df("Y,i,j");
        if (clr_threshold != 0 && tensor::detail::recompress) {
            TA::foreach_inplace(
                  dG_df,
                  [](tensor::Tile<tensor::DecomposedTensor<double>> &t_tile) {
                      auto &t = t_tile.tile();
                      auto input_norm = norm(t);

                      auto compressed_norm = input_norm;
                      if (t.cut() != 0.0) {
                          if (t.ndecomp() == 1) {
                              auto test
                                    = tensor::algebra::two_way_decomposition(t);
                              if (!test.empty()) {
                                  t = test;
                                  compressed_norm = norm(t);
                              }
                          } else {
                              tensor::algebra::recompress(t);
                              compressed_norm = norm(t);
                          }
                      }

                      // Both are always larger than or equal to the real
                      // norm.
                      return std::min(input_norm, compressed_norm);
                  });
        } else {
            dG_df.truncate();
        }

        world.gop.fence();
        tensor::detail::recompress = old_compress;

        auto storage = storage_for_array(dC_df);
        std::cout << "Dense C_df " << storage[0] << "\n"
                  << "Sparse C_df " << storage[1] << "\n"
                  << "CLR C_df " << storage[2] << "\n\n";

        storage = storage_for_array(dG_df);
        std::cout << "Dense G_df " << storage[0] << "\n"
                  << "Sparse G_df " << storage[1] << "\n"
                  << "CLR G_df " << storage[2] << "\n";

        auto F_soad
              = scf::fock_from_soad(world, clustered_mol, basis, eri_e, H);

        auto multi_pool
              = ints::make_1body_shr_pool("emultipole2", basis, clustered_mol);

        auto r_xyz = ints::sparse_xyz_integrals(world, multi_pool, bs_array);

        ThreeCenterScf scf(H, F_soad, S, M_inv_oh, r_xyz, occ / 2,
                           repulsion_energy, clr_threshold);

        scf.solve(20, 1e-11, eri3, dC_df, dG_df);
    }
#endif

    return 0;
}
