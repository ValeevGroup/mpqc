#include "../common/namespaces.h"
#include "../common/typedefs.h"

#include "../include/libint.h"
#include "../include/tiledarray.h"

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
#include "../utility/make_array.h"
#include "../utility/array_info.h"
#include "../ta_routines/array_to_eigen.h"

#include "../scf/diagonalize_for_coffs.hpp"
#include "../scf/soad.h"
#include "../scf/orbital_localization.h"
#include "../scf/clusterd_coeffs.h"

#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_nonintrusive_interface.h"
#include "../tensor/mpqc_tile.h"

#include "../ta_routines/diagonal_array.h"

#include <memory>

using namespace mpqc;
namespace ints = mpqc::integrals;

bool tensor::detail::recompress = true;

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
    double avg_cluster_time = 0;

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

            std::cout << "avg localization time, ";
            std::cout << "avg clustering time, ";

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

            std::cout << avg_loc_time << ", ";
            std::cout << avg_cluster_time << ", ";

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

class to_dtile {
    double clr_thresh_;

  public:
    using dtile = tensor::Tile<tensor::DecomposedTensor<double>>;
    to_dtile(double clr_thresh) : clr_thresh_(clr_thresh) {}

    dtile operator()(TA::TensorD const &tile) {
        auto range = tile.range();

        const auto extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];

        auto local_range = TA::Range{i, j};
        auto tensor = TA::TensorD(local_range, tile.data());

        auto new_tile = tensor::DecomposedTensor<double>(clr_thresh_,
                                                         std::move(tensor));

        return dtile(range, std::move(new_tile));
    }
};

class to_dtile_compress {
    double clr_thresh_;

  public:
    using dtile = tensor::Tile<tensor::DecomposedTensor<double>>;
    to_dtile_compress(double clr_thresh) : clr_thresh_(clr_thresh) {}

    dtile operator()(TA::TensorD const &tile) {
        auto range = tile.range();

        const auto extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];

        auto local_range = TA::Range{i, j};
        auto tensor = TA::TensorD(local_range, tile.data());

        auto new_tile = tensor::DecomposedTensor<double>(clr_thresh_,
                                                         std::move(tensor));

        if (clr_thresh_ != 0.0) {
            auto test = tensor::algebra::two_way_decomposition(new_tile);
            if (!test.empty()) {
                new_tile = std::move(test);
            }
        }

        return dtile(range, std::move(new_tile));
    }
};

class to_ta_tile {
    using dtile = tensor::Tile<tensor::DecomposedTensor<double>>;

  public:
    TA::TensorD operator()(dtile const &t) const {
        auto tensor = tensor::algebra::combine(t.tile());
        auto range = t.range();
        return TA::Tensor<double>(range, tensor.data());
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
    array_type C_;

    darray_type dV_inv_oh_;
    darray_type dH_;

    std::vector<array_type> r_xyz_ints_;
    TiledArray::DIIS<array_type> diis_;

    std::vector<double> k_times_;
    std::vector<double> occ_k_times_;
    std::vector<double> j_times_;
    std::vector<double> w_times_;
    std::vector<double> scf_times_;
    std::vector<double> w_sparse_store_;
    std::vector<double> w_sparse_clr_store_;
    std::vector<double> clustering_times_;
    std::vector<double> localization_times_;

    int64_t occ_;
    double repulsion_;
    double final_rms_error_;
    double final_ediff_error_;
    double final_energy_;
    double w_full_storage_;

    double clr_thresh_;
    double exchange_energy_;

    void compute_density(int64_t occ) {
        auto &world = F_.get_world();

        auto F_eig = array_ops::array_to_eigen(F_);
        auto S_eig = array_ops::array_to_eigen(S_);

        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);
        decltype(S_eig) C = es.eigenvectors().leftCols(occ);
        auto tr_ao = S_.trange().data()[0];

        auto tr_occ = scf::tr_occupied(occ_, occ_);
        C_ = array_ops::eigen_to_array<TA::TensorD>(H_.get_world(), C, tr_ao,
                                                    tr_occ);

        auto loc0 = mpqc_time::fenced_now(world);

        auto U = mpqc::scf::BoysLocalization{}(C_, r_xyz_ints_);
        C_("mu,i") = C_("mu,k") * U("k,i");

        auto loc1 = mpqc_time::fenced_now(world);
        localization_times_.push_back(mpqc_time::duration_in_s(loc0, loc1));

        auto cluster0 = mpqc_time::fenced_now(world);

        auto obs_ntiles = C_.trange().tiles().extent()[0];
        scf::clustered_coeffs(r_xyz_ints_, C_, obs_ntiles);

        auto cluster1 = mpqc_time::fenced_now(world);
        clustering_times_.push_back(
              mpqc_time::duration_in_s(cluster0, cluster1));

        D_("i,j") = C_("i,k") * C_("j,k");
    }

    template <typename Integral, typename Array>
    void form_fock(Integral const &eri3, Array const &B) {
        auto &world = F_.get_world();

        darray_type dH_ = TA::to_new_tile_type(H_, to_dtile(clr_thresh_));
        darray_type dS_ = TA::to_new_tile_type(S_, to_dtile(clr_thresh_));
        darray_type dC_ = TA::to_new_tile_type(C_, to_dtile(clr_thresh_));
        darray_type dD_ = TA::to_new_tile_type(D_, to_dtile(clr_thresh_));

        darray_type J;
        auto j0 = mpqc_time::fenced_now(world);
        auto old_thresh = TA::SparseShape<float>::threshold();
        TA::SparseShape<float>::threshold(std::min(1e-16f, old_thresh));

        J("mu, nu") = eri3("Y,mu,nu")
                      * (dV_inv_oh_("Y,Z")
                         * (dV_inv_oh_("Z,X") * (eri3("X,r,s") * dD_("r,s"))));

        auto j1 = mpqc_time::fenced_now(world);
        TA::SparseShape<float>::threshold(old_thresh);

        j_times_.push_back(mpqc_time::duration_in_s(j0, j1));

        auto w0 = mpqc_time::fenced_now(world);

        DArray<3, dtile, SpPolicy> W;
        W("X, mu, i") = B("X, mu, nu") * dC_("nu, i");

        if (clr_thresh_ != 0 && tensor::detail::recompress) {
            TA::foreach_inplace(
                  W,
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

                      // Both are always larger than or equal to the real norm.
                      return std::min(input_norm, compressed_norm);
                  });
        } else {
            W.truncate();
        }

        W("X,i,mu") = W("X,mu,i");

        auto w1 = mpqc_time::fenced_now(world);
        w_times_.push_back(mpqc_time::duration_in_s(w0, w1));

        auto w_store = utility::array_storage(W);
        w_sparse_store_.push_back(w_store[1]);
        w_sparse_clr_store_.push_back(w_store[2]);
        w_full_storage_ = w_store[0];


        auto occk0 = mpqc_time::fenced_now(world);

        // OCC RI
        darray_type K, Kij, Sc;
        K("mu, j") = W("X, i, mu") * ((W("X, i, nu") * dC_("nu, j")));
        Kij("i,j") = dC_("mu, i") * K("mu, j");
        Sc("mu, j") = dS_("mu, lam") * dC_("lam, j");
        K("mu, nu") = Sc("mu, j") * K("nu, j") + K("mu, j") * Sc("nu,j")
                      - (Sc("mu,i") * Kij("i,j")) * Sc("nu,j");

        auto occk1 = mpqc_time::fenced_now(world);
        occ_k_times_.push_back(mpqc_time::duration_in_s(occk0, occk1));

        auto k0 = mpqc_time::fenced_now(world);

        K("mu, nu") = W("X,i,mu") * W("X,i,nu");

        auto k1 = mpqc_time::fenced_now(world);
        k_times_.push_back(mpqc_time::duration_in_s(k0, k1));

        darray_type dF_;
        dF_("i,j") = dH_("i,j") + 2 * J("i,j") - K("i,j");

        F_ = TA::to_new_tile_type(dF_, to_ta_tile{});
        K_ = TA::to_new_tile_type(K, to_ta_tile{});
    }

  public:
    ThreeCenterScf(array_type const &H, array_type const &F_guess,
                   array_type const &S, array_type const &V_inv_oh,
                   std::vector<array_type> const &rxyz, int64_t occ, double rep,
                   double clr_thresh)
            : H_(H),
              F_(F_guess),
              S_(S),
              r_xyz_ints_(rxyz),
              occ_(occ),
              repulsion_(rep),
              clr_thresh_(clr_thresh) {

        dV_inv_oh_ = TA::to_new_tile_type(V_inv_oh, to_dtile(clr_thresh_));
        auto dl_sizes = utility::array_storage(dV_inv_oh_);

        if (dV_inv_oh_.get_world().rank() == 0) {
            std::cout << "V_inv storage:"
                      << "\n\tFull: " << dl_sizes[0]
                      << "\n\tSparse Only: " << dl_sizes[1]
                      << "\n\tSparse + CLR: " << dl_sizes[2] << std::endl;
        }

        compute_density(occ_);
    }

    template <typename Integral, typename Array>
    void solve(int64_t max_iters, double thresh, Integral const &eri3,
               Array const &B) {
        auto iter = 0;
        auto error = std::numeric_limits<double>::max();
        auto rms_error = std::numeric_limits<double>::max();
        auto old_energy = 0.0;
        const double volume = F_.trange().elements().volume();

        auto &world = F_.get_world();

        while (iter < max_iters && (thresh < error || thresh < rms_error)
               && (thresh / 100.0 < error && thresh / 100.0 < rms_error)) {
            auto s0 = mpqc_time::fenced_now(world);
            form_fock(eri3, B);

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

            auto s1 = mpqc_time::fenced_now(world);
            scf_times_.push_back(mpqc_time::duration_in_s(s0, s1));


            if (world.rank() == 0) {
                std::cout << "Iteration: " << (iter + 1)
                          << " energy: " << old_energy << " error: " << error
                          << " RMS error: " << rms_error;
                std::cout
                      << "\n\tW time: " << w_times_.back()
                      << "\n\tJ time: " << j_times_.back()
                      << "\n\tocc-RI K time: " << occ_k_times_.back()
                      << "\n\tK time: " << k_times_.back()
                      << "\n\titer time: " << scf_times_.back()
                      << "\n\tW sparse only storage: " << w_sparse_store_.back()
                      << "\n\tW sparse clr storage: "
                      << w_sparse_clr_store_.back()
                      << "\n\tClustering time: " << clustering_times_.back()
                      << "\n\tLocalization time: " << localization_times_.back()
                      << std::endl;
                if (iter == 0) {
                    std::cout << "\tExchange energy = " << exchange_energy_
                              << std::endl;
                }
            }


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
        js.avg_cluster_time = avg(clustering_times_);

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
        if (world.rank() == 0) {
            std::cout << "Not using rounded addition." << std::endl;
        }
        tensor::detail::recompress = false;
    } else {
        if (world.rank() == 0) {
            std::cout << "Using rounded addition." << std::endl;
        }
    }
    world.gop.fence();

    bool use_direct_ints;
    if (eri3_storage_method == "direct") {
        use_direct_ints = true;
        if (world.rank() == 0)
            std::cout << "Using direct 3 center integrals" << std::endl;
    } else if (eri3_storage_method == "stored") {
        use_direct_ints = false;
        if (world.rank() == 0)
            std::cout << "Using stored 3 center integrals" << std::endl;
    } else {
        if (world.rank() == 0) {
            std::cout
                  << "Integral storage method must be either direct or stored"
                  << std::endl;
        }
        return 0;
    }

    TiledArray::SparseShape<float>::threshold(threshold);

    if (world.rank() == 0)
        std::cout << "TA::Sparse Threshold: " << threshold << std::endl;

    auto clustered_mol = molecule::attach_hydrogens_and_kmeans(
          molecule::read_xyz(mol_file).clusterables(), nclusters);


    auto repulsion_energy = clustered_mol.nuclear_repulsion();
    if (world.rank() == 0) {
        std::cout << "Nuclear Repulsion Energy: " << repulsion_energy
                  << std::endl;
    }
    auto occ = clustered_mol.occupation(0);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));
    if(world.rank() == 0){
        std::cout << "Basis has " << basis.nfunctions() << " functions\n";
    }

    auto df_nclusters = std::max(nclusters / 2, 1);
    auto df_clustered_mol = molecule::attach_hydrogens_and_kmeans(
          molecule::read_xyz(mol_file).clusterables(), df_nclusters);

    if (world.rank() == 0) {
        auto cluster_num = 1;
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
    }

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

    decltype(H) L_inv, V_inv_oh;
    {
        auto Vmetric = ints::sparse_integrals(world, eri_e, dfbs_array);
        auto V_eig = array_ops::array_to_eigen(Vmetric);

        MatrixD Leig = Eig::LLT<MatrixD>(V_eig).matrixL();
        MatrixD L_inv_eig = Leig.inverse();

        Eig::SelfAdjointEigenSolver<MatrixD> es(V_eig);
        MatrixD V_inv_oh_eig = es.operatorInverseSqrt();

        auto tr_V = Vmetric.trange().data()[0];
        L_inv = array_ops::eigen_to_array<TA::TensorD>(world, L_inv_eig, tr_V,
                                                       tr_V);
        V_inv_oh = array_ops::eigen_to_array<TA::TensorD>(world, V_inv_oh_eig,
                                                          tr_V, tr_V);
    }

    auto three_c_array = utility::make_array(df_basis, basis, basis);

    { // Schwarz Screened ints
        const auto schwarz_thresh = 1e-12;
        if (world.rank() == 0){
            std::cout << "Schwarz Threshold: " << schwarz_thresh << std::endl;
        }
        auto sbuilder = ints::init_schwarz_screen(schwarz_thresh);
        auto shr_screen = std::make_shared<ints::SchwarzScreen>(
              sbuilder(world, eri_e, df_basis, basis));

        auto soad0 = mpqc_time::fenced_now(world);
        auto F_soad
              = scf::fock_from_soad(world, clustered_mol, basis, eri_e, H);

        auto soad1 = mpqc_time::fenced_now(world);
        auto soad_time = mpqc_time::duration_in_s(soad0, soad1);
        if (world.rank() == 0){
            std::cout << "Soad Time: " << soad_time << std::endl;
        }

        auto multi_pool
              = ints::make_1body_shr_pool("emultipole2", basis, clustered_mol);
        auto r_xyz = ints::sparse_xyz_integrals(world, multi_pool, bs_array);

        JobSummary job_summary;

        DArray<2, TA::TensorD, SpPolicy> Fao;

        auto decomp_3d = [&](TA::TensorD &&t) {
            auto range = t.range();

            auto const extent = range.extent();
            const auto i = extent[0];
            const auto j = extent[1];
            const auto k = extent[2];
            auto local_range = TA::Range{i, j, k};

            auto tensor = TA::TensorD(local_range, t.data());

            tensor::DecomposedTensor<double> dc_tile(clr_threshold,
                                                     std::move(tensor));

            if (clr_threshold != 0.0) {
                auto lr_tile = tensor::algebra::two_way_decomposition(dc_tile);

                if (!lr_tile.empty()) {
                    dc_tile = std::move(lr_tile);
                }
            }

            return tensor::Tile<decltype(dc_tile)>(range, std::move(dc_tile));
        };

        auto decomp_3d_no_compress = [&](TA::TensorD &&t) {
            auto range = t.range();

            auto const extent = range.extent();
            const auto i = extent[0];
            const auto j = extent[1];
            const auto k = extent[2];
            auto local_range = TA::Range{i, j, k};

            auto tensor = TA::TensorD(local_range, t.data());

            tensor::DecomposedTensor<double> dc_tile(0.0, std::move(tensor));

            return tensor::Tile<decltype(dc_tile)>(range, std::move(dc_tile));
        };

        if (!use_direct_ints) {
            std::cout << "Always use direct eri3 for now!\n";
            return 0;
        }
        if (use_direct_ints) {


            using BTile = tensor::Tile<tensor::DecomposedTensor<double>>;
            DArray<3, BTile, SpPolicy> B;
            std::array<double, 3> eri3_storage;
            std::array<double, 3> b_storage;
            {
                world.gop.fence();
                auto int0 = mpqc_time::fenced_now(world);

                auto eri3 = ints::direct_sparse_integrals(
                      world, eri_e, three_c_array, shr_screen, decomp_3d);

                auto int1 = mpqc_time::fenced_now(world);
                auto eri3_time = mpqc_time::duration_in_s(int0, int1);


                if (world.rank() == 0) {
                    std::cout << "Eri3 Integral time: " << eri3_time
                              << std::endl;
                }
                eri3_storage = utility::array_storage(eri3.array());
                if (world.rank() == 0) {
                    std::cout << "Eri3 Integral Storage:\n";
                    std::cout << "\tAll Full: " << eri3_storage[0] << "\n";
                    std::cout << "\tBlock Sparse Only: " << eri3_storage[1]
                              << "\n";
                    std::cout << "\tCLR: " << eri3_storage[2] << "\n";
                }

                auto dL_inv
                      = TA::to_new_tile_type(L_inv, to_dtile(clr_threshold));

                world.gop.fence();
                auto old_compress = tensor::detail::recompress;
                tensor::detail::recompress = true;

                auto B0 = mpqc_time::fenced_now(world);
                B("X,mu,nu") = dL_inv("X,Y") * eri3("Y,mu,nu");
                if (clr_threshold != 0) {
                    TA::foreach_inplace(
                          B, [](tensor::Tile<tensor::DecomposedTensor<double>> &
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

                              // Both are always larger than or equal to the
                              // real
                              // norm.
                              return std::min(input_norm, compressed_norm);
                          });
                } else {
                    B.truncate();
                }
                auto B1 = mpqc_time::fenced_now(world);
                auto Btime = mpqc_time::duration_in_s(B0, B1);

                if (world.rank() == 0) {
                    std::cout << "B time = " << Btime << std::endl;
                }
                b_storage = utility::array_storage(B);
                if (world.rank() == 0) {
                    std::cout << "B Storage:\n";
                    std::cout << "\tAll Full: " << b_storage[0] << "\n";
                    std::cout << "\tBlock Sparse Only: " << b_storage[1]
                              << "\n";
                    std::cout << "\tCLR: " << b_storage[2] << "\n";
                }

                // Reset the compression flag
                tensor::detail::recompress = old_compress;
            }


            ThreeCenterScf scf(H, F_soad, S, V_inv_oh, r_xyz, occ / 2,
                               repulsion_energy, clr_threshold);


            world.gop.fence();
            auto old_thresh = TA::SparseShape<float>::threshold();
            TA::SparseShape<float>::threshold(std::min(1e-16f, old_thresh));

            auto eri3 = ints::direct_sparse_integrals(world, eri_e,
                                                      three_c_array, shr_screen,
                                                      decomp_3d_no_compress);

            world.gop.fence();
            TA::SparseShape<float>::threshold(old_thresh);

            scf.solve(30, 1e-11, eri3, B);

            job_summary = scf.init_summary();

            // Prepare for Dipole and MP2
            Fao = scf.fock_matrix();

            // Add eri3 info
            job_summary.eri3_fully_dense_storage = eri3_storage[0];
            job_summary.eri3_sparse_only_storage = eri3_storage[1];
            job_summary.eri3_sparse_clr_storage = eri3_storage[2];

            job_summary.B_fully_dense_storage = b_storage[0];
            job_summary.B_sparse_only_storage = b_storage[1];
            job_summary.B_sparse_clr_storage = b_storage[2];

            job_summary.schwarz_threshold = schwarz_thresh;
            job_summary.clr_threshold = clr_threshold;

            job_summary.obs = basis_name;
            job_summary.dfbs = df_basis_name;

            job_summary.num_obs = nclusters;
            job_summary.num_dfbs = df_nclusters;

            job_summary.obs_nfuncs = basis.nfunctions();
            job_summary.dfbs_nfuncs = df_basis.nfunctions();
        }

        auto mp2_occ = occ / 2;
        auto mp2_vir = basis.nfunctions() - mp2_occ;


        auto F_eig = array_ops::array_to_eigen(Fao);
        auto S_eig = array_ops::array_to_eigen(S);

        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);

        auto vec_ptr = std::make_shared<Eig::VectorXd>(es.eigenvalues());

        decltype(S_eig) Ceigi = es.eigenvectors().leftCols(mp2_occ);
        decltype(S_eig) Ceigv = es.eigenvectors().rightCols(mp2_vir);

        auto tr_occ = scf::tr_occupied(mp2_occ, mp2_occ);
        auto tr_vir = scf::tr_occupied(clustered_mol.nclusters(), mp2_vir);

        auto Ci = array_ops::eigen_to_array<TA::TensorD>(
              world, Ceigi, S.trange().data()[0], tr_occ);

        auto Cv = array_ops::eigen_to_array<TA::TensorD>(
              world, Ceigv, S.trange().data()[0], tr_vir);

        // Compute Dipole Moment
        try {
            decltype(S) D;
            D("mu,nu") = Ci("mu, i") * Ci("nu,i");

            double ex = -2 * (D("i,j") * r_xyz[0]("i,j")).sum();
            double ey = -2 * (D("i,j") * r_xyz[1]("i,j")).sum();
            double ez = -2 * (D("i,j") * r_xyz[2]("i,j")).sum();

            double nx = 0;
            double ny = 0;
            double nz = 0;
            for (auto const &atom : clustered_mol.atoms()) {
                nx += atom.charge() * atom.center()[0];
                ny += atom.charge() * atom.center()[1];
                nz += atom.charge() * atom.center()[2];
            }

            const auto au_to_debye = 1.0 / 0.393430307;
            auto dx = (nx + ex) * au_to_debye;
            auto dy = (ny + ey) * au_to_debye;
            auto dz = (nz + ez) * au_to_debye;

            double dip_mom = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2)
                                       + std::pow(dz, 2));

            if (world.rank() == 0) {
                std::cout << "Dipole Moment:"
                          << "\n\tx: " << (nx + ex) << "\n\ty: " << (ny + ey)
                          << "\n\tz: " << (nz + ez) << "\n\tvec: " << dip_mom
                          << std::endl;
            }

            job_summary.computed_dipole = true;
            job_summary.dipole_x = dx;
            job_summary.dipole_y = dy;
            job_summary.dipole_z = dz;
            job_summary.dipole_vector = dip_mom;
        } catch (...) {
            job_summary.print_report();
        }

        // Do MP2
        auto approx_t_elems = mp2_occ * mp2_occ * mp2_vir * mp2_vir;
        auto approx_t_storage = approx_t_elems * 8 * 1e-9;

// Only Do MP2 if we can stay below 60 GB
#if 1
        try {
            if (world.rank() == 0) {
                std::cout << "Approximate T storage = " << approx_t_storage
                          << " GB" << std::endl;
                std::cout << "Max allowed MP2 storage = " << mp2_mem_thresh
                          << " GB" << std::endl;
            }
            if (approx_t_storage < mp2_mem_thresh) {
                // For MP2 don't worry about direct integrals or CLR
                auto eri3_mp2
                      = ints::sparse_integrals(world, eri_e, three_c_array,
                                               shr_screen);

                DArray<3, TA::TensorD, TA::SparsePolicy> Wiv;
                Wiv("Y,i,a") = eri3_mp2("Y,nu,mu") * Ci("mu,i") * Cv("nu,a");
                Wiv("X,i,a") = L_inv("X,Y") * Wiv("Y,i,a");

                // struct to perform the MP2 Reduction
                struct Mp2Red {
                    using result_type = double;
                    using argument_type = TA::Tensor<double>;

                    std::shared_ptr<Eig::VectorXd> vec_;
                    unsigned int n_occ_;

                    Mp2Red(std::shared_ptr<Eig::VectorXd> vec, int n_occ)
                            : vec_(std::move(vec)), n_occ_(n_occ) {}
                    Mp2Red(Mp2Red const &) = default;

                    result_type operator()() const { return 0.0; }
                    result_type operator()(result_type const &t) const {
                        return t;
                    }
                    void operator()(result_type &me,
                                    result_type const &other) const {
                        me += other;
                    }

                    void operator()(result_type &me,
                                    argument_type const &tile) const {
                        auto const &range = tile.range();
                        auto const &vec = *vec_;
                        auto const st = range.lobound();
                        auto const fn = range.upbound();
                        auto tile_idx = 0;
                        for (auto i = st[0]; i < fn[0]; ++i) {
                            const auto e_i = vec[i];
                            for (auto a = st[1]; a < fn[1]; ++a) {
                                const auto e_ia = e_i - vec[a + n_occ_];
                                for (auto j = st[2]; j < fn[2]; ++j) {
                                    const auto e_iaj = e_ia + vec[j];
                                    for (auto b = st[3]; b < fn[3];
                                         ++b, ++tile_idx) {
                                        const auto e_iajb = e_iaj
                                                            - vec[b + n_occ_];
                                        me += 1
                                              / (e_iajb)*tile.data()[tile_idx];
                                    }
                                }
                            }
                        }
                    }
                };


                DArray<4, TA::TensorD, TA::SparsePolicy> G;
                G("i,a,j,b") = Wiv("X,i,a") * Wiv("X,j,b");

                double energy_mp2
                      = (G("i,a,j,b") * (2 * G("i,a,j,b") - G("i,b,j,a")))
                              .reduce(Mp2Red(vec_ptr, mp2_occ))
                              .get();

                if (world.rank() == 0) {
                    std::cout << "Energy from MP2: " << energy_mp2 << std::endl;
                }

                job_summary.did_mp2 = true;
                job_summary.mp2_energy = energy_mp2;
                if(world.rank() == 0){
                    job_summary.print_report();
                }
                world.gop.fence();
            } else {
                if(world.rank() == 0){
                    job_summary.print_report();
                }
                world.gop.fence();
            }
        } catch (...) {
            job_summary.print_report();
        }
#endif
    }

    return 0;
}
