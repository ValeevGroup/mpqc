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
#include "../utility/array_storage.h"
#include "../ta_routines/array_to_eigen.h"

#include "../scf/diagonalize_for_coffs.hpp"
#include "../scf/soad.h"
#include "../scf/orbital_localization.h"
#include "../scf/clusterd_coeffs.h"

#include "../tensor/decomposed_tensor.h"
#include "../tensor/decomposed_tensor_unary.h"
#include "../tensor/decomposed_tensor_algebra.h"
#include "../tensor/decomposed_tensor_gemm.h"
#include "../tensor/decomposed_tensor_addition.h"
#include "../tensor/decomposed_tensor_subtraction.h"
#include "../tensor/decomposed_tensor_multiplication.h"
#include "../tensor/tcc_tile.h"

#include <memory>

using namespace mpqc;
namespace ints = mpqc::integrals;

bool tcc::tensor::detail::recompress = true;

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

            tcc::tensor::Tile<tcc::tensor::DecomposedTensor<double>>
                  wrapper_tile = a.find(ord).get();

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
    // Times
    double avg_w_time = 0;
    double avg_occ_ri_k_time = 0;
    double avg_k_time = 0;

    // Storages
    double w_fully_dense_storage = 0;
    double avg_w_sparse_only_storage = 0;
    double avg_w_sparse_clr_storage = 0;

    double eri3_fully_dense_storage = 0;
    double eri3_sparse_only_storage = 0;
    double eri3_sparse_clr_storage = 0;

    // Thresholds
    double ta_sparse_threshold = 0;
    double clr_threshold = 0;
    double schwarz_threshold = 0;

    // Computables
    double hf_energy = 0;
    double hf_energy_diff_error = 0;
    double hf_grad_div_volume = 0;

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

        std::cout << "ta sparse threshold, ";
        std::cout << "schwarz threshold, ";
        std::cout << "clr threshold, ";

        std::cout << "avg w time, ";
        std::cout << "avg k time, ";
        std::cout << "avg occ-ri k time, ";

        std::cout << "eri3 fully dense storage, ";
        std::cout << "eri3 sparse only, ";
        std::cout << "eri3 sparse + clr storage, ";

        std::cout << "w fully dense storage, ";
        std::cout << "w sparse only, ";
        std::cout << "w sparse + clr storage, ";

        std::cout << "hf energy, ";
        std::cout << "error change in energy, ";
        std::cout << "error rms div volume, ";

        std::cout << "computed dipole, ";

        std::cout << "dipole x, ";
        std::cout << "dipole y, ";
        std::cout << "dipole z, ";
        std::cout << "dipole vec, ";

        std::cout << "computed mp2, ";
        std::cout << "mp2 energy\n";

        std::cout << ta_sparse_threshold << ", ";
        std::cout << schwarz_threshold << ", ";
        std::cout << clr_threshold << ", ";


        std::cout << avg_w_time << ", ";
        std::cout << avg_k_time << ", ";
        std::cout << avg_occ_ri_k_time << ", ";

        std::cout << eri3_fully_dense_storage << ", ";
        std::cout << eri3_sparse_only_storage << ", ";
        std::cout << eri3_sparse_clr_storage << ", ";

        std::cout << w_fully_dense_storage << ", ";
        std::cout << avg_w_sparse_only_storage << ", ";
        std::cout << avg_w_sparse_clr_storage << ", ";

        std::cout << hf_energy << ", ";
        std::cout << hf_energy_diff_error << ", ";
        std::cout << hf_grad_div_volume << ", ";

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
    using dtile = tcc::tensor::Tile<tcc::tensor::DecomposedTensor<double>>;
    to_dtile(double clr_thresh) : clr_thresh_(clr_thresh) {}

    dtile operator()(TA::TensorD const &tile) {
        auto range = tile.range();

        const auto extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];

        auto local_range = TA::Range{i, j};
        auto tensor = TA::TensorD(local_range, tile.data());

        auto new_tile
              = tcc::tensor::DecomposedTensor<double>(clr_thresh_,
                                                      std::move(tensor));

        return dtile(range, std::move(new_tile));
    }
};

class to_dtile_compress {
    double clr_thresh_;

  public:
    using dtile = tcc::tensor::Tile<tcc::tensor::DecomposedTensor<double>>;
    to_dtile_compress(double clr_thresh) : clr_thresh_(clr_thresh) {}

    dtile operator()(TA::TensorD const &tile) {
        auto range = tile.range();

        const auto extent = range.extent();
        const auto i = extent[0];
        const auto j = extent[1];

        auto local_range = TA::Range{i, j};
        auto tensor = TA::TensorD(local_range, tile.data());

        auto new_tile
              = tcc::tensor::DecomposedTensor<double>(clr_thresh_,
                                                      std::move(tensor));

        if (clr_thresh_ != 0.0) {
            auto test = tcc::tensor::algebra::two_way_decomposition(new_tile);
            if (!test.empty()) {
                new_tile = std::move(test);
            }
        }

        return dtile(range, std::move(new_tile));
    }
};

class to_ta_tile {
    using dtile = tcc::tensor::Tile<tcc::tensor::DecomposedTensor<double>>;

  public:
    TA::TensorD operator()(dtile const &t) const {
        auto tensor = tcc::tensor::algebra::combine(t.tile());
        auto range = t.range();
        return TA::Tensor<double>(range, tensor.data());
    }
};

class ThreeCenterScf {
  private:
    using array_type = DArray<2, TA::TensorD, SpPolicy>;
    using dtile = tcc::tensor::Tile<tcc::tensor::DecomposedTensor<double>>;
    using darray_type = DArray<2, dtile, SpPolicy>;

    array_type H_;
    array_type S_;

    array_type F_;
    array_type D_;
    array_type C_;
    array_type L_invV_;

    darray_type dL_invV_;
    darray_type dV_inv_;
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

    int64_t occ_;
    double repulsion_;
    double final_rms_error_;
    double final_ediff_error_;
    double final_energy_;
    double w_full_storage_;

    double recompress_w_time_;
    double clr_w_no_recompress_;

    double clr_thresh_;

    void compute_density(int64_t occ) {
        auto F_eig = tcc::array_ops::array_to_eigen(F_);
        auto S_eig = tcc::array_ops::array_to_eigen(S_);

        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);
        decltype(S_eig) C = es.eigenvectors().leftCols(occ);
        auto tr_ao = S_.trange().data()[0];

        auto occ_nclusters = (occ_ < 10) ? occ_ : 10;
        auto tr_occ = tcc::scf::tr_occupied(occ_nclusters, occ_);

        C_ = tcc::array_ops::eigen_to_array<TA::TensorD>(H_.get_world(), C,
                                                         tr_ao, tr_occ);

        auto U = mpqc::scf::BoysLocalization{}(C_, r_xyz_ints_);
        C_("mu,i") = C_("mu,k") * U("k,i");
        auto obs_ntiles = C_.trange().tiles().extent()[0];
        tcc::scf::clustered_coeffs(r_xyz_ints_, C_, obs_ntiles);

        D_("i,j") = C_("i,k") * C_("j,k");
    }

    template <typename Integral>
    void form_fock(Integral const &eri3) {
        auto &world = F_.get_world();

        world.gop.fence();
        auto w0 = tcc::utility::time::now();

        darray_type dC_ = TA::to_new_tile_type(C_, to_dtile(clr_thresh_));

        DArray<3, dtile, SpPolicy> W;
        W("Y, mu, i") = eri3("Y, mu, nu") * dC_("nu, i");

        world.gop.fence();
        auto w1 = tcc::utility::time::now();
        w_times_.push_back(tcc::utility::time::duration_in_s(w0, w1));

        auto w_store_no_recomp = storage_for_array(W);
        clr_w_no_recompress_ = w_store_no_recomp[2];

        world.gop.fence();
        auto r0 = tcc::utility::time::now();
        if (tcc::tensor::detail::recompress && clr_thresh_ != 0) {
            TA::foreach_inplace(
                  W,
                  [](tcc::tensor::Tile<tcc::tensor::DecomposedTensor<double>> &
                           t_tile) {
                      auto &t = t_tile.tile();
                      auto input_norm = norm(t);

                      auto compressed_norm = input_norm;
                      if (t.cut() != 0.0) {
                          if (t.ndecomp() == 1) {
                              auto test = tcc::tensor::algebra::
                                    two_way_decomposition(t);
                              if (!test.empty()) {
                                  t = test;
                                  compressed_norm = norm(t);
                              }
                          } else {
                              tcc::tensor::algebra::recompress(t);
                              compressed_norm = norm(t);
                          }
                      }

                      // Both are always larger than or equal to the real norm.
                      return std::min(input_norm, compressed_norm);
                  });
        } else {
            W.truncate();
        }
        world.gop.fence();
        auto r1 = tcc::utility::time::now();
        recompress_w_time_ = tcc::utility::time::duration_in_s(r0, r1);

        auto w_store = storage_for_array(W);
        w_sparse_store_.push_back(w_store[1]);
        w_sparse_clr_store_.push_back(w_store[2]);
        w_full_storage_ = w_store[0];

        darray_type J;
        // J("mu, nu") = eri3("X, mu, nu")
        //             * (dV_inv_("X, Y") * (W("Y, rho, i") * dC_("rho, i")));
        J("mu, nu")
              = eri3("X, mu, nu")
                * (dL_invV_("X,Z")
                   * (
                      dL_invV_("Z, Y") * (
                          W("Y, rho, i") * dC_("rho, i"))));
        world.gop.fence();
        auto j1 = tcc::utility::time::now();
        j_times_.push_back(tcc::utility::time::duration_in_s(w1, j1));

        // Permute W
        W("X,i,nu") = W("X,nu,i");
        world.gop.fence();
        auto occk0 = tcc::utility::time::now();

        darray_type dS_ = TA::to_new_tile_type(S_, to_dtile(clr_thresh_));

        darray_type K, Kij, Sc;
        K("mu, j") = W("X, i, mu")
                     * (dV_inv_("X,Y") * (W("Y, i, nu") * dC_("nu, j")));
        Kij("i,j") = dC_("mu, i") * K("mu, j");
        Sc("mu, j") = dS_("mu, lam") * dC_("lam, j");
        K("mu, nu") = Sc("mu, j") * K("nu, j") + K("mu, j") * Sc("nu,j")
                      - (Sc("mu,i") * Kij("i,j")) * Sc("nu,j");
        world.gop.fence();
        auto occk1 = tcc::utility::time::now();
        occ_k_times_.push_back(tcc::utility::time::duration_in_s(occk0, occk1));

        K("mu, nu") = W("X, i, mu") * (dV_inv_("X,Y") * W("Y, i, nu"));
        // W("X,i,mu") = dL_invV_("X,Y") * W("Y,i,mu");
        // K("mu, nu") = W("X,i,mu") * W("X,i,nu");
        world.gop.fence();
        auto k1 = tcc::utility::time::now();
        k_times_.push_back(tcc::utility::time::duration_in_s(occk1, k1));

        darray_type dH_ = TA::to_new_tile_type(H_, to_dtile(clr_thresh_));

        darray_type dF_;
        dF_("i,j") = dH_("i,j") + 2 * J("i,j") - K("i,j");

        F_ = TA::to_new_tile_type(dF_, to_ta_tile{});
    }

  public:
    ThreeCenterScf(array_type const &H, array_type const &F_guess,
                   array_type const &S, array_type const &L_invV,
                   std::vector<array_type> const &rxyz, int64_t occ, double rep,
                   double clr_thresh)
            : H_(H),
              F_(F_guess),
              S_(S),
              L_invV_(L_invV),
              r_xyz_ints_(rxyz),
              occ_(occ),
              repulsion_(rep),
              clr_thresh_(clr_thresh) {
        dL_invV_
              = TA::to_new_tile_type(L_invV_, to_dtile_compress(clr_thresh_));
        auto ll_sizes = storage_for_array(dL_invV_);
        std::cout << "L_inv storage:"
                  << "\n\tFull: " << ll_sizes[0]
                  << "\n\tSparse Only: " << ll_sizes[1]
                  << "\n\tSparse + CLR: " << ll_sizes[2] << std::endl;

        decltype(L_invV_) V_inv;
        V_inv("X,Z") = L_invV_("Y,X") * L_invV_("Y,Z");
        dV_inv_ = TA::to_new_tile_type(V_inv, to_dtile_compress(clr_thresh_));

        auto dl_sizes = storage_for_array(dV_inv_);
        std::cout << "V_inv storage:"
                  << "\n\tFull: " << dl_sizes[0]
                  << "\n\tSparse Only: " << dl_sizes[1]
                  << "\n\tSparse + CLR: " << dl_sizes[2] << std::endl;

        compute_density(occ_);
    }

    template <typename Integral>
    void solve(int64_t max_iters, double thresh, Integral const &eri3) {
        auto iter = 0;
        auto error = std::numeric_limits<double>::max();
        auto rms_error = std::numeric_limits<double>::max();
        auto old_energy = 0.0;
        const double volume = F_.trange().elements().volume();

        while (iter < max_iters && (thresh < error || thresh < rms_error)
               && (thresh / 100.0 < error && thresh / 100.0 < rms_error)) {
            auto s0 = tcc_time::now();
            F_.get_world().gop.fence();
            form_fock(eri3);

            auto current_energy = energy();
            error = std::abs(old_energy - current_energy);
            old_energy = current_energy;

            array_type Grad;
            Grad("i,j") = F_("i,k") * D_("k,l") * S_("l,j")
                          - S_("i,k") * D_("k,l") * F_("l,j");

            rms_error = Grad("i,j").norm().get() / volume;

            diis_.extrapolate(F_, Grad);

            // Lastly update density
            compute_density(occ_);

            F_.get_world().gop.fence();
            auto s1 = tcc_time::now();
            scf_times_.push_back(tcc_time::duration_in_s(s0, s1));


            std::cout << "Iteration: " << (iter + 1)
                      << " energy: " << old_energy << " error: " << error
                      << " RMS error: " << rms_error;
            std::cout << "\n\tW time: " << w_times_.back()
                      << "\n\tW recompress time: " << recompress_w_time_
                      << "\n\tJ time: " << j_times_.back()
                      << "\n\tocc-RI K time: " << occ_k_times_.back()
                      << "\n\tK time: " << k_times_.back()
                      << "\n\titer time: " << scf_times_.back()
                      << "\n\tW sparse only storage: " << w_sparse_store_.back()
                      << "\n\tW sparse clr storage no recompress: "
                      << clr_w_no_recompress_ << "\n\tW sparse clr storage: "
                      << w_sparse_clr_store_.back() << std::endl;

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
        tcc::tensor::detail::recompress = false;
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

    auto df_nclusters = std::max(nclusters/2, 1);
    auto df_clustered_mol = molecule::attach_hydrogens_and_kmeans(
          molecule::read_xyz(mol_file).clusterables(), df_nclusters);

    auto cluster_num = 1;
    std::cout <<"\nDF clustering:\n";
    for(auto const &c : df_clustered_mol.clusterables()){
        std::cout << "Cluster " << cluster_num++ << ":\n";
        auto atoms = c.atoms();
        std::cout << "\t" << atoms.size() << "\n";
        for(auto const &atom : c.atoms()){
            std::cout << "\t" << atom.xyz_string(true) << "\n";
        }
    }
    std::cout << "\n";

    basis::BasisSet dfbs(df_basis_name);
    basis::Basis df_basis(dfbs.get_cluster_shells(df_clustered_mol));

    libint2::init();

    const auto bs_array = tcc::utility::make_array(basis, basis);

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

    const auto dfbs_array = tcc::utility::make_array(df_basis, df_basis);
    auto eri_e = ints::make_2body_shr_pool(df_basis, basis);

    decltype(H) L_inv;
    {
        auto Vmetric = ints::sparse_integrals(world, eri_e, dfbs_array);
        auto V_eig = tcc::array_ops::array_to_eigen(Vmetric);

        // MatrixD Leig = Eig::LLT<MatrixD>(V_eig).matrixL();
        // MatrixD L_inv_eig = Leig.inverse();
        Eig::SelfAdjointEigenSolver<MatrixD> es(V_eig);
        MatrixD L_inv_eig = es.operatorInverseSqrt();

        auto tr_V = Vmetric.trange().data()[0];
        L_inv = tcc::array_ops::eigen_to_array<TA::TensorD>(world, L_inv_eig,
                                                            tr_V, tr_V);
    }

    auto three_c_array = tcc::utility::make_array(df_basis, basis, basis);

    { // Schwarz Screened ints
        const auto schwarz_thresh = 1e-12;
        std::cout << "Schwarz Threshold: " << schwarz_thresh << std::endl;
        auto sbuilder = ints::init_schwarz_screen(schwarz_thresh);
        auto shr_screen = std::make_shared<ints::SchwarzScreen>(
              sbuilder(world, eri_e, df_basis, basis));

        world.gop.fence();
        auto soad0 = tcc::utility::time::now();
        auto F_soad
              = scf::fock_from_soad(world, clustered_mol, basis, eri_e, H);

        auto soad1 = tcc::utility::time::now();
        auto soad_time = tcc::utility::time::duration_in_s(soad0, soad1);
        std::cout << "Soad Time: " << soad_time << std::endl;

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

            tcc::tensor::DecomposedTensor<double> dc_tile(clr_threshold,
                                                          std::move(tensor));

            if (clr_threshold != 0.0) {
                auto lr_tile
                      = tcc::tensor::algebra::two_way_decomposition(dc_tile);

                if (!lr_tile.empty()) {
                    dc_tile = std::move(lr_tile);
                }
            }

            return tcc::tensor::Tile<decltype(dc_tile)>(range,
                                                        std::move(dc_tile));
        };
        if (!use_direct_ints) {
            world.gop.fence();
            auto int0 = tcc::utility::time::now();
            auto eri3 = ints::sparse_integrals(world, eri_e, three_c_array,
                                               shr_screen, decomp_3d);

            world.gop.fence();
            auto int1 = tcc::utility::time::now();
            auto eri3_time = tcc::utility::time::duration_in_s(int0, int1);
            std::cout << "Eri3 Integral time: " << eri3_time << std::endl;
            auto eri3_storage = storage_for_array(eri3);
            std::cout << "Eri3 Integral Storage:\n";
            std::cout << "\tAll Full: " << eri3_storage[0] << "\n";
            std::cout << "\tBlock Sparse Only: " << eri3_storage[1] << "\n";
            std::cout << "\tCLR: " << eri3_storage[2] << "\n";

            ThreeCenterScf scf(H, F_soad, S, L_inv, r_xyz, occ / 2,
                               repulsion_energy, clr_threshold);

            scf.solve(30, 1e-11, eri3);

            job_summary = scf.init_summary();

            // Prepare for Dipole and MP2
            Fao = scf.fock_matrix();

            // Add eri3 info
            job_summary.eri3_fully_dense_storage = eri3_storage[0];
            job_summary.eri3_sparse_only_storage = eri3_storage[1];
            job_summary.eri3_sparse_clr_storage = eri3_storage[2];
            job_summary.schwarz_threshold = schwarz_thresh;
            job_summary.clr_threshold = clr_threshold;
        }
        if (use_direct_ints) {
            world.gop.fence();
            auto int0 = tcc::utility::time::now();
            auto eri3
                  = ints::direct_sparse_integrals(world, eri_e, three_c_array,
                                                  shr_screen, decomp_3d);

            world.gop.fence();
            auto int1 = tcc::utility::time::now();
            auto eri3_time = tcc::utility::time::duration_in_s(int0, int1);
            std::cout << "Eri3 Integral time: " << eri3_time << std::endl;
            auto eri3_storage = storage_for_array(eri3.array());
            std::cout << "Eri3 Integral Storage:\n";
            std::cout << "\tAll Full: " << eri3_storage[0] << "\n";
            std::cout << "\tBlock Sparse Only: " << eri3_storage[1] << "\n";
            std::cout << "\tCLR: " << eri3_storage[2] << "\n";

            ThreeCenterScf scf(H, F_soad, S, L_inv, r_xyz, occ / 2,
                               repulsion_energy, clr_threshold);

            scf.solve(30, 1e-11, eri3);

            job_summary = scf.init_summary();

            // Prepare for Dipole and MP2
            Fao = scf.fock_matrix();

            // Add eri3 info
            job_summary.eri3_fully_dense_storage = eri3_storage[0];
            job_summary.eri3_sparse_only_storage = eri3_storage[1];
            job_summary.eri3_sparse_clr_storage = eri3_storage[2];
            job_summary.schwarz_threshold = schwarz_thresh;
            job_summary.clr_threshold = clr_threshold;
        }

        auto mp2_occ = occ / 2;
        auto mp2_vir = basis.nfunctions() - mp2_occ;


        auto F_eig = tcc::array_ops::array_to_eigen(Fao);
        auto S_eig = tcc::array_ops::array_to_eigen(S);

        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);

        auto vec_ptr = std::make_shared<Eig::VectorXd>(es.eigenvalues());

        decltype(S_eig) Ceigi = es.eigenvectors().leftCols(mp2_occ);
        decltype(S_eig) Ceigv = es.eigenvectors().rightCols(mp2_vir);

        auto tr_occ = tcc::scf::tr_occupied(clustered_mol.nclusters(), mp2_occ);
        auto tr_vir = tcc::scf::tr_occupied(clustered_mol.nclusters(), mp2_vir);

        auto Ci = tcc::array_ops::eigen_to_array<TA::TensorD>(
              world, Ceigi, S.trange().data()[0], tr_occ);

        auto Cv = tcc::array_ops::eigen_to_array<TA::TensorD>(
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

            std::cout << "Dipole Moment:"
                      << "\n\tx: " << (nx + ex) << "\n\ty: " << (ny + ey)
                      << "\n\tz: " << (nz + ez) << "\n\tvec: " << dip_mom
                      << std::endl;

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
            std::cout << "Approximate T storage = " << approx_t_storage << " GB"
                      << std::endl;
            std::cout << "Max allowed MP2 storage = " << mp2_mem_thresh << " GB"
                      << std::endl;
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

                std::cout << "Energy from MP2: " << energy_mp2 << std::endl;

                job_summary.did_mp2 = true;
                job_summary.mp2_energy = energy_mp2;
                job_summary.print_report();
            } else {
                job_summary.print_report();
            }
        } catch (...) {
            job_summary.print_report();
        }
#endif
    }

    return 0;
}
