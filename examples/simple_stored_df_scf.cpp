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

#include <memory>

using namespace mpqc;
namespace ints = mpqc::integrals;

template <typename Array>
std::array<double, 3> storage_for_array(Array const &a) {
    std::atomic<long> full_a(0.0);
    std::atomic<long> sparse_a(0.0);
    std::atomic<long> clr_a(0.0);

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
    double clr = clr_a * 1e-9 * 8;

    return {full, sparse, clr};
}

class ThreeCenterScf {
  private:
    using array_type = DArray<2, TA::TensorD, SpPolicy>;
    array_type H_;
    array_type S_;

    array_type F_;
    array_type D_;
    array_type C_;
    array_type L_invV_;

    std::vector<array_type> r_xyz_ints_;
    TiledArray::DIIS<array_type> diis_;

    std::vector<double> k_times_;
    std::vector<double> occ_k_times_;
    std::vector<double> j_times_;
    std::vector<double> w_times_;
    std::vector<double> scf_times_;
    std::vector<double> w_sparse_store_;
    std::vector<double> c_sparse_store_;

    int64_t occ_;
    double repulsion_;
    double final_rms_error_;
    double final_ediff_error_;
    double final_energy_;
    double w_full_storage_;
    double c_full_storage_;

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

        auto store = storage_for_array(C_);
        c_sparse_store_.push_back(store[1]);
        c_full_storage_ = store[0];


        D_("i,j") = C_("i,k") * C_("j,k");
    }

    template <typename Integral>
    void form_fock(Integral const &eri3) {
        auto &world = F_.get_world();

        world.gop.fence();
        auto w0 = tcc::utility::time::now();
        TA::Array<double, 3, TA::TensorD, TA::SparsePolicy> W;
        W("X, mu, i") = L_invV_("X,Y") * (eri3("Y, mu, nu") * C_("nu, i"));
        W.truncate();
        world.gop.fence();
        auto w1 = tcc::utility::time::now();
        w_times_.push_back(tcc::utility::time::duration_in_s(w0, w1));

        auto w_store = storage_for_array(W);
        w_sparse_store_.push_back(w_store[1]);
        w_full_storage_ = w_store[0];


        array_type J;
        J("mu, nu") = eri3("X, mu, nu")
                      * (L_invV_("Y, X") * (W("Y, rho, i") * C_("rho, i")));
        world.gop.fence();
        auto j1 = tcc::utility::time::now();
        j_times_.push_back(tcc::utility::time::duration_in_s(w1, j1));

        // Permute W
        W("X,i,nu") = W("X,nu,i");
        world.gop.fence();
        auto occk0 = tcc::utility::time::now();

        array_type K, Kij, Sc;
        K("mu, j") = W("X, i, mu") * (W("X, i, nu") * C_("nu, j"));
        Kij("i,j") = C_("mu, i") * K("mu, j");
        Sc("mu, j") = S_("mu, lam") * C_("lam, j");
        K("mu, nu") = Sc("mu, j") * K("nu, j") + K("mu, j") * Sc("nu,j")
                      - (Sc("mu,i") * Kij("i,j")) * Sc("nu,j");
        world.gop.fence();
        auto occk1 = tcc::utility::time::now();
        occ_k_times_.push_back(tcc::utility::time::duration_in_s(occk0, occk1));

        K("mu, nu") = W("X, i, mu") * W("X, i, nu");
        world.gop.fence();
        auto k1 = tcc::utility::time::now();
        k_times_.push_back(tcc::utility::time::duration_in_s(occk1, k1));

        F_("i,j") = H_("i,j") + 2 * J("i,j") - K("i,j");
    }


  public:
    ThreeCenterScf(array_type const &H, array_type const &F_guess,
                   array_type const &S, array_type const &L_invV,
                   std::vector<array_type> const &rxyz, int64_t occ, double rep)
            : H_(H),
              F_(F_guess),
              S_(S),
              L_invV_(L_invV),
              r_xyz_ints_(rxyz),
              occ_(occ),
              repulsion_(rep) {
        compute_density(occ_);
    }

    template <typename Integral>
    void solve(int64_t max_iters, double thresh, Integral const &eri3) {
        auto iter = 0;
        auto error = std::numeric_limits<double>::max();
        auto rms_error = std::numeric_limits<double>::max();
        auto old_energy = 0.0;
        const double volume = F_.trange().elements().volume();

        while (iter < max_iters && thresh < error && thresh < rms_error) {
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
            std::cout << "\n\tW time: " << w_times_.back();
            std::cout << "\n\tJ time: " << j_times_.back()
                      << "\n\tocc-RI K time: " << occ_k_times_.back()
                      << "\n\tK time: " << k_times_.back()
                      << "\n\titer time: " << scf_times_.back()
                      << "\n\tW storage: " << w_sparse_store_.back()
                      << "\n\tC storage: " << c_sparse_store_.back() << "\n";

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

    void print_summary() {
        auto avg = [](std::vector<double> const &v) {
            auto sum = 0.0;
            for (auto const &e : v) {
                sum += e;
            }
            return sum / v.size();
        };

        auto avg_w = avg(w_times_);
        auto avg_occ_ri_k = avg(occ_k_times_);
        auto avg_k = avg(k_times_);
        auto avg_ws = avg(w_sparse_store_);
        auto avg_cs = avg(c_sparse_store_);
        std::cout << "Summary:\n";
        std::cout
              << "energy, energy diff, (rms error)/volume, avg w time, avg "
                 "occ-ri k time, avg k time, all dense w storage, all dense c "
                 "stroage, avg w storage, avg c storage\n";
        std::cout << final_energy_ << ", " << final_ediff_error_ << ", "
                  << final_rms_error_ << ", " << avg_w << ", " << avg_occ_ri_k
                  << ", " << avg_k << ", " << w_full_storage_ << ", "
                  << c_full_storage_ << ", " << avg_ws << ", " << avg_cs
                  << "\n";
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
    if (argc == 6) {
        mol_file = argv[1];
        basis_name = argv[2];
        df_basis_name = argv[3];
        nclusters = std::stoi(argv[4]);
        threshold = std::stod(argv[5]);
    } else {
        std::cout << "input is $./program mol_file basis_file df_basis_file "
                     "nclusters ";
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

    basis::BasisSet dfbs(df_basis_name);
    basis::Basis df_basis(dfbs.get_cluster_shells(clustered_mol));

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

        MatrixD Leig = Eig::LLT<MatrixD>(V_eig).matrixL();
        MatrixD L_inv_eig = Leig.inverse();

        auto tr_V = Vmetric.trange().data()[0];
        L_inv = tcc::array_ops::eigen_to_array<TA::TensorD>(world, L_inv_eig,
                                                            tr_V, tr_V);
    }

    auto three_c_array = tcc::utility::make_array(df_basis, basis, basis);
    { // Schwarz Screened ints
        auto int0 = tcc::utility::time::now();
        const auto thresh = 1e-12;
        std::cout << "Schwarz Threshold: " << thresh << std::endl;
        auto sbuilder = ints::init_schwarz_screen(thresh);
        auto shr_screen = std::make_shared<ints::SchwarzScreen>(
              sbuilder(world, eri_e, df_basis, basis));

        auto eri3
              = ints::sparse_integrals(world, eri_e, three_c_array, shr_screen);
        world.gop.fence();
        auto int1 = tcc::utility::time::now();
        auto eri3_time = tcc::utility::time::duration_in_s(int0, int1);
        std::cout << "Eri3 Integral time: " << eri3_time << std::endl;
        auto eri3_storage = storage_for_array(eri3);
        std::cout << "Eri3 Integral Storage:\n";
        std::cout << "\tAll Full: " << eri3_storage[0] << "\n";
        std::cout << "\tBlock Sparse Only: " << eri3_storage[1] << "\n";
        std::cout << "\tCLR: " << eri3_storage[2] << "\n";


        auto F_soad
              = scf::fock_from_soad(world, clustered_mol, basis, eri_e, H);

        auto soad1 = tcc::utility::time::now();
        auto soad_time = tcc::utility::time::duration_in_s(int1, soad1);
        std::cout << "Soad Time: " << soad_time << std::endl;

        auto multi_pool
              = ints::make_1body_shr_pool("emultipole2", basis, clustered_mol);
        auto r_xyz = ints::sparse_xyz_integrals(world, multi_pool, bs_array);

        ThreeCenterScf scf(H, F_soad, S, L_inv, r_xyz, occ / 2,
                           repulsion_energy);
        scf.solve(30, 1e-12, eri3);
        scf.print_summary();
    }

    return 0;
}
