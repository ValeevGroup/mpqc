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
#include "../scf/orbital_localization.h"
#include "../scf/soad.h"
#include "../scf/clusterd_coeffs.h"

#include <memory>
#include <atomic>

using namespace mpqc;
namespace ints = mpqc::integrals;

template <typename Array>
void print_storage_for_array_4(Array const &a, double cut = 1e-8) {
    std::atomic<long> full_a(0.0);
    std::atomic<long> sparse_a(0.0);
    std::atomic<long> clr_a(0.0);

    auto task_f = [&](int ord) {
        auto const &trange = a.trange();
        auto size = trange.make_tile_range(ord).volume();
        full_a += size;

        if (!a.is_zero(ord)) {
            sparse_a += size;

            TA::TensorD tile = a.find(ord).get();
            auto const &ext = tile.range().extent();
            MatrixD e_tile(ext[0] * ext[1], ext[2] * ext[3]);

            for (auto i = 0; i < tile.range().volume(); ++i) {
                *(e_tile.data() + i) = *(tile.data() + i);
            }

            Eig::ColPivHouseholderQR<MatrixD> qr(e_tile);
            qr.setThreshold(cut);
            auto rank = qr.rank();

            if (rank > std::min(ext[0] * ext[1], ext[2] * ext[3]) / 2) {
                clr_a += tile.range().volume();
            } else {
                auto l_size = ext[0] * ext[1] * rank;
                auto r_size = ext[2] * ext[3] * rank;
                clr_a += l_size + r_size;
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

    std::cout << "G Full Storage = " << full << " GB" << std::endl;
    std::cout << "G Sparse Storage = " << sparse << " GB" << std::endl;
    std::cout << "G CLR Storage = " << clr << " GB" << std::endl;
}

template <typename Array>
void print_storage_for_array(Array const &a, double cut = 1e-8) {
    std::atomic<long> full_a(0.0);
    std::atomic<long> sparse_a(0.0);
    std::atomic<long> clr_a(0.0);

    auto task_f = [&](int ord) {
        auto const &trange = a.trange();
        auto size = trange.make_tile_range(ord).volume();
        full_a += size;

        if (!a.is_zero(ord)) {
            sparse_a += size;

            TA::TensorD tile = a.find(ord).get();
            tcc::tensor::DecomposedTensor<double> dc_tile(cut, std::move(tile));
            auto lr_tile = tcc::tensor::algebra::two_way_decomposition(dc_tile);

            if (!lr_tile.empty()) {
                dc_tile = std::move(lr_tile);
            }

            if (dc_tile.ndecomp() == 1) {
                clr_a += dc_tile.tensors()[0].range().volume();
            } else {
                auto l_size = dc_tile.tensors()[0].range().volume();
                auto r_size = dc_tile.tensors()[1].range().volume();

                clr_a += l_size + r_size;
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

    std::cout << "Eri3 Full Storage = " << full << " GB" << std::endl;
    std::cout << "Eri3 Sparse Storage = " << sparse << " GB" << std::endl;
    std::cout << "Eri3 CLR Storage = " << clr << " GB" << std::endl;
}


class ThreeCenterScf {
  private:
    using array_type = DArray<2, TA::TensorD, SpPolicy>;
    array_type H_;
    array_type S_;

    array_type F_;
    array_type D_;
    array_type C_;
    array_type V_inv_;
    TiledArray::DIIS<array_type> diis_;

    std::vector<double> k_times_;
    std::vector<double> j_times_;
    std::vector<double> w_times_;
    std::vector<double> scf_times_;

    int64_t occ_;
    double repulsion_;


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

        D_("i,j") = C_("i,k") * C_("j,k");
    }

    template <typename Integral>
    void form_fock(Integral const &eri3) {
        auto &world = F_.get_world();

        world.gop.fence();
        auto w0 = tcc::utility::time::now();
        TA::Array<double, 3, TA::TensorD, TA::SparsePolicy> W;
        W("Y, mu, i") = eri3("Y, mu, nu") * C_("nu, i");
        auto w1 = tcc::utility::time::now();
        w_times_.push_back(tcc::utility::time::duration_in_s(w0, w1));


        array_type J;
        J("mu, nu") = eri3("X, mu, nu")
                      * (V_inv_("X, Y") * (W("Y, rho, i") * C_("rho, i")));
        world.gop.fence();
        auto j1 = tcc::utility::time::now();
        j_times_.push_back(tcc::utility::time::duration_in_s(w1, j1));


        // Permute W
        W("Y,i,nu") = W("Y,nu,i");
        array_type K, Kij, Sc;
        K("mu, j") = W("X, i, mu")
                     * (V_inv_("X,Y") * (W("Y, i, nu") * C_("nu, j")));
        Kij("i,j") = C_("mu, i") * K("mu, j");
        Sc("mu, j") = S_("mu, lam") * C_("lam, j");
        K("mu, nu") = Sc("mu, j") * K("nu, j") + K("mu, j") * Sc("nu,j")
                      - (Sc("mu,i") * Kij("i,j")) * Sc("nu,j");
        world.gop.fence();
        auto k1 = tcc::utility::time::now();
        k_times_.push_back(tcc::utility::time::duration_in_s(j1, k1));

        F_("i,j") = H_("i,j") + 2 * J("i,j") - K("i,j");
    }


  public:
    ThreeCenterScf(array_type const &H, array_type const &S,
                   array_type const &F_guess, array_type const &V_inv,
                   int64_t occ, double rep)
            : H_(H),
              S_(S),
              F_(F_guess),
              V_inv_(V_inv),
              occ_(occ),
              repulsion_(rep) {
        compute_density(occ_);
    }

    array_type fock() const { return F_; }

    template <typename Integral>
    void solve(int64_t max_iters, double thresh, Integral const &eri3) {
        auto iter = 0;
        auto error = std::numeric_limits<double>::max();
        auto old_energy = 0.0;

        while (iter < max_iters && thresh < error) {
            auto s0 = tcc_time::now();
            F_.get_world().gop.fence();
            form_fock(eri3);

            auto current_energy = energy();
            error = std::abs(old_energy - current_energy);
            old_energy = current_energy;

            array_type Grad;
            Grad("i,j") = F_("i,k") * D_("k,l") * S_("l,j")
                          - S_("i,k") * D_("k,l") * F_("l,j");

            diis_.extrapolate(F_, Grad);

            // Lastly update density
            compute_density(occ_);

            F_.get_world().gop.fence();
            auto s1 = tcc_time::now();
            scf_times_.push_back(tcc_time::duration_in_s(s0, s1));


            std::cout << "Iteration: " << (iter + 1)
                      << " energy: " << old_energy << " error: " << error
                      << std::endl;
            std::cout << "\tW time: " << w_times_.back() << std::endl;
            std::cout << "\tJ time: " << j_times_.back()
                      << " s K time: " << k_times_.back()
                      << " s iter time: " << scf_times_.back() << std::endl;

            ++iter;
        }
    }

    double energy() {
        return repulsion_
               + D_("i,j").dot(F_("i,j") + H_("i,j"), D_.get_world()).get();
    }
};


int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    std::string df_basis_name = "";
    int obs_nclusters = 0;
    int dfbs_nclusters = 0;
    double threshold = 1e-11;
    auto rank_thresh = 1e-16;
    std::cout << std::setprecision(15);
    if (argc >= 6) {
        mol_file = argv[1];
        basis_name = argv[2];
        df_basis_name = argv[3];
        obs_nclusters = std::stoi(argv[4]);
        dfbs_nclusters = std::stoi(argv[5]);
    }
    if (argc >= 7) {
        threshold = std::stod(argv[6]);
        std::cout << "Setting sp thresh = " << threshold << std::endl;
    }
    if (argc >= 8) {
        rank_thresh = std::stod(argv[7]);
        std::cout << "Setting rank thresh = " << rank_thresh << std::endl;
    } else {
        std::cout << "input is $./program mol_file basis_file df_basis_file "
                     "nclusters ";
        return 0;
    }
    TiledArray::SparseShape<float>::threshold(threshold);

    auto clustered_mol = molecule::attach_hydrogens_and_kmeans(
          molecule::read_xyz(mol_file).clusterables(), obs_nclusters);

    auto df_clustered_mol = molecule::attach_hydrogens_and_kmeans(
          molecule::read_xyz(mol_file).clusterables(), dfbs_nclusters);

    auto repulsion_energy = clustered_mol.nuclear_repulsion();
    std::cout << "Nuclear Repulsion Energy: " << repulsion_energy << std::endl;
    auto occ = clustered_mol.occupation(0);

    basis::BasisSet bs(basis_name);
    basis::Basis basis(bs.get_cluster_shells(clustered_mol));

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

    decltype(H) V_inv;
    {
        auto Vmetric = ints::sparse_integrals(world, eri_e, dfbs_array);
        auto V_eig = tcc::array_ops::array_to_eigen(Vmetric);
        { // Silly Drew Test not meaningful
            Eig::JacobiSVD<MatrixD> svdS(V_eig, Eig::ComputeThinU | Eig::ComputeThinV);
            std::cout << std::setprecision(4);
            std::cout << "S vals of V = " << svdS.singularValues().transpose() << std::endl;
            std::cout << "S vects U of  V = \n" << svdS.matrixU() << std::endl;
            std::cout << "S PCA T  = \n" << (svdS.matrixU() * svdS.singularValues().asDiagonal()).transpose() << std::endl;

            MatrixD Qs = V_eig.transpose() * V_eig;
            Qs = Qs.inverse();
            MatrixD Rs = V_eig * V_eig.transpose();
            Rs = Rs.inverse();

            Eig::SelfAdjointEigenSolver<MatrixD> QRes(Qs);
            Qs = QRes.operatorSqrt();
            QRes.compute(Rs);
            Rs = QRes.operatorSqrt();

            MatrixD sph_S = Qs.transpose() * V_eig * Rs;
            svdS.compute(sph_S);
            std::cout << "\nS vals of Q * V * R = " << svdS.singularValues().transpose() << std::endl;
            std::cout << "S vects U of Trans V = \n" << Qs.inverse() * svdS.matrixU() << std::endl;
            std::cout << "S PCA T  = \n" << (Qs.inverse() * svdS.matrixU() * svdS.singularValues().asDiagonal()).transpose() << std::endl;
            std::cout << std::setprecision(15);
        }
        MatrixD Vinv = V_eig.inverse();

        auto tr_V = Vmetric.trange().data()[0];
        V_inv = tcc::array_ops::eigen_to_array<TA::TensorD>(world, Vinv, tr_V,
                                                            tr_V);
    }

    auto three_c_array = tcc::utility::make_array(df_basis, basis, basis);
    { // Schwarz Screened
        std::cout << "Direct Schwarz DF" << std::endl;
        auto sbuilder = ints::init_schwarz_screen(1e-8);
        auto shr_screen = std::make_shared<ints::SchwarzScreen>(
              sbuilder(world, eri_e, df_basis, basis));

        auto eri3 = ints::direct_sparse_integrals(world, eri_e, three_c_array,
                                                  shr_screen);

        auto F_soad
              = scf::fock_from_soad(world, clustered_mol, basis, eri_e, H);
        ThreeCenterScf scf(H, S, F_soad, V_inv, occ / 2, repulsion_energy);
        scf.solve(20, 1e-8, eri3);
        auto hf_energy = scf.energy();
        std::cout << "Final HF energy = " << hf_energy << std::endl;

        auto Fao = scf.fock();

        auto F_eig = tcc::array_ops::array_to_eigen(Fao);
        auto S_eig = tcc::array_ops::array_to_eigen(S);

        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);
        auto nf = basis.nfunctions();
        auto nocc = occ / 2;
        auto nvir = nf - nocc;

        decltype(S_eig) Ceigi = es.eigenvectors().leftCols(nocc);
        decltype(S_eig) Ceigv = es.eigenvectors().rightCols(nvir);


        auto occ_nclusters = obs_nclusters;
        auto tr_occ = tcc::scf::tr_occupied(occ_nclusters, nocc);
        auto tr_vir = tcc::scf::tr_occupied(occ_nclusters, F_eig.cols() - nocc);


        auto Ci = tcc::array_ops::eigen_to_array<TA::TensorD>(
              world, Ceigi, S.trange().data()[0], tr_occ);

        auto multi_pool
              = ints::make_1body_shr_pool("emultipole2", basis, clustered_mol);
        auto r_xyz = ints::sparse_xyz_integrals(world, multi_pool, bs_array);


        std::cout << "\n\nSize for Eri3 if stored" << std::endl;
        print_storage_for_array(eri3.array());
        std::cout << std::endl;


        auto Cil = scf::BoysLocalization{}(Ci, r_xyz);

        tcc::scf::clustered_coeffs(r_xyz, Cil, obs_nclusters);
        std::cout << "Sparsity in Ci local = \n" << Cil.get_shape().sparsity()
                  << std::endl;


        DArray<3, TA::TensorD, TA::SparsePolicy> W;
        W("Y, mu, i") = eri3("Y, mu, nu") * Cil("nu, i");
        world.gop.fence();

        std::cout << "Size for W(Y,mu, i) CMO " << std::endl;
        print_storage_for_array(W);
        std::cout << std::endl;


        auto Cv = tcc::array_ops::eigen_to_array<TA::TensorD>(
              world, Ceigv, S.trange().data()[0], tr_vir);
        tcc::scf::clustered_coeffs(r_xyz, Cv, obs_nclusters);

        DArray<3, TA::TensorD, TA::SparsePolicy> Wiv;
        Wiv("Y, i, a") = W("Y, mu, i") * Cv("mu, a");
        world.gop.fence();

        std::cout << "Size for Wiv(Y, i, a) CMO " << std::endl;
        print_storage_for_array(Wiv);
        std::cout << std::endl;

        DArray<2, TA::TensorD, TA::SparsePolicy> Q;
        Q("mu, v") = S("mu, nu") * Cv("nu,v");
        tcc::scf::clustered_coeffs(r_xyz, Q, obs_nclusters);

        Wiv("Y, i, a") = W("Y, mu, i") * Q("mu, a");

        std::cout << "Size for Wiv(Y, i, a) PAO " << std::endl;
        print_storage_for_array(Wiv);
        std::cout << std::endl;

        auto eval_guess = false;
        auto Cvl = scf::BoysLocalization{}(Cv, r_xyz, eval_guess);
        tcc::scf::clustered_coeffs(r_xyz, Cvl, obs_nclusters);

        Wiv("Y, i, a") = W("Y, mu, i") * Cvl("mu, a");
        world.gop.fence();

        std::cout << "Size for Wiv(Y, i, a) Localized Vir " << std::endl;
        print_storage_for_array(Wiv);
        std::cout << std::endl;

        eval_guess = true;
        Cvl = scf::BoysLocalization{}(Cv, r_xyz, eval_guess);
        tcc::scf::clustered_coeffs(r_xyz, Cvl, obs_nclusters);
        Wiv("Y, i, a") = W("Y, mu, i") * Cvl("mu, a");
        world.gop.fence();

        std::cout << "Size for Wiv(Y, i, a) Localized Vir Eval guess "
                  << std::endl;
        print_storage_for_array(Wiv);
        std::cout << std::endl;
    }

#if 0
        std::cout << "Testing G\n";
        if (basis.nfunctions() <= 200) {
            DArray<4, TA::TensorD, TA::SparsePolicy> G;
            G("i, a, j, b") = Wiv("X,i,a") * (V_inv("X,Y") * Wiv("Y,j,b"));
            world.gop.fence();

            std::cout << "Size for G" << std::endl;
            print_storage_for_array_4(G);
            std::cout << std::endl;


            // Makeing 1/e tensor
            auto vec_ptr = std::make_shared<Eig::VectorXd>(
                  std::move(es.eigenvalues()));
            MatrixD meiajb(nocc * nvir, nocc * nvir);
            auto &eval_vec = *vec_ptr;
            auto ia = 0;
            for (auto i = 0; i < nocc; ++i) {
                const auto ei = eval_vec[i];

                for (auto a = 0; a < nvir; ++a, ++ia) {
                    const auto eia = ei - eval_vec[a + nocc];

                    auto jb = 0;
                    for (auto j = 0; j < nocc; ++j) {
                        auto eiaj = eia + eval_vec[j];

                        for (auto b = 0; b < nvir; ++b, ++jb) {
                            auto eiajb = eiaj - eval_vec[b + nocc];
                            meiajb(ia, jb) = 1.0 / eiajb;
                        }
                    }
                }
            }

            struct Mp2Mat {
                using result_type = double;
                using argument_type = TA::TensorD;

                std::shared_ptr<MatrixD> mat_;
                unsigned int n_occ_;
                unsigned int n_vir_;

                Mp2Mat(std::shared_ptr<MatrixD> mat, int n_occ, int n_vir)
                        : mat_(std::move(mat)), n_occ_(n_occ), n_vir_(n_vir) {}
                Mp2Mat(Mp2Mat const &) = default;

                result_type operator()() const { return 0.0; }
                result_type operator()(result_type const &t) const { return t; }
                void
                operator()(result_type &me, result_type const &other) const {
                    me += other;
                }

                void
                operator()(result_type &me, argument_type const &tile) const {
                    auto const &range = tile.range();
                    auto const &mat = *mat_;
                    auto const st = range.lobound();
                    auto const fn = range.upbound();
                    auto tidx = 0;
                    for (auto i = st[0]; i < fn[0]; ++i) {
                        for (auto a = st[1]; a < fn[1]; ++a) {
                            auto idx_ia = a + i * n_vir_;

                            for (auto j = st[2]; j < fn[2]; ++j) {
                                for (auto b = st[3]; b < fn[3]; ++b, ++tidx) {

                                    auto idx_jb = b + j * n_vir_;
                                    me += mat(idx_ia, idx_jb)
                                          * tile.data()[tidx];
                                }
                            }
                        }
                    }
                }
            };

            auto mat_ptr = std::make_shared<MatrixD>(std::move(meiajb));
            auto energy_mp2 = (G("i,a,j,b") * (2 * G("i,a,j,b") - G("i,b,j,a")))
                                    .reduce(Mp2Mat(mat_ptr, nocc, nvir))
                                    .get();

            std::cout << "MP2 energy = " << energy_mp2
                      << " total energy = " << hf_energy + energy_mp2
                      << std::endl;

            std::cout << "Rank threshold = " << rank_thresh << std::endl;
            (*mat_ptr).resize(nocc, nvir * nocc * nvir);
            Eig::JacobiSVD<MatrixD> svd(*mat_ptr,
                                        Eig::ComputeThinU | Eig::ComputeThinV);
            svd.setThreshold(rank_thresh);
            auto rank1 = svd.rank();
            std::cout << "R1 rank = " << rank1 << std::endl;

            MatrixD G1 = svd.matrixU().leftCols(rank1);
            MatrixD R1 = MatrixD(svd.singularValues().asDiagonal())
                               .block(0, 0, rank1, rank1)
                         * svd.matrixV().leftCols(rank1).transpose();

            R1.resize(rank1 * nvir, nocc * nvir);
            svd.compute(R1);
            auto rank2 = svd.rank();
            std::cout << "R2 rank = " << rank2 << std::endl;


            MatrixD G2 = svd.matrixU().leftCols(rank2);
            MatrixD R2 = MatrixD(svd.singularValues().asDiagonal())
                               .block(0, 0, rank2, rank2)
                         * svd.matrixV().leftCols(rank2).transpose();


            R2.resize(rank2 * nocc, nvir);
            svd.compute(R2);
            auto rank3 = svd.rank();
            std::cout << "R3 rank = " << rank3 << std::endl;

            MatrixD G3 = svd.matrixU().leftCols(rank3);
            MatrixD G4 = MatrixD(svd.singularValues().asDiagonal())
                               .block(0, 0, rank3, rank3)
                         * svd.matrixV().leftCols(rank3).transpose();


            auto org_elems = mat_ptr->size();
            auto tt_elems = nocc * rank1 + rank1 * nvir * rank2
                            + rank2 * nocc * rank3 + rank3 * nvir;
            std::cout << "Num orgi elements = " << mat_ptr->size() << " "
                      << 1e-9 * org_elems * 8 << " GB." << std::endl;
            std::cout << "TT storage = " << tt_elems << " "
                      << 1e-9 * tt_elems * 8 << " GB." << std::endl;

            std::cout << "\nReconstructing Delta iajb\n";
            MatrixD AppR2 = G3 * G4;
            AppR2.resize(rank2, nocc * nvir);

            MatrixD AppR1 = G2 * AppR2;
            AppR1.resize(rank1, nvir * nocc * nvir);

            MatrixD AppD = G1 * AppR1;
            AppD.resize(nocc * nvir, nocc * nvir);
            mat_ptr->resize(nocc * nvir, nocc * nvir);

            auto diff = (AppD - *mat_ptr).lpNorm<2>();
            std::cout << "||.||_2 diff after reconstruction = " << diff
                      << std::endl;

            std::cout << "Energy with reconstructed values\n";
            *mat_ptr = AppD;
            auto app_energy_mp2
                  = (G("i,a,j,b") * (2 * G("i,a,j,b") - G("i,b,j,a")))
                          .reduce(Mp2Mat(mat_ptr, nocc, nvir))
                          .get();

            std::cout << "MP2 energy = " << app_energy_mp2
                      << " total energy = " << hf_energy + app_energy_mp2
                      << std::endl;
            auto app_mp2_diff = std::abs(app_energy_mp2 - energy_mp2);
            std::cout
                  << "MP2 error = " << app_mp2_diff << " percent of energy = "
                  << 100.0 * (1.0 - app_energy_mp2 / energy_mp2) << std::endl;
        }
    }
#endif
#if 0
        {
            auto Cl = scf::StuffFromFactorAna{}(Ci);
            auto Clv = scf::StuffFromFactorAna{}(Cv);
            tcc::scf::clustered_coeffs(r_xyz, Cl, obs_nclusters);
            tcc::scf::clustered_coeffs(r_xyz, Clv, obs_nclusters);
            DArray<3, TA::TensorD, TA::SparsePolicy> W;
            W("Y, mu, i") = eri3("Y, mu, nu") * S_oh_inv("nu, rho")
                            * Cl("rho, i");
            world.gop.fence();

            std::cout << "Size for W(Y,mu, i) SVD " << std::endl;
            print_storage_for_array(W);
            std::cout << std::endl;

            DArray<3, TA::TensorD, TA::SparsePolicy> Wiv;
            Wiv("Y, i, a") = W("Y, mu, i") * (S_oh_inv("mu, rho")
                             * Clv("rho, a"));
            world.gop.fence();

            std::cout << "Size for Wiv(Y, i, a) CMO " << std::endl;
            print_storage_for_array(Wiv);
            std::cout << std::endl;
        }

        decltype(Ci) D_orth;
        D_orth("mu,nu") = Ci("mu, i") * Ci("nu,i");

        Ci = scf::BoysLocalization{}(Ci, r_xyz);

        D_orth("mu,nu") = D_orth("mu, nu") - Ci("mu, i") * Ci("nu,i");
        auto norm_diff = D_orth("i,j").norm().get();
        std::cout << "Norm of density diff after localization = " << norm_diff
                  << std::endl;

        tcc::scf::clustered_coeffs(r_xyz, Ci, obs_nclusters);

        {
            DArray<3, TA::TensorD, TA::SparsePolicy> W;
            W("Y, mu, i") = eri3("Y, mu, nu") * S_oh_inv("nu, rho")
                            * Ci("rho, i");
            world.gop.fence();

            std::cout << "Size for W(Y,mu, i) if stored" << std::endl;
            print_storage_for_array(W);
            std::cout << std::endl;
        }
    }

    if(false){ // QVl Screened
        std::cout << "\n\nDirect Schwarz QVl" << std::endl;
        auto sbuilder = ints::init_qvl_screen{};
        ints::QVl::well_sep_threshold(1e-12);
        auto shr_screen = std::make_shared<ints::QVl>(
              sbuilder(world, eri_e, df_basis, basis, 1e-8));

        auto eri3 = ints::direct_sparse_integrals(world, eri_e, three_c_array,
                                                  shr_screen);

        ThreeCenterScf scf(H, S, L_inv, occ / 2, repulsion_energy);
        scf.solve(20, 1e-7, eri3);
    }

    if (false) { // Unscreened ints
        std::cout << "\n\nDirect Unscreened DF" << std::endl;
        auto eri3 = ints::direct_sparse_integrals(world, eri_e, three_c_array);
        ThreeCenterScf scf(H, S, L_inv, occ / 2, repulsion_energy);
        scf.solve(20, 1e-7, eri3);
    }
#endif


    return 0;
}
