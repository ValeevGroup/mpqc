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

#include <memory>

using namespace mpqc;
namespace ints = mpqc::integrals;

class ThreeCenterScf {
  private:
    using array_type = DArray<2, TA::TensorD, SpPolicy>;
    array_type H_;
    array_type S_;

    array_type F_;
    array_type D_;
    array_type C_;
    array_type L_invV_;
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

        C_ = tcc::array_ops::eigen_to_array<TA::TensorD>(H_.get_world(), C, tr_ao, tr_occ);

        D_("i,j") = C_("i,k") * C_("j,k");
    }

    template <typename Integral>
    void form_fock(Integral const &eri3) {
        auto &world = F_.get_world();

        world.gop.fence();
        auto w0 = tcc::utility::time::now();
        TA::Array<double, 3, TA::TensorD, TA::SparsePolicy> W;
        W("X, mu, i") = L_invV_("X,Y") * (eri3("Y, mu, nu") * C_("nu, i"));
        world.gop.fence();
        auto w1 = tcc::utility::time::now();
        w_times_.push_back(tcc::utility::time::duration_in_s(w0,w1));

        
        array_type J;
        J("mu, nu") = eri3("X, mu, nu")
                      * (L_invV_("Y, X") * (W("Y, rho, i") * C_("rho, i")));
        world.gop.fence();
        auto j1 = tcc::utility::time::now();
        j_times_.push_back(tcc::utility::time::duration_in_s(w1,j1));

    
        // Permute W
        W("X,i,nu") = W("X,nu,i");
        array_type K;
        K("mu, nu") = W("X, i, mu") * W("X, i, nu");
        world.gop.fence();
        auto k1 = tcc::utility::time::now();
        k_times_.push_back(tcc::utility::time::duration_in_s(j1,k1));

        F_("i,j") = H_("i,j") + 2 * J("i,j") - K("i,j");
    }


  public:
    ThreeCenterScf(array_type const &H, array_type const &S,
                  array_type const &L_invV, int64_t occ, double rep)
            : H_(H), S_(S), L_invV_(L_invV), occ_(occ), repulsion_(rep) {
        F_ = H_;
        compute_density(occ_);
    }

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
    int nclusters = 0;
    std::cout << std::setprecision(15);
    double threshold = 1e-12;
    double well_sep_threshold = 0.1;
    integrals::QQR::well_sep_threshold(well_sep_threshold);
    if (argc == 5) {
        mol_file = argv[1];
        basis_name = argv[2];
        df_basis_name = argv[3];
        nclusters = std::stoi(argv[4]);
    } else {
        std::cout << "input is $./program mol_file basis_file df_basis_file "
                     "nclusters ";
        return 0;
    }
    TiledArray::SparseShape<float>::threshold(threshold);

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
    { // Unscreened ints
        auto eri3 = ints::sparse_integrals(world, eri_e, three_c_array);
        ThreeCenterScf scf(H, S, L_inv, occ / 2, repulsion_energy);
        scf.solve(20, 1e-7, eri3);
    }

    return 0;
}
