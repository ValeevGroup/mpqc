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

#include <memory>

using namespace mpqc;
namespace ints = mpqc::integrals;

class FourCenterSCF {
  private:
    using array_type = DArray<2, TA::TensorD, SpPolicy>;
    array_type H_;
    array_type S_;

    array_type F_;
    array_type D_;
    TiledArray::DIIS<array_type> diis_;

    std::vector<double> k_times_;
    std::vector<double> j_times_;
    std::vector<double> scf_times_;

    int64_t occ_;
    double repulsion_;

    template <typename Integral>
    array_type compute_k(Integral const &eri4) {

        array_type K;
        auto &world = eri4.get_world();
        world.gop.fence();
        auto k0 = tcc_time::now();
        K("i,j") = eri4("i,k,j,l") * D_("k,l");
        world.gop.fence();
        auto k1 = tcc_time::now();
        k_times_.push_back(tcc_time::duration_in_s(k0, k1));

        return K;
    }

    template <typename Integral>
    array_type compute_j(Integral const &eri4) {

        array_type J;
        auto &world = eri4.get_world();
        world.gop.fence();
        auto j0 = tcc_time::now();
        J("i,j") = eri4("i,j,k,l") * D_("k,l");
        world.gop.fence();
        auto j1 = tcc_time::now();
        j_times_.push_back(tcc_time::duration_in_s(j0, j1));

        return J;
    }


    void compute_density(int64_t occ) {
        auto F_eig = tcc::array_ops::array_to_eigen(F_);
        auto S_eig = tcc::array_ops::array_to_eigen(S_);

        Eig::GeneralizedSelfAdjointEigenSolver<decltype(S_eig)> es(F_eig,
                                                                   S_eig);
        decltype(S_eig) C = es.eigenvectors().leftCols(occ);
        MatrixD D_eig = C * C.transpose();

        auto tr_ao = S_.trange().data()[0];

        D_ = tcc::array_ops::eigen_to_array<TA::TensorD>(F_.get_world(), D_eig,
                                                         tr_ao, tr_ao);
    }

    template <typename Integral>
    void form_fock(Integral const &eri4) {
        F_("i,j") = H_("i,j") + 2 * compute_j(eri4)("i,j")
                    - compute_k(eri4)("i,j");
    }


  public:
    FourCenterSCF(array_type const &H, array_type const &S, int64_t occ,
                  double rep)
            : H_(H), S_(S), occ_(occ), repulsion_(rep) {
        F_ = H_;
        compute_density(occ_);
    }

    template <typename Integral>
    void solve(int64_t max_iters, double thresh, Integral const &eri4) {
        auto iter = 0;
        auto error = std::numeric_limits<double>::max();
        auto old_energy = 0.0;

        while (iter < max_iters && thresh < error) {
            auto s0 = tcc_time::now();
            F_.get_world().gop.fence();
            form_fock(eri4);

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

TA::TensorD ta_pass_through(TA::TensorD &&ten) { return std::move(ten); }

int main(int argc, char *argv[]) {
    auto &world = madness::initialize(argc, argv);
    std::string mol_file = "";
    std::string basis_name = "";
    int nclusters = 0;
    std::cout << std::setprecision(15);
    double threshold = 1e-12;
    double well_sep_threshold = 0.1;
    integrals::QQR::well_sep_threshold(well_sep_threshold);
    if (argc == 4) {
        mol_file = argv[1];
        basis_name = argv[2];
        nclusters = std::stoi(argv[3]);
    } else {
        std::cout << "input is $./program mol_file basis_file nclusters ";
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

    libint2::init();

    const auto bs_array = tcc::utility::make_array(basis, basis);

    auto overlap_e = ints::make_1body_shr_pool("overlap", basis, clustered_mol);
    auto S = ints::sparse_integrals(world, overlap_e, bs_array);

    auto kinetic_e = ints::make_1body_shr_pool("kinetic", basis, clustered_mol);
    auto T = ints::sparse_integrals(world, kinetic_e, bs_array);

    auto nuclear_e = ints::make_1body_shr_pool("nuclear", basis, clustered_mol);
    auto V = ints::sparse_integrals(world, nuclear_e, bs_array);

    decltype(T) H;
    H("i,j") = T("i,j") + V("i,j");

    auto bs4_array = tcc::utility::make_array(basis, basis, basis, basis);
    auto eri_e = ints::make_2body_shr_pool(basis);
    { // Unscreened SCF
        std::cout << "Computing HF with No Screening" << std::endl;
        world.gop.fence();
        auto eri40 = tcc_time::now();
        auto eri4 = ints::direct_sparse_integrals(world, eri_e, bs4_array);
        world.gop.fence();
        auto eri41 = tcc_time::now();
        std::cout << "Took " << tcc_time::duration_in_s(eri40, eri41)
                  << " s to initialize integrals." << std::endl;

        FourCenterSCF scf(H, S, occ / 2, repulsion_energy);
        scf.solve(50, 1e-12, eri4);
    }

#if 0
    { // Do schwarz 
        std::cout << "Computing HF with Schwarz Screening" << std::endl;
        world.gop.fence();
        auto eri40 = tcc_time::now();

        auto eri4 = mpqc_ints::
              direct_sparse_integrals(
                    world, eri_e, bs4_array, );
        world.gop.fence();
        auto eri41 = tcc_time::now();
        std::cout << "Took " << tcc_time::duration_in_s(eri40, eri41)
                  << " s to initialize integrals." << std::endl;

        FourCenterSCF scf(H, S, occ / 2, repulsion_energy);
        scf.solve(20, 1e-8, eri4);
        std::cout << "\n\n";
    }
    { // Do then QQR
        std::cout << "Computing HF with QQR Based Integrals" << std::endl;
        world.gop.fence();
        auto eri40 = tcc_time::now();
        auto eri4
              = mpqc_ints::direct_sparse_integrals<mpqc_ints::init_qqr_screen>(
                    world, eri_pool, bs4_array, ta_pass_through);
        world.gop.fence();
        auto eri41 = tcc_time::now();
        std::cout << "Took " << tcc_time::duration_in_s(eri40, eri41)
                  << " s to initialize integrals." << std::endl;

        FourCenterSCF scf(H, S, occ / 2, repulsion_energy);
        scf.solve(20, 1e-8, eri4);
        std::cout << "\n\n";
    }
    { // Finally unscreened
        std::cout << "Computing HF with Unscreened Integrals" << std::endl;
        world.gop.fence();
        auto eri40 = tcc_time::now();
        auto eri4
              = mpqc_ints::direct_sparse_integrals(world, eri_pool, bs4_array,
                                                   ta_pass_through);
        world.gop.fence();
        auto eri41 = tcc_time::now();
        std::cout << "Took " << tcc_time::duration_in_s(eri40, eri41)
                  << " s to initialize integrals." << std::endl;

        FourCenterSCF scf(H, S, occ / 2, repulsion_energy);
        scf.solve(50, 1e-8, eri4);
    }
#endif

    libint2::cleanup();
    madness::finalize();
    return 0;
}
