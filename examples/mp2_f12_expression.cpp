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

#include "../scf/scf.h"

#include "../basis/atom_basisset.h"
#include "../basis/basis_set.h"
#include "../basis/basis.h"

#include "../integrals/integrals.h"
#include "../integrals/atomic_integral.h"

#include "../utility/time.h"
#include "../utility/wcout_utf8.h"
#include "../utility/array_info.h"
#include "../ta_routines/array_to_eigen.h"
#include "../f12/utility.h"

#include <memory>
#include <locale>

using namespace mpqc;
namespace ints = mpqc::integrals;


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

    basis::BasisSet dfbs("cc-pVTZ-RI");
    basis::Basis df_basis(dfbs.get_cluster_shells(clustered_mol));

    basis::BasisSet abs("aug-cc-pVDZ-CABS");
    basis::Basis abs_basis(abs.get_cluster_shells(clustered_mol));

    f12::GTGParams gtg_params(1.2,6);

    libint2::init();

    integrals::AtomicIntegral<TA::TensorD, TA::SparsePolicy> ao_int
            (world,
             ta_pass_through,
            std::make_shared<molecule::Molecule>(clustered_mol),
            std::make_shared<basis::Basis>(basis),
             std::make_shared<basis::Basis>(df_basis),
            std::make_shared<basis::Basis> (abs_basis),
             gtg_params.compute()
            );

    // Overlap ints
    auto S = ao_int.compute(L"(κ|λ)");

    auto T = ao_int.compute(L"(κ| T|λ)");

    auto V = ao_int.compute(L"(κ| V|λ)");

    decltype(T) H;
    H("i,j") = T("i,j") + V("i,j");

    // Unscreened four center stored RHF.
    auto eri4 = ao_int.compute(L"(κ λ| G|κ1 λ1)");
    world.gop.fence();

    RHF scf(H, S, occ / 2, repulsion_energy);
    scf.solve(50, 1e-8, eri4);


    // CABS fock build

    // integral

    auto S_cabs = ao_int.compute(L"(α|β)");
    auto S_ribs = ao_int.compute(L"(ρ|σ)");
    auto S_obs_ribs = ao_int.compute(L"(μ|σ)");
    auto S_obs = scf.get_overlap();


    // construct cabs
    auto S_obs_eigen = array_ops::array_to_eigen(S_obs);
    MatrixD X_obs_eigen = Eigen::LLT<MatrixD>(S_obs_eigen).matrixL();
    MatrixD X_obs_eigen_inv = X_obs_eigen.inverse();

    auto S_ribs_eigen = array_ops::array_to_eigen(S_ribs);
    MatrixD X_ribs_eigen = Eigen::LLT<MatrixD>(S_ribs_eigen).matrixL();
    MatrixD X_ribs_eigen_inv = X_ribs_eigen.inverse();


    MatrixD S_obs_ribs_eigen = array_ops::array_to_eigen(S_obs_ribs);

    auto S_obs_ribs_ortho_eigen = X_obs_eigen_inv.transpose()*S_obs_ribs_eigen * X_ribs_eigen_inv;

    Eigen::JacobiSVD<MatrixD> svd(S_obs_ribs_ortho_eigen, Eigen::ComputeFullV);
    MatrixD V_eigen = svd.matrixV();
    // but there could be more!! (homework!)
    size_t nbf_ribs = S_obs_ribs_ortho_eigen.cols();
    //size_t nbf_obs = S_obs_ribs_ortho.rows();
    auto nbf_cabs = nbf_ribs - svd.nonzeroSingularValues();
    MatrixD Vnull(nbf_ribs, nbf_cabs);

    //Populate Vnull with vectors of V that are orthogonal to AO space
    Vnull = V_eigen.block(0, svd.nonzeroSingularValues(), nbf_ribs, nbf_cabs);

    //Un-orthogonalize coefficients
    MatrixD C_cabs_eigen = X_ribs_eigen_inv * Vnull;

    auto tr_cabs = S_cabs.trange().data()[0];
    auto tr_ribs = S_ribs.trange().data()[0];

    auto C_cabs = array_ops::eigen_to_array<TA::TensorD>(world,C_cabs_eigen,tr_ribs,tr_cabs);
//    std::cout << "C_cabs" << std::endl;
//    std::cout << C_cabs << std::endl;

    // compute fock
    auto T_ribs = ao_int.compute(L"(ρ|T|σ)");
    auto V_ribs = ao_int.compute(L"(ρ|V|σ)");

    auto J_ribs_obs = ao_int.compute(L"(ρ σ|G| μ ν)");
    auto K_ribs_obs = ao_int.compute(L"(ρ μ|G| σ ν)");



    // compute r12 integral

//    auto r12 = ao_int.compute4(L"(α β |R| γ δ)");
    auto gr12 = ao_int.compute(L"( γ |GR|  δ)");
//    auto r12_2 = ao_int.compute4(L"(α β |R2| γ δ)");






    madness::finalize();
    libint2::cleanup();
    return 0;
}
