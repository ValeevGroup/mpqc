#include "mpqc/chemistry/qc/lcao/integrals/integrals.h"
#include "mpqc/chemistry/qc/lcao/scf/eigen_solve_density_builder.h"

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/linalg/cholesky_inverse.h"
#include "mpqc/math/linalg/conditioned_orthogonalizer.h"
#include "mpqc/math/linalg/sqrt_inv.h"
#include "mpqc/math/tensor/clr/minimize_storage.h"

#include "mpqc/chemistry/qc/lcao/scf/util.h"

#include "mpqc/chemistry/qc/lcao/scf/clusterd_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/diagonalize_for_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/orbital_localization.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy>
ESolveDensityBuilder<Tile, Policy>::ESolveDensityBuilder(
    typename ESolveDensityBuilder<Tile, Policy>::array_type const &S,
    std::vector<typename ESolveDensityBuilder<Tile, Policy>::array_type> r_xyz,
    int64_t nocc, int64_t nclusters, double TcutC,
    std::string const &metric_decomp_type, bool localize)
    : S_(S),
      r_xyz_ints_(r_xyz),
      TcutC_(TcutC),
      localize_(localize),
      n_coeff_clusters_(nclusters),
      metric_decomp_type_(metric_decomp_type),
      nocc_(nocc) {
  auto inv0 = mpqc::fenced_now(S_.world());
  if (metric_decomp_type_ == "cholesky_inverse") {
    M_inv_ = array_ops::cholesky_inverse(S_);
  } else if (metric_decomp_type_ == "inverse_sqrt") {
    M_inv_ = array_ops::inverse_sqrt(S_);
  } else if (metric_decomp_type_ == "conditioned") {
    M_inv_ = array_ops::conditioned_orthogonalizer(
        S_, 1.0 / std::numeric_limits<double>::epsilon());
  } else {
    throw "Error did not recognize overlap decomposition in "
              "EsolveDensityBuilder";
  }
  auto inv1 = mpqc::fenced_now(S_.world());
  inverse_time_ = mpqc::duration_in_s(inv0, inv1);
}

template <typename Tile, typename Policy>
std::pair<typename ESolveDensityBuilder<Tile, Policy>::array_type,
          typename ESolveDensityBuilder<Tile, Policy>::array_type>
ESolveDensityBuilder<Tile, Policy>::operator()(
    typename ESolveDensityBuilder<Tile, Policy>::array_type const &F) {
  typename ESolveDensityBuilder<Tile, Policy>::array_type Fp, C_occ, C_occ_ao,
      D;
  auto &world = F.world();

  auto e0 = mpqc::fenced_now(world);
  Fp("i,j") = M_inv_("i,k") * F("k,l") * M_inv_("j,l");

  auto Fp_eig = array_ops::array_to_eigen(Fp);

  // TODO must solve on rank-0 only, propagate to the rest
  Eigen::SelfAdjointEigenSolver<decltype(Fp_eig)> es(Fp_eig);
  auto eps = es.eigenvalues();
  decltype(Fp_eig) C_eig = es.eigenvectors();
  decltype(Fp_eig) C_occ_eig = C_eig.leftCols(nocc_);
  auto tr_ao = Fp.trange().data()[0];

  auto tr_occ = scf::tr_occupied(n_coeff_clusters_, nocc_);
  C_occ = array_ops::eigen_to_array<Tile, Policy>(Fp.world(), C_occ_eig, tr_ao,
                                                  tr_occ);

  // Get back to AO land
  C_occ_ao("i,j") = M_inv_("k,i") * C_occ("k,j");
  C_occ_ao.truncate();
  if (!localize_) {
    eps_ = eps;
    auto nobs = eps.rows();
    auto tr_obs = scf::tr_occupied(n_coeff_clusters_, nobs);
    auto C = array_ops::eigen_to_array<Tile, Policy>(Fp.world(), C_eig, tr_ao,
                                                     tr_obs);
    C_("i,j") = M_inv_("k,i") * C("k,j");
  }
  auto e1 = mpqc::fenced_now(world);

  // Compute D to full accuracy
  D("i,j") = C_occ_ao("i,k") * C_occ_ao("j,k");
  D.truncate();
  density_storages_.push_back(detail::array_storage(D));

  if (localize_) {
    auto l0 = mpqc::fenced_now(world);
    auto U = mpqc::scf::FosterBoysLocalization{}(C_occ_ao, r_xyz_ints_);
    C_occ_ao("mu,i") = C_occ_ao("mu,k") * U("k,i");
    auto l1 = mpqc::fenced_now(world);

    auto obs_ntiles = C_occ_ao.trange().tiles_range().extent()[0];
    scf::clustered_coeffs(r_xyz_ints_, C_occ_ao, obs_ntiles);
    auto c1 = mpqc::fenced_now(world);
    localization_times_.push_back(mpqc::duration_in_s(l0, l1));
    clustering_times_.push_back(mpqc::duration_in_s(l1, c1));
  }

#if TA_DEFAULT_POLICY == 1
  if (TcutC_ != 0) {
    minimize_storage(C_occ_ao, TcutC_);
  }
#endif

  esolve_times_.push_back(mpqc::duration_in_s(e0, e1));
  coeff_storages_.push_back(detail::array_storage(C_occ_ao));

  return std::make_pair(D, C_occ_ao);
}

}  // namespace scf
}  // namespace mpqc
