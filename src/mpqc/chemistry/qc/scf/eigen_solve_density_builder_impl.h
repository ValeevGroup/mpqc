#include "mpqc/chemistry/qc/scf/eigen_solve_density_builder.h"
#include "mpqc/chemistry/qc/integrals/integrals.h"

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/linalg/cholesky_inverse.h"
#include "mpqc/math/linalg/conditioned_orthogonalizer.h"
#include "mpqc/math/linalg/sqrt_inv.h"
#include "mpqc/math/tensor/clr/minimize_storage.h"

#include "mpqc/chemistry/qc/scf/util.h"

#include "mpqc/chemistry/qc/scf/clusterd_coeffs.h"
#include "mpqc/chemistry/qc/scf/diagonalize_for_coeffs.h"
#include "mpqc/chemistry/qc/scf/orbital_localization.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy>
ESolveDensityBuilder<Tile,Policy>::ESolveDensityBuilder(
    typename ESolveDensityBuilder<Tile,Policy>::array_type const &S, std::vector<typename ESolveDensityBuilder<Tile,Policy>::array_type> r_xyz, int64_t occ,
    int64_t nclusters, double TcutC, std::string const &metric_decomp_type,
    bool localize)
    : S_(S),
      r_xyz_ints_(r_xyz),
      occ_(occ),
      n_coeff_clusters_(nclusters),
      TcutC_(TcutC),
      metric_decomp_type_(metric_decomp_type),
      localize_(localize) {
  auto inv0 = mpqc::fenced_now(S_.world());
  if (metric_decomp_type_ == "cholesky_inverse") {
    M_inv_ = array_ops::cholesky_inverse(S_);
  } else if (metric_decomp_type_ == "inverse_sqrt") {
    M_inv_ = array_ops::inverse_sqrt(S_);
  } else if (metric_decomp_type_ == "conditioned") {
     M_inv_ = array_ops::conditioned_orthogonalizer(S_, 1.0 / std::numeric_limits<double>::epsilon());
  } else {
    throw "Error did not recognize overlap decomposition in "
              "EsolveDensityBuilder";
  }
  auto inv1 = mpqc::fenced_now(S_.world());
  inverse_time_ = mpqc::duration_in_s(inv0, inv1);
}

template <typename Tile, typename Policy>
std::pair<typename ESolveDensityBuilder<Tile,Policy>::array_type, typename ESolveDensityBuilder<Tile,Policy>::array_type>
ESolveDensityBuilder<Tile,Policy>::operator()(
    typename ESolveDensityBuilder<Tile,Policy>::array_type const &F) {
  typename ESolveDensityBuilder<Tile,Policy>::array_type Fp, C, Cao, D;
  auto &world = F.world();

  auto e0 = mpqc::fenced_now(world);
  Fp("i,j") = M_inv_("i,k") * F("k,l") * M_inv_("j,l");

  auto Fp_eig = array_ops::array_to_eigen(Fp);
  Eigen::SelfAdjointEigenSolver<decltype(Fp_eig)> es(Fp_eig);

  decltype(Fp_eig) C_eig = es.eigenvectors().leftCols(occ_);
  auto tr_ao = Fp.trange().data()[0];

  auto tr_occ = scf::tr_occupied(n_coeff_clusters_, occ_);
  C = array_ops::eigen_to_array<Tile,Policy>(Fp.world(), C_eig, tr_ao, tr_occ);

  // Get back to AO land
  Cao("i,j") = M_inv_("k,i") * C("k,j");
  Cao.truncate();
  auto e1 = mpqc::fenced_now(world);

  // Compute D to full accuracy
  D("i,j") = Cao("i,k") * Cao("j,k");
  D.truncate();
  density_storages_.push_back(detail::array_storage(D));

  if (localize_) {
    auto l0 = mpqc::fenced_now(world);
    auto U = mpqc::scf::BoysLocalization{}(Cao, r_xyz_ints_);
    Cao("mu,i") = Cao("mu,k") * U("k,i");
    auto l1 = mpqc::fenced_now(world);

    auto obs_ntiles = Cao.trange().tiles_range().extent()[0];
    scf::clustered_coeffs(r_xyz_ints_, Cao, obs_ntiles);
    auto c1 = mpqc::fenced_now(world);
    localization_times_.push_back(mpqc::duration_in_s(l0, l1));
    clustering_times_.push_back(mpqc::duration_in_s(l1, c1));
  }

#if TA_DEFAULT_POLICY == 1
  if (TcutC_ != 0) {
    minimize_storage(Cao, TcutC_);
    auto shape_c = Cao.shape();
    auto &norms = shape_c.data();
    auto norm_mat = TA::eigen_map(norms, norms.range().extent_data()[0],
                                  norms.range().extent_data()[1]);
  }
#endif

  esolve_times_.push_back(mpqc::duration_in_s(e0, e1));
  coeff_storages_.push_back(detail::array_storage(Cao));

  return std::make_pair(D, Cao);
}

}  // namespace scf
}  // namespace mpqc
