#include "mpqc/chemistry/qc/lcao/integrals/integrals.h"
#include "mpqc/chemistry/qc/lcao/scf/eigen_solve_density_builder.h"

#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/linalg/cholesky_inverse.h"
#include "mpqc/math/linalg/conditioned_orthogonalizer.h"
#include "mpqc/math/linalg/sqrt_inv.h"
#include "mpqc/math/tensor/clr/minimize_storage.h"

#include "mpqc/chemistry/qc/lcao/expression/trange1_engine.h"
#include "mpqc/chemistry/qc/lcao/scf/util.h"

#include "mpqc/chemistry/qc/lcao/scf/clusterd_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/orbital_localization.h"

namespace mpqc {
namespace lcao {
namespace scf {

template<typename Tile, typename Policy>
ESolveDensityBuilder<Tile, Policy>::ESolveDensityBuilder(
    typename ESolveDensityBuilder<Tile, Policy>::array_type const &S,
    std::vector<typename ESolveDensityBuilder<Tile, Policy>::array_type> r_xyz,
    int64_t nocc, int64_t ncore, int64_t nclusters, double TcutC,
    std::string const &metric_decomp_type, double s_tolerance,
    std::shared_ptr<OrbitalLocalizer < Tile, Policy>>
localizer,
bool localize_core,
bool clustered_coeffs
)
:
S_ (S),
r_xyz_ints_(r_xyz),
TcutC_(TcutC),
localizer_(localizer),
localize_core_(localize_core),
n_coeff_clusters_(nclusters),
clustered_coeffs_(clustered_coeffs),
metric_decomp_type_(metric_decomp_type),
nocc_(nocc),
ncore_(ncore) {
  auto inv0 = mpqc::fenced_now(S_.world());
  if (metric_decomp_type_ == "cholesky_inverse") {
    M_inv_ = math::cholesky_inverse(S_);
  } else if (metric_decomp_type_ == "inverse_sqrt") {
    M_inv_ = math::inverse_sqrt(S_);
  } else if (metric_decomp_type_ == "conditioned") {
    M_inv_ = math::conditioned_orthogonalizer(S_, s_tolerance);
  } else {
    throw ProgrammingError("invalid metric_decomp_type in EsolveDensityBuilder",
                           __FILE__, __LINE__);
  }

  MPQC_ASSERT(!(clustered_coeffs_ && localize_core_ == false));

  auto inv1 = mpqc::fenced_now(S_.world());
  inverse_time_ = mpqc::duration_in_s(inv0, inv1);
}

template<typename Tile, typename Policy>
std::pair<typename ESolveDensityBuilder<Tile, Policy>::array_type,
          typename ESolveDensityBuilder<Tile, Policy>::array_type>
ESolveDensityBuilder<Tile, Policy>::operator()(
    typename ESolveDensityBuilder<Tile, Policy>::array_type const &F) {
  using scalar_type = typename Tile::scalar_type;
  typename ESolveDensityBuilder<Tile, Policy>::array_type Fp, C_occ, C_occ_ao, D;
  auto &world = F.world();

  auto e0 = mpqc::fenced_now(world);
  Fp("i,j") = M_inv_("i,k") * F("k,l") * M_inv_("j,l");

  auto Fp_eig = math::array_to_eigen(Fp);

  // solve on rank-0 only, propagate to the rest
  EigenVector<scalar_type> eps_eig;
  decltype(Fp_eig) C_eig;
  if (world.rank() == 0) {
    Eigen::SelfAdjointEigenSolver<decltype(Fp_eig)> es(Fp_eig);
    eps_eig = es.eigenvalues();
    C_eig = es.eigenvectors();

    // warn about "exact" degeneracies among valence occupied orbitals
    bool exact_degeneracies = false;
    for(auto i=ncore_+1; i<nocc_ && !exact_degeneracies; ++i) {
      if (std::abs(eps_eig(i) - eps_eig(i - 1)) < std::numeric_limits<scalar_type>::epsilon() * 1e3) {
        exact_degeneracies = true;
        ExEnv::out0() << indent
                      << "WARNING: nearly exactly degenerate occupied orbital energies, phases are NOT fixed, so take care in interpreting orbital coefficients"
                      << std::endl;
      }
    }
  }
  world.gop.broadcast_serializable(eps_eig, 0);
  world.gop.broadcast_serializable(C_eig, 0);
  decltype(Fp_eig) C_occ_eig = C_eig.leftCols(nocc_);

  auto tr_ao = Fp.trange().data()[0];
  auto tr_occ = utility::compute_trange1(nocc_, n_coeff_clusters_);
  C_occ = math::eigen_to_array<Tile, Policy>(Fp.world(), C_occ_eig, tr_ao,
                                                  tr_occ);

  // Get back to AO land
  C_occ_ao("i,j") = M_inv_("k,i") * C_occ("k,j");
  C_occ_ao.truncate();
  if (!localizer_) {
    eps_ = eps_eig;
    auto nobs = eps_eig.rows();
    auto tr_obs = utility::compute_trange1(nobs, n_coeff_clusters_);
    auto C = math::eigen_to_array<Tile, Policy>(Fp.world(), C_eig, tr_ao,
                                                     tr_obs);
    C_("i,j") = M_inv_("k,i") * C("k,j");
  }
  auto e1 = mpqc::fenced_now(world);

  // Compute D to full accuracy
  D("i,j") = C_occ_ao("i,k") * C_occ_ao("j,k");
  D.truncate();
  density_storages_.push_back(mpqc::detail::array_storage(D));

  if (localizer_) {
    auto l0 = mpqc::fenced_now(world);
    auto U = localizer_->compute(C_occ_ao, (!localize_core_ ? ncore_ : 0));
    C_occ_ao("mu,i") = C_occ_ao("mu,k") * U("k,i");
    auto l1 = mpqc::fenced_now(world);

    auto obs_ntiles = C_occ_ao.trange().tiles_range().extent()[0];
    if (clustered_coeffs_) {
      scf::clustered_coeffs(r_xyz_ints_, C_occ_ao, obs_ntiles);
    }
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
  coeff_storages_.push_back(mpqc::detail::array_storage(C_occ_ao));

  return std::make_pair(D, C_occ_ao);
}

}  // namespace scf
}  // namespace lcao
}  // namespace mpqc
