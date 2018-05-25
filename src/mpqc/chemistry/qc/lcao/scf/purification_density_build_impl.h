#include "mpqc/chemistry/qc/lcao/integrals/integrals.h"

#include "mpqc/chemistry/qc/lcao/scf/purification_density_build.h"
#include "mpqc/math/external/eigen/eigen.h"
#include "mpqc/math/linalg/cholesky_inverse.h"
#include "mpqc/math/linalg/sqrt_inv.h"
#include "mpqc/math/tensor/clr/minimize_storage.h"

#include "mpqc/chemistry/qc/lcao/scf/clusterd_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/diagonalize_for_coeffs.h"
#include "mpqc/chemistry/qc/lcao/scf/orbital_localization.h"

namespace mpqc {
namespace scf {

template <typename Tile, typename Policy>
PurificationDensityBuilder<Tile, Policy>::PurificationDensityBuilder(
    PurificationDensityBuilder<Tile, Policy>::array_type const &S,
    std::vector<PurificationDensityBuilder<Tile, Policy>::array_type> r_xyz,
    int64_t occ, int64_t ncore, int64_t nclusters, double TcutC, bool localize,
    std::string localization_method, bool clustered_coeffs)
    : S_(S),
      r_xyz_ints_(r_xyz),
      TcutC_(TcutC),
      localize_(localize),
      localization_method_(localization_method),
      n_coeff_clusters_(nclusters),
      clustered_coeffs_(clustered_coeffs),
      occ_(occ),
      ncore_(ncore) {
  M_inv_ = array_ops::inverse_sqrt(S_);
  I_ = array_ops::create_diagonal_matrix(S_, 1.0);
}

template <typename Tile, typename Policy>
typename PurificationDensityBuilder<Tile, Policy>::array_type
PurificationDensityBuilder<Tile, Policy>::purify(
    typename PurificationDensityBuilder<Tile, Policy>::array_type const &F) {
  auto &world = F.world();

  typename PurificationDensityBuilder<Tile, Policy>::array_type Fp, D, D2;
  Fp("i,j") = M_inv_("i,k") * F("k,l") * M_inv_("j,l");

  auto eig_pair = array_ops::eval_guess(Fp);
  auto emax = eig_pair[1];
  auto emin = eig_pair[0];
  auto scale = 1.0 / (emax - emin);

  D("i,j") = scale * (emax * I_("i,j") - Fp("i,j"));

  auto trace = D("i,j").trace().get();

  auto iter = 0;
  auto error = std::abs(trace - occ_);
  while (error >= 1e-13 && iter <= 100) {
    // Compute D2
    D2("i,j") = D("i,k") * D("k,j");
    if (trace > occ_) {
      D = D2;
    } else {
      D("i,j") = 2 * D("i,j") - D2("i,j");
    }
    D.truncate();
    trace = D("i,j").trace().get();
    error = std::abs(trace - occ_);
    ++iter;
  }

  D("i,j") = M_inv_("i,k") * D("k,l") * M_inv_("l,j");
  D.truncate();
  world.gop.fence();

  return D;
}

template <typename Tile, typename Policy>
typename PurificationDensityBuilder<Tile, Policy>::array_type
PurificationDensityBuilder<Tile, Policy>::orbitals(
    typename PurificationDensityBuilder<Tile, Policy>::array_type const &D) {
  auto D_eig = array_ops::array_to_eigen(D);
  tensor::algebra::piv_cholesky(D_eig);

  auto tr_ao = D.trange().data()[0];
  auto tr_occ = scf::tr_occupied(n_coeff_clusters_, occ_);

  auto Cao =
      array_ops::eigen_to_array<Tile, Policy>(D.world(), D_eig, tr_ao, tr_occ);

  if (localize_) {
    auto U = ((localization_method_ == "rrqr") ||
              (localization_method_ == "rrqr(valence)"))
                 ? mpqc::scf::RRQRLocalization{}(
                       Cao, S_,
                       (localization_method_ == "rrqr(valence)" ? ncore_ : 0))
                 : mpqc::scf::FosterBoysLocalization{}(
                       Cao, r_xyz_ints_,
                       (localization_method_ == "boys-foster(valence)" ? ncore_
                                                                       : 0));
    Cao("mu,i") = Cao("mu,k") * U("k,i");

    auto obs_ntiles = Cao.trange().tiles_range().extent()[0];
    if (clustered_coeffs_) {
      scf::clustered_coeffs(r_xyz_ints_, Cao, obs_ntiles);
    }
  }
#if TA_DEFAULT_POLICY == 1
  if (TcutC_ != 0) {
    minimize_storage(Cao, TcutC_);
  }
#endif
  return Cao;
}

template <typename Tile, typename Policy>
std::pair<typename PurificationDensityBuilder<Tile, Policy>::array_type,
          typename PurificationDensityBuilder<Tile, Policy>::array_type>
PurificationDensityBuilder<Tile, Policy>::operator()(
    typename PurificationDensityBuilder<Tile, Policy>::array_type const &F) {
  auto D = purify(F);
  auto C = orbitals(D);
  return std::make_pair(D, C);
}

}  // namespace scf
}  // namespace mpqc
