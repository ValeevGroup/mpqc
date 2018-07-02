#include "mpqc/chemistry/qc/lcao/scf/orbital_localization.h"

#include <algorithm>
#include <limits>

#include "mpqc/util/keyval/forcelink.h"

namespace mpqc {
namespace lcao {
namespace scf {

/// Foster-Boys maximizes this function
double fb_objective_function(std::array<Mat, 3> const &xyz) {
  auto sum = 0.0;
  for (auto i = 0; i < xyz[0].cols(); ++i) {
    sum += xyz[0](i, i) * xyz[0](i, i);
    sum += xyz[1](i, i) * xyz[1](i, i);
    sum += xyz[2](i, i) * xyz[2](i, i);
  }

  return sum;
}

double compute_angle(double Aij, double Bij) {
  auto AB = std::sqrt(Aij * Aij + Bij * Bij);
  auto cos_4gamma = -Aij / AB;
  auto sin_4gamma = Bij / AB;
  auto gamma = 0.25 * std::acos(cos_4gamma) * ((sin_4gamma < 0) ? -1 : 1);
  return (std::abs(gamma) < 1e-7) ? 0.0 : gamma;
};

bool fb_jacobi_sweeps(Mat &Cm, Mat &U, std::vector<Mat> const &ao_xyz,
                      double convergence_threshold, size_t max_iter) {
  std::array<Mat, 3> mo_xyz;
  mo_xyz[0] = Cm.transpose() * ao_xyz[0] * Cm;
  mo_xyz[1] = Cm.transpose() * ao_xyz[1] * Cm;
  mo_xyz[2] = Cm.transpose() * ao_xyz[2] * Cm;

  auto &mx = mo_xyz[0];
  auto &my = mo_xyz[1];
  auto &mz = mo_xyz[2];

  //auto D = fb_objective_function(mo_xyz);
  decltype(max_iter) iter = 1;
  double max_abs_angle_prev_iter = std::numeric_limits<double>::max();
  double error = max_abs_angle_prev_iter;
  while (error > convergence_threshold && iter <= max_iter) {
    double max_abs_angle = 0.0;
    for (auto i = 0; i < Cm.cols(); ++i) {
      for (auto j = 0; j < i; ++j) {
        Vector3d vij = {mx(i, j), my(i, j), mz(i, j)};
        Vector3d vii = {mx(i, i), my(i, i), mz(i, i)};
        Vector3d vjj = {mx(j, j), my(j, j), mz(j, j)};

        double Aij = vij.squaredNorm() - 0.25 * (vii - vjj).squaredNorm();
        double Bij = (vii - vjj).dot(vij);

        double gamma;
        gamma = compute_angle(Aij, Bij);
        max_abs_angle = std::max(max_abs_angle, std::abs(gamma));
        auto cg = std::cos(gamma);
        auto sg = std::sin(gamma);

        Eigen::VectorXd col_Ui = U.col(i);
        Eigen::VectorXd col_Uj = U.col(j);

        U.col(i) = cg * col_Ui + sg * col_Uj;
        U.col(j) = -sg * col_Ui + cg * col_Uj;

        for (auto z = 0; z < 3; ++z) {
          auto &m = mo_xyz[z];
          Eigen::VectorXd z_i = m.col(i);
          Eigen::VectorXd z_j = m.col(j);

          m.col(i) = cg * z_i + sg * z_j;
          m.col(j) = -sg * z_i + cg * z_j;

          z_i = m.row(i);
          z_j = m.row(j);

          m.row(i) = cg * z_i + sg * z_j;
          m.row(j) = -sg * z_i + cg * z_j;
        }
      }
    }

    //D = fb_objective_function(mo_xyz);
    error = iter > 1 ? std::abs(max_abs_angle - max_abs_angle_prev_iter)
                     : max_abs_angle;
    ++iter;
    max_abs_angle_prev_iter = max_abs_angle;
  }

  return error <= convergence_threshold;
}

}  // namespace scf
}  // namespace lcao
}  // namespace mpqc

#if TA_DEFAULT_POLICY == 0

template class mpqc::lcao::scf::FosterBoysLocalizer<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("FosterBoysLocalizer", mpqc::lcao::scf::FosterBoysLocalizer<TA::TensorD, TA::DensePolicy>);

template class mpqc::lcao::scf::RRQRLocalizer<TA::TensorD, TA::DensePolicy>;
MPQC_CLASS_EXPORT2("RRQRLocalizer", mpqc::lcao::scf::RRQRLocalizer<TA::TensorD, TA::DensePolicy>);

#elif TA_DEFAULT_POLICY == 1

template class mpqc::lcao::scf::FosterBoysLocalizer<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("FosterBoysLocalizer", mpqc::lcao::scf::FosterBoysLocalizer<TA::TensorD, TA::SparsePolicy>);

template class mpqc::lcao::scf::RRQRLocalizer<TA::TensorD, TA::SparsePolicy>;
MPQC_CLASS_EXPORT2("RRQRLocalizer", mpqc::lcao::scf::RRQRLocalizer<TA::TensorD, TA::SparsePolicy>);

#endif
