#include "mpqc/math/quadrature/gaussian.h"

namespace mpqc {
namespace math {

void gauss_legendre(int N, Eigen::VectorXd &w, Eigen::VectorXd &x) {
    const double a = 0.0;
    const double b = 1.0;
    Eigen::MatrixXd J;
    J.setZero(N, N);
    for (auto i = 0; i < N; i++) {
      if (i < N - 1) {
        J(i, i + 1) = sqrt(1 / (4 - pow(i + 1, -2)));
      }
    }
    Eigen::MatrixXd Jfin = J + J.transpose();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(Jfin);
    x = es.eigenvalues();
    Eigen::MatrixXd V = es.eigenvectors();

    // Ensure w is the correct size
    w = Eigen::VectorXd::Zero(N);
    for (auto i = 0; i < N; i++) {
      w(i) = 0.5 * 2.0 * (b - a) * V(0, i) * V(0, i);
      x(i) = (b - a) * 0.5 * x(i) + (b + a) * 0.5;
    }
  }

}
}
