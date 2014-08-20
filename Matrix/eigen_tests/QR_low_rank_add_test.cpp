#include "../../include/eigen.h"
#include <iostream>
#include <string>
#include <chrono>

using namespace Eigen;

int main(int argc, char **argv) {
  int input = (argc > 1) ? std::stoi(argv[1]) : 500;
  int num_zero_vals = (argc > 2) ? std::stoi(argv[2]) : input / 2;

  MatrixXd mat = MatrixXd::Random(input, input);

  // Make mat symetric
  mat = mat * mat.transpose();

  SelfAdjointEigenSolver<MatrixXd> es(mat);
  MatrixXd C = es.eigenvectors();
  VectorXd V = es.eigenvalues();

  // Create Low Rank Matrix
  std::for_each(V.data(), V.data() + num_zero_vals, [](double &x) { x = 0.0; });
  mat = C * V.asDiagonal() * C.transpose();

  // Time Matrix Formation
  std::chrono::steady_clock::time_point qr_total0 =
      std::chrono::steady_clock::now();
  ColPivHouseholderQR<MatrixXd> qr(mat);
  MatrixXd Q = MatrixXd(qr.householderQ()).leftCols(qr.rank());
  MatrixXd R = qr.matrixR()
                   .topLeftCorner(qr.rank(), qr.matrixQR().cols())
                   .template triangularView<Upper>();
  R *= qr.colsPermutation().transpose();
  std::chrono::steady_clock::time_point qr_total1 =
      std::chrono::steady_clock::now();
  double qr_time = std::chrono::duration_cast<std::chrono::duration<double>>(
      qr_total1 - qr_total0).count();

  // Create 2 * rank mats for add
  MatrixXd La(mat.rows(), qr.rank() * 2);
  MatrixXd Ra(qr.rank() * 2, mat.cols());

  // Fill matrices with data
  La.leftCols(qr.rank()) = Q;
  La.rightCols(qr.rank()) = Q;
  Ra.topRows(qr.rank()) = R;
  Ra.bottomRows(qr.rank()) = R;

  // Assuming C = A + B, C currently has rank(C) = rank(A) + rank(B), we should
  // be able to reduce this at least a little.

  // Lets try using ColPivQr to reduce the rank of L^C and R^C
  qr.compute(La);
  ColPivHouseholderQR<MatrixXd> qrR(Ra);
  std::cout << "La rank = " << qr.rank() << std::endl;
  std::cout << "Ra rank = " << qrR.rank() << std::endl;

  // At this point just do the same thing we did for multiplication
  MatrixXd LaQ = MatrixXd(qr.householderQ()).leftCols(qr.rank());
  MatrixXd RaQ = MatrixXd(qrR.householderQ()).leftCols(qrR.rank());
  MatrixXd LaR = qr.matrixR()
                     .topLeftCorner(qr.rank(), qr.matrixQR().cols())
                     .template triangularView<Upper>();
  LaR *= qr.colsPermutation().transpose();
  MatrixXd RaR = qrR.matrixR()
                     .topLeftCorner(qrR.rank(), qrR.matrixQR().cols())
                     .template triangularView<Upper>();
  RaR *= qrR.colsPermutation().transpose();

  // New low rank L and low rank R.  Arbitrarrialy picked L to absorb the small
  // matrix.
  MatrixXd FinalL = LaQ * (LaR * RaQ);
  MatrixXd FinalR = RaR;

  // Calculate Full added matrix and check if correct!
  MatrixXd approx = MatrixXd(FinalL * FinalR);
  std::cout << "Is add correct ? " << (mat + mat).isApprox(approx) << "\n";

  return 0;
}
