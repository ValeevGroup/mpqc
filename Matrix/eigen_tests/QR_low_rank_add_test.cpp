#include "../../include/eigen.h"
#include <iostream>
#include <string>
#include <chrono>

using namespace Eigen;

int main(int argc, char **argv) {
  int input = (argc > 1) ? std::stoi(argv[1]) : 500;
  int num_zero_vals = (argc > 2) ? std::stoi(argv[2]) : input / 2;

  MatrixXd mat = MatrixXd::Random(input, input);
  MatrixXd mat2 = MatrixXd::Random(input, input);

  // Make mat symetric
  mat = mat.transpose() * mat;
  mat2 = mat2.transpose() * mat2;

  SelfAdjointEigenSolver<MatrixXd> es(mat);
  MatrixXd C = es.eigenvectors();
  VectorXd V = es.eigenvalues();

  // Create Low Rank Matrix
  std::for_each(V.data(), V.data() + num_zero_vals, [](double &x) { x = 0.0; });
  mat = C * V.asDiagonal() * C.transpose();

  es.compute(mat2);
  C = es.eigenvectors();
  V = es.eigenvalues();

  // Create Low Rank Matrix
  std::for_each(V.data(), V.data() + num_zero_vals, [](double &x) { x = 0.0; });
  mat2 = C * V.asDiagonal() * C.transpose();

  ColPivHouseholderQR<MatrixXd> qr(mat);
  MatrixXd Q = MatrixXd(qr.householderQ()).leftCols(qr.rank());
  MatrixXd R = qr.matrixR()
                   .topLeftCorner(qr.rank(), qr.matrixQR().cols())
                   .template triangularView<Upper>();
  R *= qr.colsPermutation().transpose();
  auto rank1 = qr.rank();
  std::cout << "Rank of mat1 = " << rank1 << std::endl;

  ColPivHouseholderQR<MatrixXd> qr2(mat2);
  MatrixXd Q2 = MatrixXd(qr2.householderQ()).leftCols(qr2.rank());
  MatrixXd R2 = qr2.matrixR()
                    .topLeftCorner(qr2.rank(), qr2.matrixQR().cols())
                    .template triangularView<Upper>();
  R2 *= qr2.colsPermutation().transpose();
  auto rank2 = qr2.rank();
  std::cout << "Rank of mat2 = " << rank2 << std::endl;

  // Create 2 * rank mats for add
  MatrixXd La(mat.rows(), rank1 + rank2);
  MatrixXd Ra(rank1 + rank2, mat.cols());

  // Fill matrices with data
  La.leftCols(rank1) = Q;
  La.rightCols(rank2) = Q2;
  Ra.topRows(rank1) = R;
  Ra.bottomRows(rank2) = R2;

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
  std::cout << "Is add correct ? " << (mat + mat2).isApprox(approx) << "\n\n";

  // Trying something new.
  SelfAdjointEigenSolver<MatrixXd> esR((Ra * Ra.transpose()));
  std::cout << "RaTRa evals = " << esR.eigenvalues().transpose()
            << std::endl;
  MatrixXd evecsR = esR.eigenvectors().rightCols(qr.rank());
  std::cout << "Ra = \n" << Ra << "\n" << std::endl;
  VectorXd evalsR = esR.eigenvalues().bottomRows(qr.rank());
  MatrixXd LR_Ra = evalsR.asDiagonal() * evecsR.transpose() * Ra;
  std::cout << "LR_Ra = \n" << LR_Ra << "\n" << std::endl;

  SelfAdjointEigenSolver<MatrixXd> esL(La.transpose() * La);
  std::cout << "LaTLa evals = " << esL.eigenvalues().transpose()
            << std::endl;
  MatrixXd evecsL = esL.eigenvectors().rightCols(qr.rank());
  VectorXd evalsL = esL.eigenvalues().bottomRows(qr.rank());
  MatrixXd LR_La = La * evecsL;
  std::cout << "La = \n" << La << "\n" << std::endl;
  std::cout << "LR_La = \n" << LR_La << "\n" << std::endl;

  std::cout << "Evecs inner = \n" << evecsL.transpose() * evecsR << std::endl;

  MatrixXd FR_approx = MatrixXd(La * Ra);
  MatrixXd LR_approx = MatrixXd(LR_La * LR_Ra);
  std::cout << "LR - FR = \n" << LR_approx - FR_approx << std::endl;
  std::cout << "Is diff add correct ? " << (FR_approx).isApprox(LR_approx)
            << "\n";


  return 0;
}
