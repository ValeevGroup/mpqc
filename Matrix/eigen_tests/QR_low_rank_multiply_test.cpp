#include "../../include/eigen.h"
#include <iostream>
#include <string>
#include <chrono>

using namespace Eigen;

int main(int argc, char** argv){
  int input = (argc > 1) ? std::stoi(argv[1]) : 500;
  int num_sig  = (argc > 2) ? std::stoi(argv[2]) : input/2;

  MatrixXd mat = MatrixXd::Random(input,input);

  // Make symetric
  mat = mat*mat.transpose();

  SelfAdjointEigenSolver<MatrixXd> es(mat);
  MatrixXd C = es.eigenvectors();
  VectorXd V = es.eigenvalues();

  // Delete half the entries
  std::for_each(V.data(), V.data() + num_sig, [](double &x){
    x = 0.0;
  });

  mat = C * V.asDiagonal() * C.transpose();

  std::chrono::steady_clock::time_point qr_total0 = std::chrono::steady_clock::now();
  ColPivHouseholderQR<MatrixXd> qr;
  qr.compute(mat);

  MatrixXd Q = MatrixXd(qr.householderQ()).leftCols(qr.rank());
  MatrixXd R = qr.matrixR().topLeftCorner(qr.rank(), qr.matrixQR().cols()).template triangularView<Upper>();
  R *= qr.colsPermutation().transpose();
  std::chrono::steady_clock::time_point qr_total1 = std::chrono::steady_clock::now();
  double qr_time = std::chrono::duration_cast<std::chrono::duration<double>>(qr_total1-qr_total0).count();

  std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
  MatrixXd final2 = Q * (R * Q) * R;
  std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
  double reduced_time = std::chrono::duration_cast<std::chrono::duration<double>>(t1-t0).count();

  std::chrono::steady_clock::time_point m0 = std::chrono::steady_clock::now();
  MatrixXd mat2 = mat*mat;
  std::chrono::steady_clock::time_point m1 = std::chrono::steady_clock::now();
  double full_time = std::chrono::duration_cast<std::chrono::duration<double>>(m1-m0).count();

  std::cout << input << ", " << qr.rank()
  << ", " << qr_time << ", " << reduced_time << ", " << qr_time+reduced_time << ", "
  << full_time << "\n";

  return 0;
}
