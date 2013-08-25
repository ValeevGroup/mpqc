#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>

typedef Eigen::MatrixXd Matrix;

template<class T>
void test() {
  Matrix m = Matrix::Random(5, 5);
  Eigen::MatrixXd b = Matrix::Random(5, 5);
  m.row(0) = b.col(0);
  Eigen::SelfAdjointEigenSolver<Matrix> eig(m+m.transpose());
  Matrix m_invsqrt = eig.operatorInverseSqrt();
  std::cout << m_invsqrt << std::endl;
}

int main() {
  test<void>();
}
